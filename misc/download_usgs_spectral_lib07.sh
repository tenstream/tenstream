#!/bin/bash
set -eu -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT="$(readlink -f $SCRIPTDIR/../)"
DST=$PROJECT_ROOT/external/USGS
mkdir -p $DST
cd $DST

USGS_ZIP=usgs_splib07.zip

if [[ ! -e $USGS_ZIP ]]; then
  echo ""
  echo " Could not find file $DST/$USGS_ZIP"
  echo ""
  echo " To download zip from the official source, go to "
  echo "  - https://dx.doi.org/10.5066/F7RR1WDJ"
  echo ""
  echo " and download the file"
  echo "  - usgs_splib07.zip"
  echo ""
  exit 1
fi

unzip -u $USGS_ZIP

cat > read_ascii_files.py << EOF
#!/usr/bin/env python3
from bs4 import BeautifulSoup as BS
import copy
import numpy as np
import random
import re

# indices in table
descr_col = 0
wvl_col   = 4
spec_col  = 2

with open('indexes/datatable_splib07a.html') as fh:
    html = BS(fh.read())

# Use all materials
trs = html.find_all("tr")

# or use only entries after Veg Chapter
veg_chapter = html.find(text=re.compile("Chapter 7: Vegetation")).find_parent("tr")
trs = veg_chapter.find_all_next("tr")

# cache entries for wvl tables
wvl_tables = {}

# aggregate everything in result
res = {}

for row, tr in enumerate(trs[:]):
    tds = tr.find_all('td')
    if len(tds) == 14:
        description = tds[descr_col].text
        wvl_file    = 'indexes/' + tds[wvl_col].a.get("href")
        spec_file   = 'indexes/' + tds[spec_col].a.get("href")
        print(f"Reading Material: {description}")

        if wvl_file not in wvl_tables:
            print(f"Parsing wvl table: {wvl_file}")
            wvl_tables[wvl_file] = np.loadtxt(wvl_file, comments=('!','#', 'splib'))

        res[description] = {
            'lambda': wvl_tables[wvl_file],
            'albedo': np.loadtxt(spec_file, comments=('!','#', 'splib')),
        }

# remove NaN values
for k, mat in res.items():
    valid_entries = mat['albedo'] > 0
    res[k]['albedo'] = mat['albedo'][valid_entries]
    res[k]['lambda'] = mat['lambda'][valid_entries]

# try to simplify results, i.e. if linear interpolation gives the same values, drop it
wgt    = lambda a0, a1, x: (x-a0) / (a1-a0)
spline = lambda x0, x1, x, y0, y1: y0 + wgt(x0, x1, x) * (y1 - y0)

for k, mat in res.items():
    print(f"Pruning Material {k}")
    mat_orig = copy.deepcopy(mat)
    while True:
        changed = False
        iterator = list(range(1, len(mat['albedo'])-1))
        random.shuffle(iterator)

        for i in iterator:
            l0, l, l1 = [ mat['lambda'][_] for _ in (i-1,i,i+1) ]
            a0, a, a1 = [ mat['albedo'][_] for _ in (i-1,i,i+1) ]
            orig_a = np.interp(l, mat_orig['lambda'], mat_orig['albedo'])
            if abs(orig_a - spline(l0, l1, l, a0, a1)) < 5e-3:
                #print(f"Dropping {i} of {len(mat['albedo'])} because original albedo @ {l} = {orig_a}",
                #      f" and we now get {spline(l0, l1, l, a0, a1)}",
                #      f" difference: {orig_a - spline(l0, l1, l, a0, a1)}")
                mat['albedo'] = np.delete(mat['albedo'], i)
                mat['lambda'] = np.delete(mat['lambda'], i)
                changed = True
                break
        if not changed:
            break


outfile = 'usgs_material_albedi.txt'
with open(outfile, 'w') as fh:
    fh.write("# 3 lines each, first the material name, then the wavelength [micrometer], i.e. [1e-6 m] and third the albedo\n")
    for k, mat in res.items():
        wvl    = ' '.join(map(str, mat['lambda'].round(decimals=3)))
        albedo = ' '.join(map(str, mat['albedo'].round(decimals=3)))
        fh.write(k+'\n')
        fh.write(wvl+'\n')
        fh.write(albedo+'\n')
EOF

python3 read_ascii_files.py
