import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
import os

outfile = sys.argv[1]
files = sys.argv[2:]

if os.path.exists(outfile):
    raise Exception(f"outfile {outfile} already exists")

res = {}

for f in sorted(files):
    print(f"Reading {f}")
    with open(f) as fh:
        lines = fh.read()

    tsw = None
    tlw = None

    try:
        tsw = re.search(": tenstr_repwvl.repwvl_solar: (.*?) ", lines).groups()[0]
    except:
        pass
    try:
        tlw = re.search(": tenstr_repwvl.repwvl_thermal: (.*?) ", lines).groups()[0]
    except:
        pass

    try:
        tsw = re.search(": pprts_rrtmg_solar: (.*?) ", lines).groups()[0]
    except:
        pass
    try:
        tlw = re.search(": pprts_rrtmg_thermal: (.*?) ", lines).groups()[0]
    except:
        pass


    label=os.path.basename(f).split(".")[0]
    print(f"{label} {tsw} {tlw}")

    if label not in res.keys():
        res[label] = {}

    if tsw:
        res[label]['sw'] = float(tsw)
    if tlw:
        res[label]['lw'] = float(tlw)

print(res)

plt.figure(figsize=(8,8))
plt.clf()
w=.2
hoff = 150
voff = w/4
for i, (l,v) in enumerate(res.items()):
    print(i,l,v)
    if 'sw' in v:
        plt.barh(i - w, v["sw"]          , height=w                           , align="edge")
        plt.text(v["sw"] +hoff         , i - w + voff, f'sw {int(v["sw"])}s'             , fontsize='small')
    if 'lw' in v:
        plt.barh(i    , v["lw"]          , height=w                           , align="edge")
        plt.text(v["lw"] + hoff        , i + voff   , f'lw {int(v["lw"])}s'             , fontsize='small')
    if ('sw' in v) and ('lw' in v):
        plt.barh(i+w  , v["sw"] + v["lw"], height=w                           , align="edge")
        plt.text(v["sw"] + v["lw"] + hoff, i+w+ voff , f'total {int(v["sw"] + v["lw"])}s', fontsize='small')

labels = [_ for _ in res.keys()]
plt.gca().set_yticks([ _ + w/2 for _ in range(len(labels))])
plt.gca().set_yticklabels(labels)
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x//3600))
plt.gca().xaxis.set_major_formatter(ticks)
plt.xlabel('runtime [hrs]')

print(f"Plotting to {outfile}")
plt.savefig(outfile, bbox_inches='tight')
