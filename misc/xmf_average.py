import numpy as np
import glob
import xml.etree.ElementTree as ET

def read_data(fname):
    print("Reading", fname)
    r = {}
    doc = ET.parse(fname)
    for a in doc.iter('Attribute'):
        r[a.attrib['Name']] = np.array(a[0].text.split()).astype(float)
    return r

def main(basepattern, output):
    files = glob.glob(basepattern)

    data = [ read_data(f) for f in files ]

    # build average of data
    avg = {}
    for k in data[0].keys():
      avg[k] = np.average([d[k] for d in data],axis=0)

    # embed avg results into structure of first file
    doc = ET.parse(files[0])
    for a in doc.iter('Attribute'):
        k = a.attrib['Name']
        a[0].text = ' '.join(map(str, avg[k]))

    # and write output
    print("Writing to", output)
    with open(output,'w') as fh:
        fh.write(ET.tostring(doc.getroot()).decode('utf8'))


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Average variables found in tenstream xmf files')
    parser.add_argument('basepattern', type=str, help='base filename pattern')
    parser.add_argument('output',default='avg.xmf', type=str, help='output file')
    args = parser.parse_args()

    main(args.basepattern, args.output)
