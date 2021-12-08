import numpy as np
import xml.etree.ElementTree as ET

rmse  = lambda a,b: np.sqrt(np.mean((a-b)**2))
rrmse = lambda a,b: rmse(a,b) / np.mean(b) * 100
bias  = lambda a,b: (np.mean(a) / np.mean(b) -1) * 100

def read_data(fname):
    r = {}
    doc = ET.parse(fname)
    for a in doc.iter('Attribute'):
        r[a.attrib['Name']] = np.array(a[0].text.split()).astype(float)
    return r

def main(f1, f2):
    d1, d2 = [ read_data(f) for f in (f1,f2) ]

    for k in d1.keys():
        a, b = [ _[k] for _ in (d1,d2) ]
        print(f"{k:20} : rmse {rmse(a,b):12.4f} {rrmse(a,b):8.1f}% bias {bias(a,b):8.2f}%")


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Compare two xmf outputs as found in tenstream xmf files')
    parser.add_argument('f1', type=str, help='first xmf file')
    parser.add_argument('f2', type=str, help='second xmf file, if suitable, should be the benchmark/truth')
    args = parser.parse_args()

    main(args.f1, args.f2)
