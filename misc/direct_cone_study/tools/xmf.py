"""read surface fluxes from pprts -pprts_xdmf output (inline XML data) into global 2D arrays [x, y]"""
import sys
import re
import xml.etree.ElementTree as ET
import numpy as np


def read(fname):
    txt = open(fname).read()
    txt = re.sub(r"<!DOCTYPE[^>]*>", "", txt)
    root = ET.fromstring(txt)
    cells = {}  # (xc, yc) -> {var: val}
    for grid in root.iter("Grid"):
        if grid.get("GridType") == "Collection":
            continue
        geo = grid.find("Geometry/DataItem")
        nver = [int(v) for v in geo.get("Dimensions").split()]
        xyz = np.array(geo.text.split(), dtype=float).reshape(nver)
        for att in grid.findall("Attribute"):
            item = att.find("DataItem")
            dims = [int(v) for v in item.get("Dimensions").split()]
            nx, ny = dims  # Fortran order, x is fastest running in the data
            data = np.array(item.text.split(), dtype=float).reshape(ny, nx)
            verts = xyz.reshape(ny + 1, nx + 1, 3)
            xc = 0.5 * (verts[:-1, :-1, 0] + verts[1:, 1:, 0])
            yc = 0.5 * (verts[:-1, :-1, 1] + verts[1:, 1:, 1])
            for j in range(ny):
                for i in range(nx):
                    cells.setdefault((xc[j, i], yc[j, i]), {})[att.get("Name")] = data[j, i]
    xs = sorted({k[0] for k in cells})
    ys = sorted({k[1] for k in cells})
    names = sorted({n for v in cells.values() for n in v})
    out = {n: np.full((len(xs), len(ys)), np.nan) for n in names}
    for (x, y), v in cells.items():
        for n, val in v.items():
            out[n][xs.index(x), ys.index(y)] = val
    return out


if __name__ == "__main__":
    a = read(sys.argv[1])
    for n, v in a.items():
        print(n, v.shape, "mean", np.nanmean(v), "min", np.nanmin(v), "max", np.nanmax(v))
