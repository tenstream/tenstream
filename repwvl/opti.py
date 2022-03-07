from scipy import optimize as O
from scipy.ndimage import gaussian_filter as gf
from pylab import *

Nwvl = 100
wvl = np.linspace(4,200,Nwvl)

rng = np.random.default_rng(1234)
func = lambda wvl: sin(wvl)**2 * exp(wvl/wvl.max())
func = lambda wvl: gf(rng.random(len(wvl)) * wvl * exp(-wvl/wvl.max()*2),2)
benchmark = func(wvl)
sum_benchmark = np.sum(benchmark)

def f0(pt):
    pts = np.sort(pt)
    idx = np.floor(pts).astype(int)
    #print(pts, idx)
    mid_idx = (idx[:-1] + idx[1:]) / 2
    lb = np.hstack(([0,],mid_idx))
    rb = np.hstack((mid_idx, [Nwvl,]))
    wgt = (pts - idx) * (rb-lb)
    y = np.sum(benchmark[idx] * wgt)
    return y

R = []

def f(pt):
    y = f0(pt)
    c = np.abs(y - sum_benchmark)
    R.append(pt)
    print(f'f @ {pt} : {y} : {c}')
    return c

N = 2
bounds = np.tile(np.array([0.0, Nwvl - 1e-6]), (N,1))
x0 = np.random.rand(N) * (Nwvl+1)

ret = O.dual_annealing(f, bounds, x0=x0)
