from pylab import *
import py_eddington as P

# Define the benchmark grid... we will use eddington coeffs on that highres grid to find suitable LUT supports
Ntau, Nw0, Ng, Nmu = 200, 50, 20, 90
tau = logspace(-10, 2, Ntau)
w0 = linspace(0, .99999, Nw0)
g = linspace(0, .5, Ng)

theta = np.linspace(0, 90, Nmu)
mu = np.cos(np.deg2rad(theta))

def compute_eddington_data(tau, w0, g, mu):
    """ compute eddington coefficients on the provided grid """
    a = np.zeros((7, tau.size, w0.size, g.size, mu.size))

    for itau, vtau in enumerate(tau):
        # print("Computing tau={:12.8e} ({:3.0f}%)".format(vtau,100.*itau/(tau.size-1.)))
        for iw0, vw0 in enumerate(w0):
            for ig, vg in enumerate(g):
                for imu, vmu in enumerate(mu):
                    a[:, itau, iw0, ig, imu] = P.m_py_eddington.py_eddington_coeff_zdun(vtau, vw0, vg, vmu)

    return a


def find_supports(c, N, dim):
    """ find support points that minimize the curvature along the remaining dimension with about N points or less """
    reduc_dims = tuple([x for x in xrange(len(c.shape)) if x != dim])

    dc = np.gradient(c)[dim]
    dc = np.mean(dc, axis=reduc_dims)

    d2c = np.abs(np.gradient(dc))**.5

    curvature = np.cumsum(d2c)/ np.sum(d2c)
    support_indices = np.interp(np.linspace(0,1,N), curvature, np.arange(c.shape[dim]))

    support_indices = np.hstack(([0,c.shape[dim]-1],support_indices))

    return np.sort(np.unique(support_indices))

def find_new_tau(a, N=20, iw0=int(Nw0//1.1), ig=Ng//2, imu=Nmu//2, lplot=True, linear_supports=False):
    """ find new support points for tau axis """
    a11, a12, a13, a23, a33, b1, b2 = a

    #new_tau_a11 = tau[find_supports(a11, 20, 0).astype(int)]
    #new_tau_a12 = tau[find_supports(a12, 20, 0).astype(int)]
    #new_tau_a33 = tau[find_supports(a33, 20, 0).astype(int)]
    #support_indices = np.hstack([find_supports(c, 20, 0) for c in (a11, a12,a13,a23,a33)]).astype(int)

    support_indices = find_supports(a11, N, 0)
    if linear_supports:
        support_indices = np.linspace(0, tau.size, N)

    new_tau = np.interp(support_indices, np.arange(tau.size), tau)

    #new_tau = tau[find_supports(a11, N, 0).astype(int)]

    new_a11 = [ P.m_py_eddington.py_eddington_coeff_zdun(vtau, w0[iw0], g[ig], mu[imu])[0] for vtau in new_tau ]
    new_a12 = [ P.m_py_eddington.py_eddington_coeff_zdun(vtau, w0[iw0], g[ig], mu[imu])[1] for vtau in new_tau ]
    new_a33 = [ P.m_py_eddington.py_eddington_coeff_zdun(vtau, w0[iw0], g[ig], mu[imu])[4] for vtau in new_tau ]

    if lplot:
        figure(1); clf()
        plot(tau, a11[:, iw0, ig, imu], 'r.-', label='a11'); plot(new_tau, new_a11, 'ro--')
        plot(tau, a12[:, iw0, ig, imu], 'g.-', label='a12'); plot(new_tau, new_a12, 'go--')
        plot(tau, a33[:, iw0, ig, imu], 'b.-', label='a33'); plot(new_tau, new_a33, 'bo--')
        xlabel('tau'); xscale('log'); tight_layout(); legend()

    return new_tau


def find_new_w0(a, N=10, itau=int(Ntau//1.2), ig=Ng//2, imu=Nmu//2, lplot=False, linear_supports=False):
    a11, a12, a13, a23, a33, b1, b2 = a

    support_indices = find_supports(a11, N, 1)
    if linear_supports:
        support_indices = np.linspace(0, w0.size, N)

    new_w0 = np.interp(support_indices, np.arange(w0.size), w0)
    #new_w0 = w0[find_supports(a11, N, 1).astype(int)]

    new_a11 = [ P.m_py_eddington.py_eddington_coeff_zdun(tau[itau], vw0, g[ig], mu[imu])[0] for vw0 in new_w0 ]
    new_a12 = [ P.m_py_eddington.py_eddington_coeff_zdun(tau[itau], vw0, g[ig], mu[imu])[1] for vw0 in new_w0 ]
    new_a33 = [ P.m_py_eddington.py_eddington_coeff_zdun(tau[itau], vw0, g[ig], mu[imu])[4] for vw0 in new_w0 ]

    if lplot:
        figure(2); clf()
        plot(w0, a11[itau, :, ig, imu], 'r.-', label='a11'); plot(new_w0, new_a11, 'ro--')
        plot(w0, a12[itau, :, ig, imu], 'g.-', label='a12'); plot(new_w0, new_a12, 'go--')
        plot(w0, a33[itau, :, ig, imu], 'b.-', label='a33'); plot(new_w0, new_a33, 'bo--')
        xlabel('w0'); tight_layout(); legend()

    return new_w0


def find_new_g(a, N=3, itau=int(Ntau//1.2), iw0=-1, imu=Nmu//2, lplot=False, linear_supports=False):
    a11, a12, a13, a23, a33, b1, b2 = a

    support_indices = find_supports(a11, N, 2)
    if linear_supports:
        support_indices = np.linspace(0, g.size, N)

    new_g = np.interp(support_indices, np.arange(g.size), g)
    #new_g = g[find_supports(a11, N, 2).astype(int)]

    new_a11 = [ P.m_py_eddington.py_eddington_coeff_zdun(tau[itau], w0[iw0], vg, mu[imu])[0] for vg in new_g ]
    new_a12 = [ P.m_py_eddington.py_eddington_coeff_zdun(tau[itau], w0[iw0], vg, mu[imu])[1] for vg in new_g ]
    new_a33 = [ P.m_py_eddington.py_eddington_coeff_zdun(tau[itau], w0[iw0], vg, mu[imu])[4] for vg in new_g ]

    if lplot:
        figure(3); clf()
        plot(g, a11[itau, iw0, :, imu], 'r.-', label='a11'); plot(new_g, new_a11, 'ro--')
        plot(g, a12[itau, iw0, :, imu], 'g.-', label='a12'); plot(new_g, new_a12, 'go--')
        plot(g, a33[itau, iw0, :, imu], 'b.-', label='a33'); plot(new_g, new_a33, 'bo--')
        xlabel('g'); tight_layout(); legend()

    return new_g


def check_results(a, tau, w0, g, mu, na, ntau, nw0, ng, nmu, ival=0):
    """ evaluate new dimensions by linearly interpolating from new grid onto the highres grid and compare """
    from scipy import interpolate
    me = np.vstack([ m.flatten() for m in meshgrid(tau, w0, g, mu, indexing='ij') ])

    ia = interpolate.interpn((ntau, nw0, ng, nmu), na[ival], me.T, bounds_error=False)

    rmse = lambda a,b: np.sqrt(np.nanmean((a-b)**2))
    rel_rmse = lambda a,b: rmse(a,b)/np.nanmean(b)*1e2
    maxerr = lambda a,b: np.nanmax(np.abs(a-b))
    bias = lambda a,b: (np.mean(a)/np.mean(b) - 1.)*1e2

    error = (rel_rmse(ia, a[ival].flatten()), maxerr(ia, a[ival].flatten()), bias(ia, a[ival].flatten()))

    return ia.reshape((tau.size, w0.size, g.size, mu.size)), error


def print_fortran_code(ntau, nw0, ng, ntheta):
    """ print copy and paste output for tenstream src code """
    print(""" ! Copy the following code into src/optprop_parameters.f90
    logical, parameter :: use_prescribed_LUT_dims=.True.
    """)
    print("integer(iintegers), parameter :: Ntau={}, Nw0={}, Ng={}, Ntheta={}".format(ntau.size, nw0.size, ng.size, ntheta.size))
    print("real(ireals), parameter :: preset_tau({}) = [{}]".format(ntau.size, ','.join([str(x) for x in ntau])))
    print("real(ireals), parameter :: preset_w0({}) = [{}]".format(nw0.size, ','.join([str(x) for x in nw0])))
    print("real(ireals), parameter :: preset_g({}) = [{}]".format(ng.size, ','.join([str(x) for x in ng])))
    print("real(ireals), parameter :: preset_theta({}) = [{}]".format(ntheta.size, ','.join([str(x) for x in ntheta])))


def _main(Ntau=20, Nw0=15, Ng=3, Ntheta=19, lplot=False):
    """ run a set of new dimensions """
    a = compute_eddington_data(tau, w0, g, mu)

    kwargs = {
            'lplot' : lplot,
            'linear_supports' : False,
            }

    ntau = find_new_tau(a, N=Ntau, **kwargs)
    nw0  = find_new_w0 (a, N=Nw0 , **kwargs)
    ng   = find_new_g  (a, N=Ng  , **kwargs)
    ntheta = np.linspace(0, 90, Ntheta)
    nmu = np.sort(np.cos(np.deg2rad(ntheta)))
    print("New Indices:",ntau.size, nw0.size, ng.size, nmu.size)

    na = compute_eddington_data(ntau, nw0, ng, nmu)
    ia, err = check_results(a, tau, w0, g, mu, na, ntau, nw0, ng, nmu)
    print(Ntau, Nw0, Ng, Ntheta, "RMSE & max-error with linear interpolation:", err)

    if lplot:
        figure(4); clf()
        imshow(ia[:,:,-1,-1], aspect='auto',interpolation='nearest')
        cb = colorbar(); cb.set_label('interpolated coeffs')
        xlabel('w0 index'); ylabel('tau index')
        title('RMSE {:.2f} % & max-error {:.4f} & bias {:.2f} %'.format(err[0], err[1], err[2]))

    print_fortran_code(ntau, nw0, ng, ntheta)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='LUT support point finder')
    parser.add_argument('Ntau'  , type=int, help='Number of supports in dimension tau')
    parser.add_argument('Nw0'   , type=int, help='Number of supports in dimension w0')
    parser.add_argument('Ng'    , type=int, help='Number of supports in dimension g')
    parser.add_argument('Ntheta', type=int, help='Number of supports in dimension zenith theta')
    parser.add_argument('-plot', action='store_const', const=True, help='create some plots to check', default=False)
    args = parser.parse_args()

    _main(args.Ntau, args.Nw0, args.Ng, args.Ntheta, args.plot)

