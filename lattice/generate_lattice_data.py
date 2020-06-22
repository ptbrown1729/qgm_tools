import numpy as np
import lattice_single_particle as latt

# lattice parameters
fname_out = "C:\\Users\\Peter\\Documents\\MATLAB\\PeterB-analysis\\lattice\\lattice_data2.txt"
theta = 91.6267
alpha = 0.
# paramters to loop over
#depths = np.concatenate((np.arange(1, 20, 2), np.arange(20, 100, 5)))
depths = np.array([1, 3, 5, 7, 10])
rs = np.arange(1, 0, -0.1)

labels = ['Depth (Er)', 'r', 'theta (deg)', 'alpha (deg)',
          'tx wannier (Er)', 'ty wannier (Er)', 'td wannier (Er)', 'wannier integral (Er)',
          'tx bs (Er)', 'ty bs (Er)', 'td bs(Er)',
          'mean 0->d1', 'gap 0->d1', 'max 0->d1', 'bw d1',
          'mean 0->d2', 'gap 0->d2', 'max 0->d2', 'bw d2',
          'mean 0->d3', 'gap 0->d3', 'max 0->d3', 'bw d3']

for ii in range(0, rs.size):
    hoppings_bs, ts_bw, hoppings_wannier, hoppings_wannier_dens_dep, wannier_integral, energies, wannier_fn, potential, xpts, ypts = \
        latt.solve_d4_lattice(depths, theta, rs[ii], alpha,
                              dx=0.1, n1_sites=12, n2_sites=12, n1_recp_max=10, n2_recp_max=10, normalize_er=True)

    gband = energies[:, :, 0, :]

    d_dat = np.zeros((depths.size, 4 * 3))
    for jj, bi in enumerate([3, 4, 5]):
        d_es = energies[:, :, bi, :]
        gd_mean = np.mean(np.mean(d_es - gband, axis=1), axis=0)
        gd_gap = np.min(np.min(d_es - gband, axis=1), axis=0)
        gd_max = np.max(np.max(d_es -gband, axis=1), axis=0)
        d_bw = np.max(np.max(d_es, axis=1), axis=0) - np.min(np.min(d_es, axis=1), axis=0)

        d_dat[:, 4 * jj : 4 * jj + 4] =  np.concatenate((gd_mean[:, None], gd_gap[:, None], gd_max[:, None], d_bw[:, None]), axis=1)

    data_curr_r = np.concatenate((depths[:, None], rs[ii] * np.ones((depths.size, 1)), theta * np.ones((depths.size, 1)), alpha * np.ones((depths.size, 1)),
                                  hoppings_wannier[:, 0:3], np.array([wannier_integral]).transpose(), hoppings_bs[:, 0:3], d_dat), axis=1)
    if ii == 0:
        data = data_curr_r
    else:
        data = np.concatenate((data, data_curr_r), axis=0)

# save data
# TODO: do this in loop ... to avoid data loss
header = labels[0]
for ii in range(1, len(labels)):
    header = header + '\t%s' % labels[ii]
np.savetxt(fname_out, data, delimiter='\t', header=header)



