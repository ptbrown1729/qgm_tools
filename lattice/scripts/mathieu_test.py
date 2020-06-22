import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import lattice_single_particle as latt

# all energies in units of Er

# set lattice depth
latt_depth = 2 # in Er

# get Wannier function

# get wannier functions and hopping for lattice of this depth
v_fourier_components = latt_depth * np.array([0, -0.25, -0.25])
# v_fourier_components = depth * np.array([0.5, -0.25, -0.25])
v_fourier_indices1 = [0, -1, 1]
v_fourier_indices2 = [0, 0, 0]
recp_basis_vect1 = np.array([[2. * np.pi], [0.]])
recp_basis_vect2 = np.array([[0.], [2. * np.pi]])

# diagonalize lattice problem
eigvals1d, eigvects1d, bz_vects1d, recp_vects1d = \
    latt.lattice_single_particle(recp_basis_vect1, recp_basis_vect2, v_fourier_components, v_fourier_indices1,
                            v_fourier_indices2, n1_sites=31, n2_sites=0, n1_recp_max=20, n2_recp_max=0)

# compute wannier functions
xpts1d = np.linspace(-5.0, 5.0, 301)
ypts1d = np.array([0.])
bloch_wfuns1d, wannier_fn_gridded1d, x_grid1d, y_grid1d = latt.get_bloch_wfns(eigvects1d[:, :, 0], bz_vects1d, recp_vects1d,
                                                                         xpts1d, ypts1d)

# computing hopping and interaction from wannier function
# hopping1d, interaction1d = latt.get_hopping_wannier(wannier_fn, x_grid1d, y_grid1d)
hopping_params, _ = latt.get_hopping_bsfit(bz_vects1d[0, :], bz_vects1d[1, :], eigvals1d[:, 0])
tunneling = hopping_params[0]

# W_exp = 0.5 * m * w^2 * d^2 / Er (i.e. in units of Er)
# Note that W_exp/Er = (a/l)^4 / np.pi^2
W = tunneling / 50.
# so we get the oscillator length in units of the lattice spacing as
osc_len = (1. / (W * np.pi ** 2) ) ** 0.25
# ratio of 1D bandwidth to tunneling
alpha = 4 * tunneling / W

orders = 20
num_js = 30
js, es, fjs = latt.get_wannier_weights_lattice_plus_harmonic(alpha, orders, num_js)

wave_fns, xs_all = latt.get_wavefn_lattice_plus_harmonic(fjs, js, wannier_fn_gridded1d, x_grid1d)
sho_wave_fns = latt.get_sho_wavefn(osc_len, xs_all, range(0, fjs.shape[1]) )

eigvals, eigvects, xpts = latt.get_wavefn_real_space(osc_len, latt_depth, num_eigs=fjs.shape[1], npts=5000)

# solve density distributions
beta = 10
n_particles = 19

n_exact = np.sum( np.abs(wave_fns[:, 0 : n_particles]) ** 2, axis=1)
n_sho = np.sum( np.abs(sho_wave_fns[:, 0 : n_particles]) ** 2, axis=1)
n_realsp = np.sum( np.abs(eigvects[:, 0 : n_particles]) ** 2, axis=1)

# _,  = scipy.signal.butter(2, 0.1)

#nr_semiclass, xx = latt.get_semiclassical_profile(osc_len_start, tunneling, beta, 0, nsites=30)
n_sc, x_sc, _, mu_sc = latt.get_semiclassical_profile_fixed_num(osc_len, tunneling, beta, n_particles, nsites=30)

plt.figure()
nrows = 2
ncols = 2

plot_index = 0
plt.subplot(nrows, ncols, 1)
# plt.plot(xs_all, wave_fns, '.')
plt.plot(xs_all, np.abs(wave_fns[:, plot_index]))
plt.plot(xs_all, np.abs(sho_wave_fns[:, plot_index]))
plt.plot(xpts, np.abs(eigvects[:, plot_index]))
plt.legend(['exact', 'sho', 'real-space'])
plt.xlabel('position (sites)')
plt.title('eigenstate')


plt.subplot(nrows, ncols, 2)
plt.plot(x_grid1d.ravel(), wannier_fn_gridded1d.ravel())
plt.plot(np.array([np.min(x_grid1d), np.max(x_grid1d)]), np.array([0,0]), 'k')
plt.plot(np.array([0.5, 0.5]), np.array([0, 2]), 'k--')
plt.plot(np.array([-0.5, -0.5]), np.array([0, 2]), 'k--')
plt.plot(np.array([1.5, 1.5]), np.array([0, 2]), 'k--')
plt.plot(np.array([-1.5, -1.5]), np.array([0, 2]), 'k--')
plt.xlabel('position (sites)')
plt.title('wannier fn')

plt.subplot(nrows, ncols, 3)
plt.plot(xs_all, n_exact)
plt.plot(xs_all, n_sho)
plt.plot(xpts, n_realsp)
plt.plot(x_sc, n_sc)
plt.legend(['exact', 'sho', 'real-space', 'semiclassical'])
plt.xlabel('position (sites)')
plt.title('density distribution')