import datetime
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import lattice_single_particle as latt
import fermi_gas as fg

# set up parameters
# ####################################
# real parameters
# ####################################
latt_sp = 750e-9 # m
hbar = 1.054e-34 # Js
h = 6.626e-34 # Js
omega_hz_init = 2 * np.pi * 160 # Hz eigvals 0-58 are in first band
omega_hz_expand = 2 * np.pi * 360 # Hz
m = 9.9883e-27 # kg
Er = hbar**2 * (np.pi / latt_sp)**2 / (2*m) # J
bmap_time = 200e-6 # s
beta = 30 # TODO: implement
mu = 2.3
depth_start = 6.5

# ####################################
# simulation parameters
# ####################################
# Distance is measured in lattice spacing, energy is measured in Er's, and time is measured in hbar / Er

# spatial sampling
n_index = 100
nsites = 30
npts = 40 * nsites
dx = float(nsites) / float(npts)

# time sampling
n_samples = 50
t_end = bmap_time / (hbar / Er) #1.
dt = t_end / n_samples # 0.25
# dt = 0.09

# lattice
depth_end = 0.
latt_depths_fn = lambda t: depth_start * (1 - t / t_end) + (t/ t_end) * depth_end

# initial harmonic trap
osc_len_init = np.sqrt(hbar / (m * omega_hz_init)) / latt_sp
omega_init = 2. / (osc_len_init ** 2 * np.pi ** 2)
period_init = np.array([2 * np.pi / omega_init])
W_init = 1. / (osc_len_init ** 4 * np.pi ** 2) # 0.5 * m * w^2 * a^2 / Er

# harmonic trap for expansion
osc_len_exp = np.sqrt(hbar / (m * omega_hz_expand)) / latt_sp # ~4.
osc_len_fn = lambda t: osc_len_exp * np.ones(t.shape)
# harmonic trap time scales
omega_exp = 2. / (osc_len_exp ** 2 * np.pi ** 2)
period_exp = np.array([2 * np.pi / omega_exp])
# lattice to harmonic trap scales
W_exp = 1. / (osc_len_exp ** 4 * np.pi ** 2) # 0.5 * m * w^2 * a^2 / Er

# ####################################
# solve for real space of initial state
# ####################################
eigvals, eigvects, xpts = latt.get_wavefn_real_space(osc_len_init, latt_depths_fn(0), num_eigs=n_index + 1, nsites=nsites, npts=npts)
thermal_weights = np.diag(fg.fermi(beta, mu, eigvals))
nr_init = np.sum(np.dot( np.abs(eigvects)**2, thermal_weights), axis=1)
N = np.trapz(nr_init, xpts)

# solve lattice problem to determine hopping
v_fourier_components = depth_start * np.array([0, -0.25, -0.25])
v_fourier_indices1 = [0, -1, 1]
v_fourier_indices2 = [0, 0, 0]
recp_basis_vect1 = np.array([[2. * np.pi], [0.]])
recp_basis_vect2 = np.array([[0.], [2. * np.pi]])

eigvals1d, eigvects1d, bz_vects, recp_vects = \
    latt.lattice_single_particle(recp_basis_vect1, recp_basis_vect2, v_fourier_components, v_fourier_indices1,
                            v_fourier_indices2, n1_sites=50, n2_sites=0, n1_recp_max=20, n2_recp_max=0)

# compute wannier functions
xpts1d = np.linspace(-5.0, 5.0, 301)
ypts1d = np.array([0.])
bloch_wfuns1d, wannier_fn, x_grid1d, y_grid1d = latt.get_bloch_wfns(eigvects1d[:, :, 0], bz_vects, recp_vects,
                                                                    xpts1d, ypts1d)

bloch_wfns, _, _, _ = latt.get_bloch_wfns(eigvects1d[:, :, 0], bz_vects, recp_vects, xpts, np.array([0.]))
# bloch wfns appropriately normalized: should integrate to one over a single site
bloch_wfns = bloch_wfns.transpose()

# project eigenstates onto bloch_wfns
kpts2 = bz_vects[0, :]
dk2 = bz_vects[0, 1] - bz_vects[0, 0]
overlaps = np.dot(bloch_wfns.conj().transpose(), eigvects)
nq_numerical = np.sum( np.dot( np.abs(overlaps)**2, thermal_weights), axis=1)
# TODO: why doesn't give me orrect normalization?
nq_numerical = nq_numerical * N / np.trapz(nq_numerical, kpts2)

# computing hopping and interaction
hopping_params, _ = latt.get_hopping_bsfit(bz_vects[0, :], bz_vects[1, :], eigvals1d[:, 0])
tunneling = hopping_params[0]

# calculate quasimomentum distribution from mathieu solution
alpha = 4 * tunneling / W_init
orders = range(0, n_index + 1)

# real space
nr_init_semiclass, xs_semicl, _, mu_sc = latt.get_semiclassical_profile_fixed_num(osc_len_init, tunneling, beta, nparticles=N, nsites=nsites)
#nr_init_semiclass, xs_semicl, _, _ = latt.get_semiclassical_profile(osc_len_init, tunneling, np.inf, 0.135, nsites=nsites)
# js, es, fjs = latt.get_wannier_weights_lattice_plus_harmonic(alpha, orders, num_js)
# wavefns_math, xs_math = latt.get_wavefn_lattice_plus_harmonic(fjs, js, wannier_fn, x_grid1d)

# q-space, Mathieu solution
qsample = np.linspace(-np.pi, np.pi, 200)
psi_q = latt.get_wavefn_kspace_lattice_plus_harmonic(alpha, orders, qsample)
if n_index == 0:
    psi_q = psi_q[:, None]
dens_q = np.sum( np.dot(np.abs(psi_q) ** 2, thermal_weights), axis=1)

# ####################################
# time evolve initial states
# ####################################
psis, times, latt_depths, osc_lens, xpts = latt.solve_time_evolve_real_space(eigvects, osc_len_fn, latt_depths_fn, dt, t_end, nsites=nsites, npts=npts)
if n_index == 0:
    psis = psis[:, :, None]
nr_bmap = np.sum( np.dot(np.abs( np.squeeze(psis[:, -1, :]) )**2, thermal_weights), axis=1)

# ####################################
# T/4 expansion
# ####################################
xpts_extra = np.arange(xpts[-1] + dx, 2 * xpts[-1] + dx, dx)
npad = xpts_extra.size
xpts_expansion = np.concatenate((np.flip(-xpts_extra, axis=0), xpts, xpts_extra), axis=0)
kpts = (1 / osc_len_exp) ** 2 * xpts_expansion
if psis.ndim == 2:
    psi_to_expand = np.concatenate((np.zeros(npad), psis[:, -1], np.zeros(npad)), axis=0)
else:
    psi_to_expand = np.concatenate(( np.zeros((npad, psis.shape[2])), psis[:, -1, :], np.zeros((npad, psis.shape[2])) ), axis=0)

# seems to be messing up the normalization...
psi_final = latt.evolve_harmonic_trap(osc_len_exp, omega_exp, xpts_expansion, psi_to_expand, 0.25 * period_exp)

if n_index == 0:
    psi_final = psi_final[:, None]
dens_final = np.sum( np.dot(np.abs(psi_final) ** 2, thermal_weights), axis=1)
Nf = np.trapz(dens_final, xpts_expansion)
dens_final = dens_final * N/Nf # should of course correct the problem instead of doing this...

# have to integrate this against dk (not dk/(2pi) )
dens_final_toq = dens_final * osc_len_exp ** 2

# also scale to get best fit value
nk_interp_fn = lambda ks: np.interp(ks, kpts, dens_final_toq)
fit_fn = lambda amp:  np.sum( (nq_numerical - amp * nk_interp_fn(kpts2) ) ** 2)

# do fitting
init_guess = np.array([1.0])
fit_handle = scipy.optimize.minimize(fit_fn, init_guess)
fit_params = fit_handle.x

dens_final_toq_scaled = dens_final_toq * fit_params

# max deviation
max_dev = np.max( np.abs( nq_numerical - nk_interp_fn(kpts2) ) )
max_dev_corrected = np.max( np.abs( nq_numerical - fit_params * nk_interp_fn(kpts2)) )

devs = np.divide( np.abs(nq_numerical - fit_params * nk_interp_fn(kpts2)), nq_numerical)

# ####################################
# plot results
# ####################################
figh = plt.figure(figsize=(14,9))

nrows = 3
ncols = 3

plt.subplot(nrows, ncols, 1)
plt.plot(xpts, np.abs(psis[:, 0, 0]), 'entries')
plt.plot(xpts, W_init * xpts ** 2, 'r')
plt.plot(xpts, W_exp * xpts ** 2 , 'r--')
plt.ylim([0, 1])
plt.xlabel('position (lattice sites)')
plt.ylabel('|psi(x)|')
plt.title('psi(t = 0)')
plt.legend(['psi_g', 'V init', 'V expand'])

plt.subplot(nrows, ncols, 2)
plt.plot(xpts, nr_init)
plt.plot(xs_semicl, nr_init_semiclass)
plt.ylim([0, 2.5])
plt.xlabel('position (lattice sites)')
plt.ylabel('n(x, t=0)')
plt.title('n(x, t=0)')
plt.legend(['n(r)', 'semiclassical'])

plt.subplot(nrows, ncols, 3)
plt.plot(times, latt_depths)
plt.xlabel('time (hbar / Er)')
plt.ylabel('lattice depth (Er)')
plt.title('lattice depth vs. time')

#
plt.subplot(nrows, ncols, 4)
plt.plot(xpts, np.abs(psis[:, -1, 0]))
plt.ylim([0, 1])
plt.xlabel('position (lattice sites)')
plt.ylabel('|psi(x)|')
plt.title('psi after bandmap')

# initial distribution vs semiclassical distribution
plt.subplot(nrows, ncols, 5)
plt.plot(xpts, nr_bmap)
plt.plot(xs_semicl, nr_init_semiclass)
plt.ylim([0, 1.2])
plt.xlabel('position (lattice sites)')
plt.ylabel('n(x, t=bmap)')
plt.title('n after bandmap')

# momentum distribution vs. bandmap
plt.subplot(nrows, ncols, 6)
plt.plot(kpts / np.pi, dens_final_toq)
plt.plot(kpts2 / np.pi, nq_numerical)
plt.plot(qsample / np.pi, dens_q)

plt.xlabel('quasimomentum (pi)')
plt.ylabel('n(k)')
plt.title('n(r) and n(k), after T/4 expansion')
plt.legend(['bandmap', 'numerical', 'mathieu thry'])

# eigenstate energies and fermi function
plt.subplot(nrows, ncols, 7)
plt.plot(eigvals, np.diag(thermal_weights), 'entries')
p1, = plt.plot(eigvals, np.diag(thermal_weights), 'bo')
p2, = plt.plot(np.array([mu, mu]), np.array([0, 1]), 'r')

plt.ylim([-0.02, 1.02])
plt.xlabel('energy (Er)')
plt.ylabel('fermi weight')
plt.title('Energy distribution')
plt.legend([p1, p2], ['n_f', 'mu'])


# momentum distributions scaled
plt.subplot(nrows, ncols, 9)
plt.plot(kpts / np.pi, dens_final_toq_scaled)
plt.plot(kpts2 / np.pi, nq_numerical)

plt.xlabel('quasimomentum (pi)')
plt.ylabel('n(k)')
plt.title('n(r) scaled to match n(k), after T/4 expansion')
plt.legend(['bandmap', 'numerical'])

plt.suptitle('Bmap = %.0f us = %0.1f hbar/Er, N = %0.1f, Nstates = %d, dt = %0.2f, dx = %0.2f, beta*Er = %0.2f, mu/Er =%0.2f \n'
             'Er/h = %0.1f KHz, hbar/Er = %.0f us, t/h = %.0f Hz, t / Er = %.3f\n'
             'osc len/a = %.1f, omega_init = 2pi %.0f Hz = %0.3f Er/hbar, T/4 = %.0f us = %.0f hbar/Er\n'
             'osc len/a = %.1f, omega_exp = 2pi %.0f Hz = %0.3f Er/hbar, T/4 = %.0f us = %.0f hbar/Er'
             % (bmap_time / 1e-6, t_end, N, n_index + 1, dt, dx, beta, mu,
                Er / h / 1e3, hbar / Er / 1e-6, tunneling * Er / h, tunneling,
                osc_len_init, omega_hz_init / (2 * np.pi), omega_init, 0.25 * 2 * np.pi / omega_hz_init / 1e-6, 0.25 * period_init,
                osc_len_exp, omega_hz_expand / (2 * np.pi), omega_exp, 0.25 * 2 * np.pi / omega_hz_expand / 1e-6, 0.25 * period_exp))

# ####################################
# save results
# ####################################
now = datetime.datetime.now()
date_str = "%04d;%02d;%02d_%02dh_%02dm" % (now.year, now.month, now.day, now.hour, now.minute)

# save figure
fname = "%s_bandmap.png" % date_str
figh.savefig(fname)

# export data
fname = "%s_bandmap_profile.txt" % date_str
data = np.concatenate((kpts[:, None], dens_final_toq[:, None], dens_final_toq_scaled[:, None]), axis=1)
np.savetxt(fname, data, header='k, density, density  scaled', delimiter=', ')

fname = "%s_momentum_dist.txt" % date_str
data = np.concatenate((kpts2[:, None], nq_numerical[:, None]), axis=1)
np.savetxt(fname, data, header='k, nk bloch waves', delimiter=', ')