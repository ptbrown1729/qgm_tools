import datetime
import numpy as np
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
omega_hz_init = 2 * np.pi * 150 # Hz eigvals 0-58 are in first band
omega_hz_expand = 2 * np.pi * 350 # Hz
m = 9.9883e-27 # kg
Er = hbar**2 * (np.pi / latt_sp)**2 / (2*m) # J
bmap_time = 50e-6 # s
beta = 1.e2 # TODO: implement
mu = 1.95

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
depth_start = 5.
depth_end = 0.
latt_depths_fn = lambda t: depth_start * (1 - t / t_end) + (t/ t_end) * depth_end

# initial harmonic trap
osc_len_init = np.sqrt(hbar / (m * omega_hz_init)) / latt_sp
omega_init = 2. / (osc_len_init ** 2 * np.pi ** 2)
period_init = np.array([2 * np.pi / omega_init])
W_init = 1. / (osc_len_init ** 4 * np.pi ** 2) # 0.5 * m * w^2 * a^2 / Er

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

# compute spectral functions of trapped system
# 58 is number of states in first band for this particular choice of parameters
lk_sqr = np.abs(overlaps[:, 0:58].conj().transpose())**2
lor = lambda w, eta: np.divide(eta / np.pi, eta ** 2 + w ** 2)
ws = np.linspace(1.4, 2.6, 1000)
afns = np.zeros( (ws.size, lk_sqr.shape[1]) )
peak_pos = np.zeros( (lk_sqr.shape[1]))
eta = 0.1
for ii in range(0, lk_sqr.shape[1]):
    for jj in range(0, lk_sqr.shape[0]):
        afns[:, ii] = afns[:, ii] + lk_sqr[jj, ii] * lor(ws - eigvals[jj], eta)
        ind_sort = np.argsort(afns[:, ii])
        peak_pos[ii] = ws[ ind_sort[-1]]

plt.figure
plt.subplot(1, 2, 1)
plt.plot(kpts2, peak_pos, 'o')

plt.subplot(1, 2, 2)
dk = kpts2[1] - kpts2[0]
dw = ws[1] - ws[0]
extent = [(kpts2[0] - dk) / np.pi, (kpts2[-1] + dk)/np.pi, ws[0] - dw, ws[-1] + dw]
plt.imshow(afns, aspect='auto', extent=extent, origin='lower')

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

# q-space
qsample = np.linspace(-np.pi, np.pi, 200)
psi_q = latt.get_wavefn_kspace_lattice_plus_harmonic(alpha, orders, qsample)
if n_index == 0:
    psi_q = psi_q[:, None]
dens_q = np.sum( np.dot(np.abs(psi_q) ** 2, thermal_weights), axis=1)

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

plt.subplot(nrows, ncols, 4)
plt.plot(xpts, np.abs(psis[:, -1, 0]))
plt.ylim([0, 1])
plt.xlabel('position (lattice sites)')
plt.ylabel('|psi(x)|')
plt.title('psi after bandmap')

plt.subplot(nrows, ncols, 5)
plt.plot(xpts, nr_bmap)
plt.plot(xs_semicl, nr_init_semiclass)
plt.ylim([0, 1.2])
plt.xlabel('position (lattice sites)')
plt.ylabel('n(x, t=bmap)')
plt.title('n after bandmap')

plt.subplot(nrows, ncols, 6)
plt.plot(kpts / np.pi, dens_final_toq * 0.644)
plt.plot(qsample / np.pi, dens_q)
plt.plot(kpts2 / np.pi, nq_numerical * np.max(dens_q) / np.max(nq_numerical))

plt.xlabel('quasimomentum (pi)')
plt.ylabel('n(k)')
plt.title('n(r) and n(k), after T/4 expansion')
plt.legend(['bandmap', 'thry', 'numerical'])

plt.subplot(nrows, ncols, 7)
plt.plot(eigvals, np.diag(thermal_weights), 'entries')
plt.plot(eigvals, np.diag(thermal_weights), 'bo')
plt.plot(np.array([mu, mu]), np.array([0, 1]), 'r')

plt.ylim([-0.02, 1.02])
plt.xlabel('energy (Er)')
plt.ylabel('fermi weight')
plt.title('Energy distribution')
plt.legend(['n_f', 'mu'])


# plt.suptitle('Bmap = %.0f us_interspecies = %0.1f hbar/Er, N = %0.1f, Nstates = %d, dt = %0.2f, dx = %0.2f, beta*Er = %0.2f, mu/Er =%0.2f \n'
#              'Er/h = %0.1f KHz, hbar/Er = %.0f us_interspecies, t/h = %.0f Hz, t / Er = %.3f\n'
#              'osc len/a = %.1f, omega_init = 2pi %.0f Hz = %0.3f Er/hbar, T/4 = %.0f us_interspecies = %.0f hbar/Er\n'
#              'osc len/a = %.1f, omega_exp = 2pi %.0f Hz = %0.3f Er/hbar, T/4 = %.0f us_interspecies = %.0f hbar/Er'
#              % (bmap_time / 1e-6, t_end, N, n_index + 1, dt, dx, beta, mu,
#                 Er / h / 1e3, hbar / Er / 1e-6, tunneling * Er / h, tunneling,
#                 osc_len_init, omega_hz_init / (2 * np.pi), omega_init, 0.25 * 2 * np.pi / omega_hz_init / 1e-6, 0.25 * period_init,
#                 osc_len_exp, omega_hz_expand / (2 * np.pi), omega_exp, 0.25 * 2 * np.pi / omega_hz_expand / 1e-6, 0.25 * period_exp))

# ####################################
# save results
# ####################################
now = datetime.datetime.now()
# fname = "%04d;%02d;%02d_%02dh_%02dm_bandmap.png" % (now.year, now.month, now.day, now.hour, now.minute)
# figh_corr.savefig(fname)
