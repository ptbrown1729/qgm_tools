import datetime
import numpy as np
import matplotlib.pyplot as plt
import lattice_single_particle as latt

# set up parameters
# real parameters
latt_sp = 750e-9
hbar = 1.054e-34
h = 6.626e-34
omega_hz = 2*np.pi * 350
m = 1.0e-26
Er = hbar**2 * (np.pi / latt_sp)**2 / (2*m)
bmap_time = 200e-6

# simulation parameters
# Distance is measured in lattice spacing, energy is measured in Er's, and time is measured in hbar / Er

# spatial sampling
n_index = 200
nsites = 100
npts = 40 * nsites
dx = float(nsites) / float(npts)
# time sampling
n_samples = 50
t_end = bmap_time / (hbar / Er) #1.
dt = t_end / n_samples # 0.25
# harmonic trap
osc_len = np.sqrt(hbar / (m * omega_hz)) / latt_sp# 4.
osc_len_fn = lambda t: osc_len * np.ones(t.shape)
# harmonic trap time scales
omega = 2. / (osc_len ** 2 * np.pi ** 2)
period = np.array([2 * np.pi / omega])
# lattice
depth_start = 5.
depth_end = 0.
latt_depths_fn = lambda t: depth_start * (1 - t / t_end) + (t/ t_end) * depth_end
# lattice to harmonic trap scales
W = 1. / (osc_len ** 4 * np.pi ** 2)


# solve for eigenstates
eigvals, eigvects, xpts = latt.get_wavefn_real_space(osc_len, latt_depths_fn(0), num_eigs=n_index+1, nsites=nsites, npts=npts)
_, eigvects_end, _ = latt.get_wavefn_real_space(osc_len, latt_depths_fn(t_end), num_eigs=n_index+1, nsites=nsites, npts=npts)
eigvects_end = eigvects_end[:, n_index]

# time evolve eigenstates
# psi_init = eigvects[:, n_index]
psis, times, latt_depths, osc_lens, xpts = latt.solve_time_evolve_real_space(eigvects, osc_len_fn, latt_depths_fn, dt, t_end, nsites=nsites, npts=npts)
if n_index == 0:
    psis = psis[:, :, None]

# T/4 expansion
xpts_extra = np.arange(xpts[-1] + dx, 2 * xpts[-1] + dx, dx)
npad = xpts_extra.size
xpts_expansion = np.concatenate((np.flip(-xpts_extra, axis=0), xpts, xpts_extra), axis=0)
kpts = (1 / osc_len)**2 * xpts_expansion
if psis.ndim == 2:
    psi_to_expand = np.concatenate((np.zeros(npad), psis[:, -1], np.zeros(npad)), axis=0)
else:
    psi_to_expand = np.concatenate(( np.zeros((npad, psis.shape[2])), psis[:, -1, :], np.zeros((npad, psis.shape[2])) ), axis=0)

psi_final = latt.evolve_harmonic_trap(osc_len, omega, xpts_expansion, psi_to_expand, 0.25 * period)
if n_index == 0:
    psi_final = psi_final[:, None]
dens_final = np.sum(np.abs(psi_final) ** 2, axis=1)
# TODO: why is this not giving me the right height?
dens_final_toq = dens_final * (2 * np.pi) * osc_len ** 2

# get expected q-space distribution
v_fourier_components = depth_start * np.array([0, -0.25, -0.25])
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
hopping_params, _ = latt.get_hopping_bsfit(bz_vects1d[0, :], bz_vects1d[1, :], eigvals1d[:, 0])
tunneling = hopping_params[0]

# calculate quasimomentum distribution
alpha = 4 * tunneling / W
qsample = np.linspace(-np.pi, np.pi, 1000)
psi_q = latt.get_wavefn_kspace_lattice_plus_harmonic(alpha, range(0, n_index + 1), qsample)
if n_index == 0:
    psi_q = psi_q[:, None]
dens_q = np.sum(np.abs(psi_q) ** 2, axis=1)

# plot results
figh = plt.figure(figsize=(14,9))

plt.subplot(2, 2, 1)
plt.plot(xpts, np.abs(psis[:, 0, :]))
plt.xlabel('position (lattice sites)')
plt.ylabel('|psi(x)|')
plt.title('psi(t = 0)')

plt.subplot(2, 2, 2)
plt.plot(times, latt_depths)
plt.xlabel('time (hbar / Er)')
plt.ylabel('lattice depth (Er)')
plt.title('lattice depth vs. time')

plt.subplot(2, 2, 3)
plt.plot(xpts, np.abs(psis[:, -1, :]))
# plt.plot(xpts, np.abs(eigvects_end), 'green')
plt.xlabel('position (lattice sites)')
plt.ylabel('|psi(x)|')
plt.title('psi after bandmap')
# plt.legend(['time evolve', 'ground state'])

plt.subplot(2, 2, 4)
# plt.plot(kpts / np.pi, dens_final)
# plt.plot(qsample / np.pi, dens_q * np.max(dens_final) / np.max(dens_q))
plt.plot(kpts/ np.pi, dens_final_toq * 0.65)
plt.plot(qsample / np.pi, dens_q)

plt.xlabel('quasimomentum (pi)')
plt.ylabel('n(k)')
plt.title('n(r) and n(k), after T/4 expansion')
plt.legend(['bandmap', 'thry'])

plt.suptitle('N = %d, dt = %0.2f, dx = %0.2f Bmap = %.0f us_interspecies = %0.1f hbar/Er\n Er/h = %0.1f KHz, hbar/Er = %.0f us_interspecies, t/h = %.0f Hz, t / Er = %.3f \n osc len/a = %.1f, omega_exp = 2pi %.0f Hz = %0.3f Er/hbar, T/4 = %.0f us_interspecies = %.0f hbar/Er'
             % (n_index+1, dt, dx, bmap_time/1e-6, t_end, Er/h/1e3, hbar/Er/1e-6, tunneling * Er / h, tunneling, osc_len, omega_hz/(2*np.pi), omega, 0.25 * 2*np.pi / omega_hz/1e-6, 0.25 * period))

now = datetime.datetime.now()
fname = "%04d;%02d;%02d_%02dh_%02dm_bandmap.png" % (now.year, now.month, now.day, now.hour, now.minute)
figh.savefig(fname)
