# test time evolution solver by comparing results to exactly solvable simple harmonic oscillator problem
import datetime
import numpy as np
import matplotlib.pyplot as plt
import lattice_single_particle as latt

# spatial sampling
nsites = 30
npts = 3000
dx = float(nsites) / float(npts)
# time sampling
dt = 5.
t_end = 500.
ts = np.arange(0, t_end + dt, dt)
# oscillator lengths
osc_len_start = 4.
osc_len_end = osc_len_start
osc_len_fn = lambda t: 4 * np.ones(t.shape)
# oscillator time scales
omega_start = 2. / (osc_len_start ** 2 * np.pi ** 2)
period_start = 2 * np.pi / omega_start
# lattice depths
latt_depths_fn = lambda t: 0. * t

# solve for eigenstate
# eigvals, eigvects, xpts = latt.get_wavefn_real_space(osc_len_start, latt_depths_fn(0), num_eigs=1, nsites=nsites, npts=npts)
eigvals, eigvects, xpts = latt.get_wavefn_real_space(osc_len_start, 10., num_eigs=1, nsites=nsites, npts=npts)
# psi_init = eigvects

# initial SHO wavefns
psi_sho = latt.get_sho_wavefn(osc_len_start, xpts, 0)
psi_sho_1 = latt.get_sho_wavefn(osc_len_start, xpts, 1)

# final SHO wavefns
psi_sho_final = latt.get_sho_wavefn(osc_len_end, xpts, 0)
psi_sho_final_1 = latt.get_sho_wavefn(osc_len_end, xpts, 1)

# initial wavefunction
psi_init = psi_sho / np.sqrt(2) + psi_sho_1 / np.sqrt(2)
psi_final = psi_sho / np.sqrt(2) + psi_sho_1 / np.sqrt(2) * np.exp(-1j * omega_start * t_end)


psis_sho = latt.evolve_harmonic_trap(osc_len_start, omega_start, xpts, psi_init, ts)

## solve time-dependent problem
psis, times, latt_depths, osc_lens, xpts = latt.solve_time_evolve_real_space(psi_init, osc_len_fn, latt_depths_fn, dt, t_end, nsites=nsites, npts=npts)

## plot results
figh = plt.figure()

nplots = 9
step = np.floor(len(times) / (nplots - 1.))
indices = range(0, len(times), int(step))
indices[-1] = len(times) - 1

for ii in range(0, nplots):

    plt.subplot(2, 5, ii + 1)
    plt.plot(xpts, np.abs(psis[:, indices[ii]]), 'blue')
    plt.plot(xpts, np.abs(psis_sho[:, indices[ii]]), 'orange')
    plt.xlabel('position (lattice sites)')
    plt.ylabel('|psi(x)|')
    plt.title('psi vs. time')
    plt.legend(['numerical', 'exact'])
    plt.title('time = %0.2fs' % times[indices[ii]])

plt.subplot(2, 5, nplots + 1)
plt.plot(times, osc_lens)
plt.xlabel('time')
plt.ylabel('oscillator length')

# plt.legend(['t=0', 'tend', 'expected end'])

plt.suptitle('dt = %0.2f, dx = %0.2f' % (dt, dx))

# save fig
now = datetime.datetime.now()
fname = "%04d;%02d;%02d_%02dh_%02dm_sho_time_evolve.png" % (now.year, now.month, now.day, now.hour, now.minute)
figh.savefig(fname)