import numpy as np
import scipy.special
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import lattice_single_particle as latt
import mathieu_test


# measure distance in units of the lattice constant, a
# measure energy in terms of the recoial energy
# For a 1D lattice potential of the form V(x) = V_o/2 * cos(kx), the lattice spacing is a = (2pi)/k. We define the
#     recoil energy as the potential energy for wavevector halfway across the Brillouin zone. i.e.
# Er = hbar^2 * (k/2)^2/(2m) = hbar^2 *k^2 / 8m

# k = 2 * pi / a
# R = 0.5 * m * w^2 * a^2
omega = 2 * np.pi * 420 # Hz
m = 9.988345515692e-27 # mass of 6-li in kg
a = 750e-9 # lattice spacing, meters
hbar = 1.054e-34 # Js
er = hbar ** 2 * ( 2 * np.pi / a) ** 2 / (8 * m ) # lattice recoil in J

# set lattice depth
latt_depth = 2 # in Er

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
                            v_fourier_indices2, n1_sites=11, n2_sites=0, n1_recp_max=6, n2_recp_max=0)

# compute wannier functions
xpts1d = np.linspace(-2.0, 2.0, 101)
ypts1d = np.array([0.])
bloch_wfuns1d, wannier_fn_gridded1d, x_grid1d, y_grid1d = latt.get_bloch_wfns(eigvects1d[:, :, 0], bz_vects1d, recp_vects1d,
                                                                         xpts1d, ypts1d)

# computing hopping and interaction from wannier function
# hopping1d, interaction1d = latt.get_hopping_wannier(wannier_fn, x_grid1d, y_grid1d)
hopping_params, _ = latt.get_hopping_bsfit(bz_vects1d[0, :], bz_vects1d[1, :], eigvals1d[:, 0])
tx = hopping_params[0]

# set harmonic trap params
R = 0.5 * m * omega ** 2 * a ** 2 / er
# alpha = 1D band width / trap energy
alpha = 4 * tx / R
v = lambda x: latt_depth * 0.5 * (1 - np.cos( 2 * np.pi * x)) + R * x ** 2

nsites = 100
npts = nsites * 20
xpts = np.linspace(-nsites, nsites, npts)
dx = xpts[1] - xpts[0]

# build hamiltonian
vpart = sp.diags( v(xpts) )
dx_sqr = (1 / dx ** 2) * ( -2. * sp.identity(npts)
                           + sp.diags( np.ones(npts-1), offsets=1, shape=(npts, npts))
                           + sp.diags( np.ones(npts-1), offsets=-1, shape=(npts, npts)) )
# hbar^2/2m / Er = a ** 2 / pi ** 2
ham = vpart - 1 / np.pi ** 2 * dx_sqr

eigvals, eigvects = eigsh(ham, 50, which='SM')

plt.figure()

plt.subplot(2, 2, 1)
plt.plot(xpts, v(xpts))
plt.xlabel('position (sites)')
plt.title('potential')

plt.subplot(2, 2, 2)
plt.plot( np.ones(eigvals.size), eigvals, 'o' )
plt.title('energies')

plt.subplot(2, 2, 3)
plt.plot(xpts, np.abs(eigvects[:, 0]) ** 2)
plt.xlabel('position (sites)')
plt.ylabel('|psi|^2')
plt.title("ground state")

plt.subplot(2, 2, 4)
plt.plot( np.ravel(x_grid1d), np.ravel(np.abs(wannier_fn_gridded1d)) ** 2)
plt.xlabel('position (sites)')
plt.ylabel('|w|^2')
plt.title('wannier function')

plt.suptitle('depth = %0.2f Er, alpha = %0.2f, t / Er = %0.2f' % (latt_depth, alpha, tx))