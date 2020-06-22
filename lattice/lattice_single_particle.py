"""
Tools for simulating various features of 1D and 2D optical lattice potentials and optical lattices in
the presence of harmonic trapping potentials, and single-particle Hamiltonians.
Useful for computing band structures, tightbinding parameters, etc.
"""

from __future__ import print_function
import math
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import numpy.linalg
import scipy.optimize
import scipy.special
import scipy.integrate
import scipy.interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import fermi_gas as fg

# solve lattice problem in momentum space
def lattice_single_particle(b1, b2, v_fourier_components, fourier_inds1, fourier_inds2,
                            n1_sites=21, n2_sites=21, n1_recp_max=10, n2_recp_max=10, print_results=1):
    """
    Solve Schrodinger equation for a single particle in a periodic potential. Energies are returned in terms of the
    lattice recoil.

    Suppose we have a lattice with vectors a_1 and a_2 and reciprocal vectors b_1 and b_2 satisfying
    a_i.dot(b_j) = 2*pi*delta_{ij}. The reciprocal lattice vectors are then rec_{n, m} = n * b_1 + m * b_2. We only
    keep a finite number of reciprocal lattice vectors for the numerical calculation, so we choose a cutoffs n_c and m_c
    and restrict -n_c <= n <= n_c, -m_c <= m <= m_c. We also work with a finite lattice of N_1 and N_2 sites in each
    direction. This implies a finite number of vectors in the Brillouin zone, which are given by b_i * n/N_i for
    n = 0, ... N_i - 1.

    We work in units of the lattice recoil energy and distance units are the same as a_1 and a_2 are given in

    For a 1D lattice potential of the form V(x) = Vo/2 * cos(kx), the lattice spacing is a = (2pi)/k, and
    the reciprocal lattice vector is entries = k
    We define the recoil energy as
    Er = hbar^2 * (entries/2)^2 / (2m)

    For a 2D lattice, the lengths of the reciprocal vectors may not be the same, therefore we use
    entries = sqrt(|b1| * |b2|)
    Er = hbar^2 * (entries/2)^2 / (2m),
    where b1 and b2 are reciprocal vectors of the lattice.

    For example,
    b1 = np.array([2 * np.pi, 0])
    b2 = np.array([0, 2 * np.pi])
    fourier_inds1 = [0, 1, -1, 0,  0]
    v_fourier_inds2 = [0, 0,  0, 1, -1]
    v_fourier_components = [0, 0.25, 0.25, 0.25, 0.25]

    :param b1: Reciprocal basis vector, b1
    :param b2: Second reciprocal basis vector, b2
    :param v_fourier_components: A list of Fourier components of the lattice potential
    :param fourier_inds1: A list of indices which helps to specify the K-vectors corresponding to the Fourier
    components. If the given indiex are n1, n2, then this has K = n1*b1 + n2*b2
    :param fourier_inds2: Indices n2 corresponding to the Fourier components and n2*b2
    :param n1_sites: Number of sites to include in calculation along the dirction of the first lattice
    basis vector a1. This is the number of wavectors in the Brillouin zone along direction a1.
    :param n2_sites: Number of sites to include in the calculation along the second lattice basis vector a2
    :param n1_recp_max: The total number of reciprocal lattice vectors which are multiples of b1 is
    (2 * n1_recp_max)
    :param n2_recp_max:

    :return eigvals: A N_bz x N_recp array of eigenvalues. eigvals[ii, jj] is an eigenvalue in the jjth band associated with
    the iith vector in the Brillouin zone.
    :return eigvects: A N_bz x N_recp x N_recp array of eigenvectors written in a momentum space basis.
    eigvects[ii, jj, kk] is the jjth component of the eigenvector in the kkth band associated with the iith vector in
    the Brillouin zone. This eigenvector has eigenvalue eigvals[ii, kk]. It is associated with the wavevector
    bz_vects[:, ii] - recp_vects[:, jj]. The bloch wavefunction can be obtained from this as
    Bloch[ii, kk](r) = \sum_{jj} exp(ij * (bz_vects[:, ii] - recp_vects[:, jj]) * r) * eigvects[ii, jj, kk]
    :return bz_vects: A 2 x N_bz array of Brillouin zone vectors
    :return recp_vects: A 2 x N_recp array of reciprocal lattice vectors
    """

    t_start = time.time()

    ####################################
    # checks inputs
    ####################################
    b1 = np.array(b1)
    if b1.size == 2:
        b1 = np.reshape(b1, [2, 1])
    if b1.shape != (2, 1):
        raise Exception('b1 was not of shape (2,1)')

    b2 = np.array(b2)
    if b2.size == 2:
        b2 = np.reshape(b2, [2, 1])
    if b2.shape != (2, 1):
        raise Exception('b2 was not of shape (2,1)')

    v_fourier_components = np.asarray(v_fourier_components)
    if v_fourier_components.ndim != 1:
        raise Exception('Expected v_fourier_components of dimension 1')

    fourier_inds1 = np.asarray(fourier_inds1)
    if fourier_inds1.ndim != 1:
        raise Exception('Expected fourier_inds1 of dimension 1.')

    fourier_inds2 = np.asarray(fourier_inds2)
    if fourier_inds2.ndim != 1:
        raise Exception('Expected fourier_inds2 of dimension 1.')

    if not (v_fourier_components.size == fourier_inds1.size and fourier_inds1.size == fourier_inds2.size):
        raise Exception('v_fourier_components, fourier_inds1, and fourier_inds2 must all have the same number of elements.')

    ####################################
    # define energy units
    ####################################
    e_recoil = np.sqrt(np.sum(np.square(b1)) * np.sum(np.square(b2))) / 4.

    # total number of reciprocal lattice sites
    n_recp1 = (2 * n1_recp_max + 1)
    n_recp2 = (2 * n2_recp_max + 1)
    n_recp_total = n_recp1 * n_recp2

    # construct all reciprocal lattice vectors and store as a list.
    b1b1_inds, b2b2_inds = np.meshgrid(range(-n1_recp_max, n1_recp_max + 1),
                                       range(-n2_recp_max, n2_recp_max + 1))
    # convert 2D to 1D ordering.
    b1_inds = np.reshape(b1b1_inds, [1, b1b1_inds.size])
    b2_inds = np.reshape(b2b2_inds, [1, b2b2_inds.size])

    # recp vectors list. A collection of column vectors
    # note that this accounts for different spacing between the lattice directions
    recp_vects = np.kron(b1_inds, b1) + np.kron(b2_inds, b2)

    ####################################
    # construct Vk, using same ordering for k as in recp_vects
    ####################################
    # first remove any components of weight zero
    inds_to_use = v_fourier_components != 0
    v_fourier_components = v_fourier_components[inds_to_use]
    fourier_inds1 = fourier_inds1[inds_to_use]
    fourier_inds2 = fourier_inds2[inds_to_use]

    vk = np.zeros((n_recp_total))
    for ii in range(0, v_fourier_components.size):
        # typically if we want to map between a[i,j] and a[k] we would say k = i + j * n (if indexing starts at zero and ends at n)
        # here, we need to runner_offset our i and j because they are symmetric, positive and negative
        index = (fourier_inds1[ii] + n1_recp_max) + (fourier_inds2[ii] + n2_recp_max) * n_recp1
        vk[index] = v_fourier_components[ii]

    # construct Vkq
    # as above, convert from fourier components to single index
    nxs, nys = np.meshgrid(b1_inds, b1_inds)
    mxs, mys = np.meshgrid(b2_inds, b2_inds)
    b1_diff_inds = nys - nxs
    b2_diff_inds = mys - mxs
    single_indices = (b1_diff_inds + n1_recp_max) + (b2_diff_inds + n2_recp_max) * n_recp1

    # assign values to vkq
    exclude_indices = np.logical_or(np.abs(b1_diff_inds) > n1_recp_max, np.abs(b2_diff_inds) > n2_recp_max)
    single_indices[exclude_indices] = 0

    vkq = vk[single_indices]
    vkq[exclude_indices] = 0

    ####################################
    # construct vectors in brillouin zone
    ####################################
    # k =  [0, ..., n_sites - 1] / n_sites
    # But we prefer to write these as symmetricly about zero as possible
    # k = [-floor((n_sites - 1) / 2), ... , ceil((nsites - 1) / 2)]
    # if n_sites is odd this reduces to:
    # k=  [-(n_sites - 1) / 2, ..., 0, .. (n_sites - 1) / 2]

    if n1_sites != 0:
        sites1_range = np.arange(-np.floor((n1_sites - 1) / 2.), np.ceil((n1_sites - 1) / 2.) + 1)
    else:
        sites1_range = np.array([0.])

    if n2_sites != 0:
        sites2_range = np.arange(-np.floor((n2_sites - 1) / 2.), np.ceil((n2_sites -1) / 2.) + 1)
    else:
        sites2_range = np.array([0.])

    nbz, mbz = np.meshgrid(sites1_range, sites2_range)
    nbz_flat = np.reshape(nbz, [1, nbz.size])
    mbz_flat = np.reshape(mbz, [1, mbz.size])

    nbz_vects = nbz_flat.size

    # handle situation where we have a 1D lattice
    if n1_sites != 0:
        bz_basis_vect1 = b1 / n1_sites
    else:
        bz_basis_vect1 = 0.

    if n2_sites != 0:
        bz_basis_vect2 = b2 / n2_sites
    else:
        bz_basis_vect2 = 0

    # allowed vectors in the Brillouin zone as a 2 x n array
    bz_vects = np.kron(nbz_flat, bz_basis_vect1) + np.kron(mbz_flat, bz_basis_vect2)

    eigvals = np.zeros((bz_vects.shape[1], n_recp_total))
    eigvects = np.zeros((bz_vects.shape[1], n_recp_total, n_recp_total))

    for ii in range(0, nbz_vects):
        bz_vect = bz_vects[:, ii][:, None]
        vects = np.kron(np.ones((1, n_recp_total)), bz_vect) - recp_vects
        epsilons = np.diag( np.sum( vects ** 2,  0) )
        # put in units of Er
        epsilons = epsilons / e_recoil

        eigvals[ii, :], eigvects[ii, :, :] = np.linalg.eigh(epsilons + vkq)

    t_end = time.time()
    if print_results:
        print("Solving bandstructure, diagonalizing %d %dx%d matrices took %0.2fs"
          % (bz_vects.shape[1], epsilons.shape[0], epsilons.shape[0], t_end - t_start))

    # TODO: maybe should reshape these to be 2D etc.

    return eigvals, eigvects, bz_vects, recp_vects

def get_bloch_wfns(eigvects_band, bz_vects, recp_vects, xpts, ypts, vmin_position=np.array([0, 0]), print_results=1):
    """
    Compute bloch wavefunctions and wannier function at real space points given the solution to the lattice single-particle
    problem. Most of the necessary arguments to this function are generated using the lattice_single_particle function.

    :param eigvects: An N_bz x N_recp x N_recp array of eigenvectors written in a momentum space basis.
    eigvects[ii, jj, kk] is the jjth component of the eigenvector in the kkth band associated with the iith vector in
    the Brillouin zone. This eigenvector has eigenvalue eigvals[ii, kk]. It is associated with the wavevector
    bz_vects[:, ii] - recp_vects[:, jj]. The bloch wavefunction can be obtained from this as
    Bloch[ii, kk](r) = \sum_{jj} exp(ij * (bz_vects[:, ii] - recp_vects[:, jj]) * r) * eigvects[ii, jj, kk]'''
    :param bz_vects:
    :param recp_vects:
    :param xpts: A NumPy array of xpts. With ypts, this is used create a grid on which the Wannier function is computed
    :param ypts: A NumPy array of ypts.
    :param band_index:
    :return bloch_wfuns:
    :return wannier_fn:
    :return x_grid:
    :return y_grid:
    """

    # TODO: could save a lot of time by using only symmetric part of function

    t_start = time.time()
    # calculate overlaps between neighboring bloch states. Note that these eigenvectors are real (in this momentum
    # space representation) by construction
    overlaps = np.sum(eigvects_band[0:-1, :] * eigvects_band[1:], 1)

    # do neighboring bloch states have a consistent phase?
    neighbor_flip = np.power(-1., overlaps < 0)

    # build global flip array
    nbz_vects = bz_vects.shape[1]
    global_flip = np.zeros([nbz_vects])
    global_flip[0] = 1.
    for ii in range(1, nbz_vects):
        global_flip[ii] = global_flip[ii - 1] * neighbor_flip[ii - 1]
    global_flip_matrix = np.diag(global_flip)

    # left multiplying by this diagonal matrix changes the sign of each row of eigvects_band
    eigvects_band_gaugefixed = global_flip_matrix.dot(eigvects_band)

    # choose real space points to compute Bloch and Wannier functions
    # change coordinates to put minimum of potential at (0,0)
    # we rely on the fact that if we shift our coordinates so the potential mininmum is at (xo, yo), then
    # the new wannier function w' is related to the old by
    #w'_{R + ro}(r-ro) = w_R(r)
    xpts = xpts - vmin_position[0]
    ypts = ypts - vmin_position[1]
    R = np.array([0., 0.]) + vmin_position

    # coordinate grids
    x_grid, y_grid = np.meshgrid(xpts, ypts)
    x_flat = np.reshape(x_grid, [x_grid.size])
    y_flat = np.reshape(y_grid, [y_grid.size])

    # for small number of real-space pts, vectorize fully. For large number of real-space pts we will run out of
    # memory, so in this case use a loop
    # TODO: having two different implementations is a little scary because errors could be very hard to debug
    max_gb_per_array = 0.5
    if recp_vects.shape[1] * bz_vects.shape[1] * x_flat.size * 8 / 1e9 < max_gb_per_array:
        # [ii, jj, r] = [bz_vects, recp_vects, r]
        recp_vects_x, bz_vects_x, x_vals = np.meshgrid(recp_vects[0, :], bz_vects[0, :], x_flat)
        recp_vects_y, bz_vects_y, y_vals = np.meshgrid(recp_vects[1, :], bz_vects[1, :], y_flat)
        exp_factors = np.exp(1j * (bz_vects_x - recp_vects_x) * x_vals +
                             1j * (bz_vects_y - recp_vects_y) * y_vals) * \
                      np.exp(-1j * R[0] * bz_vects_x
                             -1j * R[1] * bz_vects_y)

        eigvects_band_repeated = np.repeat(eigvects_band_gaugefixed[:, :, None], x_flat.size, 2)
        bloch_wfuns = np.sum(exp_factors * eigvects_band_repeated, 1)
    else:
        bloch_wfuns = np.zeros( (bz_vects.shape[1], x_flat.size), dtype=np.complex )
        for ii in range(0, x_flat.size):
            x = x_flat[ii]
            y = y_flat[ii]
            recp_vects_x, bz_vects_x = np.meshgrid(recp_vects[0, :], bz_vects[0, :])
            recp_vects_y, bz_vects_y = np.meshgrid(recp_vects[1, :], bz_vects[1, :])
            exp_factors = np.exp(1j * (bz_vects_x - recp_vects_x) * x +
                                 1j * (bz_vects_y - recp_vects_y) * y) * \
                          np.exp(-1j * R[0] * bz_vects_x
                                 - 1j * R[1] * bz_vects_y)
            bloch_wfuns[:, ii] = np.sum( eigvects_band_gaugefixed * exp_factors, 1)

    # wannier function
    wannier_fn = (np.sum(bloch_wfuns, 0)) / np.sqrt(bloch_wfuns.shape[0])

    # reshape wannier function to have 2d grid like shape
    wannier_fn = np.reshape(wannier_fn, x_grid.shape)

    # normalize wannier function
    wannier_probability = np.abs(wannier_fn) ** 2
    if ypts.size > 1:
        y_int = np.trapz(wannier_probability, ypts, axis=1)
    else:
        y_int = (wannier_probability).flatten()

    norm_sqr = np.trapz(y_int, xpts, axis=0)
    wannier_fn = wannier_fn / np.sqrt(norm_sqr)

    t_end = time.time()
    if print_results:
        print("Computing real-space Bloch functions and wannier functions at %d space points took %0.2fs"
          % (x_grid.size, t_end - t_start))

    # change coordinates back to original before output
    x_grid = x_grid + vmin_position[0]
    y_grid = y_grid + vmin_position[1]

    return bloch_wfuns, wannier_fn, x_grid, y_grid

def get_potential_realspace(recp_basis_vect1, recp_basis_vect2, v_fourier_components,
                            v_fourier_indices1, v_fourier_indices2, xpts, ypts):
    """
    Return the potential on a grid of real space points given its k-space representation.

    Given x_j x-positions and y_gi y-positions, the potential is returned as the matrix V[i,j] = V(x_j, y_i)

    :param recp_basis_vect1:
    :param recp_basis_vect2:
    :param v_fourier_components: List of Fourier components
    :param v_fourier_indices1: List of
    :param v_fourier_indices2:
    :param xpts: 1D array of xpts.
    :param ypts: 1D array of ypts.
    :return v_real_space: a ypts.size x xpts.size matrix, where V[i, j] = V(x_j, y_i)
    """
    # create grid of real space points
    x_grid, y_grid = np.meshgrid(xpts, ypts)
    x_flat = x_grid.flatten()
    y_flat = y_grid.flatten()

    # Get all relevant reciprocal vectors as a 2 x N array
    recp_vects = np.kron(v_fourier_indices1, recp_basis_vect1) + np.kron(v_fourier_indices2, recp_basis_vect2)

    # Get list of all possible pairs of reciprocal vectors with (x,y) coordinates
    x_vals, recp_vects_x = np.meshgrid(x_flat, recp_vects[0, :])
    y_vals, recp_vects_y = np.meshgrid(y_flat, recp_vects[1, :])

    # compute phase factor for Fourier transform
    exp_factors = np.exp(1j * (-recp_vects_x) * x_vals + 1j * (-recp_vects_y) * y_vals)

    # reshape v_fourier_components to match size of exp_factors
    v_fourier_components_repeated = np.repeat(v_fourier_components[:, None], x_flat.size, 1)

    # fourier transform
    v_real_space = np.sum(exp_factors * v_fourier_components_repeated, 0)
    # reshape to be 2D array
    v_real_space = np.reshape(v_real_space, [ypts.size, xpts.size])

    if np.any( np.round(np.imag(v_real_space), 12) > 0 ):
        raise Exception("Expected potential to be real, but it has non-zero complex part")

    v_real_space = np.real(v_real_space)

    return v_real_space

# get hopping, interaction, and other parameters
def get_hopping_bw(gband_2d, mode='2d'):
    """
    Estimate hopping parameters from band structure by taking differences between points in the band
    :param gband_2d: num_ky x num_kx array of energies for the ground band
    :return hoppings:
    """

    center_index = int((gband_2d.shape[0] - 1) / 2)
    t = 0.125 * (np.max(np.ravel(gband_2d)) - np.min(np.ravel(gband_2d)))

    if mode == '1d':
        hoppings = np.array([t])
    elif mode =='2d':
        tx = 0.125 * (gband_2d[center_index, -1] - gband_2d[center_index, center_index] + gband_2d[-1, -1] - gband_2d[-1, center_index])
        ty = 0.125 * (gband_2d[-1, center_index] - gband_2d[center_index, center_index] + gband_2d[-1, -1] - gband_2d[center_index, -1])
        td = 1./16 * (gband_2d[center_index, -1] + gband_2d[-1, center_index] - gband_2d[center_index, center_index] - gband_2d[-1, -1])
        hoppings = np.array([t, tx, ty, td])
    else:
        raise Exception("mode should be '1d' or '2d', but it was %s.", mode)

    return hoppings

def get_hopping_bsfit(kxs_2d, kys_2d, energies_2d, ax=1, ay=1):
    """
    Fit to tight binding band structure to determine tunn, ty, and td
    :param kxs_2d: kxs arranged in a 2d grid, as produced by oned_2_two_kvects
    :param kys_2d:
    :param energies_2d: band energies arranged in a 2d grid, as prodced by oned_2_twod_kvects
    :param ax: lattice constant along x-direction
    :param ay:
    :return fit_params = [tunn, ty, td, tx_nn, ty_nn, runner_offset]
    Here tunn and ty are nearest neighbor hopping in the x and y directions. td is diagonal hopping along the x+y and x-y
    directions. tx_nn and ty_nn are next-nearest-neighbor hopping (two lattice sites in the x or y directions).


    """

    # functions to be used in fitting
    # params = [tunn, ty, td, tx_nn, ty_nn, runner_offset]

    # handle 1D vs. 2D case. Keep parameter order the same for all, but for 1D don't use some terms
    is_1d_x = kxs_2d.shape[0] == 1 and np.all( kys_2d == 0)
    is_1d_y = kxs_2d.shape[1] == 1 and np.all( kxs_2d == 0)

    if is_1d_x:
        dispersion_fn = lambda params, kx, ky: params[0] * (2 - 2 * np.cos(kx * ax)) + \
                                               params[3] * (2 - 2 * np.cos(2 * kx * ax)) + \
                                               params[5]
    elif is_1d_y:
        dispersion_fn = lambda params, kx, ky: params[1] * (2 - 2 * np.cos(ky * ay)) + \
                                               params[4] * (2 - 2 * np.cos(2 * ky * ay)) + \
                                               params[5]
    else:
        dispersion_fn = lambda params, kx, ky: params[0] * (2 - 2 * np.cos(kx * ax)) + \
                                               params[1] * (2 - 2 * np.cos(ky * ay)) + \
                                               params[2] * (2 - 2 * np.cos(kx * ax + ky * ay)) + \
                                               params[2] * (2 - 2 * np.cos(kx * ax - ky * ay)) + \
                                               params[3] * (2 - 2 * np.cos(2 * kx * ax)) + \
                                               params[4] * (2 - 2 * np.cos(2 * ky * ay)) + \
                                               params[5]

    band_fn = lambda params: dispersion_fn(params, kxs_2d, kys_2d)
    fit_fn = lambda params: 0.5 * np.sum(np.sum((band_fn(params) - energies_2d) ** 2))

    # do fitting
    init_guess = np.array([0.1, 0.1, 0.01, 0.01, 0.01, np.min(energies_2d.ravel())])
    fit_handle = scipy.optimize.minimize(fit_fn, init_guess)
    fit_params = fit_handle.x

    if is_1d_x:
        fit_params[1] = 0
        fit_params[2] = 0
        fit_params[4] = 0
    elif is_1d_y:
        fit_params[0] = 0
        fit_params[2] = 0
        fit_params[3] = 0
    else:
        pass

    return fit_params, dispersion_fn

def get_hopping_wannier(wannier_fn, potential, xpts, ypts, dx=1., dy=0., print_results=0):
    """
    Compute hopping and density dependent hopping from wannier functions. A reference for the density dependent hopping
    expression is Phys. Rev. Lett. 113, 193003.
    :param wannier_fn: 2D complex numpy array, where W[i,j] = W(x_j, y_i)
    :param xpts: 1D array. must be equally spaced and spacing must evenly divide ax
    :param ypts: 1D array.
    :return hopping:
    :return density_dep_hopping:
    """

    #TODO: ability to shift by non-integer distances?

    t_start = time.time()

    # #####################################
    # Get wannier 2nd_derivative
    # #####################################

    # cannot get derivative at edge points, so crop these off
    if xpts.size > 1 and ypts.size > 1:
        grid_dx = xpts[1] - xpts[0]
        wannier_delx = 1/grid_dx**2 * (wannier_fn[1:-1, 2:] + wannier_fn[1:-1, 0:-2] - 2 * wannier_fn[1:-1, 1:-1])

        grid_dy = ypts[1] - ypts[0]
        wannier_dely = 1 / grid_dy ** 2 * (wannier_fn[2:, 1:-1] + wannier_fn[0:-2, 1:-1] - 2 * wannier_fn[1:-1, 1:-1])

        # get wannier function and potential on new grid
        wannier = wannier_fn[1:-1, 1:-1]
        potential = potential[1:-1, 1:-1]

    elif xpts.size == 1 and ypts.size > 1:
        grid_dx = 1.
        wannier_delx = np.zeros((wannier_fn.shape[0] - 2))

        grid_dy = ypts[1] - ypts[0]
        wannier_dely = 1 / grid_dy ** 2 * (wannier_fn[2:, :] + wannier_fn[0:-2, :] - 2 * wannier_fn[1:-1, :])

        wannier = wannier_fn[1:-1, :]
        potential = potential[1:-1, :]

    elif ypts.size == 1 and xpts.size > 1:
        grid_dx = xpts[1] - xpts[0]
        wannier_delx = 1 / grid_dx ** 2 * (wannier_fn[:, 2:] + wannier_fn[:, 0:-2]
                                           - 2 * wannier_fn[:, 1:-1])

        grid_dy = 1.
        wannier_dely = np.zeros((1, wannier_fn.shape[1] - 2))

        wannier = wannier_fn[:, 1:-1]
        potential = potential[:, 1:-1]
    else:
        raise Exception()

    # #####################################
    # get shifted wannier function
    # #####################################
    # x_shift_index = np.argsort(np.abs(xpts - (xpts[0] + dx)))[0]
    # x_shift_index = int(x_shift_index)
    # dx_real = xpts[x_shift_index] - xpts[0]
    #
    # y_shift_index = np.argsort(np.abs(ypts - (ypts[0] + dy)))[0]
    # y_shift_index = int(y_shift_index)
    # dy_real = ypts[y_shift_index] - ypts[0]
    #
    # wannier_shifted = wannier[y_shift_index:, x_shift_index:]
    #
    ######## get other functions on new grid
    # x_index_end = wannier.shape[1] - x_shift_index
    # y_index_end = wannier.shape[0] - y_shift_index
    #
    # wannier_delx = wannier_delx[:y_index_end, :x_index_end]
    # wannier_dely = wannier_dely[:y_index_end, :x_index_end]
    # wannier = wannier[:y_index_end, :x_index_end]
    # potential = potential[:y_index_end, :x_index_end]

    # interpolate wannier function at points
    # New method
    # points to evaluate normal wannier function at to get shifted wannier function at xpts
    xpts_shifted = xpts + dx
    x_ntrim_end = np.sum(xpts_shifted > np.max(xpts))
    x_ntrim_start = np.sum(xpts_shifted < np.min(xpts))
    if x_ntrim_end > 0:
        xpts_shifted = xpts_shifted[x_ntrim_start:-x_ntrim_end]
        xpts_crop = xpts[x_ntrim_start:-x_ntrim_end]
    else:
        xpts_shifted = xpts_shifted[x_ntrim_start:]
        xpts_crop = xpts[x_ntrim_start:]

    ypts_shifted = ypts + dy
    y_ntrim_end = np.sum(ypts_shifted > np.max(ypts))
    y_ntrim_start = np.sum(ypts_shifted < np.min(ypts))
    if y_ntrim_end > 0:
        ypts_shifted = ypts_shifted[y_ntrim_start:-y_ntrim_end]
        ypts_crop = ypts[y_ntrim_start:-y_ntrim_end]
    else:
        ypts_shifted = ypts_shifted[y_ntrim_start:]
        ypts_crop = ypts[y_ntrim_start:]

    # interpolate wannier function at new points
    # since these are in matrix order, should invert x and y as below
    if xpts.size > 3 and ypts.size > 3:
        wannier_real_fn = scipy.interpolate.RectBivariateSpline(ypts, xpts, np.real(wannier_fn))
        wannier_imag_fn = scipy.interpolate.RectBivariateSpline(ypts, xpts, np.imag(wannier_fn))
        wannier_shifted = wannier_real_fn(ypts_shifted, xpts_shifted) + 1j * wannier_imag_fn(ypts_shifted, xpts_shifted)
    elif xpts.size > 3 and ypts.size == 1:
        wannier_real_fn = scipy.interpolate.InterpolatedUnivariateSpline(xpts, np.real(wannier_fn))
        wannier_imag_fn = scipy.interpolate.InterpolatedUnivariateSpline(xpts, np.imag(wannier_fn))
        wannier_shifted = wannier_real_fn(xpts_shifted) + 1j * wannier_imag_fn(xpts_shifted)
        wannier_shifted = wannier_shifted[None, :]
    elif xpts.size == 1 and ypts.size > 3:
        wannier_real_fn = scipy.interpolate.InterpolatedUnivariateSpline(ypts, np.real(wannier_fn))
        wannier_imag_fn = scipy.interpolate.InterpolatedUnivariateSpline(ypts, np.imag(wannier_fn))
        wannier_shifted = wannier_real_fn(ypts_shifted) + 1j * wannier_imag_fn(ypts_shifted)
        wannier_shifted = wannier_shifted[:, None]
    else:
        raise Exception('At least one of xpts and ypts should have more than three points.')

    # #####################################
    # crop off extra piece where couldn't interpolate wannier function because went outside sampling area.
    # #####################################
    if x_ntrim_start == 0 and x_ntrim_end > 0:
        wannier_shifted = wannier_shifted[:, 1:]
        wannier_delx = wannier_delx[:, :1-x_ntrim_end]
        wannier_dely = wannier_dely[:, :1-x_ntrim_end]
        wannier = wannier[:, :1-x_ntrim_end]
        potential = potential[:, :1-x_ntrim_end]

    elif x_ntrim_start > 0 and x_ntrim_end == 0:
        wannier_shifted = wannier_shifted[:, :-1]
        wannier_delx = wannier_delx[:, x_ntrim_start-1:]
        wannier_dely = wannier_dely[:, x_ntrim_start-1:]
        wannier = wannier_fn[:, x_ntrim_start:]
        potential = potential[:, x_ntrim_start:]
    elif x_ntrim_start == 0 and x_ntrim_end == 0 and xpts.size > 1:
        wannier_shifted = wannier_shifted[:, 1:-1]
    elif xpts.size == 1:
        pass
    else:
        raise Exception("Problem cropping wannier function, x direction")

    if y_ntrim_start == 0 and y_ntrim_end > 0:
        wannier_shifted = wannier_shifted[1:, :]
        wannier_delx = wannier_delx[:1 - y_ntrim_end, :]
        wannier_dely = wannier_dely[:1 - y_ntrim_end, :]
        wannier = wannier[:1 - y_ntrim_end, :]
        potential = potential[:1 - y_ntrim_end, :]

    elif y_ntrim_start > 0 and y_ntrim_end == 0:
        wannier_shifted = wannier_shifted[:-1, :]
        wannier_delx = wannier_delx[y_ntrim_end - 1:, :]
        wannier_dely = wannier_dely[y_ntrim_end - 1:, :]
        wannier = wannier[y_ntrim_end - 1:, :]
        potential = potential[y_ntrim_end - 1:, :]

    elif y_ntrim_start == 0 and y_ntrim_end == 0 and ypts.size > 1:
        wannier_shifted = wannier_shifted[1:-1, :]
    elif ypts.size == 1:
        pass
    else:
        raise Exception("Problem cropping wannier fn, y direction")

    # #####################################
    # Do wannier integral
    # #####################################

    Er = (2 * np.pi / 2) ** 2
    # extra negative sign due to fact we typically define -t \sum_{ij} c^\dag_i c_j as Hamiltonian
    hopping_integrand = -wannier_shifted.conj() * ( -1/ Er * (wannier_delx + wannier_dely) + potential * wannier)
    density_dep_hopping_integrand = -wannier_shifted.conj() * wannier.conj() * wannier**2

    # do integrals
    if xpts.size > 1 and ypts.size > 1:
        hopping_int_x = np.trapz(hopping_integrand, dx=grid_dx, axis=1)
        hopping = np.real(np.trapz(hopping_int_x, dx=grid_dy, axis=0))

        dens_dep_int_x = np.trapz(density_dep_hopping_integrand, dx=grid_dx, axis=1)
        density_dep_hopping = np.real(np.trapz(dens_dep_int_x, dx=grid_dy, axis=0))

    elif xpts.size == 1 and ypts.size > 1:
        hopping = np.real(np.trapz(hopping_integrand, dx=grid_dy, axis=0)).flatten()
        density_dep_hopping = np.real(np.trapz(density_dep_hopping_integrand, dx=grid_dy, axis=0)).flatten()

    elif ypts.size == 1 and xpts.size > 1:
        hopping = np.real(np.trapz(hopping_integrand, dx=grid_dx, axis=1)).flatten()
        density_dep_hopping = np.real(np.trapz(density_dep_hopping_integrand, dx=grid_dx, axis=1)).flatten()

    else:
        raise Exception()

    t_end = time.time()
    if print_results:
        print("Calculating hopping from wannier function with %dx%d points took %0.2f" % (xpts.size + 2, ypts.size + 2, t_end - t_start))

    return hopping, density_dep_hopping

def get_interaction_wannier(wannier_fn, x_grid, y_grid):
    """
    Compute wannier interaction integral, \int |w(x,y)|^4 dxdy

    :param wannier_fn: Wannier function as a 2D array
    :param x_grid: grid of xpts
    :param y_grid: grid of ypts
    :return interaction:
    """
    wannier_probability = np.abs(wannier_fn) ** 2
    xpts = x_grid[0, :].flatten()
    ypts = y_grid[:, 0].flatten()

    if ypts.size > 1:
        y_int = np.trapz(wannier_probability ** 2, ypts, axis=1)
    else:
        y_int = (wannier_probability ** 2).flatten()

    interaction = np.trapz(y_int, xpts, axis=0)
    interaction = np.array([interaction])

    return interaction

#
def project_to_wannier(states, xs_full, wannier_fn, xs, print_results=0):

    # xs in units of lattice sites and equispaced

    t_start = time.time()

    # ensure correct shapes
    states = np.asarray(states)
    if states.ndim <= 1:
        states = states.reshape([states.size, 1])

    # check compatibility of grids
    dx_wn = xs[1] - xs[0]
    dx = xs_full[1] - xs_full[0]
    if not np.round(dx_wn, 14) == np.round(dx, 14):
        raise Exception('xs_full grid and xs grid have different spacing, and are therefore incompatible')

    # get width of site for wannier fn
    site_width = (xs.size - 1) / (xs.ravel()[-1] - xs.ravel()[0])
    if not np.mod(np.round(site_width, 14), 1) == 0:
        raise Exception('Both 0 and 1 must be x coordinates so that we can shift one site with a shift of our indices')
    site_width = int(site_width)
    half_size = 0.5 * (xs.size - 1)

    # get number of lattice site centers
    all_inds = np.arange(0, xs_full.size)
    inds = all_inds[np.mod(np.round(xs_full, 14), 1) == 0]
    nsites = inds.size

    # wannier fns projector
    wfns = np.zeros((xs_full.size, nsites))

    for ii in range(0, nsites):
        start_ind = inds[ii] - half_size
        wn_start_ind = 0
        if start_ind < 0:
            wn_start_ind = -start_ind
            start_ind = 0

        end_ind = inds[ii] + half_size
        wn_end_ind = wannier_fn.size
        if end_ind >= xs_full.size:
            wn_end_ind = wn_end_ind - (end_ind - xs_full.size)
            end_ind = xs_full.size

        wfns[start_ind:end_ind, ii] = wannier_fn[wn_start_ind:wn_end_ind]

    # project states
    proj = wfns.conj().transpose().dot(states)

    t_end = time.clock()
    if print_results:
        print("Projecting %d states at %d points onto wannier function evaluated at %d points took %0.2fs"
          % (states.shape[1], states.shape[0], wannier_fn.size, t_end - t_start) )

    return proj, wfns

# lattice plus harmonic trap

def get_omega_osclen(osc_len):
    """
    Get omega in units of Er given oscillator length in units of the lattice spacing
    W_exp = 0.5 * m * omega_start ** 2 * a ** 2 / Er
    Er = hbar^2 pi^2 / (2m)
    W_exp/Er = (a/l)^4 / np.pi ** 2

    :param osc_len: harmonic oscillator in units of the lattice constant. Given by sqrt(hbar /(m*w)) / a.
    :return:
    W: in units of the recoil energy Er
    """
    W = 1. / (osc_len ** 4 * np.pi ** 2)
    return W

def get_omega_realunits():
    pass

def get_wannier_weights_lattice_plus_harmonic(alpha, orders, num_js):
    """
    Compute the amplitude of the Wannier function on site j in the eigenstate of a combined lattice and harmonic
    potential. Also compute the energy of this eigenstate. The nth eigenstate real-space wavefunction is then given as
    Psi(n, x) = \sum_j f_jn * W_j(x). See PRA 79 063605 (2009) for more details.

    :param alpha:  4t/W_exp where t is the hopping energy and W_exp = 0.5 * m * w^2 * a^2 with m the mass, w the trapping
    frequency and a the lattice spacing.
    :param num_orders: Integer, the number of orders to compute
    :param num_js: Half the number of lattice sites to use
    :return:
    js: lattice site indices
    es: energies of each eigenstate
    fjs: num_sites x num_orders array. fjs[j, n] = f_jn
    """

    # for some reason, the argument of the mathieu cosine is in degrees instead of radians ...
    energy_even = lambda order: scipy.special.mathieu_a(order, alpha)
    energy_odd = lambda order: scipy.special.mathieu_b(order + 1, alpha)
    int_even = lambda x, j, order: scipy.special.mathieu_cem(order, -alpha, x * 180. / np.pi)[0] * np.cos(2 * j * x)
    int_odd = lambda x, j, order: scipy.special.mathieu_sem(order + 1, -alpha, x * 180. / np.pi)[0] * np.sin(2 * j * x)

    # eigenvectors
    #orders = np.arange(0, num_orders)
    orders = np.array([orders])
    # lattice sites
    js = np.arange(-num_js, num_js + 1)

    es = np.zeros(orders.size)
    fjs = np.zeros((len(js), len(orders)))
    for ii in range(0, len(orders)):
        for jj in range(0, len(js)):

            if np.mod(orders[ii], 2) == 0:
                integrand = lambda x: int_even(x, js[jj], orders[ii])
                es[ii] = energy_even(orders[ii])
            else:
                integrand = lambda x: int_odd(x, js[jj], orders[ii])
                es[ii] = energy_odd(orders[ii])

            fjs[jj, ii] = 1 / np.pi * scipy.integrate.quad(integrand, 0, 2 * np.pi)[0]

    return js, es, fjs

def get_wavefn_lattice_plus_harmonic(fjs, js, wannier_fn, xs):
    """
    Compute the sum of wannier states of different lattice sites each with given amplitude. This function accepts the
    output of get_wannier_weights_lattice_plus_harmonic as input. Together, they can be used to compute the
    eigenfunctions of a combined lattice and harmonic trapping potential in real space.
    See PRA 79 063605 (2009) for more details.

    :param fjs: A Numpy array of size nsites x neigenfunctions, where fjs[ii, jj] represents
    the weight of the Wannier function centered on sites ii for the eigenstate jj
    :param js: Site indices. Should run from -n to n, where 2n + 1 = nsites
    :param wannier_fn: A Numpy array giving the Wannier function at points x
    :param xs: The x-coordinates corresponding to wannier_fn. These points must be chosen so that both
    0 and 1 are included so that we can easily shift by one site
    :return:
    wave_fns: A Numpy array of size npts x number of functions

    xs_all: x-coordinates for the wave functions
    """
    t_start = time.clock()

    dx = xs.ravel()[1] - xs.ravel()[0]
    x_start = js[0] + xs.ravel()[0]
    x_end = js[-1] + xs.ravel()[-1]
    xs_all = np.arange(x_start, x_end + dx, dx)

    npts = xs_all.size
    nsites = fjs.shape[0]
    nfns = fjs.shape[1]
    wave_fns = np.zeros([npts, nfns], dtype=np.complex)

    site_width_indices = (xs.size - 1) / (xs.ravel()[-1] - xs.ravel()[0])
    if not np.mod( np.round(site_width_indices, 14), 1) == 0:
        raise Exception('Both 0 and 1 must be x coordinates so that we can shift one site with a shift of our indices')
    site_width_indices = int(site_width_indices)


    full_size = xs.size
    for jj in range(0, nfns):
        for ii in range(0, nsites):
            offset = ii * site_width_indices
            wave_fns[offset: offset + full_size, jj] = \
                wave_fns[offset: offset + full_size, jj] \
                + fjs[ii, jj] * wannier_fn.ravel()

    # ensure normalized
    # not sure why I need to do this???
    wave_fns = wave_fns / np.sqrt(2)

    t_end = time.clock()
    print("Calculating realspace wavefunctions using %d sites %0.2fs" % (nsites, t_end - t_start))

    return wave_fns, xs_all

def get_wavefn_kspace_lattice_plus_harmonic(alpha, orders, qs):
    """
    Compute the quasimomentum representation of the eigenstates for a combined lattice and harmonic potential. See PRA
    79 063605 (2009) for more details.
    :param alpha:  4t/W_exp where t is the hopping energy and W_exp = 0.5 * m * w^2 * a^2 with m the mass, w the trapping
    frequency and a the lattice spacing.
    :param num_orders: Integer, the number of orders to compute
    :param qs: array of quasimomentum values at which the wavefunctions will be computed
    :return:
    """

    t_start = time.clock()

    # ensure orders is a zero dimension vector
    orders = np.asarray([orders])
    orders = orders.reshape(orders.size)

    qb = np.pi
    int_even = lambda q, order: (-1) ** (order/2) / np.sqrt(np.pi) * scipy.special.mathieu_cem(order, alpha, np.pi/2. * (1 - q / qb) * 180. / np.pi)[0]
    int_odd = lambda q, order: (-1) ** ((order - 1) / 2) / np.sqrt(np.pi) * scipy.special.mathieu_sem(order + 1, alpha, np.pi/2. * (1 - q / qb) * 180. / np.pi)[0]

    # eigenvectors
    fqs = np.zeros((qs.size, orders.size))
    for ii in range(0, orders.size):
        for jj in range(0, len(qs)):

            if np.mod(orders[ii], 2) == 0:
                fqs[:, ii] = int_even(qs, orders[ii])
            else:
                fqs[:, ii] = int_odd(qs, orders[ii])

    t_end = time.clock()
    print("Calculating %d q-space wavefunctions at %d q-points took %0.2fs" % (orders.size, qs.size, t_end - t_start))

    return fqs

def get_wavefn_real_space(osc_len, latt_depth, num_eigs=10, nsites=100, npts=2000):
    """
    Solve for eigenfunction of combined lattice and harmonic trapping potential in real space. The lattice potential
    is given by V(r) = 0.5 * D * (1 - cos(2*pi*x) ) + W * x**2, where x is measured in lattice sites.

    :param osc_len: harmonic oscillator in units of the lattice constant. Given by sqrt(hbar /(m*w)) / a. This
           determines the harmonic trap frequency, which is W = 1/(osc_len ** 4 * pi ** 2)
    :param latt_depth: lattice depth in units of lattice recoil.
    :param num_eigs: number of eigenvalues and eigenvectors to solve for
    :param nsites: number of lattice sites to include in calculation
    :param npts: number of real space points to include in calculation. There will be approximately npts/nsites points
           used for each lattice site
    :return:
    eigvals: eigenvalues
    eigvects: matrix of eigenvectors. Each column represents one eigenvector.
    xpts: real space points corresponding to eigvects
    """

    t_start = time.clock()
    # set harmonic trap params
    W = 1. / (osc_len ** 4 * np.pi ** 2)
    # alpha = 1D band width / trap energy
    v = lambda x: latt_depth * 0.5 * (1 - np.cos(2 * np.pi * x)) + W * x ** 2

    # nsites = 100
    # npts = nsites * 20
    xpts = np.linspace(-nsites, nsites, npts)
    dx = xpts[1] - xpts[0]

    # build hamiltonian
    vpart = sp.diags(v(xpts))
    dx_sqr = (1 / dx ** 2) * (-2. * sp.identity(npts)
                              + sp.diags(np.ones(npts - 1), offsets=1, shape=(npts, npts))
                              + sp.diags(np.ones(npts - 1), offsets=-1, shape=(npts, npts)))
    # hbar^2/2m / Er = a ** 2 / pi ** 2
    ham = vpart - 1 / np.pi ** 2 * dx_sqr

    # since this is tridiagonal we should be able to use a tridiagonal matrix solver...
    eigvals, eigvects = scipy.sparse.linalg.eigsh(ham, num_eigs, which='SM')

    # normalize eigvects. They are normalized as vectors, but not if we integrate over x
    #norms = np.trapz(np.abs(eigvects) ** 2, xpts, axis=0)
    norm = np.trapz(np.abs(eigvects[:, 0]) ** 2, xpts, axis=0)
    eigvects = eigvects / np.sqrt(norm)

    t_end = time.clock()

    print("solving real space potential at %d space points for %d eigenvalues took %0.2fs" % (npts, num_eigs, t_end - t_start))

    return eigvals, eigvects, xpts

def solve_time_evolve_real_space(psi_init, osc_len_fn, latt_depths_fn, dt, t_end, nsites=100, npts=2000):
    """
    Solve for the time evolution of a quantum state in a 1D potential with a lattice and a harmonic trap. The total
    potential is V(x) = 0.5 * depth * (1 - cos(2*pi*x) ) + W*x**2.

    Distance is measured in terms of the lattice spacing, energy is measured in terms of the recoil energy,
    and time is measured in terms of hbar / Er.

    :param psi_init: Initial state psi, or an array where each column represents a different initial state
    :param osc_len_fn: Function which returns the harmonic oscillator length in units of the lattice constant as a function
           of time. The oscillator length is defined by sqrt(hbar /(m*w)) / a. This determines the harmonic trap
            frequency, which is W = 1/(osc_len ** 4 * pi ** 2).
    :param latt_depths_fn: Function which gives the lattce depth in units of Er as a function of time.
    :param dt: Time spacing in units of hbar/Er
    :param t_end: End time in units of hbar/Er
    :param nsites: Number of  lattice sites to use in the calculation
    :param npts: Number of total points to use in the calculation. Number of points per lattice site is approximately npts/nsites.
    :return:
    psis: wavefunction at given times
    times: [0, dt, ... t_end]
    latt_depths: lattice depths at given times
    osc_lens: oscillator lengths at given times
    xpts: real-space points
    """

    # ensure psi_init is a column vector
    psi_init = np.asarray(psi_init)
    if psi_init.ndim <= 1:
        psi_init = psi_init.reshape([psi_init.size, 1])

    # alpha = 1D band width / trap energy
    v_latt_fn = lambda x: 0.5 * (1 - np.cos(2 * np.pi * x))
    v_harm_fn = lambda x: x ** 2

    # TODO: maybe better to take xpts as an argument ??? Handle case where not uniformly spaced.
    xpts = np.linspace(-nsites, nsites, npts)
    dx = xpts[1] - xpts[0]

    # build hamiltonian
    v_latt_noscale = sp.diags(v_latt_fn(xpts))
    v_harm = sp.diags(v_harm_fn(xpts))
    dx_sqr = (1 / dx ** 2) * (-2. * sp.identity(npts)
                              + sp.diags(np.ones(npts - 1), offsets=1, shape=(npts, npts))
                              + sp.diags(np.ones(npts - 1), offsets=-1, shape=(npts, npts)))

    # time dependent pieces
    times = np.arange(0, t_end + dt, dt)
    latt_depths = latt_depths_fn(times)
    osc_lens = osc_len_fn(times)
    Ws = 1. / (osc_lens ** 4 * np.pi ** 2)

    if psi_init.ndim == 1:
        nstates = 1
    else:
        nstates = psi_init.shape[1]

    psis = np.zeros((xpts.size, times.size, nstates), dtype=np.complex)
    psis[:, 0, :] = psi_init
    for ii in range(1, len(times)):
        t_iter_start = time.clock()

        # hbar^2/2m / Er = a ** 2 / pi ** 2
        ham = latt_depths[ii - 1] * v_latt_noscale + Ws[ii - 1] * v_harm - 1 / np.pi ** 2 * dx_sqr
        num = sp.identity(ham.shape[0]) - dt * 0.5 * 1j * ham

        # scipy.linalg.eigh_tridiagonal does not support complex matrices, so
        # diagonalize hamiltonian and use this to construct our desired inverse
        eig_vals, eig_vects = scipy.linalg.eigh_tridiagonal(ham.diagonal(0), ham.diagonal(1))
        # get 1 / (1 + 0.5 * dt * 1j * Ham)
        mat_inv_diag = sp.diags(np.divide(1., 1. + dt * 0.5 * 1j * eig_vals)).tocsc()
        # convert this to our original basis
        denom = eig_vects.dot(mat_inv_diag.dot(eig_vects.conj().transpose()))

        for jj in range(0, nstates):
            psis[:, ii, jj] = num.dot(denom.dot(psis[:, ii - 1, jj]))

        t_iter_end = time.clock()
        print("solved step %d/%d for %d states in %0.2fs" % (ii, len(times), nstates, t_iter_end - t_iter_start))
    psis = np.squeeze(psis)

    return psis, times, latt_depths, osc_lens, xpts

# semiclassical profiles

def get_semiclassical_profile(osc_len, tunneling, beta, mu_o, nsites=30, nks=50, dim='1d'):
    """
    Solve for atom density on each site for a 1D semiclassical lattice with a tight-binding dispersion in a harmonic trap.
    Energies are in units of the recoil energy, and lengths are in units of the lattice spacing.

    :param osc_len: Oscillator length for harmonic oscillator in units of the lattice spacing. Given by sqrt(hbar /(m*w)) / a
    :param tunneling: Tunneling energy in units of the recoil energy
    :param beta: Inverse temperature in units of the recoil energy.
    :param mu_o: Chemical potential at the center of the trap in units of the recoil energy.
    :param nsites: Number of lattice sites to compute the real space density profile at.
    :param nks: number of k-vectors to include in the computation. Use more k-vectors for a more accurate computation.
    :return:
    nr_sc: lattice density profile
    xpts: real space points corresponding to lattice density profile
    jacobian: Chemical potential derivative of total atom number, dn/dmu
    """

    W = get_omega_osclen(osc_len)

    if dim == '1d':
        xpts = np.arange(-nsites, nsites + 1)
        ypts = np.zeros(xpts.shape)
        mus = mu_o - W * xpts ** 2
    elif dim == '2d':
        pts = np.arange(-nsites, nsites + 1)
        xpts, ypts = np.meshgrid(pts, pts)
        mus = mu_o - W * (xpts**2 + ypts**2)
    else:
        raise Exception("dim should be '1d' or '2d' but was %s." % dim)

    # TODO: finish implementing x and y
    nr_sc = np.zeros(xpts.shape)
    jacob_terms = np.zeros(xpts.shape)
    for ii in range(0, xpts.size):
        coord = np.unravel_index(ii, xpts.shape)
        # convert to tunneling units for fg functions
        nr_sc[coord] = fg.fg_density(beta * tunneling, mus[coord] / tunneling, nsites=nks, dim=dim)
        jacob_terms[coord] = fg.fg_compressibility(beta * tunneling, mus[coord] / tunneling, nsites=nks, dim=dim)

    # This is the jacobian of the function N(mu) = np.sum(nr_sc)
    jacobian = np.sum(np.ravel(jacob_terms))

    return nr_sc, xpts, ypts, jacobian

def get_semiclassical_profile_fixed_num(osc_len, tunneling, beta, nparticles, nsites=30, nks=50, dim='1d'):
    """
    Solve for the chemical potential giving a specified atom number for a 1D semiclassical lattice

    :param osc_len: Oscillator length for harmonic oscillator. Given by sqrt(hbar /(m*w)) / a
    :param tunneling: Tunneling energy in units of the recoil energy
    :param beta: Inverse temperature in units of the recoil energy
    :param nparticles: number of atoms
    :param nsites: number of lattice sites to compute profile at
    :param nks: number of k-vectors to include in the computation. Use more k-vectors for a more accurate computation.
    :return:
    nr_sc:
    xpts:
    mu_final:
    """
    fn = lambda mu: np.sum(get_semiclassical_profile(osc_len, tunneling, beta, mu, nsites=nsites, nks=nks, dim=dim)[0])
    jacobian = lambda mu: get_semiclassical_profile(osc_len, tunneling, beta, mu, nsites=nsites, nks=nks, dim=dim)[3]
    jacobian_fsqr = lambda mu: 2 * (nparticles - fn(mu)) * jacobian(mu)
    #min_fn = lambda mu: (nparticles - fn(mu)) ** 2

    # had better luck using root finding than trying to specify jacobian in minimization function.
    # Not sure if I'm doing it wrong or what?
    fit_handle = scipy.optimize.root_scalar(jacobian_fsqr, x0=-0.1, x1=0.1)
    # fit_handle = scipy.optimize.root_scalar(jacobian_fsqr, x0=-0.01, x1=0.01, bracket = [-0.1, 0.1])
    mu_final = fit_handle.root
    # fit_handle = scipy.optimize.minimize(min_fn, np.array([0.1]), method='BFGS', jac=jacobian_fsqr)
    # fit_handle = scipy.optimize.minimize(min_fn, np.array([0.1]), jac=jacobian_fsqr)
    # mu_final = fit_handle.x

    nr_sc, xpts, ypts, jacob_final = get_semiclassical_profile(osc_len, tunneling, beta, mu_final, nsites=nsites, nks=nks, dim=dim)

    return nr_sc, xpts, ypts, mu_final

def solve_atomic_limit(beta, mu, U):
    """
    Return thermodynamic parameters for the atomic limit at inverse temperature beta, chemical potential mu, and interaction U.
    Beta, mu, and U must be given in the same units.

    :param beta: inverse temperature
    :param mu: chemical potential
    :param U: interaction energy
    :return:
    density: density per lattice site
    energy: energy, in the same units as T, mu, and U
    entropy: Dimensionless (i.e. S/kb)
    """

    Z = 1 + 2 * np.exp(beta * mu) + np.exp(-beta * (U - 2 * mu))
    density = np.divide(2 * np.exp(beta * mu) + 2 * np.exp(-beta * (U - 2 * mu)), Z)
    energy = np.divide(U * np.exp(-beta * (U - 2 * mu)), Z)
    # grand canonical ensemble calculation, S = - d\Omega/dT|_{v,\mu}
    entropy = np.log(Z) - np.divide(2 * np.exp(beta * mu) * mu * beta + (2 * mu - U) * beta * np.exp(-beta * (U - 2*mu)), Z)

    return density, energy, entropy

# SHO functions
def get_sho_wavefn(osc_length, xs, orders):
    """
    Get normalized harmonic oscillator wavefunctions
    :param osc_length: oscillator length defined as np.sqrt(hbar / (m*w) ). Give in the same units as xs, which will
    typically be in lattice sites
    :param xs: x-coordinates to evaluate the harmonic oscillator wavefunctions at. Give in the same units as osc_length
    :param orders: order of the wavefunction to evaluate
    :return:
    vals
    """
    # length measured in units of lattice spacing

    # oscillator length is defined as np.sqrt(hbar / (m*w) ). We meausre it in units of the lattice constant, so we
    # we should actually use np.sqrt(hbar/(m*w)) / a, with a the lattice constant

    xs = np.array(xs)
    orders = np.array(orders)

    if orders.size > 1:
        vals = np.zeros((xs.size, orders.size))
        for ii in range(0, orders.size):
            vals[:, ii] = get_sho_wavefn(osc_length, xs, orders[ii])
    else:
        factor = 1 / math.sqrt( 2 ** orders) / math.sqrt(math.factorial(orders) ) * 1 / np.sqrt(osc_length) * np.pi ** (-0.25)
        hpoly = scipy.special.hermite(orders)
        vals = factor * np.exp(- 0.5 * xs ** 2 / osc_length ** 2) * hpoly(xs / osc_length)

    return vals

def evolve_harmonic_trap(osc_length, omega, xs, psi_init, times):
    """
    TODO: seems like in some cases this is not preserving the normalization of the wavefunctions!? Am I missing a factor? Need to check carefully.
    Evolve wavefunction in harmonic trap using propogator. See doi: 10.1119/1.4960479 for a detailed discussion of the
    harmonic oscillator propagator, as well as a discussion of the propagator at half-period points.
    :param osc_length: Harmonic oscillator length defined as np.sqrt(hbar / (m*omega_start)) given in the same units as xs
    :param w: Harmonic oscillator angular frequency omegas. If times are measured in units of A, then omega_start should be
    measured in units of 1 / A
    :param xs: spatial sample points, given in same units as harmonic oscillator length. Measured from the center of
    the oscillator potential
    :param psi_init: Initial wavefunction. Either a single wavefunction, or an array where each column represents a
    different wavefunction
    :param times: times to calculate evolved wavefunction at
    :return:
    psis: a psi_init.size x len(times) array where each column gives the initial wavefunction evolved to the supplied time
    """
    # ensure psi_init is a column vector
    t_start = time.clock()

    psi_init = np.asarray(psi_init)
    if psi_init.ndim <= 1:
        psi_init = psi_init.reshape([psi_init.size, 1])

    times = np.asarray(times)

    dx = xs[1] - xs[0]
    xi, xf = np.meshgrid(xs, xs)
    # the following expression holds except when
    prop = lambda t: 1. / osc_length * 1. / np.sqrt( 2 * np. pi * 1j * np.sin(omega * t) ) * \
                     np.exp(1j / (2 * osc_length ** 2 * np.sin( omega * t)) * ( (xi**2 + xf**2) * np.cos( omega * t ) - 2 * xi * xf) )

    if psi_init.ndim == 1:
        nstates = 1
    else:
        nstates = psi_init.shape[1]

    psis = np.zeros((xs.size, times.size, nstates), dtype=np.complex)
    for ii in range(0, times.size):

        if np.round( np.mod( omega * times[ii] / (2 * np.pi), 1), 14) == 0:
            psis[:, ii, :] = psi_init
        elif np.round( np.mod( omega * times[ii] / (np.pi), 1), 14) == 0:
            psis[:, ii, :] = -psi_init
            # what do do in this case
        else:
            # Psi[x, t] = \int dx' K[x, t; x', 0] Psi[x', 0]
            for jj in range(0, nstates):
                psis[:, ii, jj] = prop(times[ii]).dot(psi_init[:, jj]) * dx

    psis = np.squeeze(psis)

    t_end = time.clock()

    print("Evolving %d points and %d states at %d times in harmonic potential took %0.2fs" % (xs.size, nstates, times.size, t_end - t_start))

    return psis

# reshaping function

def oned_2_twod_kvects(bz_vects, eigvals, n1_sites, n2_sites, eigvects=None):
    """
    Reshapes the 1D representation of the brillouin zone vectors to a 2d/grid representation.
    TODO: write this function in a more general way to handle generic arrays

    :param bz_vects:
    :param eigvals:
    :param n1_sites:
    :param n2_sites:
    :return: kxs_2d, kys_2d, bands_2d
    """

    # reshape Brillouin zone vectors
    kxs_2d = np.reshape(bz_vects[0, :], [n2_sites, n1_sites])
    kys_2d = np.reshape(bz_vects[1, :], [n2_sites, n1_sites])

    # reshape eigenvalues
    if eigvals.ndim == 1:
        # if e.g. only reshaping one band
        bands_2d = np.reshape(eigvals, [n2_sites, n1_sites])
    elif eigvals.ndim == 2:
        bands_2d = np.reshape(eigvals, [n2_sites, n1_sites, eigvals.shape[1]])
    else:
        raise Exception()

    # reshape eigenvectors
    if eigvects != None:
        if eigvects.ndim == 2:
            eigvects_2d = np.reshape(eigvects, [n2_sites, n1_sites, eigvects.shape[2]])
        elif eigvects.ndim == 3:
            eigvects_2d = np.reshape(eigvects, [n2_sites, n1_sites, eigvects.shape[2], eigvects.shape[3]])
    else:
        eigvects_2d = 0

    return kxs_2d, kys_2d, bands_2d, eigvects_2d

def twod_2_highsymm(xx, yy, array, preserve_size=False):
    """
    Given a 2d array representing values at a grid of kx and ky points, return the values at the array along the high
    symmetry lines of Brillouin zone, namely the Gamma-X-M-Gamma and Gamma-Y-M-Gamma lines
    :param xx:
    :param yy:
    :param array:
    :return: x_gxmg, y_gxmg, array_gxmg, x_gymg, y_gymg, array_gymg, linear_coord

    x_gxmg are the kx points along the Gamma-X-M-Gamma line as a 1d NumPy array
    array_gxmg is the array along the Gamma-X-M-Gamma line, with the first dimension representing points on this line.

    similarly for the gymg arguments, except along the Gamma-Y-M-Gamma line.

    linear_coord gives linear coordinates which can be used as an x-value for plotting quantities along these lines of
    high symmetry. 0 = Gamma, 1 = X (or Y), 2 = M, and 3 = Gamma. These coordinates are useful for plotting quantities
    which are sampled at different k-points (for example the computed band structure and a fit to the band structure
    sampled on a much finer k-grid).
    """
    nx_sites = xx.shape[1]
    ny_sites = xx.shape[0]

    if nx_sites != ny_sites:
        raise Exception('twod_2_highsymm not implemented for different number of sites in x and y')

    if nx_sites % 2 == 0:
        print("Warning number of coordinates was not odd.")

    center_x = (nx_sites - 1) / 2
    center_y = (ny_sites - 1) / 2

    # Gamma to X
    x_gx = xx[center_y, center_x:]
    y_gx = yy[center_y, center_x:]
    array_gx = array[center_y, center_x:, ...]

    # X to M
    # add 1 to exclude point exactly at X, which we already have on the GX line
    x_xm = xx[center_y + 1:, -1]
    y_xm = yy[center_y + 1:, -1]
    array_xm = array[center_y + 1:, -1, ...]

    # M to Gamma
    #
    yindex = range(center_y, xx.shape[1])
    yindex.reverse()
    xindex = range(center_x, xx.shape[0])
    xindex.reverse()
    x_mg = xx[yindex, xindex]
    y_mg = yy[yindex, xindex]
    array_mg = array[yindex, xindex, ...]

    # alternatively, path to y
    # Gamma to Y
    x_gy = xx[center_y:, center_x]
    y_gy = yy[center_y:, center_x]
    array_gy = array[center_y:, center_x, ...]

    # Y to M
    # add 1 to exclude point exactly at X, which we already covered
    x_ym = xx[-1, center_x + 1:]
    y_ym = yy[-1, center_x + 1:]
    array_ym = array[-1, center_x + 1:, ...]

    # collect results
    x_gxmg = np.concatenate((x_gx, x_xm, x_mg), 0)
    y_gxmg = np.concatenate((y_gx, y_xm, y_mg), 0)
    array_gxmg = np.concatenate((array_gx, array_xm, array_mg), 0)

    x_gymg = np.concatenate((x_gy, x_ym, x_mg), 0)
    y_gymg = np.concatenate((y_gy, y_ym, y_mg), 0)
    array_gymg = np.concatenate((array_gy, array_ym, array_mg), 0)

    # linear coordinates, with 0 = gamma, 1 = X, 2 = M, 3 = Gamma
    # TODO: correct for possibility k vectors may be different in different directions
    # TODO: how best to deal with that? in that case may not have k-pts exactly on this line
    # TODO: I think this current method works
    gx_coord = x_gx / x_gx[-1]
    xm_coord = y_xm / y_xm[-1] + 1.
    mg_coord = (1. - x_mg / x_mg[0]) + 2.
    if preserve_size is True:
        mg_coord = np.sqrt(2) * mg_coord

    dist_gxmg = np.concatenate((x_gx, y_xm + x_gx[-1], np.sqrt(x_mg[0]**2 + y_mg[0]**2) - np.sqrt(x_mg ** 2 + y_mg ** 2) + x_gx[-1] + y_xm[-1]))
    dist_gymg = np.concatenate((y_gy, x_ym + y_gy[-1], np.sqrt(x_mg[0]**2 + y_mg[0]**2) - np.sqrt(x_mg ** 2 + y_mg ** 2) + y_gy[-1] + x_ym[-1]))

    linear_coord = np.concatenate((gx_coord, xm_coord, mg_coord))

    return x_gxmg, y_gxmg, dist_gxmg, array_gxmg,\
           x_gymg, y_gymg, dist_gymg, array_gymg, linear_coord

# helper functions
def get_interaction_real_units(latt_const, scattering_length, omega_z, wannier_integral):
    """
    Get interaction energy in real units given the wannier integral \int |w(x,y)|^4 dxdy.

    :param latt_const: lattice constant in meters
    :param scattering_length: scattering length in bohr radii
    :param omega_z: axial trapping frequency
    :param wannier_integral: wannier integral
    :return:
    """
    abohr = 5.2917721067e-11
    hbar = 1.054e-34
    h = 2 * np.pi * hbar
    # constants for our system
    mli = 9.9883e-27  # kg
    k = 2 * np.pi / latt_const
    Er = hbar ** 2 * (k / 2) ** 2 / (2 * mli)

    # third direction factor. Units of 1/length
    zfactor = np.sqrt(mli * omega_z / h)

    # interaction factor. Units of Energy * volume
    g = 4 * np.pi * scattering_length * abohr * hbar ** 2 / mli

    # interaction in Er
    # wannier_integral is dimensionless. To get correct units of 1/length^2, divide by lattice constant squared
    interaction_real_units =  g / Er * zfactor * wannier_integral / latt_const ** 2

    return interaction_real_units

# exact results for 1D
def get_1d_bandwidth_exact(latt_depth, band_index=0):
    """
    Get exact bandwidth for 1D lattice problem using Mathieu solution.

    We use the fact that at q=0 the solutions are periodic and even, hence their eigenvalues are the Mathieu characteristic
    a values. At q=pi the solutions are periodic (with different peirod) and odd, hence their eigenvalues are the Mathieu
    entries characteristic valus

    :param latt_depth: lattice depth in Er
    :param band_index: The ground band is counted as zero. First excited band is 1, etc.
    :return:
    """
    q = latt_depth * 0.25

    bw_mathieu = scipy.special.mathieu_b(band_index+1, q) - scipy.special.mathieu_a(band_index, q)

    return bw_mathieu

def hopping_1d_asymptotic_mathieu(latt_depth):
    """
    Exact solution for a deep 1D lattice.

    See for example arXiv:1710.00657v1
    Other good references are "Mott-Hubbard transition of cold atoms in optical lattices" by Wilhelm Zwerger https://doi.org/10.1088/1464-4266/5/2/352
    where he gives this expression as equation (8).
    and the thesis "One-Dimensional Mathematical Model for Quantum Particles in Weakly Curved Optical Lattices"
    by Sandro Godtel
    :param latt_depth:
    :return:
    """
    return 4 / np.sqrt(np.pi) * latt_depth ** 0.75 * np.exp(-2 * np.sqrt(latt_depth))

def get_1d_lattice_zero_momentum_wavefn_exact(latt_depth, band_index, xpts, detuning="red"):
    """
    Compute mathieu result for q=0 wavefunctions in 1D lattice
    :param latt_depth:
    :param band_index:
    :param xpts:
    :param detuning:
    :return:
    """
    if detuning == "red":
        mathieu_q = 0.25 * latt_depth
    elif detuning == "blue":
        mathieu_q = -0.25 * latt_depth
    else:
        raise Exception("detuning argument should have been either 'red' or 'blue', but instead it was %s", detuning)

    if band_index % 2 == 0:
        mathieu_wavefn = scipy.special.mathieu_cem(band_index, mathieu_q, xpts * np.pi * (180. / np.pi))[0]
    else:
        mathieu_wavefn = scipy.special.mathieu_sem(band_index + 1, mathieu_q, xpts * np.pi * (180. / np.pi))[0]

    dx = xpts[1] - xpts[0]
    norm = np.sum(mathieu_wavefn ** 2 * dx)
    mathieu_wavefn = mathieu_wavefn / np.sqrt(norm)

    return mathieu_wavefn

# lattices
def get_fourfold_latt(theta=90., r=1., alpha=0.):
    """
    Return Fourier components for the four-fold interfering lattice geometry used in the experiment at unit
    depth. Fourier components for a given depth (in Er) can be obtained by multiplying these by the depth.
    The components are chosen with the correct phase to put the potential minimum at the origin.

    TODO: should I use length of units of wavelength of light instead of always setting a_mean = 1?
    :param theta: angle between the lattice beams in degrees
    :param r: retroreflected electric field attenuation factor
    :param alpha: polarization angle in degrees, with an angle of zero corresponding to vertical polarization
    :param a_mean: geometric mean of lattice constants a_x and a_y

    :return:
    a_x:
    a_y:
    b1:
    b2:
    fourier_indices1:
    fourier_indices2:
    fourier_comp:
    """
    # theta = 91.6267
    # r = 0.5
    theta = theta * np.pi / 180.
    alpha = alpha * np.pi / 180.

    # put distances in units of the ideal lattice, theta = 90. i.e. in units of lambda/sqrt(2)
    a_mean = 1 / np.sqrt(np.sin(theta))

    # lattice spacing quantities
    spacing_ratio = np.tan(0.5 * theta)
    a_x = a_mean * np.sqrt(spacing_ratio)
    a_y = a_mean / np.sqrt(spacing_ratio)

    # define fourier components
    pol_factor = np.cos(alpha) ** 2 - np.cos(theta) * np.sin(alpha) ** 2
    depth_denom = 4 * ((1 + r ** 2) * pol_factor + 2 * r)
    fact_dc = - 2 * (1 + r ** 2) / depth_denom
    fact_a = -2 * r * pol_factor / depth_denom
    fact_b = -(1 + r ** 2) * pol_factor / depth_denom
    fact_c = -r / depth_denom

    fourier_comp = np.array([fact_dc, fact_a, fact_a, fact_b, fact_b, fact_c, fact_c, fact_c, fact_c])

    # reciprocal basis vectors
    b1 = np.array([[2. * np.pi / a_x], [0.]])
    b2 = np.array([[0.], [2. * np.pi / a_y]])
    fourier_indices1 = [0, -1, 1, 0, 0, 1, 1, -1, -1]
    fourier_indices2 = [0, 0, 0, -1, 1, 1, -1, 1, -1]

    return a_x, a_y, b1, b2, fourier_indices1, fourier_indices2, fourier_comp

def solve_d4_lattice(depth, theta=90., r=1., alpha=0., dx=0.05, n1_sites=12, n2_sites=12, n1_recp_max=10, n2_recp_max=10, normalize_er=False):
    """
    solve d4 lattice for most commonly used quantities including hoppings, wannier integral, etc.

    :param depth: Total depth in Er
    :param theta: angle between beams. Ideal angles is 90 degrees
    :param r: retrorefelection e-field attenuation
    :param alpha: polarization angle in degrees. 0 corresponds to vertical polarization
    :param dx: Step size to use in calculating wannier functions
    :param n1_sites: number of k-points in Brillouin zone in direction 1
    :param n2_sites: number of k-points in Brillouin zone in direction 2
    :param n1_recp_max: number of reciprocal lattice vectors in direction 1
    :param n2_recp_max: number of reciprocal lattice vectors in direction 2
    :param normalize_er: If false, use the true Er of the lattice, which is based on the geometric mean of the reciprocal
    lattice vectors. If true, use the Er of the ideal lattice (i.e. alpha = 0, theta = 90, r = 1).
    :return:
    """

    depth = np.asarray(depth)
    if depth.ndim == 0:
        depth = depth[None]

    if depth.size > 1:
        # solve other depths
        for ii in range(0, depth.size):
            print("solving depth %d/%d" % (ii + 1, depth.size))
            tbs, tbw, tw, twn, wi, ev, wf, v_real, xpts, ypts = \
                solve_d4_lattice(depth[ii], theta=theta, r=r, alpha=alpha, dx=dx, n1_sites=n1_sites, n2_sites=n2_sites,
                                 n1_recp_max=n1_recp_max, n2_recp_max=n2_recp_max, normalize_er=normalize_er)

            if ii == 0:
                # variables to store outputs
                ts_bs = np.zeros((depth.size, tbs.size))
                ts_bw = np.zeros((depth.size, tbw.size))
                ts_wannier = np.zeros((depth.size, tw.size))
                ts_wannier_n = np.zeros((depth.size, twn.size))
                wannier_integral = np.zeros(depth.size)
                eigvals = np.zeros((ev.shape[0], ev.shape[1], ev.shape[2], depth.size))
                wannier_fn = np.zeros((wf.shape[0], wf.shape[1], depth.size), dtype=np.complex)

            # assign to variables
            ts_bs[ii, :] = tbs
            ts_bw[ii, :] = tbw
            ts_wannier[ii, :] = tw
            ts_wannier_n[ii, :] = twn
            wannier_integral[ii] = wi
            eigvals[:, :, :, ii] = ev
            wannier_fn[:, :, ii] = wf
    else:
        # get lattice fourier components
        ax, ay, b1, b2, fourier_inds1, fourier_inds2, fourier_comp = get_fourfold_latt(theta=theta, r=r, alpha=alpha)

        # define recoil energies
        er_ideal = (2 * np.pi / 2.) ** 2
        a_rel = 1 / np.sqrt(np.sin(theta * np.pi / 180.))
        er = (2 * np.pi / a_rel / 2.) **2

        if normalize_er == True:
            # interpret depths as being in terms of er_ideal, and translate them to being in units of er
            depth = depth * er_ideal / er

        #
        v_fourier_components = depth * fourier_comp

        # solve lattice problem
        ts_bs, ts_bw, ts_wannier, ts_wannier_n, wannier_integral, eigvals, wannier_fn, v_real, xpts, ypts = \
            solve_lattice(ax, ay, b1, b2, fourier_inds1, fourier_inds2, v_fourier_components,
                          dx=dx, n1_sites=n1_sites, n2_sites=n2_sites, n1_recp_max=n1_recp_max, n2_recp_max=n2_recp_max)

        if normalize_er:
            # calculation spits these out in units of er, convert to er_ideal
            ts_bs = ts_bs * er_ideal / er
            ts_bw = ts_bw * er_ideal / er
            ts_wannier = ts_wannier * er_ideal / er
            ts_wannier_n = ts_wannier_n * er_ideal / er
            eigvals = eigvals * er_ideal / er
            v_real = v_real * er_ideal / er

    return ts_bs, ts_bw, ts_wannier, ts_wannier_n, wannier_integral, eigvals, wannier_fn, v_real, xpts, ypts

def get_2d_separable_latt(a=1.):
    """
    Return Fourier components for the 2D separable lattice geometry at unit
    depth. Fourier components for a given depth (in Er) can be obtained by multiplying these by the depth.
    The components are chosen with the correct phase to put the potential minimum at the origin.

    # TODO: could expand this function to deal with lattices not perfectly at 90deg, with different depths, etc
    :param theta: angle between the lattice beams in degrees
    :param r: retroreflected electric field attenuation factor
    :param alpha: polarization angle in degrees, with an angle of zero corresponding to vertical polarization
    :return:
    b1:
    b2:
    fourier_inds1:
    fourier_inds2:
    v_fourier_components:
    """

    v_fourier_components = 0.5 * np.array([-1., -0.25, -0.25, -0.25, -0.25])
    fourier_inds1 = [0, -1, 1, 0, 0]
    fourier_inds2 = [0, 0, 0, -1, 1]
    b1 = np.array([[2. * np.pi / a], [0.]])
    b2 = np.array([[0.], [2. * np.pi / a]])

    return b1, b2, fourier_inds1, fourier_inds2, v_fourier_components

def solve_separable_latt(depth, a=1., dx=0.05, n1_sites=12, n2_sites=12, n1_recp_max=10, n2_recp_max=10):

    depth = np.asarray(depth)
    if depth.ndim == 0:
        depth = depth[None]

    if depth.size > 1:
        # solve first problem, get sizes
        ts_bs1, ts_bw1, ts_w1, ts_wn1, wi1, ev1, wf1, v_real, xpts, ypts = solve_separable_latt(depth[0], a, dx, n1_sites, n2_sites, n1_recp_max, n2_recp_max)

        # variables to store outputs
        ts_bs = np.zeros((depth.size, ts_bs1.size))
        ts_bw = np.zeros((depth.size, ts_bw1.size))
        ts_wannier = np.zeros((depth.size, ts_w1.size))
        ts_wannier_n = np.zeros((depth.size, ts_wn1.size))
        wannier_integral = np.zeros(depth.size)
        eigvals2d = np.zeros((ev1.shape[0], ev1.shape[1], ev1.shape[2], depth.size))
        wannier_fn = np.zeros((wf1.shape[0], wf1.shape[1], depth.size), dtype=np.complex)

        # assign first depth
        ts_bs[0, :] = ts_bs1
        ts_bw[0, :] = ts_bw1
        ts_wannier[0, :] = ts_w1
        ts_wannier_n[0, :] = ts_wn1
        wannier_integral[0] = wi1
        eigvals2d[:, :, :, 0] = ev1
        wannier_fn[:, :, 0] = wf1

        # solve other depths
        for ii in range(1, depth.size):
            print("solving depth %d/%d" % (ii + 1, depth.size))
            tbs, tbw, tw, twn, wi, ev, wf, _, _, _ = \
                solve_separable_latt(depth[ii], a=a, dx=dx, n1_sites=n1_sites, n2_sites=n2_sites, n1_recp_max=n1_recp_max, n2_recp_max=n2_recp_max)

            ts_bs[ii, :] = tbs
            ts_bw[ii, :] = tbw
            ts_wannier[ii, :] = tw
            ts_wannier_n[ii, :] = twn
            wannier_integral[ii] = wi
            eigvals2d[:, :, :, ii] = ev
            wannier_fn[:, :, ii] = wf

    else:
        # get lattice fourier components
        b1, b2, fourier_inds1, fourier_inds2, fourier_comp = get_2d_separable_latt(a=a)
        v_fourier_components = depth * fourier_comp

        # solve problem
        ts_bs, ts_bw, ts_wannier, ts_wannier_n, wannier_integral, eigvals2d, wannier_fn, v_real, xpts, ypts = \
            solve_lattice(a, a, b1, b2, fourier_inds1, fourier_inds2, v_fourier_components,
                          dx=dx, n1_sites=n1_sites, n2_sites=n2_sites, n1_recp_max=n1_recp_max, n2_recp_max=n2_recp_max)

    return ts_bs, ts_bw, ts_wannier, ts_wannier_n, wannier_integral, eigvals2d, wannier_fn, v_real, xpts, ypts

def get_1d_latt(a=1.):
    """
    Return Fourier components for the 1D lattice geometry at unit
    depth. Fourier components for a given depth (in Er) can be obtained by multiplying these by the depth.
    The components are chosen with the correct phase to put the potential minimum at the origin.
    :param a:
    :return:
    """
    v_fourier_components = np.array([-0.5, -0.25, -0.25])
    fourier_inds1 = [0, -1, 1]
    fourier_inds2 = [0, 0, 0]
    b1 = np.array([[2. * np.pi / a], [0.]])
    b2 = np.array([[0.], [2. * np.pi / a]])

    return b1, b2, fourier_inds1, fourier_inds2, v_fourier_components

def solve_1d_latt(depth, a=1., dx=0.05, n1_sites=12, n1_recp_max=10):
    """
    solve 1d lattice problem at any depth

    :param depth:
    :param a:
    :param dx:
    :param n1_sites:
    :param n1_recp_max:
    :return:
    """
    depth = np.asarray(depth)
    if depth.ndim == 0:
        depth = depth[None]

    if depth.size > 1:
        # solve first problem, get sizes
        ts_bs1, ts_bw1, ts_w1, ts_wn1, wi1, ev1, wf1, v_real, xpts, ypts = solve_1d_latt(depth[0], a, dx, n1_sites, n1_recp_max)

        # variables to store outputs
        ts_bs = np.zeros((depth.size, ts_bs1.size))
        ts_bw = np.zeros((depth.size, ts_bw1.size))
        ts_wannier = np.zeros((depth.size, ts_w1.size))
        ts_wannier_n = np.zeros((depth.size, ts_wn1.size))
        wannier_integral = np.zeros(depth.size)
        eigvals = np.zeros((ev1.shape[0], ev1.shape[1], ev1.shape[2], depth.size))
        wannier_fn = np.zeros((wf1.shape[0], wf1.shape[1], depth.size), dtype=np.complex)

        # assign first depth
        ts_bs[0, :] = ts_bs1
        ts_bw[0, :] = ts_bw1
        ts_wannier[0, :] = ts_w1
        ts_wannier_n[0, :] = ts_wn1
        wannier_integral[0] = wi1
        eigvals[:, :, :, 0] = ev1
        wannier_fn[:, :, 0] = wf1

        # solve other depths
        for ii in range(1, depth.size):
            print("solving depth %d/%d" % (ii + 1, depth.size))
            tbs, tbw, tw, twn, wi, ev, wf, _, _, _ = \
                solve_1d_latt(depth[ii], a=a, dx=dx, n1_sites=n1_sites, n1_recp_max=n1_recp_max)

            ts_bs[ii, :] = tbs
            ts_bw[ii, :] = tbw
            ts_wannier[ii, :] = tw
            ts_wannier_n[ii, :] = twn
            wannier_integral[ii] = wi
            eigvals[:, :, :, ii] = ev
            wannier_fn[:, :, ii] = wf

    else:
        # get lattice fourier components
        b1, b2, fourier_inds1, fourier_inds2, fourier_comp = get_1d_latt(a=a)
        v_fourier_components = depth * fourier_comp

        # solve problem
        ts_bs, ts_bw, ts_wannier, ts_wannier_n, wannier_integral, eigvals, wannier_fn, v_real, xpts, ypts = \
            solve_lattice(a, a, b1, b2, fourier_inds1, fourier_inds2, v_fourier_components,
                          dx=dx, n1_sites=n1_sites, n2_sites=1, n1_recp_max=n1_recp_max, n2_recp_max=0)

    return ts_bs, ts_bw, ts_wannier, ts_wannier_n, wannier_integral, eigvals, wannier_fn, v_real, xpts, ypts

def solve_lattice(ax, ay, b1, b2, fourier_inds1, fourier_inds2, v_fourier_components,
                  dx=0.05, n1_sites=12, n2_sites=12, n1_recp_max=10, n2_recp_max=10):
    """
    Solve lattice potential for most commonly used quantities including hoppings, wannier integral, etc.

    :param dx: Step size to use in calculating wannier functions
    :param n1_sites: number of k-points in Brillouin zone in direction 1
    :param n2_sites: number of k-points in Brillouin zone in direction 2
    :param n1_recp_max: number of reciprocal lattice vectors in direction 1
    :param n2_recp_max: number of reciprocal lattice vectors in direction 2
    :return:
    ts_bs, ts_wannier, ts_wannier_n, wannier_integral, eigvals2d, wannier_fn
    """

    # get lattice potential in real space
    xpts = np.arange(-3, 3 + dx, dx)
    ypts = xpts

    v_real = get_potential_realspace(b1, b2, v_fourier_components, fourier_inds1, fourier_inds2, xpts, ypts)

    # diagonalize lattice problem
    eigvals, eigvects, bz_vects, recp_vects = \
        lattice_single_particle(b1, b2, v_fourier_components, fourier_inds1, fourier_inds2,
                                n1_sites=n1_sites, n2_sites=n2_sites,
                                n1_recp_max=n1_recp_max, n2_recp_max=n2_recp_max)

    # compute wannier functions for lowest band
    bloch_wfuns, wannier_fn, x_grid, y_grid = get_bloch_wfns(eigvects[:, :, 0], bz_vects, recp_vects, xpts, ypts)

    # compute interaction from wannier function
    wannier_integral = get_interaction_wannier(wannier_fn, x_grid, y_grid)

    # compute hopping from wannier function
    tx, tx_n = get_hopping_wannier(wannier_fn, v_real, xpts, ypts, dx=ax, dy=0)
    ty, ty_n = get_hopping_wannier(wannier_fn, v_real, xpts, ypts, dx=0, dy=ay)
    td, td_n = get_hopping_wannier(wannier_fn, v_real, xpts, ypts, dx=ax, dy=ay)
    t2x, t2x_n = get_hopping_wannier(wannier_fn, v_real, xpts, ypts, dx=2*ax, dy=0)
    t2y, t2y_n = get_hopping_wannier(wannier_fn, v_real, xpts, ypts, dx=0, dy=2*ay)

    ts_wannier = np.array([tx, ty, td, t2x, t2y])
    ts_wannier_n = np.array([tx_n, ty_n, td_n, t2x_n, t2y_n])

    # compute hopping from fit band structure to tight-binding
    kxs2d, kys2d, eigvals2d, _ = oned_2_twod_kvects(bz_vects, eigvals, n1_sites, n2_sites)
    ts_bs, _ = get_hopping_bsfit(kxs2d, kys2d, eigvals2d[:, :, 0], ax=ax, ay=ay)

    # compute hopping from bandwidths
    ts_bw = get_hopping_bw(eigvals2d[:, :, 0], mode='2d')

    return ts_bs, ts_bw, ts_wannier, ts_wannier_n, wannier_integral, eigvals2d, wannier_fn, v_real, xpts, ypts

# plotting functions
def plot_wannier_fns(x_grid, bloch_wfuns, wannier_fn, lattice_depth):
    """
    Plot Wannier functions and bloch_wfuns
    :param x_grid:
    :param bloch_wfuns:
    :param wannier_fn:
    :param lattice_depth:
    :return:
    """

    # TODO: distinguish between 1D and 2D case
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(x_grid.flatten(), bloch_wfuns.real.transpose())
    plt.grid()
    plt.title('Bloch wavefn, real part')

    plt.subplot(2, 2, 2)
    plt.plot(x_grid.flatten(), bloch_wfuns.imag.transpose())
    plt.grid()
    plt.title('Bloch wavefn, imaginary part')

    plt.subplot(2, 2, 3)
    # plt.plot(x_grid.flatten(), wannier_fn.real.flatten())
    # plt.plot([x_grid.flatten()[0], x_grid.flatten()[-1]], [0, 0], 'k--')
    plt.imshow(wannier_fn.real)
    plt.grid()
    plt.title('wannier fn, real part')

    plt.subplot(2, 2, 4)
    # plt.plot(x_grid.flatten(), wannier_fn.imag.flatten())
    # plt.plot([x_grid.flatten()[0], x_grid.flatten()[-1]], [0, 0], 'k--')
    plt.imshow(wannier_fn.imag)
    plt.grid()
    plt.title('wannier fn, imaginary part')

    plt.suptitle('Lattice depth = %0.2f Er' % lattice_depth)

def plot_bands(kxs2d, kys2d, eigvals2d, lattice_depth, bands_to_plot=1):
    """
    Plot band structure with 2D wireframe plot

    :param kxs2d: 2d grid of kxs
    :param kys2d: 2d grid of kys
    :param eigvals2d: 2d matrix of eigvals, E[i,j] = E(kx_j, ky_i)
    :param lattice_depth:
    :param bands_to_plot:
    :return:
    """

    # TODO: detect 1D vs. 2D bands and handle appropriately

    # plot bands
    fig_handle = plt.figure()
    axis = Axes3D(fig_handle)
    for ii in range(0, bands_to_plot):
        axis.plot_wireframe(kxs2d, kys2d, eigvals2d[:, :, ii])
        # ax.scatter(bz_vects[0, :], bz_vects[1, :], eigvals[:, ii])
    plt.xlabel('kx')
    plt.ylabel('ky')
    plt.suptitle('Lattice depth = %0.2f Er' % lattice_depth)

    return fig_handle

def plot_fit_band(fit_params, dispersion_fn, kxs_2d, kys_2d, band_energies_2d):
    """
    Plot fit to band structure, i.e. plot the results from get_hopping_bsfit

    :param fit_params:
    :param dispersion_fn:
    :param kxs_2d:
    :param kys_2d:
    :param band_energies_2d:
    :return:
    """

    # sample fit result
    ks_sample = np.linspace(-np.pi, np.pi, 41)

    # handle 1D or 2D lattices
    if kxs_2d.ndim == 1 and np.all( kys_2d == 0 ):
        kxs_sample, kys_sample = np.meshgrid(ks_sample, np.array([0]))
    elif kxs_2d.ndim == 1 and  np.all( kxs_2d == 0):
        kxs_sample, kys_sample = np.meshgrid(np.array([0]), ks_sample)
    else:
        kxs_sample, kys_sample = np.meshgrid(ks_sample, ks_sample)

    # sample dispersion points
    dispersion_sample = dispersion_fn(fit_params, kxs_sample, kys_sample)

    # plot results
    fig_handle = plt.figure()

    # 1D case
    if kxs_2d.ndim == 1 and np.all( kys_2d == 0 ):
        plt.plot(kxs_2d, band_energies_2d, 'o')
        plt.plot( np.transpose(kxs_sample) , np.transpose(dispersion_sample) )

    elif kxs_2d.ndim == 1 and  np.all( kxs_2d == 0):
        plt.plot(kys_2d, band_energies_2d, 'o')
        plt.plot(np.transpose(kys_sample), np.transpose(dispersion_sample))

    else:
        kxs_sample_gxmg, kys_sample_gxmg, _, dispersion_sample_gxmg, \
        kxs_sample_gymg, kys_sample_gymg, _, dispersion_sample_gymg, linear_coord_sample = \
            twod_2_highsymm(kxs_sample, kys_sample, dispersion_sample)

        # reshape data
        # kxs2d, kys2d, band_energies_2d = oned_2_twod_kvects(bz_vects, band_energies, 11, 11)
        kxs_gxmg, kys_gxmg, _, band_energies_gxmg, kxs_gymg, kys_gymg, _, band_energies_gymg, linear_coord = \
            twod_2_highsymm(kxs_2d, kys_2d, band_energies_2d)

        axis = fig_handle.add_subplot(1, 2, 1, projection='3d')
        # axis = Axes3D(fig_handle)
        # ax.plot_wireframe(kxs2d, kys2d, band_energies_2d)
        axis.scatter(kxs_2d, kys_2d, band_energies_2d)
        axis.plot_wireframe(kxs_sample, kys_sample, dispersion_sample)

        plt.subplot(1, 2, 2)
        plt.plot(linear_coord, band_energies_gxmg, 'bo')
        plt.plot(linear_coord_sample, dispersion_sample_gxmg, 'entries-')
        plt.plot(linear_coord, band_energies_gymg, 'ro')
        plt.plot(linear_coord_sample, dispersion_sample_gymg, 'r-')

        plt.ylabel('Energy (Er)')
        plt.xlabel('Position')
        plt.grid()

if __name__=="__main__":
    pass

