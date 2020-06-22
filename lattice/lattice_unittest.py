import unittest

# run unit tests from the command prompt via python -m unittest hubbard_test

import numpy as np
import lattice_single_particle as latt

class TestLattice(unittest.TestCase):

    def setUp(self):
        pass

    def test_1d_latt_bw(self):
        """
        Compare the 1d lattice bandwidth with the exact result from the Mathieu equation.

        A very good reference for understanding this is "Approximate solutions to Mathieu's equation" arXiv:1710.006571
        :return:
        """
        depth = 3
        n_qpts = 501
        n_recp_pts = 201

        # #############################################
        # solve 1D lattices
        # #############################################
        v_fourier_components = depth * np.array([0, -0.25, -0.25])

        v_fourier_indices1 = [0, -1, 1]
        v_fourier_indices2 = [0, 0, 0]
        recp_basis_vect1 = np.array([[2. * np.pi], [0.]])
        recp_basis_vect2 = np.array([[0.], [2. * np.pi]])

        # diagonalize lattice problems
        eigvals, eigvects, bz_vects, recp_vects = \
            latt.lattice_single_particle(recp_basis_vect1, recp_basis_vect2, v_fourier_components,
                                         v_fourier_indices1, v_fourier_indices2,
                                         n1_sites=n_qpts, n2_sites=0, n1_recp_max=n_recp_pts, n2_recp_max=0)

        # bandwidths from diagonalization
        bws = np.max(eigvals, axis=0) - np.min(eigvals, axis=0)

        # exact solution
        nbands = bws.size
        bws_mathieu = np.zeros(nbands)
        for ii in range(0, nbands):
            bws_mathieu[ii] = latt.get_1d_bandwidth_exact(depth, ii)

        fractional_dev = np.divide(np.abs(bws - bws_mathieu), bws)

        # first two bands 0.01%
        self.assertTrue( np.round(np.abs(bws[0] - bws_mathieu[0]), 4) == 0)
        self.assertTrue( np.round(np.abs(bws[1] - bws_mathieu[1]), 4) == 0)
        # test all considered bands to better than 1%
        self.assertTrue( np.round(np.max(fractional_dev), 2) == 0)

    def test_1d_latt_wavefn(self):
        """
        Test the wavefunctions of the 1D lattice at k=0 point using Mathieu solution
        :return:
        """
        depth = 3
        n_qpts = 501
        n_recp_pts = 201

        xpts = np.linspace(-1, 1, 150)
        dx = xpts[1] - xpts[0]
        ypts = np.array([0])

        # #############################################
        # solve 1D lattices
        # #############################################
        # for red detuned case, where potential V(x) = -s * sin^2(pi*x) = - s/2 * ( 1 - cos(2*pi*x) )
        # i.e. Potential has maximum at x = 0 and Fourier components are [-0.5, 0.25, 0.25]
        # mathieu parameters are q = s / 4 and a = E
        # for blue detune case V(x) = s * sin^2(pi*x) = s/2 * (1 - cos(2*pi*x) )
        # Potential has minimum at x=0 and Fourier components are [0.5, -0.25, -0.25]
        # mathieu parameters are q = -s/4 and a = E

        # v_fourier_components = depth * np.array([-0.5, 0.25, 0.25])
        #detune = "red"

        v_fourier_components = depth * np.array([0.5, -0.25, -0.25])
        detune = "blue"


        v_fourier_indices1 = [0, -1, 1]
        v_fourier_indices2 = [0, 0, 0]
        recp_basis_vect1 = np.array([[2. * np.pi], [0.]])
        recp_basis_vect2 = np.array([[0.], [2. * np.pi]])

        # diagonalize lattice problems
        eigvals, eigvects, bz_vects, recp_vects = \
            latt.lattice_single_particle(recp_basis_vect1, recp_basis_vect2, v_fourier_components,
                                         v_fourier_indices1, v_fourier_indices2,
                                         n1_sites=n_qpts, n2_sites=0, n1_recp_max=n_recp_pts, n2_recp_max=0)

        # get bloch wavefunctions for q=0 and several bands
        # Above band six seems order of bands can get easily confused.
        nbands = 6
        center_wavefns = np.zeros((xpts.size, nbands))
        wannier_fns_mag = np.zeros((xpts.size, nbands))
        i_center = int((n_qpts - 1) / 2)
        for ii in range(0, nbands):
            bloch_wavefns, wannier_fn, _, _ = latt.get_bloch_wfns(eigvects[:, :, ii], bz_vects, recp_vects, xpts, ypts)

            # bloch function should be purely real or imaginary
            wavefn = bloch_wavefns[i_center, :]
            if np.sum( np.abs(np.imag(wavefn)) ) < np.sum( np.abs(np.real(wavefn))):
                wavefn = np.real(wavefn)
            else:
                wavefn = np.imag(wavefn)

            center_wavefns[:, ii] = wavefn
            # normalize bloch fn
            norm = np.sum(center_wavefns[:, ii]**2) * dx
            center_wavefns[:, ii] = center_wavefns[:, ii] / np.sqrt(norm)

            # wannier fn
            wannier_fns_mag[:, ii] = np.abs(wannier_fn)**2

        # get mathieu wavefunctions for q=0
        mathieu_wavefns = np.zeros((xpts.size, nbands))
        for ii in range(0, nbands):
            mathieu_wavefns[:, ii] = latt.get_1d_lattice_zero_momentum_wavefn_exact(depth, ii, xpts, detuning=detune)
            if np.sum(mathieu_wavefns[:, ii] * center_wavefns[:, ii]) < 0:
                mathieu_wavefns[:, ii] = -1 * mathieu_wavefns[:, ii]

        # plot results for testing purposes
        # nrows = np.floor( np.sqrt(nbands) )
        # ncols = np.ceil(nbands / nrows)
        #
        # v, _, _ = latt.get_potential_realspace(b1, b2, v_fourier_components, v_fourier_inds1, v_fourier_inds2, xpts, ypts)
        # v = v[0, :]
        #
        # import matplotlib.pyplot as plt
        # plt.figure()
        # for ii in range(0, nbands):
        #     plt.subplot(nrows, ncols, ii + 1)
        #     plt.plot(xpts, center_wavefns[:, ii])
        #     plt.plot(xpts, mathieu_wavefns[:, ii])
        #     plt.plot(xpts, wannier_fns_mag[:, ii])
        #     plt.plot(xpts, v)
        #     if ii == 0:
        #         plt.legend(['numerical', 'mathieu', '|wannier fn|^2', 'v(x)'])
        #     plt.title('band %d' % ii)


        max_diff = np.max( np.abs(center_wavefns - mathieu_wavefns) )
        self.assertTrue( np.round(max_diff, 10) == 0)

    def test_1d_latt_different_spacing(self):
        """
        Test that working with different length units does not effect the results of the bandstructure
        calculation. Test that changing the lattice constant produces the same results.
        :return:
        """
        depth = 3
        n_qpts = 101
        n_recp_pts = 51

        # #############################################
        # solve 1D lattices
        # #############################################
        v_fourier_components = depth * np.array([0, -0.25, -0.25])
        v_fourier_indices1 = [0, -1, 1]
        v_fourier_indices2 = [0, 0, 0]


        # diagonalize lattice problems
        a = 1.
        b1 = np.array([[2. * np.pi / a], [0.]])
        b2 = np.array([[0.], [2. * np.pi / a]])
        eigvals, eigvects, bz_vects, recp_vects = \
            latt.lattice_single_particle(b1, b2, v_fourier_components,
                                         v_fourier_indices1, v_fourier_indices2,
                                         n1_sites=n_qpts, n2_sites=0, n1_recp_max=n_recp_pts, n2_recp_max=0)

        # use different spacing
        a = 2.
        b1 = np.array([[2. * np.pi / a], [0.]])
        b2 = np.array([[0.], [2. * np.pi / a]])
        eigvals2, eigvects2, bz_vects2, recp_vects2 = \
            latt.lattice_single_particle(b1, b2, v_fourier_components,
                                         v_fourier_indices1, v_fourier_indices2,
                                         n1_sites=n_qpts, n2_sites=0, n1_recp_max=n_recp_pts, n2_recp_max=0)

        self.assertTrue(np.array_equal(eigvals, eigvals2))
        self.assertTrue(np.array_equal(eigvects, eigvects2))

    def test_1d_latt_hopping(self):
        """
        Compare hopping calculated using band structure versus using wannier function

        :return:
        """

        for depth in range(5, 15):
            n_qpts = 30
            n_recp_pts = 31

            dx = 0.005
            xpts = np.arange(-3, 3 + dx, dx)
            ypts = np.array([0])

            # #############################################
            # solve 1D lattices
            # #############################################
            v_fourier_comps = depth * np.array([0, -0.25, -0.25])

            fourier_inds1 = [0, -1, 1]
            fourier_inds2 = [0, 0, 0]
            b1 = np.array([[2. * np.pi], [0.]])
            b2 = np.array([[0.], [2. * np.pi]])

            # diagonalize lattice problems
            eigvals, eigvects, bz_vects, recp_vects = \
                latt.lattice_single_particle(b1, b2, v_fourier_comps, fourier_inds1, fourier_inds2,
                                             n1_sites=n_qpts, n2_sites=0, n1_recp_max=n_recp_pts, n2_recp_max=0)

            bloch_wavefns, wannier_fn, _, _ = latt.get_bloch_wfns(eigvects[:, :, 0], bz_vects, recp_vects, xpts, ypts)

            potential = latt.get_potential_realspace(b1, b2, v_fourier_comps, fourier_inds1, fourier_inds2, xpts, ypts)
            ts_wann, _ = latt.get_hopping_wannier(wannier_fn, potential, xpts, ypts, dx=1., dy=0)

            kxs2d, kys2d, eigvals2d, _ = latt.oned_2_twod_kvects(bz_vects, eigvals[:, 0], n_qpts, 1)
            ts_bs, _ = latt.get_hopping_bsfit(kxs2d, kys2d, eigvals2d)

            # print(np.round(np.abs(ts_bs[0] - ts_wann[0]), 5))
            self.assertTrue( np.round(np.abs(ts_bs[0] - ts_wann[0]), 5) == 0, 'hopping calculated from band structure vs. wannier function differed by more than 0.5e-5 for depth %d' % depth)



    # 2D lattice tests
    def test_2d_separable_latt(self):
        """
        Compare 2D separable square lattice with 1D lattice.
        """

        # depth of the 1d lattices in Er
        depth_a = 3
        depth_b = 2

        # set number of reciprocal lattice vectors, points in the brillouin zone, and real space points
        n_qpts = 11
        n_recp_pts = 6

        # #############################################
        # solve 1D lattices
        # #############################################
        v_fourier_components_latt_a = depth_a * np.array([0, -0.25, -0.25])
        v_fourier_components_latt_b = depth_b * np.array([0, -0.25, -0.25])

        v_fourier_indices1 = [0, -1, 1]
        v_fourier_indices2 = [0, 0, 0]
        recp_basis_vect1 = np.array([[2. * np.pi], [0.]])
        recp_basis_vect2 = np.array([[0.], [2. * np.pi]])

        # diagonalize lattice problems
        eigvals1d_a, eigvects1d_a, bz_vects1d, recp_vects1d = \
            latt.lattice_single_particle(recp_basis_vect1, recp_basis_vect2, v_fourier_components_latt_a, v_fourier_indices1,
                                    v_fourier_indices2, n1_sites=n_qpts, n2_sites=0, n1_recp_max=n_recp_pts, n2_recp_max=0)

        eigvals1d_b, eigvects1d_b, bz_vects1d, recp_vects1d = \
            latt.lattice_single_particle(recp_basis_vect1, recp_basis_vect2, v_fourier_components_latt_b, v_fourier_indices1,
                                         v_fourier_indices2, n1_sites=n_qpts, n2_sites=0, n1_recp_max=n_recp_pts, n2_recp_max=0)

        # #############################################
        # separable square 2D lattice
        # #############################################
        v_fourier_components = np.array([0, -0.25 * depth_a, -0.25 * depth_a, -0.25 * depth_b, -0.25 * depth_b])
        v_fourier_indices1 = [0, -1, 1, 0, 0]
        v_fourier_indices2 = [0, 0, 0, -1, 1]
        recp_basis_vect1 = np.array([[2. * np.pi], [0.]])
        recp_basis_vect2 = np.array([[0.], [2. * np.pi]])

        # diagonalize lattice problem
        eigvals2d, eigvects2d, bz_vects2d, recp_vects2d = \
            latt.lattice_single_particle(recp_basis_vect1, recp_basis_vect2, v_fourier_components, v_fourier_indices1,
                                    v_fourier_indices2, n1_sites=n_qpts, n2_sites=n_qpts, n1_recp_max=n_recp_pts,
                                    n2_recp_max=n_recp_pts)

        # #############################################
        # compare energy bands
        # #############################################
        # E_{sums, ny}(kx, ky) = E_nx(kx) + E_ny(ky)
        exex, eyey = np.meshgrid(np.ravel(eigvals1d_a), np.ravel(eigvals1d_b))
        E_from_1d = np.sort(np.ravel(exex + eyey))
        E_2d = np.sort(np.ravel(eigvals2d))

        max_diff = np.max(np.abs(E_from_1d - E_2d))
        self.assertTrue( np.round(max_diff, 10) == 0)

    def test_4fold_noretro(self):
        """
        Test 4-fold lattice potential with no retroreflection vs. 1D lattice calculation. This is mostly
        an exercise in changing units to get the appropriate Er
        :return:
        """
        import numpy as np
        import lattice_single_particle as latt
        import matplotlib.pyplot as plt

        depth = 3
        n_qpts = 21
        n_recp_pts = 10

        ax, ay, b1, b2, fourier_inds1, fourier_inds2, v_fourier = latt.get_fourfold_latt(theta=75., r=0.)

        er = np.sqrt(np.sum(np.square(b1)) * np.sum(np.square(b2))) / 4
        erx = np.sum(np.square(b1)) / 4
        ery = np.sum(np.square(b2)) / 4
        # want lattice to have depth ery, but to write it in units of er
        v_fourier = v_fourier * ery / er

        eigvals2d, eigvects2d, bz_vects2d, recp_vects2d = \
            latt.lattice_single_particle(b1, b2, depth * v_fourier, fourier_inds1, fourier_inds2,
                                         n1_sites=0, n2_sites=n_qpts, n1_recp_max=0, n2_recp_max=n_recp_pts)
        # want output in units of ery
        eigvals2d = eigvals2d * er / ery

        b1_1d = np.array([2 * np.pi / ay, 0])
        b2_1d = np.array([0, 2 * np.pi / ay])
        v_fourier1d = depth * np.array([-0.5, -0.25, -0.25])
        eigvals1d, eigvects1d, bz_vects1d, recp_vects1d = \
            latt.lattice_single_particle(b1_1d, b2_1d, v_fourier1d, [0, 0, 0], [0, 1, -1],
                                         n1_sites=0, n2_sites=n_qpts, n1_recp_max=0, n2_recp_max=n_recp_pts)

        max_diff = np.max(np.abs(eigvals1d - eigvals2d))
        self.assertTrue( np.round(max_diff, 11) == 0)

        # why don't eigenvectors match
        # max_diff_eigvect = np.max(np.abs(eigvects1d) - np.abs(eigvects2d))
        # self.assertTrue( np.round(max_diff_eigvect, 11) == 0)

    def test_4fold_hopping(self):
        """
        Test hopping calculation on fourfold lattice with some non-idealities to break x-y symmetry. Ensure results agree
        to within 1% for nearest-neighbor and diagonal neighbor
        :return:
        """

        theta = 91.6267
        r = 0.6
        # theta = 90
        # r = 1
        alpha = 0
        depth = 4

        hoppings_bs, hoppings_bw, hoppings_wannier, hoppings_wannier_dens_dep, wannier_integral, energies, wannier_fn, potential, xpts, ypts = \
            latt.solve_d4_lattice(depth, theta, r, alpha, dx=0.05, n1_sites=13, n2_sites=13, n1_recp_max=10, n2_recp_max=10, normalize_er=True)

        tx_fr = np.abs(hoppings_bs[0] - hoppings_wannier[0]) / hoppings_bs[0]
        ty_fr = np.abs(hoppings_bs[1] - hoppings_wannier[1]) / hoppings_bs[1]
        td_fr = np.abs(hoppings_bs[2] - hoppings_wannier[2]) / hoppings_bs[2]
        t2x_fr = np.abs(hoppings_bs[3] - hoppings_wannier[3]) / np.abs(hoppings_bs[3])
        t2y_fr = np.abs(hoppings_bs[4] - hoppings_wannier[4]) / np.abs(hoppings_bs[4])

        self.assertTrue(np.round(tx_fr, 2) == 0)
        self.assertTrue(np.round(ty_fr, 2) == 0)
        self.assertTrue(np.round(td_fr, 2) == 0)



    # other tests
    def test_real_space_potential(self):
        """
        Test function that returns real space potential from the Fourier components
        :return:
        """
        depth = 3

        xpts = np.linspace(-0.5, 0.5, 100)
        ypts = np.array([0])

        # #############################################
        # real space
        # #############################################
        vfn = lambda x: - 0.5 * depth * ( 1 - np.cos(2 * np.pi * x) )
        vfn_eval = vfn(xpts)

        # #############################################
        # fourier space
        # #############################################
        v_fourier_components = depth * np.array([-0.5, 0.25, 0.25])

        v_fourier_indices1 = [0, -1, 1]
        v_fourier_indices2 = [0, 0, 0]
        recp_basis_vect1 = np.array([[2. * np.pi], [0.]])
        recp_basis_vect2 = np.array([[0.], [2. * np.pi]])
        v = latt.get_potential_realspace(recp_basis_vect1, recp_basis_vect2, v_fourier_components,
                                               v_fourier_indices1, v_fourier_indices2, xpts, ypts)
        v = v[0, :]

        max_diff = np.max( np.abs(v - vfn_eval))
        self.assertTrue( np.round(max_diff, 14) == 0)

if __name__ == '__main__':
    unittest.main()