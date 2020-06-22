# ##################################
# Correlators
# ##################################

def getCorrExpVals(self, states, op, species, projector=None, format="fermion", print_results=0):
        """
        TODO: need nstates as argument??? or correct it in here?
        Compute correlators of form <Op Op>_c for all combinations of sites for an arbitrary number of states.
        Returns results in an NSites x NSites x num_states matrix.
        :param states:
        :param op:
        :param species:
        :param projector:
        :param print_results:
        :return:
        """
        # TODO: write this to work for a general geometry?
        nsites = self.geometry.nsites

        if print_results:
            tstart = time.clock()
        if states.ndim is 1:
            states = states[:, None]

        num_states = states.shape[1]
        xs, ys = np.meshgrid(range(0, nsites), range(0, nsites))
        i_corr = xs[xs > ys]
        j_corr = ys[xs > ys]
        nstates = 2 ** nsites
        if projector is None:
            projector = sp.eye(nstates)

        # get expectation values for each site.
        op_exp = self.getSingleSiteOpExpVals(states, op, projector=projector)
        n_corrs = len(i_corr)
        # correlators stored as matrices such that corr[ii,jj] is the correlation between the iith and jjth sites
        corr_exp_mat = np.zeros([nsites, nsites, num_states])

        # two different site correlators
        for ii in range(0, n_corrs):
            full_op = projector * self.get_two_site_op(i_corr[ii], species, j_corr[ii], species, op, op, format=format) * projector.transpose().conj()
            corr_exp_mat[i_corr[ii], j_corr[ii], :] = self.get_exp_vals(states, full_op) - np.multiply(op_exp[i_corr[ii], :], op_exp[j_corr[ii], :])

                    # this matrix has double the amount of info we actually need...but we've only filled in half of it so far
        corr_exp_mat = corr_exp_mat + corr_exp_mat.transpose((1, 0, 2))
        # same site correlator
        for ii in range(0, nsites):
            full_op = projector * self.get_two_site_op(ii, species, ii, species, op, op) * projector.transpose().conj()
            corr_exp_mat[ii, ii, :] = self.get_exp_vals(states, full_op) - op_exp[ii, :] * op_exp[ii, :]

        if print_results:
            tend = time.clock()
            print("getCorrExpVals for %d states took %0.2f s" % (num_states, tend - tstart))

        return np.squeeze(corr_exp_mat)

def computeSiteCorrMatrix(self, states, site, species, op, xsites, ysites, projector=None, cutoff_dist=1.5, format="fermion", print_results=0):
        """
        For faster computation, compute only the correlation matrix at a given site, as opposed to combining
        getCorrExpVals() to get all correlators, and then calling getSiteCorrMatrix() to organize only the ones you
        want. Returns the correlation matrix for the given site and the site index in the correlation matrix (i.e.
        which entry represents <Op(Site) Op(Site)>_c. Can also handle multiple states.
        :param states:
        :param site:
        :param op:
        :param xsites:
        :param ysites:
        :param xlocs:
        :param ylocs:
        :param print_results:
        :return:
        """
        #TODO: redo this function...doesn't make much sense with new class
        xlocs = self.geometry.xlocs
        ylocs = self.geometry.ylocs

        if print_results:
            tstart = time.clock()

        if states.ndim is 1:
            states = states[:, None]

        if isinstance(site, float):
            site = int(site)

        num_states = states.shape[1]

        if projector is None:
            projector = sp.eye(states.shape[0])

        # first, figure out what size the matrix should be.
        Mat = np.zeros([ysites, xsites, num_states])
        # y-matrix position and real space position axes are opposite. Matrix axis is positive going down, but real
        # space axis is positive going up.

        SiteIndexInMat = [int(ylocs.max() - ylocs[site]), int(xlocs[site] - xlocs.min())]

        XDists = xlocs - xlocs[site]
        YDists = ylocs - ylocs[site]
        XDists = XDists.astype(int)
        YDists = YDists.astype(int)
        # loop over states.
        for aa in range(0, num_states):
            # get Sz's, etc.
            Exp_Sites = self.getSingleSiteOpExpVals(states[:, aa], op, projector=projector)

            for ii in range(0, self.geometry.nsites):
                if np.sqrt(XDists[ii] ** 2 + YDists[ii] ** 2) < cutoff_dist:
                    Mat[YDists[ii] - SiteIndexInMat[0], SiteIndexInMat[1] + XDists[ii], aa] = \
                                   self.get_exp_vals(states[:, aa],
                                                     projector * self.get_two_site_op(ii, species, site, species, op, op, format=format) *
                                                     projector.transpose().conj()) - Exp_Sites[ii] * Exp_Sites[site]


        if print_results:
            tend = time.clock()
            print("computeSiteCorrMatrix with site %d for %d states took %0.2f s" % (site, num_states, (tend - tstart)))
        return Mat, SiteIndexInMat

def getSiteCorrMatrix(self, site, full_corr_mat, xsites, ysites):
        """
        For a given site, produces its correlation matrix, where M[ii,jj] = <Sz_Site Sz_(Sitex+ii,Sitey+jj)>_c. Does
        this by reshaping the full correlation matrix produced using getCorrExpVals, which has form
        FullCorrelations[ii,jj] = <Sz_Siteii Sz_Sitejj>_c. Can think of this as transforming from site labelling to 2D
        spatial labelling. We assume that the geometry is rectangular (including chains)
        :param site:
        :param full_corr_mat:
        :param xsites:
        :param ysites:
        :param xlocs:
        :param ylocs:
        :return:
        """
        xlocs = self.geometry.xlocs
        ylocs = self.geometry.ylocs

        # first, figure out what size the matrix should be.
        Mat = np.zeros([ysites, xsites])
        # y-matrix position and real space position axes are opposite. Matrix axis is positive going down, but real
        # space axis is positive going up.
        SiteIndexInMat = [ylocs.max() - ylocs[site], xlocs[site] - xlocs.min()]

        XDists = xlocs - xlocs[site]
        YDists = ylocs - ylocs[site]

        for ii in range(0, self.geometry.nsites):
            Mat[YDists[ii] - SiteIndexInMat[0], SiteIndexInMat[1] + XDists[ii]] = full_corr_mat[ii, site]

        return Mat, SiteIndexInMat

def plotQuench(self, times, evolved_states, central_site, xsites, ysites, rabi, detune):
    """
    Plot results of quench
    :param times:
    :param evolved_states:
    :param central_site:
    :param xsites:
    :param ysites:
    :param rabi:
    :param detune:
    :return:
    """
    # TODO: this doesn't belong in this class..
    nsites = xsites * ysites

    # find special states that we want to compare to
    SpStates, Desc = self.find_special_states(xsites, ysites)
    # TODO: find these states using their descriptions
    # ['AllUp','AllDn','AFM1','AFM2','AllPlus','AllMinus']
    AllUpState = SpStates[:, 0][:, None]
    AllDnState = SpStates[:, 1][:, None]
    AFMState1BasisState = SpStates[:, 2][:, None]
    AFMState2BasisState = SpStates[:, 3][:, None]
    AllPlusState = SpStates[:, 4][:, None]
    AllMinusState = SpStates[:, 5][:, None]

    f, axarray = plt.subplots(5, 2)
    # non interacting system results...see e.g. Laser Cooling and Trapping, Metcalf and Van der
    # Straaten (chapter 1.2)
    EffRabi = np.sqrt(detune ** 2 + rabi ** 2)
    gamp = lambda t: np.multiply(np.cos(EffRabi * t / 2.) + 1j * detune / EffRabi * np.sin(EffRabi * t / 2.),
                                 np.exp(-1j * detune * t / 2.))
    ramp = lambda t: 1j * rabi / EffRabi * np.multiply(np.sin(EffRabi * t / 2.),
                                                       np.exp(1j * detune * t / 2.))
    # Probabilities in each state
    # gfrac = lambda t: 1.-(self.Rabi)**2/((self.Detune)**2+(self.Rabi)**2)* \
    # np.square(np.sin(np.sqrt((self.Rabi)**2+(self.Detune)**2)*t/2.))
    # rfrac = lambda t: 1. - gfrac(t)
    rfrac = lambda t: np.square(np.abs(ramp(t)))
    gfrac = lambda t: np.square(np.abs(gamp(t)))  # 1. - rfrac(t)
    plusfrac = lambda t: 0.5 * np.square(np.abs(gamp(t) + ramp(t)))
    minusfrac = lambda t: 0.5 * np.square(np.abs(gamp(t) - ramp(t)))
    interp_times = np.linspace(min(times), max(times), 300)
    # overlap with |gggg> and |rrrr>
    NonInt_gggfrac = np.power(gfrac(interp_times), float(nsites))
    NonInt_rrrfrac = np.power(rfrac(interp_times), nsites)
    GGFrac = np.square(np.abs(self.get_overlaps(evolved_states, AllDnState)))
    RRFrac = np.square(np.abs(self.get_overlaps(evolved_states, AllUpState)))

    axarray[0, 0].plot(times, GGFrac, 'bx')
    axarray[0, 0].plot(interp_times, NonInt_gggfrac, 'entries')
    axarray[0, 0].set_ylabel('|ggg...>')
    axarray[0, 0].set_ylim([-0.05, 1.05])
    axarray[0, 1].plot(times, RRFrac, 'rx')
    axarray[0, 1].plot(interp_times, NonInt_rrrfrac, 'r')
    axarray[0, 1].set_ylabel('|rrr...>')
    axarray[0, 1].set_ylim([-0.05, 1.05])

    # overlap of central site with |g> and |r>
    CSite_gfrac = self.get_exp_vals(evolved_states, self.get_single_site_op(central_site, 0, self.ng, format="boson"))
    CSite_rfrac = self.get_exp_vals(evolved_states, self.get_single_site_op(central_site, 0, self.nr, format="boson"))
    axarray[1, 0].plot(times, CSite_gfrac, 'bx')
    axarray[1, 0].plot(interp_times, gfrac(interp_times), 'entries')
    axarray[1, 0].set_ylabel('|g_c>')
    axarray[1, 0].set_ylim([-0.05, 1.05])
    axarray[1, 1].plot(times, CSite_rfrac, 'rx')
    axarray[1, 1].plot(interp_times, rfrac(interp_times), 'r')
    axarray[1, 1].set_ylabel('|r_c>')
    axarray[1, 1].set_ylim([-0.05, 1.05])

    # overlap with |AFM1> and |AFM2>
    AFM1_frac = np.square(np.abs(self.get_overlaps(evolved_states, AFMState1BasisState)))
    AFM2_frac = np.square(np.abs(self.get_overlaps(evolved_states, AFMState2BasisState)))
    axarray[2, 0].plot(times, AFM1_frac, 'bx')
    axarray[2, 0].set_ylabel('|AFM1>')
    axarray[2, 0].set_ylim([-0.05, 1.05])
    axarray[2, 1].plot(times, AFM2_frac, 'rx')
    axarray[2, 1].set_ylabel('|AFM2>')
    axarray[2, 1].set_ylim([-0.05, 1.05])

    # overlap with |+++> and |--->
    NonInt_pppfrac = np.power(plusfrac(interp_times), nsites)
    NonInt_mmmfrac = np.power(minusfrac(interp_times), nsites)
    PlusPlus_frac = np.square(np.abs(self.get_overlaps(evolved_states, AllPlusState)))
    MinusMinus_frac = np.square(np.abs(self.get_overlaps(evolved_states, AllMinusState)))
    axarray[3, 0].plot(times, PlusPlus_frac, 'bx')
    axarray[3, 0].plot(interp_times, NonInt_pppfrac, 'entries')
    axarray[3, 0].set_ylabel('|+++>')
    axarray[3, 0].set_ylim([-0.05, 1.05])
    axarray[3, 1].plot(times, MinusMinus_frac, 'rx')
    axarray[3, 1].plot(interp_times, NonInt_mmmfrac, 'r')
    axarray[3, 1].set_ylabel('|--->')
    axarray[3, 1].set_ylim([-0.05, 1.05])

    # overlap of central site with |+> and |->
    CSite_plusfrac = self.get_exp_vals(evolved_states,
                                       self.get_single_site_op(central_site, 0, self.nplus, format="boson"))
    CSite_minusfrac = self.get_exp_vals(evolved_states,
                                        self.get_single_site_op(central_site, 0, self.nminus, format="boson"))
    axarray[4, 0].plot(times, CSite_plusfrac, 'bx')
    axarray[4, 0].plot(interp_times, plusfrac(interp_times), 'entries')
    axarray[4, 0].set_ylabel('|+_c>')
    axarray[4, 0].set_ylim([-0.05, 1.05])
    axarray[4, 1].plot(times, CSite_minusfrac, 'rx')
    axarray[4, 1].plot(interp_times, minusfrac(interp_times), 'r')
    axarray[4, 1].set_ylabel('|-_c>')
    axarray[4, 1].set_ylim([-0.05, 1.05])

    # final touches
    axarray[4, 0].set_xlabel('Time')
    axarray[4, 1].set_xlabel('Time')
    plt.draw()
    return f
