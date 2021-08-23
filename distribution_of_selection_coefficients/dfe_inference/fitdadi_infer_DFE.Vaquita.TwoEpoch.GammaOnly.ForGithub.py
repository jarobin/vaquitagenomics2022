"""
Uses dadi and fitdadi to infer the DFE of a given synonynmous sfs.

JCM 201907011 -- code originally by Jon Mah
ACB 20201104 -- modified by Annabel Beichman
"""

import sys

sys.path.append(
    "/Users/annabelbeichman/Documents/UW/general_software/fitdadi/dadi")  # need to add path to import Selection.py, wherever it is
import os
import logging
import time
import argparse
import warnings
import numpy
import dadi
import Selection
import scipy.stats.distributions
import scipy.integrate
import scipy.optimize


# class ArgumentParserNoArgHelp(argparse.ArgumentParser):
#     """Like *argparse.ArgumentParser*, but prints help when no arguments."""
#
#     def error(self, message):
#         """Print error message, then help."""
#         sys.stderr.write('error: %s\n\n' % message)
#         self.print_help()
#         sys.exit(2)


class DemographicAndDFEInference():
    """Wrapper class to allow functions to reference each other."""

    def ExistingFile(self, fname):
        """Return *fname* if existing file, otherwise raise ValueError."""
        if os.path.isfile(fname):
            return fname
        else:
            raise ValueError("%s must specify a valid file name" % fname)

    # def inferDFEParser(self):
    #     """Return *argparse.ArgumentParser* for ``fitdadi_infer_DFE.py``."""
    #     parser=argparse.ArgumentParser(description='Infer demography and the DFE under a 3 epoch model')
    #     parser.add_argument(
    #         'syn_input_sfs', type=self.ExistingFile,
    #         help=('Synonynomous site-frequency spectrum from which the '
    #               'demographic parameters should be inferred.'))
    #     parser.add_argument(
    #         'Lcds', type=float,
    #         help=(
    #             'The number of called cds (not "exon" because that includes UTRs) sites. This will get scaled by the supplied NS_S_ScalingFactor to get Lsyn and Lns and be used to scale thetasyn.'))
    #     parser.add_argument(
    #         'nonsyn_input_sfs', type=self.ExistingFile,
    #         help=('Nonsynonynomous site-frequency spectrum from which the '
    #               'distribution of fitness effects should be inferred.'))
    #     parser.add_argument(
    #         'outprefix', type=str,
    #         help='The file prefix for the output `*inferred_demography.txt`.')
    #     parser.add_argument(
    #         'mutationRate', type=float,
    #         help='The genome-wide mutation rate (mut/bp/gen)')
    #     parser.add_argument(
    #         'exonicMutationRateScalingFactorexonicMutationRateScalingFactor', type=float,
    #         help='The scaling factor to scale up the exonic mutation rate relative to the genome-wide rate (e.g. 1.25)')
    #     parser.add_argument(
    #         'NS_S_ScalingFactor', type=float,
    #         help='The scaling factor to go from thetasyn to thetanonsyn (e.g. 2.31 for humans) to account for greater seq len of NS sites and differing GC content between SYN and NS sites.')
    #     return parser

    def snm(self, notused, ns, pts):
        """Return a standard neutral model.

        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def two_epoch(self, params, ns, pts):
        """Define a two-epoch demography, i.e., an instantaneous size change.

        params = (nu, T)
            nu: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nu, T = params  # Define given parameters.
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        phi = dadi.Integration.one_pop(phi, xx, T, nu)  # Integrate.

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def two_epoch_sel(self, params, ns, pts):
        """Define a two-epoch demography, i.e., an instantaneous size change.

        This method incorporates a gamma parameter.

        params = (nu, T, gamma)
            nu: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
            gamma: Parameter tuple describing a gamma distribution.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nu, T, gamma = params  # Define given parameters.
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)  # Define initial phi

        phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def growth(self, params, ns, pts):
        """Exponential growth beginning some time ago.

        params = (nu, T)
            nu: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nu, T = params  # Define given parameters.
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        def nu_func(t): return numpy.exp(numpy.log(nu) * t / T)  # Exp growth.

        phi = dadi.Integration.one_pop(phi, xx, T, nu_func)  # Integrate.

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def bottlegrowth(self, params, ns, pts):
        """Instantaneous size change followed by exponential growth.

        params = (nuB, nuF, T)
            nuB: Ratio of population size after instantaneous change to ancient
                population size.
            nuF: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nuB, nuF, T = params  # Define given parameters.

        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi

        # Exponential growth function
        def nu_func(t): return nuB * numpy.exp(numpy.log(nuF / nuB) * t / T)

        phi = dadi.Integration.one_pop(phi, xx, T, nu_func)  # Integrate.

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def three_epoch(self, params, ns, pts):
        """Define a three-epoch demography.

        params = (nuB, nuF, TB, TF)
            nuB: Ratio of bottleneck population size to ancient
                population size.
            nuF: Ratio of contemporary to ancient population size.
            TB: Length of bottleneck, in units of 2 * N_a.
            TF: Time since bottleneck recovery, in units of 2 * N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nuB, nuF, TB, TF = params  # Define given parameters.

        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        phi = dadi.Integration.one_pop(phi, xx, TB, nuB)  # Integrate 1 to 2.
        phi = dadi.Integration.one_pop(phi, xx, TF, nuF)  # Integrate 2 to 3.

        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    # Annabel added three_epoch_sel:
    def three_epoch_sel(self, params, ns, pts):
        """Define a three-epoch demography. This function includes selection with a GAMMA distribution


        params = (nuB, nuF, TB, TF,gamma)
            nuB: Ratio of bottleneck population size to ancient
                population size.
            nuF: Ratio of contemporary to ancient population size.
            TB: Length of bottleneck, in units of 2 * N_a.
            TF: Time since bottleneck recovery, in units of 2 * N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
            gamma: Parameter tuple describing a gamma distribution.
        pts: Number of grid points to use in integration.
        """
        nuB, nuF, TB, TF, gamma = params  # Define given parameters.

        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)  # Define initial phi.

        phi = dadi.Integration.one_pop(phi, xx, TB, nuB, gamma=gamma)  # Integrate 1 to 2.
        phi = dadi.Integration.one_pop(phi, xx, TF, nuF, gamma=gamma)  # Integrate 2 to 3.

        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def four_epoch(self, params, ns, pts):
        """Define a four-epoch demography.

        params = (Na, Nb, Nc, Ta, Tb, Tc)
            Na: ratio of population size between epoch 1 and 2.
            Nb: ratio of population size between epoch 2 and 3.
            Nc: ratio of population size between epoch 3 and 4.
            Ta: Bottleneck length between epoch 1 and 2, in units of 2 * N_a.
            Tb: Length of bottleneck between epoch 2 and 3,
                in units of 2 * N_a.
            Tc: Length of bottleneck between epoch 3 and 4,
                in units of 2 * N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        Na, Nb, Nc, Ta, Tb, Tc = params  # Define given parameters.

        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        # Integrate epochs.
        phi = dadi.Integration.one_pop(phi, xx, Ta, Na)  # Integrate 1 to 2.
        phi = dadi.Integration.one_pop(phi, xx, Tb, Nb)  # Integrate 2 to 3.
        phi = dadi.Integration.one_pop(phi, xx, Tc, Nc)  # Integrate 3 to 4.

        # Construct spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def gamma_dist(self, mgamma, alpha, beta):
        """Define a gamma distribution.

        self: reference to this instance of a gamma distribution.
        mgamma: float which describes the mean value of this gamma
             distribution.
        alpha: shape parameter of the gamma distribution.
        beta: scale parameter of the gamma distribution.
        """
        return scipy.stats.distributions.gamma.pdf(-mgamma, alpha, scale=beta)

    def neugamma(self, mgamma, pneu, alpha, beta):
        """Define a neutral-gamma distribution.

        self: reference to this instance of a neutral-gamma distribution.
        mgamma: float which describes the mean value of this neutral-gamma
            distribution.
        pneu: proportion of elements which are assumed to be neutral, i.e.,
            equal to 0.
        alpha: shape parameter of the non-neutral elements of the
            neutral-gamma distribution.
        beta: scale parameter of the non-neutral elements of the
            neutral-gamma distribution.
        """
        mgamma = -mgamma
        # Assume anything with gamma < 1e-4 is neutral
        if (0 <= mgamma) and (mgamma < 1e-4):
            return pneu / (1e-4) + (1 - pneu) * self.gamma_dist(
                -mgamma, alpha, beta)
        else:
            return self.gamma_dist(-mgamma, alpha, beta) * (1 - pneu)

    def main(self):
        """Execute main function."""
        # Parse command line arguments
        parser = argparse.ArgumentParser(description='Infer demography and the DFE under a 3 epoch model')
        parser.add_argument(
            '--syn_input_sfs', type=self.ExistingFile,
            help=('Synonynomous site-frequency spectrum from which the '
                  'demographic parameters should be inferred.'))
        parser.add_argument(
            '--Lcds', type=float,
            help=(
                'The number of called cds (not "exon" because that includes UTRs) sites. This will get scaled by the supplied NS_S_ScalingFactor to get Lsyn and Lns and be used to scale thetasyn.'))
        parser.add_argument(
            '--nonsyn_input_sfs', type=self.ExistingFile,
            help=('Nonsynonynomous site-frequency spectrum from which the '
                  'distribution of fitness effects should be inferred.'))
        parser.add_argument(
            '--outprefix', type=str,
            help='The file prefix for the output `*inferred_demography.txt`.')
        parser.add_argument(
            '--mutationRateList', type=str,
            help='List of mutation rates you want to test (mut/bp/gen); must be "," delimited')
        parser.add_argument(
            '--exonicMutationRateScalingFactorexonicMutationRateScalingFactor', type=float,
            help='The scaling factor to scale up the exonic mutation rate relative to the genome-wide rate (e.g. 1.25)')
        parser.add_argument(
            '--NS_S_ScalingFactor', type=float,
            help='The scaling factor to go from thetasyn to thetanonsyn (e.g. 2.31 for humans) to account for greater seq len of NS sites and differing GC content between SYN and NS sites.')
        # parser = self.inferDFEParser
        args = parser.parse_args()
        prog = parser.prog

        # Assign arguments
        syn_input_sfs = args.syn_input_sfs
        nonsyn_input_sfs = args.nonsyn_input_sfs
        outprefix = args.outprefix
        mutationRate_list = map(float, args.mutationRateList.split(",")) # makes a list of floats from str input ; must be , delim
        exonicMutationRateScalingFactor = float(args.exonicMutationRateScalingFactorexonicMutationRateScalingFactor)
        NS_S_ScalingFactor = float(args.NS_S_ScalingFactor)
        Lcds = float(args.Lcds)

        # Numpy options
        numpy.set_printoptions(linewidth=numpy.inf)

        # create output directory if needed
        outdir = os.path.dirname(args.outprefix)
        if outdir:
            if not os.path.isdir(outdir):
                if os.path.isfile(outdir):
                    os.remove(outdir)
                os.mkdir(outdir)

        # Output files: logfile
        # Remove output files if they already exist
        underscore = '' if outprefix[-1] == '/' else '_'
        inferred_demography = \
            '{0}{1}inferred_demography.txt'.format(
                args.outprefix, underscore)
        inferred_DFE = \
            '{0}{1}inferred_DFE.txt'.format(outprefix, underscore)
        most_deleterious_spectrum = \
            '{0}{1}mostDeleteriousExpectedSpectrum.txt'.format(outprefix, underscore)
        logfile = '{0}{1}log.log'.format(outprefix, underscore)
        to_remove = [logfile, inferred_demography, inferred_DFE]
        for f in to_remove:
            if os.path.isfile(f):
                os.remove(f)

        # Set up to log everything to logfile.
        logging.shutdown()
        logging.captureWarnings(True)
        logging.basicConfig(
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.INFO)
        logger = logging.getLogger(prog)
        warning_logger = logging.getLogger("py.warnings")
        logfile_handler = logging.FileHandler(logfile)
        logger.addHandler(logfile_handler)
        warning_logger.addHandler(logfile_handler)
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
        logfile_handler.setFormatter(formatter)
        logger.setLevel(logging.INFO)

        # print some basic information
        logger.info('Beginning execution of {0} in directory {1}\n'.format(
            prog, os.getcwd()))
        logger.info('Progress is being logged to {0}\n'.format(logfile))
        logger.info('Parsed the following arguments:\n{0}\n'.format(
            '\n'.join(['\t{0} = {1}'.format(*tup) for tup in vars(args).items()])))

        # Construct initial Spectrum object from input synonymous sfs.
        syn_data = dadi.Spectrum.from_file(syn_input_sfs)
        # TODO: check to make sure its folded
        syn_ns = syn_data.sample_sizes  # Number of samples.
        pts_l = [80, 100, 120]  # what pts do I want to use ? Keep as this for now?
        #pts_l = [40,60,80] # this led to failure
        # Optimize parameters for this model. # TODO: which model?
        # First set parameter bounds for optimization
        # TODO SET BOUNDS AND INITIAL GUESS FOR VAQUITA
        ####### parameters for 3 epoch model #########
        param_names = ("nu", "T")
        upper_bound_demography = [10, 0.5]
        lower_bound_demography = [1e-4, 1e-5]
        initial_guess_demography = [1, 0.01]  # p0
        with open(inferred_demography, 'w') as f:
            f.write('Beginning with demographic inference.\n')
            max_likelihood = -1e25
            # TODO put range back to 100
            for i in range(100):
                # Start at initial guess
                p0 = initial_guess_demography
                # Randomly perturb parameters before optimization.
                p0 = dadi.Misc.perturb_params(
                    p0, fold=1, upper_bound=upper_bound_demography,
                    lower_bound=lower_bound_demography)
                # Make the extrapolating version of demographic model function.
                # TODO pick demographic model that works for vaquita
                func_ex = dadi.Numerics.make_extrap_log_func(self.two_epoch)  ### for vaquita it is three epoch!
                logger.info(
                    'Beginning optimization with guess, {0}.'.format(p0))
                popt_demography = dadi.Inference.optimize_log_lbfgsb(
                    p0=p0, data=syn_data, model_func=func_ex, pts=pts_l,
                    lower_bound=lower_bound_demography, upper_bound=upper_bound_demography,
                    verbose=len(p0), maxiter=100)
                logger.info(
                    'Finished optimization with guess, ' + str(p0) + '.')
                logger.info('Best fit parameters: {0}.'.format(popt_demography))
                # Calculate the best-fit model allele-frequency spectrum.
                # Note, this spectrum needs to be multiplied by "theta".
                non_scaled_spectrum = func_ex(popt_demography, syn_ns, pts_l)
                # Likelihood of the data given the model AFS.
                multinomial_ll_non_scaled_spectrum = \
                    dadi.Inference.ll_multinom(
                        model=non_scaled_spectrum, data=syn_data)
                logger.info(
                    'Maximum log composite likelihood: {0}.'.format(
                        multinomial_ll_non_scaled_spectrum))
                theta = dadi.Inference.optimal_sfs_scaling(
                    non_scaled_spectrum, syn_data)
                logger.info(
                    'Optimal value of theta: {0}.'.format(theta))
                if multinomial_ll_non_scaled_spectrum > max_likelihood:
                    best_params = popt_demography
                    best_non_scaled_spectrum = non_scaled_spectrum
                    max_likelihood = multinomial_ll_non_scaled_spectrum
                    theta_syn = theta
            best_scaled_spectrum = theta_syn * best_non_scaled_spectrum
            theta_nonsyn = theta_syn * NS_S_ScalingFactor  # factor e.g. 2.31 to scale up thetasyn to theta nonsyn (accounting for diff seq length and diff gc content)
            # can get poisson LL from *scaled* spectrum (don't do it from the theta=1 spectrum)
            poisson_ll = dadi.Inference.ll(
                model=best_scaled_spectrum, data=syn_data)

            # estimate the best possible poisson LL (data compared to itself) as a benchmark (Not used for anything but as a reference)
            syn_data_to_data_poisson_ll=dadi.Inference.ll(
                model=syn_data, data=syn_data) # not used for anything but as a reference to see how good it could get
            f.write('Best fit parameters: {0}.\n'.format(best_params))
            f.write(
                'Maximum multinomial log composite likelihood: {0}.\n'.format(
                    max_likelihood))
            f.write(
                'Maximum poisson log composite likelihood: {0}.\n'.format(
                    poisson_ll))
            f.write('Non-scaled best-fit model spectrum: {0}.\n'.format(
                best_non_scaled_spectrum.fold())) # AB added fold() here
            f.write('Optimal value of theta_syn: {0}.\n'.format(theta_syn))
            f.write('Optimal value of theta_nonsyn: {0}.\n'.format(
                theta_nonsyn))
            f.write('Scaled best-fit model spectrum: {0}.\n'.format(
                best_scaled_spectrum.fold())) # AB added fold() here
            f.write('Best possible poisson log composite likelihood (empirical syn sfs compared to itself): {0}.\n'.format(syn_data_to_data_poisson_ll))

        logger.info('Finished demographic inference.')
        logger.info('Beginning DFE inference.')
        nonsyn_data = dadi.Spectrum.from_file(nonsyn_input_sfs)
        nonsyn_ns = nonsyn_data.sample_sizes

        demog_params = best_params

        Lsyn = (1 / (
                    1 + NS_S_ScalingFactor)) * Lcds  # Length of synonymous sites. is 1/(1+NS:S Scaling factor) * Lcds (total cds sites)
        # for reference, get best possible poisson LL of nonsyn sfs to itself for reference (not used as part of inference)
        nonsyn_data_to_data_poisson_ll = dadi.Inference.ll(
            model=nonsyn_data, data=nonsyn_data)
        with open(inferred_DFE,"w") as f:
            header="exon_mutationRate\tmodel\trunNum\tPoissonLL_ABCalc\tPoissonLL_FromNormSFS_ABCalc\tPoissonLL\tPoissonLL_Data_to_Data\tNancestral\tParameters_ScaleScaledby2Na\tParameters_NotScaledBy2Na\tExpectedSFS\tExpectedSFS_integrate_norm\tlower_bound\tupper_bound\tmax_s\n"
            f.write(header)
            for mutationRate in mutationRate_list:
                u = float(mutationRate)
                u_exon = u * exonicMutationRateScalingFactor  # I'm using 1; was set at 1.25 by Eduardo and Jon to increase mutation rate in exons.
                print("mutation rate = " + str(u) + "exon scaling factor is: "+ str(exonicMutationRateScalingFactor) + " exon mutation rate is: " + str(u_exon))

                Na = theta_syn / (4 * u_exon * Lsyn)

                max_s = 0.5 # on 20201125 tried it with max_s = 1; didn't make diff to fit or DFE parameters. Bernard suggests multiplying the singleton bin of the most deleterious spectrum (spectra.spectra[0] by nonsyntheta and by the prob mass of your DFE with s>0.5 to see how many sites you may lose from the expected SFS due to setting max_s at 0.5. It was <1 in my case.
                max_gam = max_s * 2 * Na

                #pts_l = [1200, 1400, 1600] # AB trying with smaller grid
                #pts_l = [40, 60, 80] this led to failure
                pts_l = [600,800,1200] # AB trying with smaller grid
                ############ THIS MUST USE THE APPOPRIATE DEMOGRAPHIC MODEL ! (trying two epoch for vaq.)
                spectra = Selection.spectra(demog_params, nonsyn_ns,
                                            self.two_epoch_sel,
                                            pts_l=pts_l, int_bounds=(1e-5, max_gam),
                                            Npts=300, echo=True, mp=True)
                # try to output most deleterious expected spectrum: (Bernard tip)
                spectra_0=spectra.spectra[0] # most deleterious is this one
                #spectra_minus1=spectra.spectra[-1]
                with open(most_deleterious_spectrum,"w") as d:
                    d.write("spectra.spectra[0]: "+str(spectra_0)+"\n")
                    #d.write("spectra.spectra[-1]: "+str(spectra_minus1)+"\n")

                BETAinit = max_gam / 3
                initial_guess_gamma = [0.09, BETAinit]
                upper_beta = 10 * max_gam
                lower_bound_gamma = [1e-3, 1e-2] # AB trying to the scale lower bound non-zero. 1e-2 instead instead of 0 (note it's in units of 2Na so that's quite small)
                upper_bound_gamma = [1, upper_beta]

                gamma_max_likelihoods = []
                gamma_guesses = dict()
                for i in range(25):
                    p0 = initial_guess_gamma
                    p0 = dadi.Misc.perturb_params(p0, lower_bound=lower_bound_gamma,
                                                  upper_bound=upper_bound_gamma,fold=2) # AB adding fold=2
                    logger.info('Beginning optimization with guess, {0}.'.format(p0))
                    popt_gamma = numpy.copy(Selection.optimize_log(p0, nonsyn_data,
                                                             spectra.integrate,
                                                             self.gamma_dist,
                                                             theta_nonsyn,
                                                             lower_bound=lower_bound_gamma,
                                                             upper_bound=upper_bound_gamma,
                                                             verbose=len(p0),
                                                             maxiter=100)) # trying maxiter 100
                    logger.info('Finished optimization, results are {0}.'.format(popt_gamma))

                    gamma_max_likelihoods.append(popt_gamma[0])
                    #gamma_guesses[popt_gamma[0]] = popt_gamma
                    gamma_guesses[i] = popt_gamma # AB indexing by runnumber not likelihood



                logger.info('Finished DFE inference.')
                # trying to fix issue reported by Jon Mah that best LL is not always on top in sort by adding key=float bc perhaps it was sorting as a str (according to google) https://stackoverflow.com/questions/1513727/python-sort-not-working-as-expected
                # AB: don't need to sort these because I don't care what order they are output
                #gamma_max_likelihoods.sort(key=float)


                logger.info('Integrating expected site-frequency spectrum.')

                logger.info('Outputing results.')

                # f.write('Assuming a gamma-distributed DFE...\n')
                # f.write('Outputting best 25 MLE estimates.\n')s
                # %% gamma
                # note that this sorts things by increasing LL so the "run num" doesn't correspond to the runNums above.
                for i in range(25):
                    #best_popt_gamma = gamma_guesses[gamma_max_likelihoods[i]] # indexing by runNum not Lhood AB
                    best_popt_gamma = gamma_guesses[i]
                    expected_sfs_gamma = spectra.integrate(
                        best_popt_gamma[1], self.gamma_dist, theta_nonsyn).fold() # AB adding "fold"
                    # just to check, am recalculating the LL to make sure it matches what's output by the inference
                    poisson_ll_gamma_ABcalc = dadi.Inference.ll(
                        model=expected_sfs_gamma, data=nonsyn_data)
                    # experiment: try integrate_norm (?)
                    expected_sfs_gamma_normalized = spectra.integrate_norm(
                        best_popt_gamma[1], self.gamma_dist, theta_nonsyn).fold() # AB adding "fold"
                    poisson_ll_gamma_norm_ABcalc = dadi.Inference.ll(
                        model=expected_sfs_gamma_normalized, data=nonsyn_data)
                    # f.write(
                    #    'The population-scaled best-fit parameters: {0}.\n'.format(
                    #        best_popt))
                    # popt[0] is the LL
                    # popt[1] is an array of parameters
                    # don't want to divide shape by anything (so divide by 1), but need scale to be divided by 2Na
                    best_popt_gamma_unscaled = numpy.divide(best_popt_gamma[1], numpy.array([1, 2 * Na]))
                    f.write("\t".join([str(u_exon), "gamma", str(i), str(poisson_ll_gamma_ABcalc), str(poisson_ll_gamma_norm_ABcalc), str(best_popt_gamma[0]), str(nonsyn_data_to_data_poisson_ll), str(Na), str(best_popt_gamma[1]), str(best_popt_gamma_unscaled), str(expected_sfs_gamma),str(expected_sfs_gamma_normalized),str(lower_bound_gamma),str(upper_bound_gamma),str(max_s)]))
                    f.write("\n")
                    # Divide output scale parameter by 2 * N_a
                    # f.write(
                    #    'The non-scaled best-fit parameters: '
                    #    '[{0}, array({1})].\n'.format(
                    #        best_popt[0],
                    #        numpy.divide(best_popt[1], numpy.array([1, 2 * Na]))))
                    # f.write('The expected SFS is: {0}.\n\n'.format(expected_sfs))
                # f.write('Assuming a neutral-gamma-distributed DFE...\n')
                # f.write('Outputting best 25 MLE estimates.\n')

        f.close() # ab added this to flush out final output
        logger.info('Pipeline executed successfully.')


if __name__ == '__main__':
    DemographicAndDFEInference().main()
