# implements a phylogenetic tree emissions for hmmlearn based hmms.

import numpy as ny

from hmmlearn import hmm
from phyloTree import phyloTree

class phyloHMM(hmm._BaseHMM):
    """phyloHMM is a subclass of hmmlearn hmms for which the emissions are phylogenetic models"""
    def __init__(self, treeSchema, n_components=1, covariance_type='diag',
                 min_covar=1e-3,
                 startprob_prior=1.0, transmat_prior=1.0,
                 means_prior=0, means_weight=0,
                 covars_prior=1e-2, covars_weight=1,
                 algorithm="viterbi", random_state=None,
                 n_iter=10, tol=1e-2, verbose=False,
                 params="stmc", init_params="stmc"):

        # this function initializes the hmmm object itself.

        hmm._BaseHMM.__init__(self, n_components,
                          startprob_prior=startprob_prior,
                          transmat_prior=transmat_prior,
                          algorithm=algorithm,
                          random_state=random_state,
                          n_iter=n_iter, tol=tol, verbose=verbose,
                          params=params, init_params=init_params)


        # initialize the the phyloHMM models
        self.phyloTrees_ = []
        for index in range(self.n_components):
            self.phyloTrees_.append( phyloTree(treeSchema))



    # functions that we must overwrite from the _BaseHMM class
    def _init(self, X, length):
        # this function initalizes the model parameters for HMM prior to doing the fitting
        pass

    def _check(self):
        # validatese model parameters prior to fitting
        pass

    def _generate_sample_from_state(self, state, random_state=None):
        # Generates a random sample from a given component.
        pass

    def _compute_log_likelihood(self, X):
        # Computes per-component log probability under the model.
        pass

    def _initialize_sufficient_statistics():
        #Initializes sufficient statistics required for M-step.
        pass

    def _accumulate_sufficient_statistics(stats, X, framelogprob, posteriors, fwdlattice, bwdlattice):
        #Updates sufficient statistics from a given sample.
        pass

    def _do_mstep(stats):
        # Performs the M-step of EM algorithm.
        pass













