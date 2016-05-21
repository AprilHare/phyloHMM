# implements a phylogenetic tree emissions for hmmlearn based hmms.

import numpy as ny
from hmmlearn import hmm

class phyloHMM(object):
    """phyloHMM is a subclass of hmmlearn hmms for which the emissions are phylogenetic models"""
    def __init__(self, arg):
        super(phyloHMM, self).__init__()
        self.arg = arg
        



