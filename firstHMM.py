# a first run through implementing a phylo HMM
import numpy as np

import Bio
from Bio import HMM
from Bio import Phylo

from hmmlearn import hmm


# set up an example HMM

testModel = hmm.GaussianHMM(n_components=3, covariance_type="full")
testModel.startprob_ = np.array([0.6, 0.3, 0.1])
testModel.transmat_ = np.array([[0.7, 0.2, 0.1],
                                 [0.3, 0.5, 0.2],
                                 [0.3, 0.3, 0.4]])

testModel.means_ = np.array([[0.0, 0.0], [3.0, -3.0], [5.0, 10.0]])
testModel.covars_ = np.tile(np.identity(2), (3, 1, 1))
X, Z = model.sample(100)  # X is the sequence of emissions, Z is the sequence of states


# etc. 



