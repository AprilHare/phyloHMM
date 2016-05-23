# test for the likelihood function
# the total likelihood should be 1.

from importScience import *
from phyloTree import phyloTree

runningTotal = 0.
fullSpectrum = []

treeSchema = ( ( ('A', 'B'), ('C','D') ), ('E','F','G') )
a = phyloTree( treeSchema)
a.edges *= 1

for first in ['A', 'C', 'G', 'T']:
    for second in ['A', 'C', 'G', 'T']:
        for third in ['A', 'C', 'G', 'T']:
            print third
            for fourth in ['A', 'C', 'G', 'T']:
                for fifth in ['A', 'C', 'G', 'T']:
                    for sixth in ['A', 'C', 'G', 'T']:
                        for seventh in ['A', 'C', 'G', 'T']:
                            currString = first + second + third + fourth + fifth + sixth + seventh
                            out =  a.findLikelihood(currString)[0][0]
                            runningTotal += out
                            fullSpectrum.append( out )
                            if out > 0.00039:
                                print currString


print runningTotal
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.hist(fullSpectrum)
ax.set_yscale( 'log')
plt.show()

