# test for the likelihood function
# the total likelihood should be 1.

from phyloTree import phyloTree

runningTotal = 0.

treeSchema = ( ( ('A', 'B'), ('C','D') ), ('E','F','G') )
a = phyloTree( treeSchema)

for first in ['A', 'C', 'G', 'T']:
    for second in ['A', 'C', 'G', 'T']:
        for third in ['A', 'C', 'G', 'T']:
            print third
            for fourth in ['A', 'C', 'G', 'T']:
                for fifth in ['A', 'C', 'G', 'T']:
                    for sixth in ['A', 'C', 'G', 'T']:
                        for seventh in ['A', 'C', 'G', 'T']:
                            currString = first + second + third + fourth + fifth + sixth + seventh
                            runningTotal += a.findLikelihood(currString)[0][0]

print runningTotal


