#implements a phylogenetic tree class

import numpy as np


class phyloTree(object):
    """phyloTree is a phylogenetic tree"""
    def __init__(self, treeSchema):
        # phylo-tree is a matrix of connection lengths with a collection of other useful functions.
        # Everything disconnected is set to zero, and only updated by multiplication
        super(phyloTree, self).__init__()

        def parseSchema(treeSchema):
            """parses the tree depth first. The schema is assumed to be a tuple of tuples"""

            def count(current):
                if type(current) is tuple:
                    # we are still recursing.
                    currCount = 0
                    for entry in current:
                        currCount += count(entry)
                    return currCount + 1

                else:
                    #this is taken to be the condition for the end of the tree.
                    return 1

            numNodes = count(treeSchema)
            out = np.zeros( (numNodes, numNodes ) )


            def parser(current, bigIndex, out):
                #recurser
                if type(current) is tuple:
                    # we are still recursing.
                    myIndex = bigIndex
                    running = myIndex + 1
                    for entry in current:
                        out[myIndex, running] = 1
                        out[running, myIndex] = 1

                        biggest = parser( entry, running, out)
                        running = biggest + 1

                    return running - 1

                else:
                    #this is taken to be the condition for the end of the tree.
                    myIndex = bigIndex
                    return myIndex

            parser( treeSchema, 0, out)

            return out

        # the edges tell us the branch lengths of the tree and the structure of the tree
        self.edges = parseSchema( treeSchema)

        # nucleotide rate matrices
        self.background = 0.25 * np.ones( (4,1) )   # background rates
        self.transition = 0.25 * np.ones((4,4))     # transition rates




    def checkParams(self):
        """ checks to make sure that the paramter values are legal"""
        if abs( sum(self.background) - 1.0 ) > 1e-3:
            raise Exception( "Background rates must sum to 1")
        for row in range( self.transition.shape[0]):
            if abs( sum( self.transition[row, :] ) - 1.0 ) > 1e-3:
                raise Exception( "Transistion rates must sum to 1. Row " + str(row) )


    def getConnections(self, index):
        (out,) = np.where( self.edges[index, :] != 0 )
        return out


    def setLength(self, location, value):
        # only update if the entry is non-zero
        i1, i2 = location

        prev1 = self.edges[i1,i2];
        prev2 = self.edges[i2,i1];

        if prev1 == 0 or prev2 == 0:
            raise Exception( 'Incompatible with tree schema')
        else:
            self.edges[i1,i2] = value
            self.edges[i2,i1] = value












