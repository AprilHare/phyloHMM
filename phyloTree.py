#implements a phylogenetic tree class

import numpy as np

from scipy import linalg

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

        # keep the tree schema around so that we have something to recurse on
        self.treeSchema = treeSchema
        # the edges tell us the branch lengths of the tree and the structure of the tree
        self.edges = parseSchema( treeSchema)

        # default nucleotide rate matrices. These are just a stand in for the priors.
        self.background = 0.25 * np.ones( (4,1) )   # background rates
        self.transition = 0.25 * np.ones( (4,4) )   # transition rates
        for i in range(4):
            self.transition[i,i] = -0.75


    def findLikelihood(self, alignmentCol):
        """ Felsensteins pruning algorithm to compute likelihood of a column of an alignment given the tree structure"""

        numNodes = self.edges.shape[0]

        def findMessage(nodeIndex, content):
            # the message is f(parent) = P( index | parent, branch length)

            if nodeIndex == 0:
                # check for the root node
                probMat = self.background                       # just the background rates
            else:
                branchLength = max( self.edges[nodeIndex, :nodeIndex ] ) 
                probMat = self.ratesToProbs( branchLength )             # implicitely uses the transition matrices

            XtoContent = np.dot( content, probMat)      # probability of the content given the parental node.

            return XtoContent

        memo = {}
        alignIndex = 0              # we might have to reverse the order of the alignment
        for nodeIndex in range(numNodes-1, -1, -1):
            # iterate backwards through the matrix of edges
            if max( self.edges[nodeIndex, nodeIndex:]) < 1E-4:
                # do the leaf thing
                leafValue = self.letterToVector( alignmentCol[ alignIndex ] )
                alignIndex += 1

                memo[nodeIndex] = findMessage( nodeIndex, leafValue)

            else:
                # do the node thing

                # product the messages
                runningProd = 1
                for child in range(nodeIndex, numNodes):
                    if self.edges[ nodeIndex, child ] != 0:
                        runningProd *= memo[child]

                # marginalize out the current node.
                fullMessage = findMessage(nodeIndex, runningProd )
                memo[nodeIndex] = fullMessage

        return  memo




    # utilities
    def ratesToProbs(self, time):
        """ returns nucleotide probabilities from lengths"""
        return linalg.expm(time*self.transition)


    def letterToVector(self, letter):
        """ returns the vector corresponding to a element of the alphabet"""
        if letter.lower() == 'a':
            return np.array( [1, 0, 0, 0])
        elif letter.lower() == 'c':
            return np.array( [0, 1, 0, 0])
        elif letter.lower() == 'g':
            return np.array( [0, 0, 1, 0])
        elif letter.lower() == 't':
            return np.array( [0, 0, 0, 1])


    def checkParams(self):
        """ checks to make sure that the paramter values are legal"""
        if abs( sum(self.background) - 1.0 ) > 1e-3:
            raise Exception( "Background rates must sum to 1")
        for row in range( self.transition.shape[0]):
            if abs( sum( self.transition[row, :] ) ) > 1e-3:
                raise Exception( "Transistion rates must sum to 0. Row " + str(row) )


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












