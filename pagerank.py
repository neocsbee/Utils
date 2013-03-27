

"""
 This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/license
s/>.

PageRank calculation - optional convergence plot, full list of PR vals, personalized
PR computation 
"""
def print_mat(inputmat):
    (row,col)=inputmat.shape
    for ii in range(row):
        for jj in range(col):
            print inputmat[ii,jj],
        print
def print_lotup(ilist):
    length = len(ilist)
    for ii in range(length):
        print ilist[ii][0],ilist[ii][1]
####################################################################
def calc_pr(bigccmat, numnodes, pprloc=-99):
    """
    function calc_pr calculates PageRank based on the input transition matrix
    """
    #convert to transition matrix 
    rowsum = bigccmat.sum(axis=1)
    for ii in range(numnodes):
        if rowsum[ii,0] !=0:
            bigccmat[ii,:] = bigccmat[ii,:]/rowsum[ii,0]
        else:
            #case with no outgoing links
            bigccmat[ii,ii] = 1.0

    #convert sparse matrix format
    sp_transmat_first = scisp.csr_matrix(bigccmat) 
    oldprvec = np.matrix(np.ones((numnodes,1)))/float(numnodes)
    convergevec = [1000] #some large value
    if pprloc > 0:
        onevec = np.matrix(np.zeros((numnodes, 1)))
        onevec[pprloc,0] = 0.15
    else:
        onevec = (0.15/float(numnodes))*np.matrix(np.ones((numnodes,1)))
    ii = 0
    while convergevec[-1] > 1e-5:
        newprvec = 0.85*(sp_transmat_first.T* oldprvec)
        newprvec = newprvec + onevec
        newnorm = np.linalg.norm(newprvec, 1)
        convergevec.append(sum(abs(newprvec-oldprvec))[0,0])
        oldprvec = newprvec
        ii = ii + 1
    print 'Norm of PR vector:', newnorm
    print 'Number of iterations for convergence:', ii 
    convergevec.remove(1000)
    return (newprvec, convergevec)

#################################################################################
import cPickle as cp
import gzip #used this based on suggestion from python cookbook
import networkx as nx
import numpy as np
import scipy.sparse as scisp
import argparse

parser = argparse.ArgumentParser(description="Calculation of PageRank and optional \
        display of PR vector, output of convergence information, personalized page rank\
        for a given user")
parser.add_argument("-i", "--inputfile", default="fullgraph.dat", \
        help="input gzipped, pickled full graph file")
group = parser.add_mutually_exclusive_group()
group.add_argument("-c", "--convergence", help="plot the convergence of PR \
        computation graph", action="store_true")
group.add_argument("-p", "--prvec", help="output sorted PR vec",\
        action = "store_true")
group.add_argument("-r", "--ppr", \
        help="personalized pagerank for a particular node (all letters in small case)", default="itsvalence")
args = parser.parse_args()

oyginput = gzip.open(args.inputfile, "rb")
tweetG = cp.load(oyginput)
oyginput.close()
UtweetG = tweetG.to_undirected()
ccnodelist = nx.connected_components(UtweetG)

if args.convergence:
    """
    call calc_pr twice.
    """
    bigcc = tweetG.subgraph(ccnodelist[0]).copy()
    bigcc_numedges = bigcc.size()
    sum_edges = 0
    ii = 0
    cc_coll = ccnodelist[0][:]
    for ii in range(1, len(ccnodelist)):
        cc_coll.extend(ccnodelist[ii])
        x = tweetG.subgraph(ccnodelist[ii]).copy()
        sum_edges = sum_edges + x.size()
        if sum_edges >= bigcc_numedges:
            break
    fingraph = tweetG.subgraph(cc_coll).copy()
    #first call with half the number of edges
    bigccmat = nx.to_numpy_matrix(bigcc, ccnodelist[0], order='C', weight=None) 
    newprvec1, convergevec1 = calc_pr(bigccmat, len(ccnodelist[0]))
    #sec call with all edges
    ccmat = nx.to_numpy_matrix(fingraph, cc_coll, order='C', weight=None)
    assert len(cc_coll) > len(ccnodelist[0])
    newprvec2, convergevec2 = calc_pr(ccmat, len(cc_coll))
    import matplotlib.pyplot as plt
    import math
    plt.plot(map(math.log10, convergevec1), 'g+-', \
            label="{0} links".format(bigcc_numedges))
    plt.plot(map(math.log10, convergevec2), 'ro-', \
            label="{0} links".format(fingraph.size()))
    plt.xlabel('Number of iterations')
    plt.ylabel('log10(Total absolute difference from previous iteration)')
    plt.title('Convergence of PR')
    plt.legend(loc="upper right")
    plt.savefig("prconvergence.pdf")
    plt.show()
elif args.prvec: 
    #biggest connected component
    cc_coll = [y for x in ccnodelist[0:1] for y in x]
    bigcc = tweetG.subgraph(cc_coll)
    bigccmat = nx.to_numpy_matrix(bigcc,nodelist=cc_coll, order = 'C', weight=None)
    newprvec, convergevec = calc_pr(bigccmat, len(cc_coll))
    myzip = zip(cc_coll, list(np.asarray(newprvec).flatten()))
    from operator import itemgetter
    myzip = sorted(myzip, key=itemgetter(1),reverse=True)
    print 'Sorted PR Vals:'
    print_lotup(myzip)
else:
    """
    personalized pr computation
    """
    pprnode = args.ppr
    if pprnode not in UtweetG:
        import sys
        sys.exit('{0} not in graph'.format(pprnode))
    pprccnodelist = nx.node_connected_component(UtweetG, pprnode)
    pprcc = tweetG.subgraph(pprccnodelist).copy()
    pprccmat = nx.to_numpy_matrix(pprcc, nodelist=pprccnodelist, order='C', weight=None)
    newprvec, covergevec = calc_pr(pprccmat, len(pprccnodelist),\
            pprccnodelist.index(pprnode))
    myzip = zip(pprccnodelist, list(np.asarray(newprvec).flatten()))
    from operator import itemgetter
    myzip = sorted(myzip, key=itemgetter(1),reverse=True)
    percentiles = [100*(ii-0.5)/len(myzip) for ii in range(len(myzip),0,-1)]
    lol = [list(x) for x in zip(*myzip)]
    labels = lol[0][:11]
    vals = lol[1][:11]
    myzip = zip(myzip, percentiles)
    print 'Sorted PR Vals w.r.t {0}: Top 10 values'.format(pprnode)
    print_lotup(myzip[:11])
    import matplotlib.pyplot as plt
    import math
    plt.plot(map(math.log10, vals), percentiles[:11], '*-')
    plt.xlabel('PR vals log10 scale')
    plt.ylabel('Percentiles')
    plt.title("{0}'s view PR percentile - Top 10 values".format(pprnode))
#    locs, xlabs = plt.xticks()
#    plt.xticks (locs, labels[::-1], rotation=90)
    plt.show()
    plt.savefig("ppr.pdf")
