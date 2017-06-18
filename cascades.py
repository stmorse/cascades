##########################

# Implementation of Persistent Cascades algorithm 
#  described in [](https://stmorse.github.io/docs/BigD348.pdf)
# For usage see README
# For license see LICENSE
# Author: Steven Morse
# Email: steventmorse@gmail.com
# License: MIT License (see LICENSE in top folder)

##########################

import os
import numpy as np
import multiprocessing as mp
import networkx as nx
import time as T
import copy

from utils import load_calls
from zss import Node, simple_distance

##################
# METHODS USED IN MP, CAN'T BE CLASS INSTANCE OBJECTS

def jaccard(im1, im2):
    # |A\cap B|/|A\cup B| = |A\cap B|/(|A|+|B|-|A\cap B|)
    im1 = np.asarray(im1)
    im2 = np.asarray(im2)
    intersection = len(np.intersect1d(im1, im2))
    return intersection / float(len(im1) + len(im2) - intersection)


def build_sim_mx(trees, nodelists=None, method='zss', similarity=True, speedup=-1):
    '''Returns a distance matrix of `trees` using metric `method`
    
    trees : list of trees (Node)
    method : `zss` is normalized TED, `jaccard` is on node lists, `both` gives both 
                   and can take the speedup parameter
    similarity : if False, returns the dissimilarity (1-sim)
    speedup : [0,1]. Will set TED=0 when jaccard < speedup. (Zero = no speedup)
    
    Returns sim_mx (or if `both`, returns TED then jaccard)
    '''
    
    numt = len(trees)
    lists = [[c.label for c in t.iter()] for t in trees] if nodelists==None else nodelists
    sizes = [len(li) for li in lists]
    sim_mx = np.zeros((numt, numt))
    sim_mx2 = np.zeros((numt, numt)) if method=='both' else None
    for i in xrange(numt):
        for j in xrange(i, numt):
            if i==j: 
                sim_mx[i,j] = 1
                if method=='both':
                    sim_mx2[i,j] = 1
                continue
                
            if method=='both':
                jacc = jaccard(lists[i], lists[j])
                nted = 0
                if jacc >= speedup:
                    ted = simple_distance(trees[i], trees[j])
                    nted = 1 - (2*ted / float(sizes[i]+sizes[j]+ted))
                sim_mx[i,j], sim_mx[j,i] = nted, nted
                sim_mx2[i,j], sim_mx2[j,i] = jacc, jacc
            else:
                val = 0
                if method=='zss':
                    ted = simple_distance(trees[i], trees[j])
                    val = 1 - (2*ted / float(sizes[i]+sizes[j]+ted))
                elif method=='jaccard':
                    val = jaccard(lists[i], lists[j])
                else:
                    print 'Method unrecognized.'
                    return
                sim_mx[i,j], sim_mx[j,i] = val, val
    
    if method=='both':
        if not similarity:
            sim_mx = 1 - sim_mx
            sim_mx2 = 1 - sim_mx2
        return sim_mx, sim_mx2
    else:
        return sim_mx if similarity else (1 - sim_mx)


# mp method used by build_distmx_tgraph
def do_build(root, alltrees, mintree, minsize):
    # get indices of periods with a tree
    treex = [i for i, t in enumerate(alltrees) if type(t)==Node]
    if len(treex) < mintree:
        return None

    # build lists of nodes in each tree
    nlists = []
    bigtreex = []
    for x in treex:
        li = [c.label for c in alltrees[x].iter()]
        if len(li) < minsize:
            continue 
        nlists.append(li)
        bigtreex.append(x)
    if len(nlists) < mintree:   # if num big enough trees is too small
        return None

    nted_mx, jacc_mx = build_sim_mx([alltrees[v] for v in bigtreex],
                                    nodelists=nlists,   # send lists of nodes
                                    method='both',      # retrieve jaccard and nted
                                    speedup=0.3)        # only compute TED if jaccard >= 0.3

    # build tree graph with pairwise comparisons
    g = nx.Graph()
    g.add_nodes_from(bigtreex)
    for i, x in enumerate(bigtreex):
        for j, y in enumerate(bigtreex):
            if j <= i: continue
            if nted_mx[i,j] > 0 or jacc_mx[i,j] > 0:
                g.add_edge(x, y, nted=nted_mx[i,j], jacc=jacc_mx[i,j])

    return (root, g)


###########################


class Cascades:
    '''Loads call data, extracts cascades, and finds persistence classes.

    Dependencies: numpy, multiprocessing, networkx, zss, utils

    Callable methods:
    build : build call matrix, extracts all cascades, creates similarity matrices
    build_persistence_classes : takes tgraphs and creates persistence classes

    Internal methods:
    load_call_data : loads call data
    extract_cascade : given root(s), extracts all cascades
    build_sim_mx : given set of cascades, builds similarity matrix using Jaccard/TED
    build_all_cascades : extracts all cascades using `extract_cascade`
    do_build : multiprocessing worker method
    build_distmx_tgraph : takes trees and creates tgraphs for each root node
    '''

    def __init__(self, calls=None, UTC=False, path='', city='', nMonths=1, startx=0, moyr=[]):
        '''Initialize Cascades object.  
        `calls` is an optional numpy array of the raw call data.
        `moyr` is list of tuples giving month/year combos. (Default is for the Portugal data.)'''
        
        self.moyr = [(4,2006), (5,2006), (6,2006), (7,2006), (8,2006), (11,2006), (12,2006),
                     (1,2007), (2,2007), (3,2007), (4,2007), (5,2007), (6,2007)]
        if len(moyr) > 0:
            self.moyr = moyr

        if calls is None:
            self.calls = self.load_call_data(path+city+'/', nMonths=nMonths, startx=startx)
        else:
            self.calls = calls
        self.UTC = UTC
        

    def build(self, nsample=3000, mincalls=50, period=1, sample=[], exclude=[], mintree=10, minsize=3, 
              daybegins=4*60*60, multi=True, batchsize=1000, numproc=7, verbose=True):
        '''Load call data, extract all call data, and build the similarity matrices for each 
        possible root.  Store in vars `calls`, `allC`, and `tgraphs`.
        '''

        if len(sample) == 0:
            self.allC  = self.build_all_cascades(nsample=nsample, mincalls=mincalls, period=period, 
                                                 exclude=exclude, daybegins=daybegins, verbose=verbose)
        else:
            self.allC = self.build_all_cascades(sample=sample, period=period, verbose=verbose)
        
        self.tgraphs = self.build_distmx_tgraph(mintree=mintree, minsize=minsize, 
                                                multi=multi, batchsize=batchsize, numproc=numproc)

    def build_persistence_classes(self, ell=0.8, dayrange=[], verbose=True):
        '''Create (and return) list of lists for each root with maximal persistence classes.
        Note: dayrange should be in terms of the periods, not the actual days.
        Default is to do the entire dataset.'''

        if len(dayrange)==0:
            dayrange = range(np.amax(calls[:,6]))

        # by root, list of lists
        self.pers_class_nted = {}
        self.pers_class_jacc = {}

        START = T.time()
        numg = len(self.tgraphs)
        for k, root in enumerate(self.tgraphs):
            if k % 10000 == 0 and verbose:
                print k, numg, (T.time() - START)

            # create a temporary version of tgraph[root] with low-wt edges removed
            # find nted persistence classes
            gt = self.tgraphs[root].copy()
            for n in gt.nodes():
                if n not in dayrange:
                    gt.remove_node(n)
            for e in gt.edges():
                if gt.edge[e[0]][e[1]]['nted'] < ell:
                    gt.remove_edge(e[0], e[1])
            self.pers_class_nted[root] = list(nx.find_cliques(gt))

            # find jaccard persistence classes
            gj = self.tgraphs[root].copy()
            for n in gj.nodes():
                if n not in dayrange:
                    gj.remove_node(n)
            for e in gj.edges():
                if gj.edge[e[0]][e[1]]['jacc'] < ell:
                    gj.remove_edge(e[0], e[1])
            self.pers_class_jacc[root] = list(nx.find_cliques(gj))
        
        return self.pers_class_nted, self.pers_class_jacc

    
    ##################


    def load_call_data(self, path, nMonths=1, startx=8):
        '''Load raw call data (possibly multiple months) into single object.
        `calls` must be structured caller, tower, callee, tower, time stamp, 
        duration, period (day).  
        Time stamp is structured with the last 5 digits giving the second past
        midnight, and the leading 2-3 digits giving the 24-period past Jan 2006.
        '''

        calls = load_calls(path, sort=True, 
                           month=self.moyr[startx][0], year=self.moyr[startx][1], withDays=True)
        
        for m in range(1, nMonths):  # if nMonths==1, will skip
            calls2 = load_calls(path, sort=True, month=self.moyr[startx+m][0], 
                                year=self.moyr[startx+m][1], withDays=True)
            calls2[:,6] += np.amax(calls[:,6]) + 1
            calls = np.vstack((calls, calls2))
       
        return calls

    #####

    def extract_cascade(self, calls, roots, maxdepth=-1, maxhr=24, verbose=False):
        '''Extract time-respecting, non-repeating cascades for single or list of nodes 
        Uses a recursive internal method 'subtree()'
        
        roots    : desired root (user id)
        maxdepth : cutoff number of tiers (-1 for unlimited, default)
        maxhr    : cutoff time, starting from first base root call (in hrs)
        verbose  : debug output
        '''
        
        # sort calls[] by outgoing user
        O = np.argsort( calls[:,0] )
        calls = calls[O,:]
        indx = np.cumsum( np.bincount( calls[:,0] ) )
        indx = np.append([0],indx)
        # now range(indx[u], indx[u+1]) gives the indices in calls[] of user u's calls
        
        # Initialize the trees (indexed by root)
        Cs = {}
        
        # Initialize dict of nodes and their first call in the tree
        # structure: {node: [parentnode, time]}
        firsts = {}
        
        # time cutoff, set as first call by root + maxsec
        maxt = 0
        
        def subtree(rootnode, curdepth):        
            root = rootnode.label
            tau  = rootnode.tau
            
            # stop if we have reached maxdepth
            if curdepth == maxdepth:
                return
            
            # if root didn't make any calls, return
            try:
                if indx[root] == indx[root+1]:
                    return
            except IndexError:
                # print 'Index out of bounds for root', root
                # this will occur if root didn't make any outgoing calls
                # and has an id larger than any of the outgoing users
                return

            # get slice of all this root's outgoing calls, 
            # sort by time, and remove all before tau or after maxt
            rcalls = calls[indx[root]:indx[root+1],:]
            rcalls = rcalls[np.argsort(rcalls[:,4]),:]
            rcalls = rcalls[(rcalls[:,4] > tau) & (rcalls[:,4] < maxt)]
            
            # if no valid calls, return
            if len(rcalls) == 0:
                return

            # loop thru root's calls and start new subtree
            for i, c in enumerate(rcalls[:,2]):
                # if child already reached, check if we remove old or skip
                if c in firsts:
                    if (rcalls[i,4] < firsts[c][1]):
                        # existing spot was a later event, remove it
                        oc = firsts[c][0].get(int(c))
                        firsts[c][0].children.remove(oc)
                    else:
                        # existing spot was an earlier event, skip this child
                        continue

                # update firsts
                firsts[c] = (rootnode, rcalls[i,4])
                
                # add child
                child = Node(int(c), tau=rcalls[i,4])
                rootnode.addkid(child)
                subtree(child, curdepth+1)
        ## END def subtree
        
        def order_children(rootnode):
            if len(rootnode.children) == 0:
                return
            # order children of the rootnode
            rootnode.children = sorted(Node.get_children(rootnode), key=lambda c: c.label)
            for cnode in Node.get_children(rootnode):
                order_children(cnode)
        ## END def order_children
                
        # do build...
        maxsec = maxhr * 60 * 60
        for i, r in enumerate(roots):
            try:
                if indx[r] == indx[r+1]:
                    continue
            except IndexError:
                # if verbose: print 'Index out of bounds for root', r
                continue
            r = int(r)
            Cs[r] = Node(r, tau=-1)
            t1 = np.amin(calls[indx[r]:indx[r+1],4])
            maxt = t1 + maxsec
            firsts = {r: (None, t1)}
            subtree(Cs[r], 0)
            
        # order children by label so tree edit distance is unbiased
        # if verbose: print 'Ordering by label...',
        for r in Cs:
            order_children(Cs[r])
        
        if verbose: print '(%d)' % (len(Cs)),
        
        return Cs 

    def build_all_cascades(self, nsample=50000, mincalls=10, period=1, sample=[], 
        exclude=[], daybegins=4*60*60, verbose=True):
        '''Extract all cascades for a random sample of users.
        
        calls -- raw call data, form: caller, tow1, callee, tow2, time, duration, day
        nsample -- size of random sample (-1 does all users with mincalls)
        mincalls -- min calls in the dataset
        period -- num days for a cascade.  Start of day is 4am by default.
        '''

        calls = self.calls

        #daybegins = 4 * 60 * 60  # 4am is ``start of the day''

        print '\nTotal %d total unique users.' % (len(np.unique(calls[:,0])))

        if len(sample) == 0:
            nCalls = np.bincount(calls[:,0])
            users = np.where(nCalls >= mincalls)[0]
            users = np.setdiff1d(users, exclude)
            
            print 'Found %d non-excluded users with enough outgoing calls.' % (len(users))
            
            if nsample > 0:
                print 'Sampling %d...' % (nsample)
                users = np.random.choice(users, nsample, replace=False)
        else:
            users = sample
            print 'Sampling %d users (specified)' % (len(users))

        numdays = np.amax(calls[:,6])
        numperiods = numdays / period
        print 'Num days: %d, Num periods (%d hrs): %d\n' % (numdays, period*24, numperiods)

        # all cascades, keyed by root. format: {root: [Node, 0, ...], ...}
        allC = {}

        print 'Extracting trees... (verbose format: period (trees))'
        START = T.time()
        for idx, p in enumerate(xrange(0, numdays-(period-1)-1, period)):
            if self.UTC:
                tcalls = calls[calls[:,6]==p]
                if len(tcalls)==0:
                    if verbose: print 'Skip %d (no data)' % (idx),
                    continue
                
                if verbose: print '%d' % (idx),

                C = self.extract_cascade(tcalls,
                                         users,
                                         maxdepth=-1,
                                         maxhr=(period*24),
                                         verbose=True)
            else:
                day1 = calls[calls[:,6]==p]
                if len(day1)==0:
                    if verbose: print 'Skip %d (no data)' % (idx),
                    continue
                
                # the nasty string manipulation is to deal with the time format ...
                startsec = int(''.join([str(day1[0,4])[:-5], '{0:0>5}'.format(str(int(daybegins)))]))
                endsec   = startsec + 86400
                
                if verbose: print '%d' % (idx),

                C = self.extract_cascade(calls[(calls[:,4] > startsec) & (calls[:,4] < endsec)],
                                         users,
                                         maxdepth=-1,
                                         maxhr=(period*24),
                                         verbose=verbose)

            for tree in C:
                try:
                    allC[tree][idx] = C[tree]
                except KeyError:
                    allC[tree] = [0 for _ in xrange(numperiods)]
                    allC[tree][idx] = C[tree]
            
            if verbose: print ' .. ',       
            if idx % 5 == 0 and idx > 0 and verbose:
                print 'Time %1.3f' % (T.time() - START)

        print 'Complete.\n'

        self.users = users
        
        return allC


    #########

    def build_distmx_tgraph(self, mintree=3, minsize=3, 
                        multi=True, batchsize=1000, numproc=4,
                        verbose=True):
        ''' Build distance matrices, store in a graph'''

        # all NTED and JACC similarity matrices, indexed by root
        # all_nted = {}
        # all_jacc = {}

        # all tree graphs, indexed by root
        all_tgraphs = {}
        allC = self.allC

        numt = len(self.allC)
        START = T.time()
        numbatch = (numt / batchsize) + 1   # a little sloppy?
        output = []
        for k in range(numbatch):
            if verbose:
                print '[%d to %d). (Batch %d / %d). Total: %d.  Time: %1.3f' % \
                    (k*batchsize, (k+1)*batchsize, k, numbatch, numt, (T.time() - START))

            if multi:
                results = []            
                pool = mp.Pool(processes=numproc)

                results = [pool.apply_async(do_build, 
                    args=(root, self.allC[root], mintree, minsize, )) for \
                    i, root in enumerate(self.allC) \
                        if i >= (k*batchsize) and i < ((k+1)*batchsize)]

                pool.close()
                pool.join()

                temp = [p.get() for p in results]
            else:
                temp = [do_build(root, self.allC[root], mintree, minsize) for \
                        i, root in enumerate(self.allC) \
                            if i >= (k*batchsize) and i < ((k+1)*batchsize)]
            output.extend([t for t in temp if t != None])

        if verbose: print 'Writing to tgraphs...'
        for g in output:
            all_tgraphs[g[0]] = g[1]

        if verbose:
            print 'Complete. Total roots:', len(all_tgraphs)
            print 'Total time:', T.time() - START
            print ''
            
        # return all_tgraphs
        return all_tgraphs

    ################


    

