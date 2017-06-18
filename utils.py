import numpy as np
import networkx as nx
import time as T
import os
from collections import defaultdict
from zss import Node

#############
## NODE
#############

def get_depth(tree):
	# input: Node (tree)
	# output: depth of tree

    maxdepth = [0]
    
    def plunge(node, depth):
        if len(node.children) > 0:
            if depth+1 > maxdepth[0]:
                maxdepth[0] = depth+1
                for child in node.children:
                    plunge(child, depth+1)
            else:
                return
        else:
            return
    
    plunge(tree, 0)
    
    return maxdepth[0]


#############
## CALL DATA
#############

def load_calls(path, sort=False, nDays=-1, month=-1, year=-1, withDays=False):
    files = os.listdir(path)
    files = [f for f in files if f[-3:] == 'dat']
    files.sort(key=lambda x: int(x.split("_")[4][:-8]))
    
    if month > 0 and year > 0:
        files = [f for f in files if 
                   (int(f.split('_')[4][4:6]) == month and \
                    int(f.split('_')[4][0:4]) == year)]
        print 'Reading from month', month, 'year', year
    
    if nDays > 0:
        files = files[:np.min([len(files), nDays])]
        
    print 'LOADING CALLS (%d files)...' % (len(files))
    days = []
    for i, fname in enumerate(files):
        if os.path.getsize(path + fname) > 50:
            with open(path + fname, 'rb') as f:
                #f.read(51)    # data starts at 51 bytes in
                day = np.fromfile(f, 
                                  dtype=np.dtype('i4'),
                                  count=-1,
                                  sep='')
                day = day.reshape((len(day)/6,6))
                if day.size != 0:
                    days.append(day)
        else:
            print ' Skipping (empty)', fname

    if withDays:
        print 'ADDING DAY COLUMN...'
        for j in xrange(len(days)):
            days[j] = np.c_[days[j], 
                            np.ones(len(days[j]), dtype='i4') * j]
            
    print 'STACKING DAILY CALL LOGS...'
    calls = np.vstack(days)

    if sort:
        print 'SORTING CHRONOLOGICALLY...'
        order = np.argsort(calls[:,4])
        calls = calls[order,:]
        
    print ' Num calls (lines):', len(calls)

    return calls

###########
## CASCADE - AUGMENT/ENRICH
###########

def enrich_calls(calls, allC, pers_class, minsize=2, verbose=False):
	# input: all calls, all trees, min class size
	# output: augmented calls with column indicating whether call 
	# part of pers cascade

	# construct list of cascade calls  (10% of full calls)
	casccalls = np.zeros((700000,3))
	idx = [0]    

	def add_cascade_call(subtree):
		for c in subtree.children:
			casccalls[idx[0]] = [subtree.label, c.label, c.tau]
			idx[0] += 1
			add_cascade_call(c)
		return

	START = T.time()
	for k, root in enumerate(pers_class):
		if k % 10000 == 0 and verbose:
			print k, len(pers_class), T.time() - START
		# loop on this guys' persistent cascades
		cs = [c for c in pers_class[root] if len(c) >= minsize]
		cs = np.unique([item for sublist in cs for item in sublist])
		for cx in cs:
			tree = allC[root][cx]
			add_cascade_call(tree)         
	casccalls = casccalls[:idx[0]]

	# ensure `casccalls' and `calls' are ordered by call time, then caller
	print 'Sorting ...'
	O = np.lexsort((casccalls[:,0], casccalls[:,2]))
	casccalls = casccalls[O]
	O = np.lexsort((calls[:,0],calls[:,4]))
	calls = calls[O,:]

	# augment calls with a cascade column (index 7)
	acalls = np.c_[calls, np.zeros(len(calls), dtype='i4')]
	print 'Building...'
	pk = 0
	START = T.time()
	for k, call in enumerate(acalls):
		if k % 1000000 == 0 and verbose:
			print k, len(acalls), T.time() - START
		# HACKY (checking long range because of duplicates...)
		for p, casc in enumerate(casccalls[pk:pk+20]):
			if call[4] == casc[2] and call[0] == casc[0] and call[2] == casc[1]:
				acalls[k,7] = 1
				pk += p+1
				break

	return acalls
