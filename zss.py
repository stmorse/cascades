# Zhang-Shasha implementation (https://github.com/timtadh/zhang-shasha)
# This file consolidates the entire library, and consists of:
# classes: Node, AnnotatedTree
# methods: distance, simple_distance
# 
# Soft dependencies: numpy and editdist
# Special character comments removed, see original GitHub for full
#
# Modified Node class to include `tau` attribute (time received call)


##########################

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Author: Tim Henderson
#Email: tim.tadh@gmail.com

##########################

# ORIGINAL LiCENSE:

# Zhang-Shasha Tree Edit Distance Implementation is licensed under a BSD style
# license 

# Copyright (c) 2012, Tim Henderson and Stephen Johnson
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright notice,
#       this list of conditions and the following disclaimer in the documentation
#       and/or other materials provided with the distribution.
#     * Neither the name of this software nor the names of its contributors may
#       be used to endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

###########################

from __future__ import absolute_import

import collections


class Node(object):
    """
    A simple node object that can be used to construct trees to be used with
    :py:func:`zss.distance`.
    Example: ::
        Node("f")
            .addkid(Node("a")
                .addkid(Node("h"))
                .addkid(Node("c")
                    .addkid(Node("l"))))
            .addkid(Node("e"))
    """

    def __init__(self, label, children=None, tau=-1):
        self.label = label
        self.children = children or list()
        self.tau = tau

    @staticmethod
    def get_children(node):
        """
        Default value of ``get_children`` argument of :py:func:`zss.distance`.
        :returns: ``self.children``.
        """
        return node.children

    @staticmethod
    def get_label(node):
        """
        Default value of ``get_label`` argument of :py:func:`zss.distance`.
        :returns: ``self.label``.
        """
        return node.label

    @staticmethod
    def get_tau(node):
        return node.tau

    def addkid(self, node, before=False):
        """
        Add the given node as a child of this node.
        """
        if before:  self.children.insert(0, node)
        else:   self.children.append(node)
        return self

    def get(self, label):
        """:returns: Child with the given label."""
        if self.label == label: return self
        for c in self.children:
            if label in c: return c.get(label)

    def iter(self):
        """Iterate over this node and its children in a preorder traversal."""
        queue = collections.deque()
        queue.append(self)
        while len(queue) > 0:
            n = queue.popleft()
            for c in n.children: queue.append(c)
            yield n

    def prettyprint(self):
        s = ["0: %s (%d)" % (self.label, len(self.children))]
        def addchildprint(node, level):
            for cn in node.children:
                sp = ' ' * (level*2)
                s[0] += '\n' + sp + '%d: %s (%d) (%d)' % \
                    (level, cn.label, len(cn.children), cn.tau)
                if len(cn.children) > 0: 
                    addchildprint(cn, level+1)
        addchildprint(self, 1)
        return s[0]


    def __contains__(self, b):
        if (isinstance(b, str) or isinstance(b, int)):
            if self.label == b: 
                return 1
            else:
                return sum(b in c for c in self.children)
        else: 
            if self.label == b.label: 
                return 1
            else:
                return sum(b in c for c in self.children)
        raise TypeError("Object %s is not of type str, int or Node" % repr(b))

    def __eq__(self, b):
        if b is None: return False
        if not isinstance(b, Node):
            raise TypeError("Must compare against type Node")
        return self.label == b.label

    def __ne__(self, b):
        return not self.__eq__(b)

    def __repr__(self):
        return super(Node, self).__repr__()[:-1] + " %s>" % self.label

    def __str__(self):
        s = "%d:%s" % (len(self.children), self.label)
        s = '\n'.join([s]+[str(c) for c in self.children])
        return s


################
################
################

# from __future__ import absolute_import
from six.moves import range

# import collections

try:
    import numpy as np
    zeros = np.zeros
except ImportError:
    def py_zeros(dim, pytype):
        assert len(dim) == 2
        return [[pytype() for y in range(dim[1])]
                for x in range(dim[0])]
    zeros = py_zeros

try:
    from editdist import distance as strdist
except ImportError:
    def strdist(a, b):
        if a == b:
            return 0
        else:
            return 1

# from zss.simple_tree import Node
# from simple_tree import Node


class AnnotatedTree(object):

    def __init__(self, root, get_children):
        self.get_children = get_children

        self.root = root
        self.nodes = list()  # a pre-order enumeration of the nodes in the tree
        self.ids = list()    # a matching list of ids
        self.lmds = list()   # left most descendents
        self.keyroots = None
            # k and k' are nodes specified in the pre-order enumeration.
            # keyroots = {k | there exists no k'>k such that lmd(k) == lmd(k')}
            # see paper for more on keyroots

        stack = list()
        pstack = list()
        stack.append((root, collections.deque()))
        j = 0
        while len(stack) > 0:
            n, anc = stack.pop()
            nid = j
            for c in self.get_children(n):
                a = collections.deque(anc)
                a.appendleft(nid)
                stack.append((c, a))
            pstack.append(((n, nid), anc))
            j += 1
        lmds = dict()
        keyroots = dict()
        i = 0
        while len(pstack) > 0:
            (n, nid), anc = pstack.pop()
            #print list(anc)
            self.nodes.append(n)
            self.ids.append(nid)
            #print n.label, [a.label for a in anc]
            if not self.get_children(n):
                lmd = i
                for a in anc:
                    if a not in lmds: lmds[a] = i
                    else: break
            else:
                try: lmd = lmds[nid]
                except:
                    import pdb
                    pdb.set_trace()
            self.lmds.append(lmd)
            keyroots[lmd] = i
            i += 1
        self.keyroots = sorted(keyroots.values())


def simple_distance(A, B, get_children=Node.get_children,
        get_label=Node.get_label, label_dist=strdist):
    """Computes the exact tree edit distance between trees A and B.
    Use this function if both of these things are true:
    * The cost to insert a node is equivalent to ``label_dist('', new_label)``
    * The cost to remove a node is equivalent to ``label_dist(new_label, '')``
    Otherwise, use :py:func:`zss.distance` instead.
    :param A: The root of a tree.
    :param B: The root of a tree.
    :param get_children:
        A function ``get_children(node) == [node children]``.  Defaults to
        :py:func:`zss.Node.get_children`.
    :param get_label:
        A function ``get_label(node) == 'node label'``.All labels are assumed
        to be strings at this time. Defaults to :py:func:`zss.Node.get_label`.
    :param label_distance:
        A function
        ``label_distance((get_label(node1), get_label(node2)) >= 0``.
        This function should take the output of ``get_label(node)`` and return
        an integer greater or equal to 0 representing how many edits to
        transform the label of ``node1`` into the label of ``node2``. By
        default, this is string edit distance (if available). 0 indicates that
        the labels are the same. A number N represent it takes N changes to
        transform one label into the other.
    :return: An integer distance [0, inf+)
    """
    return distance(
        A, B, get_children,
        insert_cost=lambda node: label_dist('', get_label(node)),
        remove_cost=lambda node: label_dist(get_label(node), ''),
        update_cost=lambda a, b: label_dist(get_label(a), get_label(b)),
    )


def distance(A, B, get_children, insert_cost, remove_cost, update_cost):
    '''Computes the exact tree edit distance between trees A and B with a
    richer API than :py:func:`zss.simple_distance`.
    Use this function if either of these things are true:
    * The cost to insert a node is **not** equivalent to the cost of changing
      an empty node to have the new node's label
    * The cost to remove a node is **not** equivalent to the cost of changing
      it to a node with an empty label
    Otherwise, use :py:func:`zss.simple_distance`.
    :param A: The root of a tree.
    :param B: The root of a tree.
    :param get_children:
        A function ``get_children(node) == [node children]``.  Defaults to
        :py:func:`zss.Node.get_children`.
    :param insert_cost:
        A function ``insert_cost(node) == cost to insert node >= 0``.
    :param remove_cost:
        A function ``remove_cost(node) == cost to remove node >= 0``.
    :param update_cost:
        A function ``update_cost(a, b) == cost to change a into b >= 0``.
    :return: An integer distance [0, inf+)
    '''
    A, B = AnnotatedTree(A, get_children), AnnotatedTree(B, get_children)
    treedists = zeros((len(A.nodes), len(B.nodes)), int)

    def treedist(i, j):
        Al = A.lmds
        Bl = B.lmds
        An = A.nodes
        Bn = B.nodes

        m = i - Al[i] + 2
        n = j - Bl[j] + 2
        fd = zeros((m,n), int)

        ioff = Al[i] - 1
        joff = Bl[j] - 1

        for x in range(1, m): 
            fd[x][0] = fd[x-1][0] + remove_cost(An[x+ioff])
        for y in range(1, n): 
            fd[0][y] = fd[0][y-1] + insert_cost(Bn[y+joff])

        for x in range(1, m): ## the plus one is for the xrange impl
            for y in range(1, n):
                # only need to check if x is an ancestor of i
                # and y is an ancestor of j
                if Al[i] == Al[x+ioff] and Bl[j] == Bl[y+joff]:
                    fd[x][y] = min(
                        fd[x-1][y] + remove_cost(An[x+ioff]),
                        fd[x][y-1] + insert_cost(Bn[y+joff]),
                        fd[x-1][y-1] + update_cost(An[x+ioff], Bn[y+joff]),
                    )
                    treedists[x+ioff][y+joff] = fd[x][y]
                else:
                    p = Al[x+ioff]-1-ioff
                    q = Bl[y+joff]-1-joff
                    #print (p, q), (len(fd), len(fd[0]))
                    fd[x][y] = min(
                        fd[x-1][y] + remove_cost(An[x+ioff]),
                        fd[x][y-1] + insert_cost(Bn[y+joff]),
                        fd[p][q] + treedists[x+ioff][y+joff]
                    )

    for i in A.keyroots:
        for j in B.keyroots:
            treedist(i,j)

    return treedists[-1][-1]

