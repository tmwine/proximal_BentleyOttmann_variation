'''

inherited tree structure from base class event_tree

the status tree needs to handle insertions based on node data (vs value),
data being a segment object, and the status tree needs to handle internal
node sorting based on node value

node data amounts to segment objects, containing segment point pairs,
[[x1,y1],[x2,y2]], with the segment object inequality dunder (__gt__)
overloaded to detect if a just-inserted segment is above or below 
another

running insert will emplace the node based on the comparator for
segment objects, and will decide on an internally-assigned value for that node;
checks ensure node values don't get too close, and node value rebalancing will
occur if node values are too close (having no topological effect on the
tree structure)

MIT license; T. M. Wine

'''

import copy
import math
import matplotlib.pyplot as plt
from event_tree import AVLTree, TOL_ACC
from event_tree import AVLNode as Node


class Segment():

    # class for value objects in status tree nodes
    def __init__(self,seg_vls):
        # expects a list of 2 points, amounting to segment endpoints; each point
        # is of form [x,y]; always stores points internally so that p1 has x
        # value <= p2's x value
        if seg_vls[0][0] > seg_vls[1][0]:
            self.p1 = seg_vls[1]
            self.p2 = seg_vls[0]
        else:
            self.p1 = seg_vls[0]
            self.p2 = seg_vls[1]
        self.vec = [self.p2[0]-self.p1[0],self.p2[1]-self.p1[1]]
        tmp_nrm = math.sqrt(self.vec[0]**2+self.vec[1]**2)
        self.nrm_vec = [self.vec[0]/tmp_nrm,self.vec[1]/tmp_nrm]

    def __eq__(self,other):
        # this is unlikely; 2 segments overlap w/ measure > 0
        if self.p1[0] > other.p1[0]:
            sn = -1.0
            v0 = other.nrm_vec
            v1 = [self.p1[0]-other.p2[0],self.p1[1]-other.p2[1]]
        else:
            sn = 1.0
            v0 = self.nrm_vec
            v1 = [other.p1[0]-self.p2[0],other.p1[1]-self.p2[1]]
        tmp = -v0[1]*v1[0] + v0[0]*v1[1]
        if abs(tmp) < TOL_ACC:    # segments intersect in at least one pt
            tmp_2 = -other.vec[0] * self.vec[1] + other.vec[1] * self.vec[0]
            if tmp_2 == 0:  # segment angles are equal
                if sn==1.0 and other.p1[0] < self.p2[0]:
                    return True
                elif sn==-1.0 and self.p1[0] < other.p2[0]:
                    return True
        else:
            return False

    def __gt__(self,other):
        # expects a segment object for comparison to self's object;
        # check for truth of self (caller) is > other (received object);
        # ">" means above in the y direction;
        # NB: this expects at least some x-direction overlap between the 2
        # segments
        # for the segment with the rightmost left endpoint (rle),
        # this compares the vector of the other segment (vec) with the vector
        # from right end of vec to rle; so this looks only at where rle is
        # relative to vec (above or below)

        if self.p1[0] > other.p1[0]:
            sn = -1.0
            v0 = other.nrm_vec
            v1 = [self.p1[0]-other.p2[0],self.p1[1]-other.p2[1]]
        else:
            sn = 1.0
            v0 = self.nrm_vec
            v1 = [other.p1[0]-self.p2[0],other.p1[1]-self.p2[1]]
        tmp = -v0[1]*v1[0] + v0[0]*v1[1]    # projection of v1 onto unit
        # perpendicular to v0
        if sn*tmp > TOL_ACC:
            return False
        elif abs(tmp) < TOL_ACC:   # T intersection, or common endpoints; this
            # then
            # deals with "ties," and these are sub-sorted according to CCW /
            # CW orientation of the two segments as vectors; this does not
            # really return meaningful comparison w/ 2 segments chained
            # end-to-end (right endpt of 1st = left endpt of 2nd)
            tmp_2 = -other.vec[0]*self.vec[1] + other.vec[1]*self.vec[0]
            return tmp_2 < 0
        else:
            return True

    def get_points(self):
        return [self.p1,self.p2]


class NodeData:

    def __init__(self,seg_idx,seg_obj):
        # container for "data" portion of nodes;
        # expects the segment index # and an object of type Segment
        self.sg_ix = seg_idx
        self.seg = seg_obj

    def __repr__(self):
        return "%s: %s" % (self.sg_ix,[self.seg.p1,self.seg.p2])


class AVLStatusTree(AVLTree):

    VALUE_SPREAD = 512.0    # default for initializing node values
    VALUE_MIN = 1e-8    # minimal value separation between nodes (before
    # triggering warning eg)

    def __init__(self):
        self.node_index = {}    # dictionary for storing node value to
        # segment index lookups; key,value is segment index #, node value
        AVLTree.__init__(self)

    def insert(self,data):
        # this is a "special" insert method--this relies on an object with a
        # valid comparator (< > ==) in the "data" attribute of the nodes;
        # expects an object of type NodeData(seg_index,Segment object);
        # places the node in the tree according to the contents of the
        # NodeData object;
        # returns True if successfully inserted node with that data;
        # does nothing and returns False if a node with "equal" data
        # already exists in the tree;
        # note, the object comparator in use here is "smart" about
        # ~degeneracies--if several segments have the same start point for
        # example, the comparator sorts by slope/angle
        in_ok = False  # if node w/ value "value" does not already exist
        if self.root is None:
            nod_val = 0.0
            self.setRoot(nod_val, data)
            self.node_index[data.sg_ix] = nod_val
            in_ok = True
        else:
            in_ok = self._insertNode(self.root,data)
        if in_ok:
            self.nodes_count += 1
        return in_ok

    def _insertNode(self, currentNode, data):
        # recursive helper function to insert a value into AVLTree
        node_to_rebalance = None
        reb_vls = False # flag for whether node value rebalancing is necessary
        if currentNode.data.seg > data.seg:
            if currentNode.left_child is not None:
                return self._insertNode(currentNode.left_child, data)
            else:
                currentNode.left_child = Node(None, data)
                currentNode.left_child.parent = currentNode
                node_up = currentNode
                node_dn = self.next_dn(currentNode.left_child)
                if node_dn is None:
                    nod_val = node_up.value-type(self).VALUE_SPREAD
                else:
                    nod_val = (node_up.value+node_dn.value)/2.0
                    if abs(node_up.value-node_dn.value) < type(self).VALUE_MIN:
                        reb_vls = True
                        print("node values getting close together in "
                              "AVLStatusTree; rebalancing") # DEBUG
                currentNode.left_child.value = nod_val
                self.node_index[data.sg_ix] = nod_val   # update segment
                # index to node value dictionary
                if currentNode.height == 0:
                    # ie we've changed the height of the tree along this
                    # branch, so rebalance
                    self._recomputeHeights(currentNode) # children node
                    # heights will be OK; have to change heights of parents
                    # on up
                    node = currentNode
                    while node:
                        if node.balanceFactor() in [-2, 2]:
                            node_to_rebalance = node
                            break  # we need the one that is furthest from the
                            # root
                        node = node.parent
                    if node_to_rebalance:
                        self._rebalance(node_to_rebalance)
                if reb_vls:
                    self.rebalance_values()
                return True
        elif currentNode.data.seg < data.seg:
            if currentNode.right_child is not None:
                return self._insertNode(currentNode.right_child, data)
            else:
                currentNode.right_child = Node(None, data)
                currentNode.right_child.parent = currentNode
                node_up = self.next_up(currentNode.right_child)
                node_dn = currentNode
                if node_up is None:
                    nod_val = node_dn.value + type(self).VALUE_SPREAD
                else:
                    nod_val = (node_up.value + node_dn.value) / 2.0
                    if abs(node_up.value - node_dn.value) < type(
                            self).VALUE_MIN:
                        reb_vls = True
                        print("node values getting close together in "
                              "AVLStatusTree; rebalancing")  # DEBUG
                currentNode.right_child.value = nod_val
                self.node_index[data.sg_ix] = nod_val  # update segment
                # index to node value dictionary
                if currentNode.height == 0:
                    # height of tree has changed locally
                    self._recomputeHeights(currentNode) # children node
                    # heights will be OK; have to change heights of parents
                    # on up
                    node = currentNode
                    while node:
                        if node.balanceFactor() in [-2, 2]:
                            node_to_rebalance = node
                            break  # we need the one that is furthest from the
                            # root
                        node = node.parent
                    if node_to_rebalance:
                        self._rebalance(node_to_rebalance)
                if reb_vls:
                    self.rebalance_values()
                return True
        else:
            return False    # value already in tree

    def delete_segment(self,sg_ix):
        # removes segment w/ index sg_ix from tree, and from node_index
        # dictionary
        nd_vl = self.get_node_value(sg_ix)
        res = self.delete_value(nd_vl)    # call base class delete-by-value
        try:
            del(self.node_index[sg_ix])
        except KeyError:
            print("problem with deleting segment from tree; may not exist?")
            raise
        return res

    def next_up(self,node):
        # find node w/ value next highest in the tree, via topological search
        if self.root is None:
            return None
        if node.right_child is not None:
            node = node.right_child
            while node.left_child is not None:
                node = node.left_child
            return node
        else:
            while node.parent is not None:
                if node.parent.left_child == node:
                    return node.parent
                node = node.parent
        return None

    def next_dn(self,node):
        # find node w/ value next lowest in the tree, via topological search
        if self.root is None:
            return None
        if node.left_child is not None:
            node = node.left_child
            while node.right_child is not None:
                node = node.right_child
            return node
        else:
            while node.parent is not None:
                if node.parent.right_child == node:
                    return node.parent
                node = node.parent
        return None

    def swap_segs(self,ix_1,ix_2):
        # swaps the data between segments w/ indices ix_1 and ix_2
        val_1 = self.node_index[ix_1]
        val_2 = self.node_index[ix_2]
        node_1 = self.search(val_1)
        node_2 = self.search(val_2)
        tmp_dat = node_1.data
        node_1.data = node_2.data
        node_2.data = tmp_dat
        self.node_index[ix_1] = val_2
        self.node_index[ix_2] = val_1

    def rebalance_values(self):
        # rebalances node values to VALUE_SPREAD distance apart (does not
        # change node ordering)
        tmp_dct = {}
        nm_nd = self.node_count()
        cr_vl = -type(self).VALUE_SPREAD*math.floor(nm_nd/2.0)
        cur_nod = self.get_min(val_dat=False)   # start at minimum node (
        # topologically)
        if cur_nod is None:
            return  # nothing to do
        while True:
            cur_nod.value = cr_vl
            tmp_dct[cur_nod.data.sg_ix] = cr_vl
            cr_vl += type(self).VALUE_SPREAD
            cur_nod = self.next_up(cur_nod)
            if cur_nod is None:
                break
        self.node_index = tmp_dct

    def get_node_value(self,sg_ix):
        # given a segment index, return its current node's value from segment
        # lookup dictionary
        return self.node_index[sg_ix]

# non-essential; used for debugging
def plt_lin(p1,p2,clr='blue',mrk=False):
    # p1,p2 of form [x,y]; expects an open plot window
    x1,y1 = p1
    x2,y2 = p2
    if not mrk:
        plt.plot([x1,x2],[y1,y2],color=clr)
    else:
        plt.plot([x1,x2],[y1,y2],
                 color=clr,marker='o')

