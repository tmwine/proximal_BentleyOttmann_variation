'''

AVL tree module, based on https://github.com/cool-pot/pytrees,
under the MIT license


AVL Tree.
Balanced Binary Search Tree. Guaranteed for balance.
Convention:
- "key" and "value" are almost the same in this implementation. use term "key"
for search and delete a particular node. use term "value" for other cases
API:
- insert(self, value)
- delete_value(self, key)
- search(self, key)
- getDepth(self)
- preOrder(self)
- inOrder(self)
- postOrder(self)
- node_count(self)
- buildFromList(cls, l)
Author: Yi Zhou
Date: May 19, 2018
Reference: https://en.wikipedia.org/wiki/AVL_tree
Reference: https://github.com/pgrafov/python-avl-tree/blob/master/pyavltree.py
'''

import copy
from collections import deque
import random

TOL_ACC = 1e-9 # this is synched with its use in associated
# modules, for easier updating

class AVLNode:
    def __init__(self, value=None, data=None):
        # data type of "value" has to be consistent w/ eg Node's overloaded
        # inequality operator(s)
        self.value = value
        self.data = data
        self.parent = None
        self.left_child = None
        self.right_child = None
        self.height = 0 # each node has its height from the lowest of all
        # children below it internally recorded; eg a single root node would
        # have height of 0, and so would a leaf (childless node)

    def isLeaf(self):
        return (self.height == 0)

    def maxChildrenHeight(self):
        # returns maximum height of one or the other children of the node;
        # if it's one or two childless children, then it returns 0; if no
        # children, it returns -1
        if self.left_child and self.right_child:
            return max(self.left_child.height, self.right_child.height)
        elif self.left_child and not self.right_child:
            return self.left_child.height
        elif not self.left_child and self.right_child:
            return self.right_child.height
        else:
            return -1

    def balanceFactor(self):
        # returns disparity in heights between child nodes;
        # if there is no child, the "height" is by default -1 (and if the
        # child were only a leaf, its height would be 0)
        return (self.left_child.height if self.left_child else -1) - (
            self.right_child.height if self.right_child else -1)

    def __str__(self):
        return ("AVLNode(" + str(self.value) + ", " + str(self.data) + ", "
                "Height: %d )" % self.height)


class AVLTree:
    def __init__(self):
        self.root = None
        self.rebalance_count = 0
        self.nodes_count = 0
        self.range_find = []

    def setRoot(self, value, data):
        # set root value
        self.root = AVLNode(value, data)

    def node_count(self):
        return self.nodes_count

    def getDepth(self):
        # get max depth (height from root to lowest child node); returns -1
        # if tree is empty
        if self.root is not None:
            return self.root.height
        else:
            return -1

    def _findSmallest(self, start_node):
        # returns node object at smallest topographical location in subtree
        # with head at input node start_node
        node = start_node
        while node.left_child:
            node = node.left_child
        return node

    def _findBiggest(self, start_node):
        # returns node object at largest topographical location in subtree
        # with head at input node start_node
        node = start_node
        while node.right_child:
            node = node.right_child
        return node

    def get_min(self,val_dat=True):
        # finds minimum node on tree, using topological search;
        # if val_dat is True, returns [value,data] pair corresponding to the
        # node; otherwise, returns node object itself;
        # note, this has ~mild implementation overlap with _findSmallest
        if self.root is None:
            if val_dat:
                return [None,None]
            else:
                return None
        else:
            current = self.root
            while current.left_child is not None:
                current = current.left_child
            if val_dat:
                return [current.value,current.data]
            else:
                return current

    def get_max(self):
        # returns [value,data] pair corresponding to maximum node in tree
        if self.root is None:
            return [None,None]
        else:
            current = self.root
            while current.right_child is not None:
                current = current.right_child
            return [current.value,current.data]

    def get_range(self,lo,hi):
        # fetches nodes (value,data) in range [lo,hi];
        # range is determined by node.value objects (and whatever built-in
        # __lt__, __gt__ they use; eg for coordinate tuples, .value=[x,y],
        # the ordering is lexicographic;
        # returns an empty list if no nodes fall within range
        if self.root is not None:
            self.range_find = []
            self._get_range(lo,hi,self.root)
            return copy.deepcopy(self.range_find)
        else:
            return []

    def _get_range(self,lo,hi,node):
        # stores (value,data) tuples in instance mutable self.range_find list
        if lo < node.value and node.left_child is not None: # if node is <= low
            # value in range, can ignore left child subtree
            self._get_range(lo,hi,node.left_child)
        if lo <= node.value <= hi:
            self.range_find.append((node.value,node.data))
        if node.value < hi and node.right_child is not None:
            self._get_range(lo, hi, node.right_child)

    def insert(self, value, data=None):
        # creates a node of value value and inserts it in correct point in tree;
        # returns True if made a new node--ie if node w/ value does not
        # already exist, False otherwise
        in_ok = False   # if node w/ value "value" does not already exist
        if self.root is None:
            self.setRoot(value,data)
            in_ok = True
        else:
            in_ok = self._insertNode(self.root, value, data)
        if in_ok:
            self.nodes_count += 1
        return in_ok

    def _insertNode(self, currentNode, value, data):
        # recursive helper function to insert a value into AVLTree
        node_to_rebalance = None
        if currentNode.value > value:
            if currentNode.left_child is not None:
                return self._insertNode(currentNode.left_child, value, data)
            else:
                child_node = AVLNode(value, data)
                currentNode.left_child = child_node
                child_node.parent = currentNode
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
                return True
        elif currentNode.value < value:
            if currentNode.right_child is not None:
                return self._insertNode(currentNode.right_child, value, data)
            else:
                child_node = AVLNode(value, data)
                currentNode.right_child = child_node
                child_node.parent = currentNode
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
                return True
        else:
            return False    # value already in tree

    def _rebalance(self, node_to_rebalance):
        A = node_to_rebalance
        F = A.parent  # allowed to be NULL
        if A.balanceFactor() == -2:
            if A.right_child.balanceFactor() <= 0:
                """Rebalance, case RRC 
                [Original]:                   
                        F                         
                      /  \
                 SubTree  A
                           \                
                            B
                             \
                              C
                [After Rotation]:
                        F                         
                      /  \
                 SubTree  B
                         / \  
                        A   C
                """
                B = A.right_child
                C = B.right_child
                assert (not A is None and not B is None and not C is None)
                A.right_child = B.left_child
                if A.right_child:
                    A.right_child.parent = A
                B.left_child = A
                A.parent = B
                if F is None:
                    self.root = B
                    self.root.parent = None
                else:
                    if F.right_child == A:
                        F.right_child = B
                    else:
                        F.left_child = B
                    B.parent = F
                self._recomputeHeights(A)
                self._recomputeHeights(B.parent)
            else:
                """Rebalance, case RLC 
                [Original]:                   
                        F                         
                      /  \
                 SubTree  A
                           \                
                            B
                           /
                          C
                [After Rotation]:
                        F                         
                      /  \
                 SubTree  C
                         / \  
                        A   B
                """
                B = A.right_child
                C = B.left_child
                assert (not A is None and not B is None and not C is None)
                B.left_child = C.right_child
                if B.left_child:
                    B.left_child.parent = B
                A.right_child = C.left_child
                if A.right_child:
                    A.right_child.parent = A
                C.right_child = B
                B.parent = C
                C.left_child = A
                A.parent = C
                if F is None:
                    self.root = C
                    self.root.parent = None
                else:
                    if F.right_child == A:
                        F.right_child = C
                    else:
                        F.left_child = C
                    C.parent = F
                self._recomputeHeights(A)
                self._recomputeHeights(B)
        else:   # A.balanceFactor() == +2
            if node_to_rebalance.left_child.balanceFactor() >= 0:
                """Rebalance, case LLC 
                [Original]:                   
                        F                         
                      /  \
                     A   SubTree
                    /
                   B
                  /
                 C   
                [After Rotation]:
                        F                         
                       / \  
                      B  SubTree
                     / \  
                    C   A
                """
                B = A.left_child
                C = B.left_child
                assert (not A is None and not B is None and not C is None)
                A.left_child = B.right_child
                if A.left_child:
                    A.left_child.parent = A
                B.right_child = A
                A.parent = B
                if F is None:
                    self.root = B
                    self.root.parent = None
                else:
                    if F.right_child == A:
                        F.right_child = B
                    else:
                        F.left_child = B
                    B.parent = F
                self._recomputeHeights(A)
                self._recomputeHeights(B.parent)
            else:
                """Rebalance, case LRC 
                [Original]:                   
                        F                         
                      /  \
                     A   SubTree
                    /
                   B
                    \
                     C

                [After Rotation]:
                        F                         
                       / \  
                      C  SubTree
                     / \  
                    B   A
                """
                B = A.left_child
                C = B.right_child
                assert (not A is None and not B is None and not C is None)
                A.left_child = C.right_child
                if A.left_child:
                    A.left_child.parent = A
                B.right_child = C.left_child
                if B.right_child:
                    B.right_child.parent = B
                C.left_child = B
                B.parent = C
                C.right_child = A
                A.parent = C
                if F is None:
                    self.root = C
                    self.root.parent = None
                else:
                    if (F.right_child == A):
                        F.right_child = C
                    else:
                        F.left_child = C
                    C.parent = F
                self._recomputeHeights(A)
                self._recomputeHeights(B)
        self.rebalance_count += 1

    def _recomputeHeights(self, start_from_node):
        # this recomputes all node .height attributes upward from
        # start_from_node
        changed = True
        node = start_from_node
        while node and changed:
            old_height = node.height
            node.height = (node.maxChildrenHeight() + 1 if (
                        node.right_child or node.left_child) else 0)
            changed = node.height != old_height
            node = node.parent

    def search(self, key):
        """
        Search a AVLNode satisfies AVLNode.value = key.
        if found return AVLNode, else return None.
        """
        return self._dfsSearch(self.root, key)

    def _dfsSearch(self, currentNode, key):
        """
        Helper function to search a key in AVLTree.
        """
        if currentNode is None:
            return None
        elif currentNode.value == key:
            return currentNode
        elif currentNode.value > key:
            return self._dfsSearch(currentNode.left_child, key)
        else:
            return self._dfsSearch(currentNode.right_child, key)

    def lo_neighbor(self,value):
        if self.root is not None:
            nod_fnd = self._lo_neighbor(value)
            if nod_fnd is None:
                return [None,None]
            else:
                return [nod_fnd.value, nod_fnd.data]
        else:
            return [None,None]

    def _lo_neighbor(self,value):
        # value-based (vs topological-based) search for next-lowest value to
        # input value
        child = True
        node = self.root
        best_node = None
        while child:
            if node.value < value:  # search right subtree for improvements
                if best_node is None:
                    best_node = node
                elif node.value > best_node.value:
                    best_node = node
                if node.right_child is not None:
                    node = node.right_child
                else:
                    child = False
            else:
                if node.left_child is not None:
                    node = node.left_child
                else:
                    child = False
        return best_node

    def hi_neighbor(self,value):
        if self.root is not None:
            nod_fnd = self._hi_neighbor(value)
            if nod_fnd is None:
                return [None,None]
            else:
                return [nod_fnd.value, nod_fnd.data]
        else:
            return [None,None]

    def _hi_neighbor(self,value):
        # value-based (vs topological-based) search for next-highest value to
        # input value
        child = True
        node = self.root
        best_node = None
        while child:
            if node.value > value:  # search left subtree for improvements
                if best_node is None:
                    best_node = node
                elif node.value < best_node.value:
                    best_node = node
                if node.left_child is not None:
                    node = node.left_child
                else:
                    child = False
            else:
                if node.right_child is not None:
                    node = node.right_child
                else:
                    child = False
        return best_node

    def get_data(self,value):
        # variation on search(), pared for just [value,data] (not node object)
        # retrieval
        if self.root is not None:
            nod_fnd = self._dfsSearch(self.root, value)
            if nod_fnd is not None:
                return [nod_fnd.value,nod_fnd.data]
            else:
                return [None,None]
        else:
            return [None,None]

    def traverse(self):
        # generator; traverses tree elements in order, starting at minimum
        # returns (value,data) pairs (not node objects themselves)
        if self.root is not None:
            yield from self._traverse(self.root)
        else:
            return None

    def _traverse(self,on_node):
        # traverse function based on topological (not value-based) traverse
        if on_node is not None:
            if on_node.left_child is not None:
                yield from self._traverse(on_node.left_child)
            yield [on_node.value,on_node.data]
            if on_node.right_child is not None:
                yield from self._traverse(on_node.right_child)

    def delete_value(self, key):
        # deletes a node by value lookup;
        # returns False if no node with value "key" is found, True otherwise

        # first find the node attached to the value, "key"
        node = self.search(key)

        if node is not None:    # a node with value "key" has been found
            self.nodes_count -= 1
            #     There are three cases:
            #     1) The node is a leaf.  Remove it and return.
            #     2) The node is a branch (has only 1 child). Make the pointer
            #        to this node point to the child of this node.
            #     3) The node has two children. Swap items with the successor
            #        of the node (the smallest item in its right subtree) and
            #        delete the successor from the right subtree of the node.
            if node.isLeaf():
                self._removeLeaf(node)
            elif (bool(node.left_child)) ^ (bool(node.right_child)):
                self._removeBranch(node)
            else:
                self._swapWithSuccessorAndRemove(node)
            return True
        else:
            return False

    def _removeLeaf(self, node):
        # this removes a node that has no children
        parent = node.parent
        if parent is not None:
            if parent.left_child == node:
                parent.left_child = None
            else:   # node must be parent's right child
                parent.right_child = None
            self._recomputeHeights(parent)
        else:
            self.root = None
        del node
        # rebalance
        node = parent
        while node is not None:
            if not node.balanceFactor() in [-1, 0, 1]:
                self._rebalance(node)
            node = node.parent

    def _removeBranch(self, node):
        # this removes a node that has only one child
        parent = node.parent
        if parent is not None:
            if parent.left_child == node:
                parent.left_child = (node.right_child if node.right_child
                    is not None else node.left_child)
            else: # node must be parent's right child
                parent.right_child = (node.right_child if node.right_child
                    is not None else node.left_child)
            if node.left_child is not None:
                node.left_child.parent = parent
            else:
                node.right_child.parent = parent
            self._recomputeHeights(parent)
        # rebalance
        node = parent
        while node is not None:
            if not node.balanceFactor() in [-1, 0, 1]:
                self._rebalance(node)
            node = node.parent

    def _swapWithSuccessorAndRemove(self, node):
        # this removes a node with 2 children
        successor = self._findSmallest(node.right_child)
        self._swapNodes(node, successor)
        assert (node.left_child is None)
        if node.height == 0:
            self._removeLeaf(node)
        else:
            self._removeBranch(node)

    def _swapNodes(self, node1, node2):
        assert (node1.height > node2.height)
        parent1 = node1.parent
        leftChild1 = node1.left_child
        rightChild1 = node1.right_child
        parent2 = node2.parent
        assert (not parent2 is None)
        assert (parent2.left_child == node2 or parent2 == node1)
        leftChild2 = node2.left_child
        assert (leftChild2 is None)
        rightChild2 = node2.right_child

        # swap heights
        tmp = node1.height
        node1.height = node2.height
        node2.height = tmp

        if parent1:
            if parent1.left_child == node1:
                parent1.left_child = node2
            else:
                assert (parent1.right_child == node1)
                parent1.right_child = node2
            node2.parent = parent1
        else:
            self.root = node2
            node2.parent = None

        node2.left_child = leftChild1
        leftChild1.parent = node2
        node1.left_child = leftChild2  # None
        node1.right_child = rightChild2
        if rightChild2:
            rightChild2.parent = node1
        if not (parent2 == node1):
            node2.right_child = rightChild1
            rightChild1.parent = node2
            parent2.left_child = node1
            node1.parent = parent2
        else:
            node2.right_child = node1
            node1.parent = node2


    #~~~below are optional functions to call externally~~~

    def inOrder(self):
        res = []

        def _dfs_in_order(node, res):
            if not node:
                return
            _dfs_in_order(node.left_child, res)
            res.append(node.value)
            _dfs_in_order(node.right_child, res)

        _dfs_in_order(self.root, res)
        return res

    def preOrder(self):
        res = []

        def _dfs_pre_order(node, res):
            if not node:
                return
            res.append(node.value)
            _dfs_pre_order(node.left_child, res)
            _dfs_pre_order(node.right_child, res)

        _dfs_pre_order(self.root, res)
        return res

    def postOrder(self):
        res = []

        def _dfs_post_order(node, res):
            if not node:
                return
            _dfs_post_order(node.left_child, res)
            _dfs_post_order(node.right_child, res)
            res.append(node.value)

        _dfs_post_order(self.root, res)
        return res

    @classmethod
    def buildFromList(cls, l, shuffle=True):
        """
        return a AVLTree object from l.
        suffle the list first for better balance.
        """
        if shuffle:
            random.seed()
            random.shuffle(l)
        AVL = AVLTree()
        for item in l:
            AVL.insert(item)
        return AVL

    def visualize(self):
        """
        Naive Visualization.
        Warning: Only for simple test usage.
        """
        if self.root is None:
            print("EMPTY TREE.")
        else:
            print("-----------------Visualize Tree----------------------")
            layer = deque([self.root])
            layer_count = self.getDepth()
            while len(list(filter(lambda x: x is not None, layer))):
                new_layer = deque([])
                vd_list = []
                while len(layer):
                    node = layer.popleft()
                    if node is not None:
                        vd_list.append([node.value,node.data])
                    else:
                        vd_list.append(" ")
                    if node is None:
                        new_layer.append(None)
                        new_layer.append(None)
                    else:
                        new_layer.append(node.left_child)
                        new_layer.append(node.right_child)
                vd_list = [" "] * layer_count + vd_list
                print(*vd_list, sep="  ", end="\n")
                layer = new_layer
                layer_count -= 1
            print("-----------------End Visualization-------------------")

