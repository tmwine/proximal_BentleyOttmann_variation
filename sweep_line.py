'''
variation of Bentley-Ottmann sweep line algorithm, with proximal event detection

usage:
    out_segment_list, event_list = run_sweep_line(in_segment_list),
where,
    in_segment_list is a list of segments, [seg_1,seg_2,...], where
seg_i is of form [[x1,y1],[x2,y2]]
    out_segment_list is a list of the segments after endpoint glomming /
snapping; this list is in the same order as in_segment_list
    event_list has elements of form [[x,y],{(segment_index,L/R/I),...}],
where [x,y] is the coordinate of the event, and the set contains all segments
involved in that event, along with whether the segment hit that point with
its left endpoint ('left'), its right endpoint ('right'), or the segment 
interior ('internal'); event_list records lone segment endpoints too--so those
must be filtered out for recovering true inter-segment intersections

the routine has a plotting feature; matplotlib must be available, and the
DEBUG_PLOT flag must be set to True in the header

MIT license; T. M. Wine

'''

import math
import copy
import event_tree    # for balanced binary search tree class for events
import status_tree    # for status AVL tree class

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib not found; install matplotlib or ammend code for your "
          "plotting modules for full functionality")
    matpltlib_enabled = False
else:
    matpltlib_enabled = True

DEBUG_PLOT = True

TOL_ACC = event_tree.TOL_ACC # tolerance for eg point equality


class EventTree(event_tree.AVLTree):

    # just use superclass constructor (which may not need arguments)

    def proc_segs(self):
        # ~patch AVLTree with a specialized method that returns all segments
        # in a segment list; this is useful for just after "uploading" all
        # initial segments into the EventStorage data structure, when glomming
        # may have occurred, to get a segment list with revised final
        # endpoints
        if self.root is None:
            return []
        seg_dct = {}
        found_nodes = [self.root]
        while len(found_nodes) > 0:
            new_nodes = []
            for nd in found_nodes:
                for tup in nd.data:
                    if tup[1]=='internal':  # only want segment endpoints (L/R)
                        continue
                    if tup[0] not in seg_dct:
                        seg_dct[tup[0]] = [[],[]]
                    elif len(seg_dct[tup[0]][tup[1]=='right']) > 0:
                        print("problem with segment fetch from EventTree")
                        breakpoint()
                        raise ValueError
                    seg_dct[tup[0]][tup[1]=='right'] = nd.value
                if nd.left_child is not None:
                    new_nodes.append(nd.left_child)
                if nd.right_child is not None:
                    new_nodes.append(nd.right_child)
            found_nodes = new_nodes
        return seg_dct

    def get_nearby_range(self,lo,hi):
        # based on base class get_range, which worked by alphabetical order
        # on [x,y] point pairs, this returns node (value,data) tuples with
        # proximal point values in R2 (taxicab distance);
        # expects [lo,hi] where lo=[x1,y1] and hi=[x2,y2] and
        # the two points are opposite corners of a bounding box in R2;
        # returns an empty list if no nodes fall within range
        if self.root is not None:
            self.range_find = []
            self._get_nearby_range(lo,hi,self.root)
            return copy.deepcopy(self.range_find)
        else:
            return []

    def _get_nearby_range(self,lo,hi,node):
        # stores (value,data) tuples in instance mutable self.range_find list
        if lo < node.value and node.left_child is not None: # if node is <= low
            # value in range, can ignore left child subtree
            self._get_nearby_range(lo,hi,node.left_child)
        # only record x,y pairs that fit in the bounding box defined by lo,hi
        if all([lo[ii]<=node.value[ii]<=hi[ii] for ii in range(2)]):
            self.range_find.append((node.value,node.data))
        if node.value < hi and node.right_child is not None:
            self._get_nearby_range(lo, hi, node.right_child)


class EventStorage:

    # class for storing events

    def __init__(self):
        # Q_events is an AVL tree with 2-lists as values ([x,y] coordinates);
        # node data is a list of (segment index,event status) tuples,
        # where event status is 'left', 'right', 'internal' (endpoints or
        # intersection)
        self.Q_tree = EventTree()

    def add_point(self,x,y,seg_dat):
        # updates event tree with point x,y and associated segment
        # index tuples list--(index #, 'left'/'right'/'interior) w/ index #
        # from seg_lst;
        # this intelligently updates set of segments associated with the
        # point (unioning segs w/ any existing segments if point (or close
        # prox) already exists in the tree);
        # find tree entries within x proximity range;
        # nvl_lst is list w/ tuples ([x,y],{seg ix set})
        nvl_lst = self.Q_tree.get_nearby_range([x-TOL_ACC,y-TOL_ACC],[x+TOL_ACC,
                                                            y+TOL_ACC])
        if len(nvl_lst)==0:
            # no near-points found; add this new point to tree
            self.Q_tree.insert([x,y],set(seg_dat))
        else:   # so have a list of tuples from the tree of type ([x,y],
            # {(ix#,L/R/I),...}) with all points x,y "close"
            glm_set = set(seg_dat)
            for nd in nvl_lst:
                # delete all existing in-range nodes from tree
                self.Q_tree.delete_value(nd[0])
                # glom all segment events under a single set
                glm_set.update(nd[1])
            # update Q_tree with a single "glommed" event, using [x,
            # y] of first element in nvl_lst
            self.Q_tree.insert(nvl_lst[0][0],glm_set)

    def glom_to_seg(self,seg):
        # glom segment endpoints near true vertical segments to those segments
        if seg[0][1]>seg[1][1]:
            lo_pt = seg[1]
            hi_pt = seg[0]
        else:
            lo_pt = seg[0]
            hi_pt = seg[1]
        x = seg[0][0]   # segment should be exactly vertical
        ymn,ymx = lo_pt[1],hi_pt[1]
        nvl_lst = self.Q_tree.get_nearby_range([x - TOL_ACC, ymn],
                                               [x + TOL_ACC,
                                                ymx])
        if len(nvl_lst)==0:
            return
        for nd_fo in nvl_lst:   # nd_fo will be a nodal (value,data) tuple
            ptx,pty = nd_fo[0]
            if ymn+TOL_ACC<=pty<=ymx-TOL_ACC:
                if abs(ptx-x)<TOL_ACC:
                    # shift point ptx,pty to on-segment and update tree
                    self.Q_tree.delete_value(nd_fo[0])
                    self.Q_tree.insert([x,pty],nd_fo[1])

    def get_sgs(self):
        # gets current [[x1,y1],[x2,y2]] lists associated with segments
        return self.Q_tree.proc_segs()

    def traverse_events(self):
        # generator to run through all events in the tree, low to high
        yield from self.Q_tree.traverse()

    def get_min(self):
        return self.Q_tree.get_min()

    def hi_neighbor(self,value):
        # expects a [x,y] value list (may or may not be an actual node value)
        return self.Q_tree.hi_neighbor(value)


def seg_int_pts(seg_1,seg_2):
    # checks for intersection point between segments
    # seg_1 and seg_2 of form [[x1,y1],[x2,y2]], with x1<=x2
    # five main possibilities: () no intersection; () endpoint-endpoint
    # intersection; () T intersection; () X intersection (single
    # "clean" intersection point); () "infinite" (lines overlap w/
    # measure > 0) intersection
    # returns (None,None) if no intersection; (None,float('inf')) if overlap
    # w/ measure > 0; otherwise returns ((segment 1 L/R/I, segment 2 L/R/I),
    # [x,y]);
    # note, for an "imprecise" T intersection, this returns the stem's
    # endpoint (and not the ~near internal point of the T's top)

    lr_arr = ['left','right']
    p11 = seg_1[0]
    p12 = seg_1[1]
    p21 = seg_2[0]
    p22 = seg_2[1]
    if abs(p11[0]-p21[0])<TOL_ACC and abs(p11[1]-p21[1])<TOL_ACC:
        # p11 and p21 are ~coincident
        return (('left','left'),p11)
    elif abs(p11[0]-p22[0])<TOL_ACC and abs(p11[1]-p22[1])<TOL_ACC:
        # p11 and p22 are ~coincident
        return (('left','right'),p11)
    elif abs(p12[0]-p21[0])<TOL_ACC and abs(p12[1]-p21[1])<TOL_ACC:
        # p12 and p21 are ~coincident
        return (('right','left'),p12)
    elif abs(p12[0]-p22[0])<TOL_ACC and abs(p12[1]-p22[1])<TOL_ACC:
        # p12 and p22 are ~coincident
        return (('right','right'),p22)
    tot_pts = [p11,p12,p21,p22]
    mn_vs = [list(map(lambda a,b:b-a,tot_pts[0+2*ii],tot_pts[1+2*ii])) for ii in
            range(2)]   # these are the p11-to-p12 and p21-to-p22 segment
    # vectors
    nrms = [math.sqrt(tup[0]**2+tup[1]**2) for tup in mn_vs]
    sgns = [0.0,0.0]
    axs_tst = [[False,tuple()],[False,tuple()]] # tests for endpoints on axis
    ept_int = [False,False] # tests if endpoint internal to segment
    for ii in range(2): # each pass is just a dot product
        # dot product with perpendicular:
        tmp = list(map(lambda a,b: a*b,list(map(lambda a,b:b-a,tot_pts[0],
                tot_pts[ii+2])),[mn_vs[0][1]/nrms[0],-mn_vs[0][0]/nrms[0]]))
        sgns[ii] = tmp[0]+tmp[1]
        if abs(sgns[ii])<TOL_ACC:   # point is effectively on axis of mn_vs[0]:
            # projection onto unit perpendicular gives distance away from
            # segment; if distance is within the similar "tube" as in
            # TriangleSegment's in_interior() method, then consider the point
            # on axis
            axs_tst[ii][0] = True
            # dot product with parallel:
            tmp_int = list(map(lambda a,b: a*b,list(map(lambda a,b:b-a,
                tot_pts[0],tot_pts[ii+2])),[mn_vs[0][0]/nrms[0],
                                            mn_vs[0][1]/nrms[0]]))
            tmp_pos = tmp_int[0]+tmp_int[1]
            axs_tst[ii][1] = (tmp_pos<-TOL_ACC,tmp_pos>nrms[0]+TOL_ACC)
            if not axs_tst[ii][1][0] and not axs_tst[ii][1][1]:
                ept_int[ii] = True  # endpoint ii of segment 2 is internal to
                # segment 1
    if axs_tst[0][0] or axs_tst[1][0]:  # either one or both endpoints of
        # segment p21-p22 are co-axial with segment p11-p12
        if axs_tst[0][0] != axs_tst[1][0]:  # only one endpoint is on segment
            # 1's axis
            for jj in range(2):
                if ept_int[jj]: # is on-axis endpoint internal to segment 1?
                    return (('internal',lr_arr[jj]),seg_2[jj])    # T-
                    # intersection
            return (None,None) # one is on the axis, and not an interior point,
            # and the other is not on the axis--no way to intesect
        else:   # all endpoints are co-axial; check if they straddle the segment
            if (axs_tst[0][1][0] and axs_tst[1][1][0]) or (axs_tst[
                    0][1][1] and axs_tst[1][1][1]):
                return (None,None) # segments are co-axial but disjoint
            else:
                return (None,float('inf')) # segments coincide over measure > 0
    elif (sgns[0]>0.0 and sgns[1]>0.0) or (sgns[0]<0.0 and sgns[1]<0.0):
        return (None,None) # check endpoints p21 and p22 lie either side of
        # segment 1; if not, return no intersection
    axs_tst = [[False, tuple()], [False, tuple()]]
    ept_int = [False, False]
    for ii in range(2):
        tmp = list(map(lambda a,b: a*b,list(map(lambda a,b:b-a,tot_pts[2],
                tot_pts[ii])),[mn_vs[1][1]/nrms[1],-mn_vs[1][0]/nrms[1]]))
        sgns[ii] = tmp[0]+tmp[1]
        if abs(sgns[ii]) < TOL_ACC:  # seg 1 point ii is effectively on axis of
            # mn_vs[1]
            axs_tst[ii][0] = True
            tmp_int = list(map(lambda a, b: a * b, list(map(lambda a,b: b-a,
                    tot_pts[2],tot_pts[ii])), [mn_vs[1][0]/nrms[1],
                                               mn_vs[1][1]/nrms[1]]))
            tmp_pos = tmp_int[0] + tmp_int[1]
            axs_tst[ii][1] = (tmp_pos < -TOL_ACC, tmp_pos > nrms[1] + TOL_ACC)
            if not axs_tst[ii][1][0] and not axs_tst[ii][1][1]:
                ept_int[ii] = True  # endpoint ii of segment 1 is internal to
                # segment 2
    if axs_tst[0][0] or axs_tst[1][0]:
        if axs_tst[0][0] != axs_tst[1][0]:  # only p11 or p12 are on-axis of 2
            for jj in range(2):
                if ept_int[jj]:  # is on-axis endpoint internal to segment 2?
                    return ((lr_arr[jj],'internal'), seg_1[jj])  # T
                    # intersection
            return (None, None)
        # (already checked for all co-axial)
    elif (sgns[0] > 0.0 and sgns[1] > 0.0) or (sgns[0] < 0.0 and sgns[1] < 0.0):
        return (None,None)
    # both pairs of endpoints straddle complementary segments; find actual
    # point of intersection
    rto = abs(sgns[0])/(abs(sgns[0])+abs(sgns[1]))  # proportional distance
    # from p21 to p22 to get to point of intersection
    return (('internal','internal'),[p11[0]+rto*mn_vs[0][0],p11[1]+rto*mn_vs[0][
        1]])


def check_segs(seg_lst):
    # ensures segments in seg_lst are oriented correctly--[p1,p2] in right
    # order:
    #   if segment is non-vertical, p1 is leftmost endpoint
    #   if segment is vertical, p1 is bottommost endpoint
    # changes mutable seg_lst in place
    for seg in seg_lst:
        if seg[0][1] != seg[1][1]:  # segment is non-vertical
            if seg[0][0]>seg[1][0]:
                tmp = seg[1]
                seg[1] = seg[0]
                seg[0] = tmp
        else:   # segment is vertical
            if seg[0][1]>seg[1][1]:
                tmp = seg[1]
                seg[1] = seg[0]
                seg[0] = tmp


def update_events(val,seg_lst,ix_1,ix_2,events):
    # update event tree; return [x,y] point if a "future" intersection has
    # been found
    tup,pt = seg_int_pts(seg_lst[ix_1],seg_lst[ix_2])
    if tup is None and pt==float('inf'):
        print("segments seem to overlap on infinitely many points; this could "
              "happen, for example, between two segments with a very small "
              "angle of intersection, or because of precision issues--eg try "
              "setting TOL_ACC larger")
        breakpoint()    # DEBUG
        raise ValueError
    elif tup is not None and pt > [val[0]-TOL_ACC,val[1]-TOL_ACC]:
        events.add_point(pt[0],pt[1],[(ix_1,tup[0]),
                                              (ix_2,tup[1])])
        return pt
    return None


def plt_sgs(seg_lst):
    # plots segments in seg_lst; requires open plot window in "plt"
    for jj,seg in enumerate(seg_lst):
        plt.plot([seg[0][0],seg[1][0]],[seg[0][1],seg[1][1]],color='blue')
        plt.text(seg[0][0]+(-seg[0][0]+seg[1][0])/2,seg[0][1]+
            (-seg[0][1]+seg[1][1])/2, str(jj),
                 fontsize=9, color='blue')


def run_sweep_line(inp_lst):

    # expects an input list inp_lst of type [[pt11,pt12],[pt21,pt22],
    # ...] amounting to a list of segments; returns a new updated segments list,
    # where the ordering of segments is preserved but points may be adjusted
    # slightly for glomming, and returns a segment event list;
    # each entry in returned segment event list is of form,
    #   [x,y],{(sg_ix,L/R/I),...},
    #   ie the point the event occurred at, [x,y], followed by a set;
    #   the set is one or more tuples, each tuple being a segment that goes
    #   through the point [x,y], and how it goes through it--left endpoint of
    #   segment, right endpoint of segment, or interior of segment;
    #   note, a lone segment endpoint is considered an event too, and will be
    #   returned in the event list--so to pull true inter-segment
    #   intersection points from the returned event list, have to filter out
    #   any entries that have an [x,y] associated with an event set with only
    #   one element in it (which necessarily has to be a left or right endpoint)
    seg_lst = copy.deepcopy(inp_lst)

    # pre-process for near-vertical segments
    vt_sg_ix = [] # indices of near-vertical segments (made truly vertical)
    tp_ls_ix = []
    ver_seg_chk = []
    for ii,seg in enumerate(seg_lst):
        if abs(seg[0][0] - seg[1][0]) <= TOL_ACC:   # if segment is
            # near-vertical,
            seg[0][0] = seg[1][0]   # make the segment truly vertical
            vt_sg_ix.append(ii)
            ver_seg_chk.append(seg) # for post-glomming check
        else:
            tp_ls_ix.append(ii)

    # plot initialization
    if DEBUG_PLOT and matpltlib_enabled:
        plt.ion()  # attempt to continually update plot window (not remove
        # and start over)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt_sgs(seg_lst)
        plt.show()

    # initial segment processing via the event class;
    # note an event amounts to [x,y] (event tree node value) and a set of
    # tuples, {(segment index,L/R/I),(segment index,L/R/I),...};
    check_segs(seg_lst)  # ensures segment endpoints are in proper order for
    # processing
    events = EventStorage()
    lr_label = ['left', 'right']
    seg_ixs = vt_sg_ix+tp_ls_ix   # indices ordered with vertical
    # segments (if any) first, then the rest; when placing segments in event
    # list tree, this ensures vertical segment endpoints don't get ~skewed
    for ix in seg_ixs:
        tup = seg_lst[ix]
        for jj, pt in enumerate(tup):  # tup always has lowest-x point first
            events.add_point(pt[0], pt[1], [(ix, lr_label[jj])])

    # update segment list from event tree, for cases where glomming may have
    # changed segment endpoints
    print("revised segment list (post initial glomming):")
    seg_dct = events.get_sgs()  # returns a dictionary
    seg_lst = [seg_dct[ii] for ii in range(len(seg_dct))]
    check_segs(seg_lst)  # because points can shift during glomming, ensure
    # each segment of seg_lst still has endpoints in correct order
    if ver_seg_chk != [seg_lst[ix] for ix in vt_sg_ix]:
        print("problem with initial near-vertical segments; more than 2 "
              "vertical segments may be too close together horizontally?")
    print(seg_lst)

    # handle vertical segments--glom nearby endpoints to on-segment in event
    # tree
    for ix in vt_sg_ix:
        events.glom_to_seg(seg_lst[ix])
    # update segment list, and check (again):
    seg_dct = events.get_sgs()  # returns a dictionary
    seg_lst = [seg_dct[ii] for ii in range(len(seg_dct))]
    check_segs(seg_lst)
    if ver_seg_chk != [seg_lst[ix] for ix in vt_sg_ix]:
        print("problem with initial near-vertical segments; more than 2 "
              "vertical segments may be too close together horizontally?")

    # @@@
    # sweep line method; determine segment intersections (not including
    # rectangular boundary window)
    # this follows DeBerg
    # @@@

    # initialize (empty) status tree
    T_status = status_tree.AVLStatusTree()

    val, dat = events.get_min()  # val=[x,y]; dat={(seg ix,L/R/I),()...}

    while val is not None:

        # DEBUG
        if DEBUG_PLOT and matpltlib_enabled:
            plt.scatter(*val, marker='o', s=400, color='red', facecolors='none')

        # classify segments at this event point
        # these can be:
        #   () isolated segment endpoints
        #   () intersections of segment interiors (C(p) by DeBerg)
        #   () intersections of segment right endpoints (L(p) analogs; before
        #   sweepline / interiors lie to left of sweepline)
        #   () intersections of segment left endpoints (U(p) analogs; after
        #   sweepline / interiors lie to right of sweepline)
        C_segs = [x[0] for x in dat if
                  x[1] == 'internal']  # indices of segments w/
        # interior intersections
        L_segs = [x[0] for x in dat if x[1] == 'right']  # segments ending
        U_segs = [x[0] for x in dat if x[1] == 'left']  # segments beginning

        # prepare for L's deletion
        if len(L_segs) > 0:

            # get associated ~list of node values in T_status associated with L
            L_nd_vs = [T_status.get_node_value(x) for x in L_segs]

            if len(U_segs) + len(C_segs) == 0:
                mx_Lv = max(L_nd_vs)
                L_hn, lh_dat = T_status.hi_neighbor(mx_Lv)  # [node value,
                # data] of T-neighbor just above highest segment in L;
                # None if none
                mn_Lv = min(L_nd_vs)
                L_ln, ll_dat = T_status.lo_neighbor(mn_Lv)  # [node value,
                # data] of T-neighbor just below lowest segment in L;
                # None if none

            # delete L's segments from status tree T
            for sg_ix in L_segs:
                if not (T_status.delete_segment(sg_ix)):
                    print("Error deleting an L node from T tree")  # DEBUG

        if len(C_segs) > 1:  # reverse order of segments that have this event
            # point as a mutually interior intersection point
            if len(C_segs) % 2 == 1:  # odd number of segments to reverse
                swp_lst = [(C_segs[ii], C_segs[-(ii + 1)]) for ii in
                           range(int((len(
                               C_segs) - 1) / 2))]
            else:
                swp_lst = [(C_segs[ii], C_segs[-(ii + 1)]) for ii in
                           range(int(len(
                               C_segs) / 2))]
            for tup in swp_lst:
                T_status.swap_segs(*tup)

        if len(U_segs) > 0:  # insert newly beginning segments into status
            # tree T
            for sg_ix in U_segs:
                sg_in = status_tree.Segment(seg_lst[sg_ix])
                T_status.insert(status_tree.NodeData(sg_ix, sg_in))

        if len(U_segs) + len(C_segs) == 0:
            if L_hn is not None and L_ln is not None:
                tmp = update_events(val, seg_lst, lh_dat.sg_ix, ll_dat.sg_ix,
                                    events)
                # DEBUG
                if DEBUG_PLOT and matpltlib_enabled:
                    if tmp is not None:
                        plt.scatter(*tmp, marker='o', s=400, color='green',
                                    facecolors='none')
        else:
            UC_nd_vs = [(x, T_status.get_node_value(x)) for x in
                        U_segs + C_segs]
            mx_ix, mx_UCv = max(UC_nd_vs, key=lambda x: x[1])
            UC_hn, uch_dat = T_status.hi_neighbor(mx_UCv)  # [node value,
            # data] of T-neighbor just above highest segment in U+C; None if
            # none
            if UC_hn is not None:
                tmp = update_events(val, seg_lst, uch_dat.sg_ix, mx_ix, events)
                # DEBUG
                if DEBUG_PLOT and matpltlib_enabled:
                    if tmp is not None:
                        plt.scatter(*tmp, marker='o', s=400, color='green',
                                    facecolors='none')
            mn_ix, mn_UCv = min(UC_nd_vs, key=lambda x: x[1])
            UC_ln, ucl_dat = T_status.lo_neighbor(mn_UCv)  # [node value,
            # data] of T-neighbor just below lowest segment in U+C; None if
            # none
            if UC_ln is not None:
                tmp = update_events(val, seg_lst, ucl_dat.sg_ix, mn_ix, events)
                # DEBUG
                if DEBUG_PLOT and matpltlib_enabled:
                    if tmp is not None:
                        plt.scatter(*tmp, marker='o', s=400, color='green',
                                    facecolors='none')

        if DEBUG_PLOT and matpltlib_enabled:
            plt.show()  # shows event point of focus (red) and future found
            # event point (green)
            plt.pause(0.5)
            # input("press a key: ")

        # go to next event in the event tree
        val, dat = events.hi_neighbor(val)

    # update segment list again from event tree, over chances segment
    # endpoints changed from event glomming
    seg_dct = events.get_sgs()
    seg_lst = [seg_dct[ii] for ii in range(len(seg_dct))]

    evt_lst = [ev for ev in events.traverse_events()]

    if DEBUG_PLOT:
        input("press a key: ")

    return seg_lst,evt_lst




#@@@
# MAIN
#@@@

# as a test:
if __name__=='__main__':
    sg_ls = [([0.07238034386144987, 0.48429578403727436], [0.11695825026706552, 0.05568015662335868]), ([0.3642893208466236, 0.8696977907292295], [0.2833464246350046, 0.9476079043289422]), ([0.05256128315104491, 0.8675140732794955], [0.20252569740199833, 0.2253640502605897]), ([0.6898278928312884, 0.3375948417411545], [0.49188423870733655, 0.43897942033171833]), ([0.4760620509632222, 0.027563786315106142], [0.9731694340269933, 0.38131978204717176]), ([0.4374697914472523, 0.9792353808130763], [0.6359690490261156, 0.32889197536191295]), ([0.8167642798473997, 0.6347460798979859], [0.3641090128969011, 0.7697377586324394]), ([0.5022086803049917, 0.32320352346634706], [0.9953309823763495, 0.7811686870687637]), ([0.6647216507793059, 0.5777750183980395], [0.4532196554908011, 0.21637346895703402]), ([0.10583012083111287, 0.2955320955812374], [0.5176170845893932, 0.48092853456115836])]
    run_sweep_line([list(x) for x in sg_ls])

