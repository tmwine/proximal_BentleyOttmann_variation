	This is a repository for a relatively low-dependency version of the 2-dimensional [Bentley-Ottmann](https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm) sweep line algorithm with proximal event detection ability. It employs a vertical sweep line, moving left to right. Features include:
- accepts segments not in "general position"--that is, segments can intersect endpoint-to-endpoint, endpoint-to-interior, and there can be > 2 segments through any given intersection point, and segments can be purely vertical or purely horizontal (following the procedure of de Berg [^1])
-  handles approximate intersections as best it can, via a tunable precision parameter--that is, each segment may be imagined to have a proximal cylinder around it, where if another segment is within this proximal cylinder it is counted as an intersection

	As an example of where proximal event detection may be useful, consider two segments forming a T intersection: segment 1 from (0,2) to (2,2), and segment 2 from (1,2) to (1,0). The intersection point and the upper endpoint of the stem of the T are in this case mathematically "perfect," and will be regarded as the same point by any simple intersection routine. However, if we rotate coordinates by some non-trivial angle via matrix transformation, there is no longer any guarantee the T intersection point and the upper endpoint of the T's stem are exactly the same point--there may not in fact be any intersection point at all. Depending on precision issues, the upper endpoint of the stem of the T may be slightly above or slightly below the upper segment of the T.

## Interface
	Clone or download the modules, then import sweep_line into your program. To run,
```
segment_list, event_list = sweep_line.run_sweep_line(input_segments)
```
where input_segments is a list of segments for intersection detection, of type [seg_1,seg_2,...], where each seg_i is of type [[x1,y1],[x2,y2]]. The routine returns an updated segment list, where the initial segment list points may change because of proximal glomming, and an event list. The event list is of type [event_1, event_2,...], where each event_i is of type ([x,y],{(sg_ix,L/R/I),...}). The coordinate of the event is recorded in the first element, [x,y], while the subsequent set has all segment elements associated with that event point. Each entry in the set is a tuple of type (segment index number, segment intersection details--'left,' 'right,' or 'interior'). Note that event_list will include lone segment endpoints, if any exist. For true intersections, lone segment endpoints in event_list (i.e. event [x,y]'s where the associated set has only one element) will have to be filtered out.
	The routine also offers a plotting feature, which allows watching the sweep line routine find intersections as it sweeps. To enable this feature, matplotlib needs to be installed, and the DEBUG_PLOT flag in the sweep_line module header needs to be set to True.

## Implementation notes
	The degree of tolerance for approximate intersections can be adjusted through the TOL_ACC parameter in the event_tree module. Each segment has a blur or buffer around it of radius TOL_ACC. That is, the proximal boundary involves a rectangle around each segment, where the sides of the rectangle are at distance TOL_ACC from the segment (with the ends squared, not rounded). Proximal endpoint-endpoint intersections, proximal T intersections, and mutual proximal intersections in general are handled by glomming or snapping (grouping) the involved points into a single point.
	Also note, as a result of glomming as we go, the segments can get broken into piecewise subsegments (between event points--intersections or endpoints), where the resulting subsegments are only approximately collinear (collinear under the tube condition, but not strictly collinear).

## Snapping / glomming details
For details on glomming:
- in pre-processing, near-vertical segments, within TOL_ACC, are nudged to true vertical
- vertical segments are loaded into the event tree first; they have priority with respect to endpoint glomming (a worst case, two or more vertical segments are so close they are considered the same--this should throw an error)
- next, all other segments are loaded into the event tree; as each segment is added, if its endpoint is within tolerance of any other endpoint, that endpoint gets glommed / snapped to endpoint of segment already in event tree
- a pass is made to glom segment endpoints within tolerance TOL_ACC of any vertical segment to be exactly on that vertical segment (again, worst case, 2 or more vert segs are so close they have mutually interior intersections--this again should throw an error)
- sweep line is run; for proximal T intersections, the endpoints of the stem of the T get priority--the intersection point is glommed to the upper endpoint of the vertical stem of the T; in general, if a newly-found intersection point is proximal to an existing event point in the tree (segment endpoints count), that point gets glommed / snapped to the point already in the tree

## Efficiency notes
	The routine follows that in DeBerg [^1], however because of some ambiguity in how the status tree was implemented, the routine may not reach the true optimal efficiency of pure Bentley-Ottmann. The status tree in this implementation is an AVL tree, with a separate procedure for adding a segment to it (by node data / metadata)--separate from its internal sorting mechanism (by node value). The status tree also does self-audits of node values, and redistributes node values if they get too close (such redistributing does not affect node topology). My understanding is that this may not be an optimal implementation for Bentley-Ottmann.
	Another efficiency drawback for this implementation is the case where many vertical or near-vertical segments lie within a narrow TOL_ACC horizontal band. There is a bottleneck in the add_point() routine in the event tree, where it checks a newly-added point against "neighboring" existing points in the tree, for the purposes of glomming proximal events together. A range search is conducted for points "near" the input point, but the range search is done lexicographically under node values of type [x,y]. In the example of event x,y's, [[0.0,0.0],[0.0,1.0],[0.0,2.0],[0.1,0.0],[0.1,0.5]], if we want the lexicographic range between [-0.3,-0.3] and [0.3,0.3], the search will return [[0.0,0.0],[0.0,1.0],[0.0,2.0],[0.1,0.0]], and each of these need to be checked for being within the desired taxicab distance of 0.3 from [0.0,0.0]. This potential slowdown could be remedied by rotating the segement space so that such vertical groupings don't occur, or (better) converting the event AVL tree to a k-d tree.
	Under the assumption there is little such vertical stacking, the lexicographic range searches become more like the desired taxicab distance searches, and the overall efficiency should approach O((N+k)log(N+k)), where k is the number of actual intersection points.

## Final notes
	This routine was written to help with a specific segment intersection task, for which it worked adequately. No promises are made on proximal event detection for arbitrary collections of line segments. For instance, it is possible the glomming can suffer from corner cases that can skew the desired proximal event detection, but we don't consider or try to imagine those here.

## Licenses
	The code generated by this author is licensed under the MIT license. Support for AVL trees was provided by this [repository](https://github.com/cool-pot/pytrees), under the MIT license, which in turn references this [repo](https://github.com/pgrafov/python-avl-tree/blob/master/pyavltree.py), also under an MIT license.




[^1]: de Berg, Mark; van Kreveld, Marc; Overmars, Mark; Schwarzkopf, Otfried (2000), "Chapter 2: Line segment intersection," Computational Geometry (2nd ed.), Springer-Verlag, pp. 19-44.
