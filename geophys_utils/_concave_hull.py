#!/bin/env python
"""
Calculate the concave hull of a set of points.

Adapted from Adriano Moreira and Maribel Yasmina Santos 2007.
"""

import numpy as np
import scipy.spatial as spt
from matplotlib.path import Path
import logging

logger = logging.getLogger(__name__)
logger.level = logging.INFO


def bbox(a, b):
    return {
        'll_x': min(a[0], b[0]),
        'll_y': min(a[1], b[1]),
        'ur_x': max(a[0], b[0]),
        'ur_y': max(a[1], b[1])
    }


def doBoundingBoxesIntersect(b1, b2):
    """
    Check if bounding boxes do intersect. If one bounding box touches
    the other, they do intersect.
    """
    if b1['ll_x'] > b2['ur_x']:
        return False

    if b1['ur_x'] < b2['ll_x']:
        return False

    if b1['ll_y'] > b2['ur_y']:
        return False

    if b1['ur_y'] < b2['ll_y']:
        return False

    return True


def isPointOnLine(a, b, c):
    """ Check if a point is on a line. """
    # move to origin
    bTmp = (b[0] - a[0], b[1] - a[1])
    cTmp = (c[0] - a[0], c[1] - a[1])
    r = np.cross(bTmp, cTmp)
    return np.abs(r) < 0.0000000001


def isPointRightOfLine(a, b, c):
    """
    Check if a point (c) is right of a line (a-b).
    If (c) is on the line, it is not right it.
    """
    # move to origin
    bTmp = (b[0] - a[0], b[1] - a[1])
    cTmp = (c[0] - a[0], c[1] - a[1])
    return np.cross(bTmp, cTmp) < 0


def lineSegmentTouchesOrCrossesLine(a, b, c, d):
    """
    Check if line segment (a-b) touches or crosses
    line segment (c-d).
    """
    return isPointOnLine(a, b, c) or \
           isPointOnLine(a, b, d) or \
          (isPointRightOfLine(a, b, c) ^
           isPointRightOfLine(a, b, d))


def doLinesIntersect(a, b, c, d):
    """ Check if line segments (a-b) and (c-d) intersect. """
    if not lineSegmentTouchesOrCrossesLine(a, b, c, d):
        return False

    return lineSegmentTouchesOrCrossesLine(c, d, a, b)


class PointSet:
    """ Book-keeping for an array of points. """

    def __init__(self, points):
        """ Create a kD-tree from the points. 
        @param points: n x 2 array of point coordinates
        """
        npoints, ndim = points.shape
        assert ndim == 2, 'Coordinates must be 2D, not {}D'.format(ndim)
        
        self.points = points
        
        # keep track of which points are currently still under consideration
        self.registry = np.full((npoints, ), fill_value=False, dtype='bool')
        
        # Only use first occurrence of unique values
        _unique_values, unique_value_indices = np.unique(points, return_index=True, axis=0)
        self.registry[unique_value_indices] = True
        
        # Ignore any coordinates containing NaN
        self.registry[np.any(np.isnan(points), axis=1)] = False 
        
        self.npoints = np.count_nonzero(self.registry)
        logger.debug('Computing shape for {} valid points'.format(self.npoints))
        
        self.tree = spt.cKDTree(self.points, leafsize=10)

        #assert False, 'ABORT'


    def __getitem__(self, index):
        return self.points[index]

    def start(self):
        """ The index of the starting point (the point with the lowest y-value). """
        for i in np.argsort(self[:, 1]):
            if self.registry[i]:
                return i

        raise ValueError("no starting point to pick")

    def remove(self, index):
        self.registry[index] = False
        self.npoints = self.npoints - 1

    def restore(self, index):
        self.registry[index] = True
        self.npoints = self.npoints + 1

    def nearest_neighbors(self, current_point, k):
        """
        Return the indices of the k nearest neighbors (or fewer if not enough points left).
        """
        kk = k

        while True:
            distances, indices = self.tree.query(current_point, kk)
            indices = [index for index in indices if self.registry[index]]
            assert len(indices) <= self.npoints
            if len(indices) == k or len(indices) == self.npoints:
                return indices
            kk = kk + 1


def TurningAngle(NearestPoint, currentPoint, previousAngle):
    angle = np.arctan2(NearestPoint[1] - currentPoint[1],
                       NearestPoint[0] - currentPoint[0]) - previousAngle
    angle = np.rad2deg(angle)
    if angle <= 0.:
        angle += 360.
    return angle


def GetNearestNeighbors(dataset, point, k):
    return dataset[PointSet(dataset).nearest_neighbors(point, k), :]


def sort_by_angle_indices(kNearestPoints, currentPoint, prevPoint):
    """ Sorts the k nearest points given by angle. """
    angles = np.zeros(kNearestPoints.shape[0])
    previousAngle = np.arctan2(prevPoint[1] - currentPoint[1],
                               prevPoint[0] - currentPoint[0])

    i = 0
    for NearestPoint in kNearestPoints:
        # calculate the angle
        # only positive angles
        angles[i] = TurningAngle(NearestPoint, currentPoint, previousAngle)
        i = i + 1
    return np.argsort(angles)


def SortByAngle(kNearestPoints, currentPoint, prevPoint):
    return kNearestPoints[sort_by_angle_indices(kNearestPoints, currentPoint,
                                                prevPoint), :]


def first_valid_candidate(point_set, cPoints, hull, first, step):
    # avoid intersections: select first candidate that does not intersect any
    # polygon edge
    def intersects(candidate, lastPoint):
        a = point_set[candidate]
        b = point_set[hull[step - 2]]
        b1 = bbox(a, b)

        for j in range(2, len(hull) - lastPoint):
            c = point_set[hull[step - j - 2]]
            d = point_set[hull[step - j - 1]]
            b2 = bbox(c, d)

            if doBoundingBoxesIntersect(b1, b2) and \
                    doLinesIntersect(a, b, c, d):
                return True

        return False

    for candidate in cPoints:
        if candidate == first:
            lastPoint = 1
        else:
            lastPoint = 0

        if not intersects(candidate, lastPoint):
            return candidate

    return None


def concave_hull_indices(dataset, k):
    '''\
    '''
    logger.debug('k in concave_hull_indices: {}'.format(k))
    assert k >= 3, 'k has to be greater or equal to 3.'
    assert k <= len(dataset), 'k has to be less than or equal to {}.'.format(len(dataset))

    point_set = PointSet(dataset)
    # todo: make sure that enough points for a given k can be found

    first = point_set.start()
    # init hull as list to easily append stuff
    hull = [first]
    # and remove it from dataset
    point_set.remove(first)

    current = first
    # set prevPoint to a Point righ of currentpoint (angle=0)
    prev = np.array([point_set[current, 0] - 1, point_set[current, 1]])
    step = 2

    while (first != current or (step == 2)) and point_set.npoints > 0:
        if step == 5:
            # we're far enough to close too early, so put first back again
            point_set.restore(first)

        k_nearest = point_set.nearest_neighbors(point_set[current], k)
        cPoints = sort_by_angle_indices(point_set[k_nearest, :],
                                        point_set[current], prev)
        cPoints = np.array(k_nearest)[cPoints]
        prev = point_set[current]

        current = first_valid_candidate(point_set, cPoints, hull, first, step)
        if current is None:
            return concave_hull_indices(dataset, k * 2)

        # add current point to hull
        hull.append(current)
        point_set.remove(current)

        step = step + 1

    # check if all points are inside the hull
    p = Path(point_set.points[hull, :])

    # Check whether all valid points are contained
    pContained = p.contains_points(dataset[np.all(~np.isnan(dataset), axis=1)], radius=0.0000000001) # Check filtered points with no NaNs
    logger.debug('{}/{} valid points contained'.format(np.count_nonzero(pContained), np.count_nonzero(np.all(~np.isnan(dataset), axis=1))))
    if not pContained.all():
        return concave_hull_indices(dataset, k * 2)

    return hull


def concaveHull(dataset, k=3):
    '''\
    Generate n x 2 array of coordinates for vertices of concave hull
    '''
    assert k >= 3, 'k has to be greater or equal to 3.'
    logger.debug('dataset length in concaveHull: {}'.format(len(dataset)))
    
    points = np.unique(dataset[~np.any(np.isnan(dataset), axis=1)], axis=0) # Purge duplicates and NaNs
    logger.debug('{} valid points in concaveHull'.format(len(points)))
    
    point_indices = concave_hull_indices(points, k)
    return points[point_indices, :]

