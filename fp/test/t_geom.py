import unittest
from ve.fp.geom import *

import math

from common import NumericTestCase

class ArrayTest(NumericTestCase):
    def setUp(self):
        self.a1 = Array([1,2,3])
        self.a2 = Array([2,3,4])
        
    def test_add(self):

        self.assertArrayEqual(Array([3, 5, 7]), self.a1 + self.a2)

    def test_sub(self):
        self.assertArrayEqual(Array([-1, -1, -1]), self.a1 - self.a2)

    def test_mult_array(self):
        self.assertArrayEqual(Array([2, 6, 12]), self.a1 * self.a2)

    def test_mult_number(self):
        self.assertArrayEqual(Array([2, 4, 6]), self.a1 * 2)

    def test_pow(self):
        self.assertArrayEqual(Array([1, 4, 9]), self.a1 ** 2)

class PointTest(NumericTestCase):
    def setUp(self):
        self.p1 = Point([1,2,3])
        self.p2 = Point([2,3,4])
        
    def test_add(self):
        actual = Point([3, 5, 7]) == (self.p1 + self.p2)
        self.assertTrue(actual.all())

    def test_sub(self):
        actual = Point([-1, -1, -1]) == (self.p1 - self.p2)
        self.assertTrue(actual.all())
        
    def test_dis2point(self):
        self.assertArrayEqual(self.p1.dist2point(self.p2), math.sqrt(3.0))

class GetPerpPlaneTest(NumericTestCase):
    """Test case for the get_perp_plane function"""

    def test_general_case(self):
        vec = Vector([1,3,-2])
        pt = Point([1,2,3])
        actual = get_perp_plane(vec, pt)
        expected = Plane([1,3,-2,1])
        self.assertArrayEqual(actual, expected)

class PlaneTest(NumericTestCase):
    def setUp(self):
        self.plane = Plane([1,-2,3,5])
        self.point = Point([2,3,1])

    def test_dist2point_general_case(self):
        self.assertEqual(self.plane.dist2point(self.point), -6 / math.sqrt(14))

    def test_dist2point_negative_dist(self):
        plane = Plane([1,-2,3,100])
        self.assertAlmostEqual(plane.dist2point(self.point), -26.9933854332)


    def test_plane_perp_point(self):
        actual = self.plane.get_perp_point(self.point)
        expected = Point([2.4285714,2.1428571,2.2857143])
        self.assertArrayAlmostEqual(expected, actual)


class LineTest(NumericTestCase):
    def test_line_point_dist(self):
        p1 = Point([1,0,1])
        p2 = Point([0,1,1])
        line = Line(p1,p2)
        point = Point([1,1,1])
        
        actual = line.dist2point(point)
        expected = 0.707106781187

        self.assertAlmostEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
