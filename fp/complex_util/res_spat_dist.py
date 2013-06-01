from ve.fp.complex_util.split_cylinder import Cylinder
from collections import defaultdict

class CylinderWithResPos(Cylinder):
    def get_count_distribution(self, residues):
        """
        (CylinderWithResPos, list of Residue) -> ResiduePositionDistribution
        
        get the residue spatial distribution
        """

        return ResiduePositionDistribution(self.ring_count, [self.in_which(r.get_center()) for r in residues if self.in_which(r.get_center())])
        
    def get_residue_distribution(self, residues):
        """
        (CylinderWithResPos, list of Residue) -> dict(int -> list of Residue)

        return a dict of cylinder ring index to a list of subordinating residues
        """
        d = defaultdict(list)
        for r in residues:
            idx = self.in_which(r.get_center())
            if idx:
                d[idx].append(r)
        return d
        
from collections import Counter

class ResiduePositionDistribution(Counter):
    """
    A list recording the maximum residue number in each cylinder segment

    >>> t = ResiduePositionDistribution(4, [2,3])
    >>> t[3]
    1
    >>> t += ResiduePositionDistribution(4, [1,2,2,2])
    >>> t.ring_count
    4
    >>> t[1]
    1
    >>> t[2]
    3
    >>> len(t)
    5
    >>> print t
    idx:  0, 1, 2, 3
    cnt:  0, 1, 3, 1
    """
    def __init__(self, ring_count, *args, **kws):
        self.ring_count = ring_count
        super(ResiduePositionDistribution,self).__init__(*args, **kws)

    def __getitem__(self, idx):
        if idx >= self.ring_count:
            raise IndexError("%d exceeds %d" %(idx, self.ring_count))
        return super(ResiduePositionDistribution,self).__getitem__(idx)
            
    def __add__(self, other):
        """
        merge two position distributions
        
        (ResiduePositionDistribution, ResiduePositionDistribution) -> ResiduePositionDistribution
        """
        d = self.__class__(self.ring_count, self)
        for pos, val in other.items():
            if self[pos] < val:
                d[pos] = val
            else:
                d[pos] = self[pos]
        return d
        
    def __str__(self):
        row_head = ",".join(map(lambda x: "%2d" %x,xrange(self.ring_count)))
        row_content = ",".join(map(lambda x: "%2d" %x, (self[i] for i in xrange(self.ring_count))))
        return "idx: %s\ncnt: %s" %(row_head, row_content)

    def __repr__(self):
        return str(self)
        
    def __len__(self):
        return sum((self[i] for i in xrange(self.ring_count)))

from ve.fp.complex_util.geom import GeometryTrait
from ve.fp.geom import Line, get_perp_plane

from ve.fp.complex_util.axial_plane import HasAxialPlaneTrait
from ve.fp.complex_util.split_cylinder import SplitCylinderUtility
from ve.fp.complex_util.geom import GeometryTrait

class ResidueSpatialDistributionTrait(SplitCylinderUtility, HasAxialPlaneTrait, GeometryTrait):
    """
    Distribution on residue spatial position with respect to the splitted cylinder
    """
    
    def make_cylinder(self, center, plane):
        return CylinderWithResPos(center, plane, self.radius, self.radius_step, self.height, self.height_step)

    def sort_res_dist(self, res_dist):
        """
        (CylinderSpatialDistribution, dict(int, list of Residue)) -> dict(int, list of Residue)
        """
        return sorted(res_dist, key = lambda r: self.c.get_paraepi_center().dist2point(r.get_center()))
        
    def get_atg_res_spat_dist(self):
        """
        (CylinderSpatialDistritbution) -> Counter

        return the count of antigen residues that fall in the respective rings in the cylinder
        """
        if not hasattr(self, "atg_res_spat_dist"):
            plane = get_perp_plane(self.get_atg_center() - self.get_epitope_center(), self.get_epitope_center())
            cylinder = self.make_cylinder(self.get_epitope_center(), plane)
            self.atg_res_cnt_spat_dist = cylinder.get_count_distribution(self.atg.residues)
            self.atg_res_spat_dist = cylinder.get_residue_distribution(self.atg.residues)
            
        return self.atg_res_cnt_spat_dist, self.atg_res_spat_dist

    def get_atb_res_spat_dist(self):
        """
        (CylinderSpatialDistritbution) -> Counter

        return the count of antibody residues that fall in the respective rings in the cylinder
        """
        if not hasattr(self, "atb_res_spat_dist"):
            plane = get_perp_plane(self.get_atb_center() - self.get_paratope_center(), self.get_paratope_center())
            cylinder = self.make_cylinder(self.get_paratope_center(), plane)
            self.atb_res_cnt_spat_dist = cylinder.get_count_distribution(self.atb.residues)
            self.atb_res_spat_dist = cylinder.get_residue_distribution(self.atb.residues)

        return self.atb_res_cnt_spat_dist, self.atb_res_spat_dist

if __name__ == '__main__':
    import doctest
    doctest.testmod()