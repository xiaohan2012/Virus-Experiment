from ve.fp.residue_util.res_triangle import ResTriangle

from ve.fp.complex_util.paraepi import FindParaEpiTrait
from types import MethodType

from common import make_logger
logger = make_logger("Residue Triangle Detection")

from ve.fp.complex_util.cache import TriangleCache as cache

class ResidueTriangleTrait(FindParaEpiTrait):
    """
    find triangles in the epitope part of the complex.

    Dependency:
    1, epitope finding trait
    """
    def get_triangles(self, refresh = False):
        #cached in memory already
        if hasattr(self, "triangles"):
            return self.triangles
            
        #cached in local file already
        if not refresh and cache.has_cache(self.c_id):
            self.triangles = cache.load(self.c_id, self)
            return self.triangles
            
        #stores residue tuples that are possible to form triangle with external residue
        still_possible_pairs = []

        from itertools import combinations

        #residue triple that are within triangle formation range with each other
        initial_triangles = [ResTriangle([r1,r2,r3]) for r1,r2,r3 in combinations(self.find_epitope(), 3)
                     if r1.may_form_triangle_with(r2) and r2.may_form_triangle_with(r3) and r1.may_form_triangle_with(r3)]

        #determine if a residue pair are potential to form a triangle with an non-epitope residue
        is_potential = lambda (r1,r2): len([1 for r3 in (set(self.find_epitope()) - {r1,r2})
                                            if r1.may_form_triangle_with(r2) and (not r1.may_form_triangle_with(r3) or not r2.may_form_triangle_with(r3))]) > 0

        #pairs that are close to each other but no other can form triangle with the two
        still_possible_pairs = filter(is_potential, combinations(self.find_epitope(), 2))

        #antigen residues except the epitope ones
        marginal_residues = set(self.atg.residues)-set(self.epitope)

        #find residues that may form triangles with the still_possible pairs
        newly_detected_triangles = [ResTriangle([r1,r2,r3]) for r1,r2 in still_possible_pairs for r3 in (marginal_residues - {r1,r2})
                                    if r1.may_form_triangle_with(r3) and r2.may_form_triangle_with(r3)]

        self.triangles = initial_triangles + newly_detected_triangles

        logger.info("""
                    %d triangles initially
                    %d more triangles later
                    %d in total""" %(len(initial_triangles),
                                     len(newly_detected_triangles),
                                     len(set(self.triangles))))

        #cache it in local file
        cache.dump(self.c_id, self.triangles)
        
        return self.triangles


def main():
    from ve.util.complex import BaseComplex
    from ve.fp.fp_80 import Residue as Residue80

    class Complex(BaseComplex, ResidueTriangleTrait):
        pass

    from ve.util.load_pdb import load_complexes, complex_ids
    for c in load_complexes(complex_ids(), complex_cls = Complex, residue_cls = Residue80):
        c.get_triangles(refresh = True)

if __name__ == '__main__':
    main()
        