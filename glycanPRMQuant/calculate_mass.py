"""
Module for calculating the mass of glycans
"""
import logging

from glypy.io import iupac
from glypy.io.iupac import IUPACError
from glypy.structure import ReducedEnd
from glypy.composition.composition_transform import derivatize

logger = logging.getLogger(__name__)


def calculate_mass(glycan_str, derivatization="methyl", reduced_end=True, verbose=True):
    """
    Given an IUPAC string,
calculate the mass of the glycan, optionally applying derivatization and setting a reduced end.
    """
    try:
        glycan = iupac.loads(glycan_str, dialect="simple")
    except IUPACError as e:
        if verbose:
            logger.error("Error parsing IUPAC string: %s", e)
        return None

    try:
        if reduced_end:
            glycan.set_reducing_end(ReducedEnd())
        if derivatization:
            derivatize(glycan, derivatization)
        return glycan.mass()
    except (KeyError, ValueError) as e:
        if verbose:
            logger.error("Error modifying glycan: %s", e)
        return None

if __name__ == "__main__":
    # Example usage
    glycan_str = "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    mass = calculate_mass(glycan_str, derivatization=None, reduced_end=False)
    if mass is not None:
        print(f"The mass of {glycan_str} is: {mass:.4f} Da")
