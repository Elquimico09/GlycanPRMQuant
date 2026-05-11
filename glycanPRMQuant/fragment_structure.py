"""Module for testing fragmentation of a specific glycan"""
from glypy.io import iupac
import pandas as pd
from glypy.io.iupac import IUPACError
from glypy.structure import ReducedEnd
from glypy.composition.composition_transform import derivatize

def fragment_glycan(iupac_str: str,
                    ion_series: str = "ABCXYZ",
                    max_cleavages: int = 2,
                    reduced_end: bool = True,
                    derivatization: str | None = "methyl") -> pd.DataFrame:
    """
    Fragments a glycan structure and returns the fragment information.

    Parameters
    ----------
    iupac_str : str
        The IUPAC name of the glycan.
    ion_series : str, optional
        The series of ions to fragment, by default "ABCXYZ".
    max_cleavages : int, optional
        Maximum number of glycosidic/cross-ring cleavages per fragment, by default 2.
    reduced_end : bool, optional
        Whether to use a reduced end, by default True.
    derivatization : str | None, optional
        The type of derivatization to apply, by default "methyl".

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the fragment information.
    """
    try:
        glycan = iupac.loads(iupac_str, dialect="simple")
    except IUPACError as e:
        raise ValueError(f"Error parsing IUPAC string: {e}")
    
    if reduced_end:
        glycan.set_reducing_end(ReducedEnd())
    if derivatization:
        derivatize(glycan, derivatization)
    
    PROTON_MASS = 1.007276466812
    AMMONIUM_MASS = 18.033826
    fragments = []

    for frag in glycan.fragments(kind = ion_series, max_cleavages=max_cleavages):
        neutral_mass = frag.mass() if callable(frag.mass) else frag.mass
        series = frag.series
        fragment_name = getattr(frag, "name", None) or str(frag)

        fragments.append({
            "name": fragment_name,
            "series": series,
            "neutral_mass": neutral_mass,
            "[M+H]+": neutral_mass + PROTON_MASS,
            "[M+2H]2+": (neutral_mass + 2 * PROTON_MASS) / 2,
            "[M+NH4]+": neutral_mass + AMMONIUM_MASS,
            "[M+NH4+H]2+": (neutral_mass + AMMONIUM_MASS + PROTON_MASS) / 2,
            "[M+2NH4]2+": (neutral_mass + 2 * AMMONIUM_MASS) / 2,
        })

    df = pd.DataFrame(fragments)
    return df

if __name__ == "__main__":
    iupac_str = "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    df = fragment_glycan(iupac_str)
    print(df)
