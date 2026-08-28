import pandas as pd
import pytest

from glycanPRMQuant.constants import METHANOL_MASS
from glycanPRMQuant.fragment_structure import fragment_glycan
from glycanPRMQuant import matchMS2 as ms2_module


IUPAC = "Neu5Ac(a2-6)Gal(b1-4)GlcNAc"
COMPOSITION = "11010"


def _fragment_database():
    fragments = fragment_glycan(IUPAC, ion_series="BY", max_cleavages=1)
    fragments = fragments.copy()
    fragments["Fragment"] = fragments["name"]
    fragments["FragmentType"] = fragments["series"]
    fragments["IUPAC"] = IUPAC
    fragments["NumericalComposition"] = COMPOSITION
    fragments["Composition"] = "{'Hex': 1, 'HexNAc': 1, 'Neu5Ac': 1}"
    return fragments


def test_neuac_fragments_include_charge_aware_methanol_loss():
    ions = ms2_module._build_adduct_table(_fragment_database(), COMPOSITION)
    losses = ions.loc[ions["neutral_loss"].eq("CH3OH")]
    intact = ions.loc[ions["neutral_loss"].eq("")]

    assert not losses.empty
    assert losses["contains_neuac"].all()
    assert losses["fragment_annotation"].str.endswith("-CH3OH").all()
    assert (~ions["contains_neuac"]).any()
    assert not ions.loc[~ions["contains_neuac"], "neutral_loss"].eq("CH3OH").any()

    for loss in losses.itertuples(index=False):
        base = intact.loc[
            intact["Fragment"].eq(loss.Fragment)
            & intact["Adduct"].eq(loss.Adduct)
        ].iloc[0]
        assert base["Theo_mz"] - loss.Theo_mz == pytest.approx(
            METHANOL_MASS / loss.Charge
        )


def test_methanol_loss_annotation_is_retained_in_match(tmp_path):
    ms2_module._FRAG_DB_CACHE.clear()
    db_path = tmp_path / "glycans.csv"
    pd.DataFrame(
        {
            "Numerical Composition": [COMPOSITION],
            "Composition": ["{'Hex': 1, 'HexNAc': 1, 'Neu5Ac': 1}"],
            "Condensed IUPAC": [IUPAC],
        }
    ).to_csv(db_path, index=False)

    fragment_db = ms2_module._build_fragment_db(
        str(db_path),
        COMPOSITION,
        ion_series="BY",
        max_cleavages=1,
    )
    ions = ms2_module._build_adduct_table(fragment_db, COMPOSITION)
    loss = ions.loc[
        ions["neutral_loss"].eq("CH3OH") & ions["Charge"].eq(1)
    ].iloc[0]

    matched = ms2_module.matchMS2(
        pd.DataFrame(
            {
                "scan_number": [1],
                "rt": [5.0],
                "precursor_mz": [500.0],
                "fragment_mz": [loss["Theo_mz"]],
                "fragment_intensity": [1000.0],
            }
        ),
        pd.DataFrame(
            {
                "precursor_mz": [500.0],
                "Glycan": [COMPOSITION],
                "Adduct": ["2H"],
            }
        ),
        precursor_composition=COMPOSITION,
        fragment_mass_tol=0.001,
        intensity_threshold=100.0,
        ppm_tol=10,
        db_path=str(db_path),
        ion_series="BY",
        max_cleavages=1,
    )

    assert matched.iloc[0]["neutral_loss"] == "CH3OH"
    assert matched.iloc[0]["Fragment"].endswith("-CH3OH")
    assert matched.iloc[0]["fragment_annotation"].endswith("-CH3OH")
    assert "Neu" in matched.iloc[0]["fragment_iupac"]
