import polars as pl

from realms.cluster import run_infomap
from realms.network import build_network
from realms.score import v_measure
from realms.synth import planted_network


def test_recovers_planted_realms():
    presence, _, truth = planted_network()
    partition = run_infomap(build_network(presence))
    assert v_measure(partition.leaf_labels(), truth) >= 0.95
    assert len(set(partition.leaf_labels().values())) == 4


def test_v_measure_detects_wrong_labels():
    # guards the check above: a single module must not score as a recovery
    _, _, truth = planted_network()
    assert v_measure({c: 0 for c in truth}, truth) == 0.0


def single_cell_species(families: dict[str, list[str]]) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Every cell has 3 species of its own, so no two cells share a species and only family ties them."""
    p, t = [], []
    for family, cells in families.items():
        for cell in cells:
            for i in range(3):
                species = f"{cell}-{i}"
                p.append((cell, species))
                t.append((species, f"G{cell}", family))
    return (
        pl.DataFrame(p, schema={"cell": pl.String, "species_id": pl.String}, orient="row"),
        pl.DataFrame(t, schema={"species_id": pl.String, "genus_id": pl.String, "family_id": pl.String}, orient="row"),
    )


def test_family_nodes_connect_cells_sharing_no_species():
    cells_a = [f"a{i}" for i in range(10)]
    cells_b = [f"b{i}" for i in range(10)]
    presence, taxa = single_cell_species({"FA": cells_a, "FB": cells_b})
    truth = {c: 0 for c in cells_a} | {c: 1 for c in cells_b}

    species_only = run_infomap(build_network(presence, taxa, ranks=("species",)))
    assert len(set(species_only.leaf_labels().values())) == 20

    with_family = run_infomap(build_network(presence, taxa, ranks=("species", "family")))
    assert v_measure(with_family.top_labels(), truth) == 1.0


def test_cosmopolitan_family_does_not_merge_realms():
    presence, taxa, truth = planted_network()
    in_every_cell = pl.DataFrame({"cell": list(truth), "species_id": [f"cosmo-{c}" for c in truth]})
    cosmo = pl.DataFrame(
        {"species_id": [f"cosmo-{c}" for c in truth], "genus_id": [f"Gc-{c}" for c in truth], "family_id": ["Cosmo"] * len(truth)}
    )
    net = build_network(pl.concat([presence, in_every_cell]), pl.concat([taxa, cosmo]))
    assert "family:Cosmo" not in net.taxa
    assert v_measure(run_infomap(net).leaf_labels(), truth) >= 0.95


def test_same_seed_same_result():
    presence, _, _ = planted_network()
    net = build_network(presence)
    assert run_infomap(net, seed=7).paths == run_infomap(net, seed=7).paths


def test_isolated_cell_is_reported_unassigned():
    presence, _, truth = planted_network(n_realms=2, cells_per_realm=5, endemic_species=10, widespread_species=0)
    with_empty = pl.concat([presence, pl.DataFrame({"cell": ["empty"], "species_id": [None]}, schema=presence.schema)])
    partition = run_infomap(build_network(with_empty))
    assert partition.unassigned == ["empty"]
    assert set(partition.paths) == set(truth)
