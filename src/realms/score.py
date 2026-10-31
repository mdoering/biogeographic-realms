"""Agreement between two maps of cells."""

from sklearn.metrics import v_measure_score


def v_measure(a: dict[str, int], b: dict[str, int]) -> float:
    """V-measure over the cells both maps contain. 1.0 is the same grouping, whatever the label numbers."""
    shared = sorted(a.keys() & b.keys())
    return float(v_measure_score([a[c] for c in shared], [b[c] for c in shared]))
