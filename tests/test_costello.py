import pytest

from realms.costello import csquare_adjacency, csquare_area_km2, csquare_bounds, load_matrix
from realms.maps import plot_cells


def test_csquare_se_quadrant():
    assert csquare_bounds("3603:3") == (-70.0, -65.0, 30.0, 35.0)


def test_csquare_nw_quadrant():
    assert csquare_bounds("7010:4") == (5.0, 10.0, -110.0, -105.0)


def test_csquare_ne_low():
    assert csquare_bounds("1403:1") == (40.0, 45.0, 30.0, 35.0)


def test_csquare_sw_quadrant():
    assert csquare_bounds("5117:2") == (-15.0, -10.0, -180.0, -175.0)


def test_area_equator_near_308k():
    assert csquare_area_km2("1000:1") == pytest.approx(308_700, rel=0.01)


def test_area_shrinks_towards_the_pole():
    assert csquare_area_km2("1800:3") < csquare_area_km2("1000:1") / 10


def test_adjacency_within_one_ten_degree_square():
    cells = ["1000:1", "1000:2", "1000:3", "1000:4"]
    adjacency = csquare_adjacency(cells)
    assert adjacency["1000:1"] == {"1000:2", "1000:3"}
    assert adjacency["1000:4"] == {"1000:2", "1000:3"}


def test_adjacency_crosses_the_equator():
    # 0..5N, 0..5E and 0..5S, 0..5E
    assert csquare_adjacency(["1000:1", "3000:1"])["1000:1"] == {"3000:1"}


def test_adjacency_wraps_dateline():
    # 0..5N, 175..180E and 0..5N, 175..180W
    assert csquare_adjacency(["1017:2", "7017:2"])["1017:2"] == {"7017:2"}


def test_cell_without_neighbours_has_empty_adjacency():
    assert csquare_adjacency(["1000:1", "1403:1"])["1403:1"] == set()


@pytest.fixture
def matrix(tmp_path):
    path = tmp_path / "m.csv"
    path.write_text("species,3603:3,5604:4,7109:4\nAus bus,1,0,1\nAus cus var,0,0,1\nDus eus,0,0,0\n")
    return load_matrix(path)


def test_load_matrix_small_fixture(matrix):
    assert matrix.species == ["Aus bus", "Aus cus var", "Dus eus"]
    assert matrix.cells == ["3603:3", "5604:4", "7109:4"]
    assert matrix.X.toarray().tolist() == [[1, 0, 1], [0, 0, 1], [0, 0, 0]]


def test_presence_lists_occupied_pairs(matrix):
    rows = sorted(matrix.presence().iter_rows())
    assert rows == [("3603:3", "Aus bus"), ("7109:4", "Aus bus"), ("7109:4", "Aus cus var")]


def test_genera_takes_first_word(matrix):
    genera = matrix.genera()
    assert dict(zip(genera["species_id"], genera["genus_id"])) == {"Aus bus": "Aus", "Aus cus var": "Aus", "Dus eus": "Dus"}
    assert genera["family_id"].null_count() == 3


def test_plot_cells_writes_an_image(tmp_path):
    out = tmp_path / "sub" / "map.png"
    plot_cells(["1000:1", "1000:2"], [0, 1], csquare_bounds, out, "two cells")
    assert out.read_bytes().startswith(b"\x89PNG")
