import pytest

from realms.patterns import HUMAN_PATTERNS, HumanPattern, is_human_held, like_match


def test_like_match_is_case_insensitive():
    assert like_match("Jardin BOTANIQUE de Paris", "jardin botanique")


def test_like_match_null_is_false():
    assert not like_match(None, "botanic")


def test_like_match_needs_the_term():
    assert not like_match("Forest edge", "botanic")


def test_zoo_needs_spaces():
    assert like_match("near zoo here", " zoo ")
    assert not like_match("Department of Zoology", " zoo ")


# real text from the sample of GBIF records, as (kingdom, column, text, held by people?)
RECORDS = [
    ("Plantae", "locality", "Oxford Botanic Garden, VC23 Oxfordshire", True),
    ("Animalia", "locality", "Santa Cruz Botanic Garden", False),
    ("Plantae", "locality", "Jardín Botánico de Bogotá", True),
    ("Animalia", "occurrenceremarks", "Jardin botanique de Montréal", False),
    ("Animalia", "locality", "Philadelphia Zoological Garden", True),
    ("Plantae", "locality", "Boulevard du Jardin Zoologique", False),
    ("Plantae", "locality", "Niabi Zoo, adjacent to Niabi Zoo along two track running along zoo fence", False),
    ("Animalia", "occurrenceremarks", "Sabrina - zoo animal|received at NBM 7 November 2011", True),
    ("Animalia", "occurrenceremarks", "Zoo specimen - in captivity", True),
    ("Animalia", "occurrenceremarks", "Not a Zoo animal", False),
    ("Animalia", "occurrenceremarks", "the tag is in the zoo but it was not a captive zoo animal; it was just a wild snake", False),
    ("Animalia", "occurrenceremarks", "Faecal sample from Zoo animal", True),
    ("Animalia", "locality", "Berlin Zoo", True),
    ("Animalia", "locality", "Zoo Berlin", True),
    ("Animalia", "locality", "Dudley Zoo, Staffordshire", True),
    ("Animalia", "locality", "the old zoo.", True),
    ("Animalia", "locality", "Department of Zoology", False),
    ("Animalia", "locality", "Berlin, Tierpark", True),
    ("Animalia", "locality", "East London, Aquarium Rocks", False),
    ("Animalia", "occurrenceremarks", "From tropical greenhouse - Found by Sean Birk Bek", True),
    ("Plantae", "habitat", "Cultivated, Glasshouse", True),
    ("Plantae", "habitat", "Weed on bare soil on cultivated land", False),
    ("Animalia", "occurrenceremarks", "On cultivated Quercus coccifera", False),
    ("Plantae", "occurrenceremarks", "Cultivated tree dominating manuka and cultivated shrubs.", True),
    ("Plantae", "habitat", "Champ cultivé.", False),
    ("Plantae", "habitat", "Arbre cultivé.", True),
    ("Plantae", "occurrenceremarks", "No cultivada, creciendo naturalmente en el patio de casa.", False),
    ("Plantae", "occurrenceremarks", "Cultivada por campesinos", True),
    ("Plantae", "occurrenceremarks", "nicht kultiviert, eingewandert", False),
    ("Plantae", "habitat", "Kultiviert [cultivated]", True),
    ("Plantae", "occurrenceremarks", "Planted. Planted. Planted.", True),
    ("Animalia", "locality", "Manitoba, Captive, Captive Born, 19940700,", True),
    ("Animalia", "occurrenceremarks", "Captive Birth/Hatch", True),
    ("Animalia", "occurrenceremarks", "please do not mark this observation as “captive” or “dead”", False),
    ("Animalia", "occurrenceremarks", "2 WILD FRECKLED DUCK SEEN - NOT PART OF CAPTIVE BREEDING SCHEME", False),
    ("Animalia", "locality", "In captivity, San Jose", True),
    ("Animalia", "occurrenceremarks", "spotted perching on bouldering wall (not in captivity!)", False),
    ("Animalia", "habitat", "Garden compost", False),
    ("Plantae", "habitat", "Skogsbeite 100 m N for garden", False),
]


@pytest.mark.parametrize("kingdom,column,text,expected", RECORDS)
def test_real_record_text(kingdom, column, text, expected):
    assert is_human_held({column: text}, kingdom) is expected


def test_kingdom_is_case_insensitive():
    assert is_human_held({"locality": "Oxford Botanic Garden"}, "PLANTAE")


def test_unknown_kingdom_only_gets_the_terms_for_every_kingdom():
    assert not is_human_held({"locality": "Oxford Botanic Garden"}, None)
    assert is_human_held({"locality": "greenhouse"}, None)


def test_missing_text_is_not_held():
    assert not is_human_held({"locality": None, "habitat": None}, "Plantae")


def test_dropped_terms_stay_dropped():
    assert not {p.term for p in HUMAN_PATTERNS} & {"garden", "hortus", "aquarium"}


def test_a_pattern_can_be_limited_to_one_column():
    assert not is_human_held({"occurrenceremarks": "captive"}, "Animalia")
    assert is_human_held({"locality": "captive"}, "Animalia")


def test_human_pattern_defaults_to_every_kingdom():
    assert HumanPattern("locality", "x").kingdom is None and HumanPattern("locality", "x").unless == ()
