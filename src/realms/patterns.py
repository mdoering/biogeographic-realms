"""Text in a record that gives away a garden, zoo or other human-held organism.

These are counted, not filtered, so the effect of removing them can be measured. Matching is case insensitive and substring
based, the same as the SQL: `LOWER(COALESCE(column, '')) LIKE '%term%'`.

The terms were chosen by reading samples of real records (see DECISIONS.md). A place only says something about one kind of
organism: a wild insect in a botanic garden isn't cultivated, and a weed by a zoo fence isn't captive. Hence the kingdoms.
Being held indoors in a greenhouse says the same for anything. A false hit costs little, because a species only loses its
presence in a cell if every record of it there is hit.
"""

from collections.abc import Sequence
from dataclasses import dataclass

COLUMNS = ("locality", "habitat", "occurrenceremarks")
PLACE_COLUMNS = ("locality", "habitat")


@dataclass(frozen=True)
class HumanPattern:
    column: str
    term: str
    kingdom: str | None = None
    unless: tuple[str, ...] = ()


def _patterns(terms: Sequence[str], kingdom: str | None = None, columns: Sequence[str] = COLUMNS, unless: Sequence[str] = ()) -> list[HumanPattern]:
    return [HumanPattern(c, t, kingdom, tuple(unless)) for c in columns for t in terms]


HUMAN_PATTERNS: tuple[HumanPattern, ...] = tuple(
    _patterns(
        ("botanic", "jardin botanique", "jardín botánico", "jardim botânico", "botanischer garten", "orto botanico", "arboretum"),
        kingdom="Plantae",
    )
    + _patterns((" zoo ", " zoo,", " zoo.", "zoo animal"), kingdom="Animalia", unless=("not a zoo", "not a captive"))
    + _patterns(("zoo specimen", "zoological garden", "zoologischer garten", "tierpark", "jardin zoologique"), kingdom="Animalia")
    + _patterns(("greenhouse", "glasshouse", "gewächshaus"))
    + _patterns(
        ("cultivated",),
        kingdom="Plantae",
        unless=("not cultivated", "non-cultivated", "uncultivated", "cultivated land", "cultivated field", "cultivated ground", "cultivated soil", "cultivated area"),
    )
    + _patterns(("cultivé",), kingdom="Plantae", unless=("non cultivé", "champ cultivé", "terrain cultivé"))
    + _patterns(("cultivada",), kingdom="Plantae", unless=("no cultivada", "não cultivada"))
    + _patterns(("kultiviert",), kingdom="Plantae", unless=("nicht kultiviert",))
    + _patterns(("planted",), kingdom="Plantae", unless=("not planted", "unplanted"))
    + _patterns(("captive", "in captivity"), kingdom="Animalia", columns=PLACE_COLUMNS)
    # free text remarks discuss captivity a lot without meaning it, so they need a statement of fact
    + _patterns(("captive born", "captive-born", "captive bred", "captive-bred", "captive birth", "captive reared", "captive-reared"), kingdom="Animalia", columns=("occurrenceremarks",))
)


def like_match(value: str | None, term: str) -> bool:
    """The padding lets a term with outer spaces match at the very start and end of the text, as the query does."""
    return value is not None and term in f" {value.lower()} "


def _hits(pattern: HumanPattern, text: str | None, kingdom: str | None) -> bool:
    if pattern.kingdom is not None and (kingdom is None or kingdom.lower() != pattern.kingdom.lower()):
        return False
    return like_match(text, pattern.term) and not any(like_match(text, phrase) for phrase in pattern.unless)


def is_human_held(fields: dict[str, str | None], kingdom: str | None) -> bool:
    """Would the query count this record as held by people? `fields` holds the text per column."""
    return any(_hits(p, fields.get(p.column), kingdom) for p in HUMAN_PATTERNS)
