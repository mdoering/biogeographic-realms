import polars as pl


def presence(rows):
    return pl.DataFrame(rows, schema={"cell": pl.String, "species_id": pl.String}, orient="row")


def grid(width, height):
    """Cells named "x,y", each touching the four around it."""
    cells = [f"{x},{y}" for x in range(width) for y in range(height)]
    adjacency = {c: set() for c in cells}
    for x in range(width):
        for y in range(height):
            for dx, dy in ((1, 0), (0, 1)):
                if x + dx < width and y + dy < height:
                    adjacency[f"{x},{y}"].add(f"{x + dx},{y + dy}")
                    adjacency[f"{x + dx},{y + dy}"].add(f"{x},{y}")
    return cells, adjacency
