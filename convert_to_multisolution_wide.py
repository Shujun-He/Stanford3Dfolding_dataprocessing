import pandas as pd


def to_wide(result_solution, fill_coordinates=None, select=None):
    """
    Convert long format solution dataframe to wide format and save to CSV, keeping only specified atoms.

    Keeps all atoms if select is None, otherwise filters to atoms starting with the given prefixes.
    """
    # Columns to keep that must be present
    keep = ["group_id", "resname", "resid", "solution_id"]
    columns = set(result_solution.columns)
    missing = set(keep) - columns
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Create selected only version
    coordinate_mask = result_solution.columns.str.endswith(("_x", "_y", "_z"))
    if select is None:
        selected_columns = result_solution.columns[coordinate_mask]
    else:
        selected_columns = result_solution.columns[
            result_solution.columns.str.startswith(tuple(select)) & coordinate_mask
        ]
    if len(selected_columns) == 0:
        raise ValueError(f"No columns found for selected atoms: {select}")

    # # Non-coordinate columns
    # other_columns = list(result_solution.columns[~coordinate_mask])
    # result_solution_selected = result_solution[other_columns + selected_columns]

    # Use atom coordinates and solution_id as multiindex for columns
    result_solution_selected = result_solution[keep + list(selected_columns)]

    result_solution_selected = result_solution_selected.set_index(keep).unstack(
        "solution_id"
    )
    result_solution_selected = result_solution_selected.sort_index(level=[1, 0], axis=1)
    # Fill missing coordinates if requested
    if fill_coordinates is not None:
        # We need to convert to string first to apply float formatting, doing it here
        # so the columns are already aligned and filled with nans
        coordinate_mask = result_solution_selected.columns.get_level_values(
            0
        ).str.endswith(("_x", "_y", "_z"))
        coordinate_columns = result_solution_selected.columns[coordinate_mask]
        result_solution_selected[coordinate_columns] = result_solution_selected[
            coordinate_columns
        ].map(lambda x: f"{x:.3f}" if pd.notnull(x) else fill_coordinates)

    # Map multilevel column names to {x,y,z}_index
    result_solution_selected.columns = [
        f"{i}_{j}" if j != "" else f"{i}" for i, j in result_solution_selected.columns
    ]


    # Sort by resid within each group
    result_solution_selected = result_solution_selected.sort_values(
        by=["group_id", "resid"]
    )

    # Reconstruct ID
    result_solution_selected = result_solution_selected.reset_index()
    result_solution_selected["ID"] = (
        result_solution_selected["group_id"].astype(str)
        + "_"
        + (result_solution_selected["resid"]).astype(str)
    )

    # Remove C1'_ prefix from column names
    result_solution_selected.columns = [
        col.replace("C1'_", "") for col in result_solution_selected.columns
    ]

    # Move ID to front
    cols = list(result_solution_selected.columns)
    cols.remove("ID")
    cols = ["ID"] + cols
    result_solution_selected = result_solution_selected[cols]
    return result_solution_selected


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file", type=str, help="Path to the input CSV file in long format"
    )
    parser.add_argument("output_file", type=str, help="Output CSV file in wide format")

    parser.add_argument(
        "--select",
        nargs="+",
        type=str,
        default=None,
        help="Select atoms starting with this prefix",
    )
    parser.add_argument(
        "--fill-coordinates",
        type=str,
        default=None,
        help="Fill missing values with this string",
    )

    parser.add_argument(
        "--add-column",
        nargs=2,
        action="append",
        metavar=("column_name", "value"),
        help="Add a new column with the given name and value",
    )
    parser.add_argument(
        "--drop-column", nargs="+", type=str, help="Drop these columns before saving"
    )

    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_file = Path(args.output_file)
    result_solution = pd.read_csv(input_file)
    result_solution_selected = to_wide(
        result_solution, fill_coordinates=args.fill_coordinates, select=args.select
    )

    if args.add_column:
        for col_name, value in args.add_column:
            result_solution_selected[col_name] = value
    if args.drop_column:
        result_solution_selected = result_solution_selected.drop(
            columns=args.drop_column, errors="ignore"
        )
    result_solution_selected.to_csv(
        output_file,
        float_format="%.3f",
        index=False,
    )
