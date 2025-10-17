import pandas as pd


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Get single target per group, merge metadata by taking first occurrence after sorting by specified columns or text merging if requested"
    )
    parser.add_argument(
        "input_file", type=str, help="Path to the CSV file with sequences"
    )
    parser.add_argument(
        "output_file", type=str, help="Output CSV file with single sequence per group"
    )
    parser.add_argument(
        "--sort_by",
        nargs="+",
        type=str,
        default=["temporal_cutoff, pdb_chain_id"],
        help="Columns to sort by before merging",
    )
    parser.add_argument(
        "--merge",
        nargs="+",
        type=str,
        default=None,
        help="Columns to merge after sorting",
    )

    parser.add_argument(
        "--drop",
        nargs="+",
        type=str,
        default=None,
        help="Columns to drop before saving",
    )

    parser.add_argument(
        "--group-key", type=str, default="group_id", help="Column defining groups"
    )

    parser.add_argument(
        "--rename",
        nargs=2,
        action="append",
        metavar=("old_name", "new_name"),
        help="Rename columns after processing, applied in order",
    )

    parser.add_argument(
        "--sort_by_after",
        nargs="+",
        type=str,
        default=["temporal_cutoff", "target_id"],
        help="Columns to sort by after merging",
    )
    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_file = Path(args.output_file)
    targets = pd.read_csv(input_file)

    # Sort and merge targets
    if args.sort_by:
        targets = targets.sort_values(by=args.sort_by)

    targets = targets.groupby(args.group_key).agg(list)
    merge = args.merge if args.merge else []
    for col in targets.columns:
        if col in merge:
            targets[col] = targets[col].apply(lambda x: ";".join(map(str, x)))
        else:
            targets[col] = targets[col].apply(lambda x: x[0] if len(x) > 0 else None)

    targets = targets.reset_index()

    if args.drop:
        targets = targets.drop(columns=args.drop)

    if args.rename:
        # Do the sequential rename instead of one shot to avoid conflicts, renames are applied in order
        for old_name, new_name in args.rename:
            print(f"Renaming column: {old_name} -> {new_name}")
            # drop exising columns to avoid conflicts
            targets = targets.drop(columns=[new_name], errors="ignore")
            targets = targets.rename(columns={old_name: new_name})
    if args.sort_by_after:
        targets = targets.sort_values(by=args.sort_by_after)

    targets.to_csv(
        output_file,
        float_format="%.3f",
        index=False,
    )
