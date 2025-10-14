from collections import defaultdict


class PartialSequenceMapping(dict):
    """A dictionary that allows partial key matching.

    Retrieves the value for a key if the key is a substring of any existing key or contains existing keys
    with a minimum match ratio.
    """

    def __init__(self, *args, min_match_ratio=0.8, **kwargs):
        super().__init__(*args, **kwargs)
        self.min_match_ratio = min_match_ratio
        if self.min_match_ratio >= 1:
            self.get_match_item = self._get_exact_item
            self.__contains__ = super().__contains__

    def _get_exact_item(self, key):
        value = super().__getitem__(key)
        return key, value

    def get_match_item(self, key):
        # First try exact match
        try:
            return self._get_exact_item(key)
        except KeyError:
            # Then try partial match with minimum ratio

            for k in self.keys():
                if (key in k and len(key) / len(k) >= self.min_match_ratio) or (
                    k in key and len(k) / len(key) >= self.min_match_ratio
                ):
                    return k, super().__getitem__(k)
            raise KeyError(key)

    def __contains__(self, key):
        if super().__contains__(key):
            return True
        for k in self.keys():
            if (key in k and len(key) / len(k) >= self.min_match_ratio) or (
                k in key and len(k) / len(key) >= self.min_match_ratio
            ):
                return True
        return False


def groupby_sequence(data, min_match_ratio=0.9):
    grouped_data = {}
    seq_to_id = PartialSequenceMapping(min_match_ratio=min_match_ratio)
    for i in range(len(data["sequence"])):
        seq = data["sequence"][i]
        xyz = data["xyz"][i]
        cif_file = data["data_cif_files"][i]
        pdb_id = cif_file.split("/")[-1].split(".")[0]
        if seq not in seq_to_id:
            seq_to_id[seq] = pdb_id
            ref_id = pdb_id
            grouped_data[ref_id] = defaultdict(list)
        else:
            matched_seq, ref_id = seq_to_id.get_match_item(seq)
            # We are using partial matching, make sure the new seq is recorded, keep the original ref_id
            seq_to_id[matched_seq] = ref_id

        grouped_data[ref_id]["sequence"].append(seq)
        grouped_data[ref_id]["data_cif_files"].append(cif_file)
        grouped_data[ref_id]["pdb_ids"].append(pdb_id)
        grouped_data[ref_id]["xyz"].append(xyz)

    return grouped_data


if __name__ == "__main__":
    import argparse
    import pickle
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        type=str,
        help="Pickle file with data to be grouped by sequence",
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Output directory for the grouped data to be saved",
    )
    parser.add_argument(
        "--min_match_ratio",
        type=float,
        default=1,
        help="Minimum match ratio for partial sequence matching",
    )
    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)
    min_match_ratio = args.min_match_ratio
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_file, "rb") as f:
        data = pickle.load(f)

    grouped_data = groupby_sequence(data, min_match_ratio=min_match_ratio)

    for group_id, entries in grouped_data.items():
        print(f"Group {group_id} has {len(entries['pdb_ids'])} entries.")
        output_file = output_dir / f"{group_id}_grouped.pkl"
        with open(output_file, "wb") as f:
            pickle.dump(entries, f)
