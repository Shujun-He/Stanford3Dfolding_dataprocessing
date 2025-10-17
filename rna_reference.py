from collections import UserDict
import json
from pathlib import Path

# Manually curated mapping of noncanonical residue names to RNA.
# TODO: Consider assigning based on ability to pair
# TODO: review based on parent nucleotide mapping from NKDB; https://nakb.org/modifiednt.html
# fmt: off
modified_to_unmodified = {   'A'  :'  A',   'C':'  C',   'G':'  G',   'U':'  U',\
                '5BU':'  U', 'OMC':'  C', '5MC':'  C', 'CCC':'  C', ' DC':'  C', \
                'CBR':'  C', 'CBV':'  C', 'CB2':'  C', '2MG':'  G', \
                'H2U':'  U', 'PSU':'  U', '  U':'  U', '5MU':'  U', '2MU':'  U', 'OMU':'  U', \
                'OMG':'  G', '7MG':'  G', '1MG':'  G', 'GTP':'  G', 'AMP':'  A', 'MIA':'  A', ' YG':'  G', \
                'M2G':'  G', 'YYG':'  G', ' DG':'  G', 'G46':'  G', ' IC':'  C', ' IG':'  G',  \
                'ZMP':'ZMP', 'YYG':'  G', '2MG':'  G', 'H2U':'  U', 'AG9':'  C', ' IU':'  U', '3TD':'  U', \
                'A2M':'  A', '1MA':'  A', 'MA6':'  A', 'QUO':'  G', '6MZ':'  A', \
                'C4J':'  C', '4OC':'  C', 'G7M':'  G', 'T6A':'  A', 'AET':'  A', 'I4U':'  U', 'UR3':'  U', \
                'P7G':'  G', 'B9B':'  G', 'B8H':'  U', 'E6G':'  G', 'B8W':'  G', 'B8N':'  U', '4SU':'  U', \
                'LV2':'  C', '4AC':'  C', 'UY4':'  A', 'I2T':'  C', '7SN':'  G', 'SUR':'  U', '7S3':'  G', \
                'LHH':'  C', 'FHU':'  U', 'B9H':'  C', 'M1Y':'  U', 'B8Q':'  C',\
                'M7A':'  A', 'B8K':'  G', '2PR':'  G', 'LCG':'  G', 'UFT':'  U', 'CFZ':'  C', \
                    '3AU':'  U', '9QV':'  U', 
                'CFL':'  C', 'T2T':'  T', 'N'  :'  A', 'I'  :'  G', 'GRB':'  G', 'E3C':'  C', \
                'MMX':'  C', '1W5':'  C', '8AZ':'  G', 'B8T':'  C', 'UY1':'  U', '75B':'  U', \
                '4DU':'  A', '5HM':'  C', '6FC':'  C', 'E7G':'  G', 'MHG':'  G', 'DU' :'  U', \
                '56B':'  G', 'P5P':'  A', 'UMS':'  U', 'PYO':'  U', 'JMC':'  C', 'ZJS':'  A', \
                '6IA':'  A', 'CM0':'  U', '2MA':'  A', 'RSP':'  U', 'UD5':'  U', 'MUM':'  U', \
                'IU' :'  U', '12A':'  A', '70U':'  U', 'U8U':'  U',  'YG':'  G', 'BRU':'  U', \
                'ATP':'  A', 'CTP':'  C', 'UTP':'  U', '5IU':'  I', 'GDP':'  G', '5IC':'  C', \
                # These are DNA residues that we map to RNA, so it will effectively extract DNA too if these are included
                # 'DA' :'  A',  'DC':'  C',  'DG':'  G',  'DT':'  U', ' DU':'  U'
            }
# fmt: on

## Atom names for RNA residues
phosphate = ["P", "OP1", "OP2", "OP3"]  # OP3 will be missing in polymers
sugar = ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'", "O5'", "O2'", "O3'"]

# Get reference from the NAKB modified nucleotide database
DATA_DIR = Path(__file__).parent / "data"


def _load_nakb_mapping():
    """Load the NAKB modified nucleotide data and create a mapping of modified residue names to standard bases."""
    reference_file = DATA_DIR / "nakb_modified_nt_20251014.json"
    with open(reference_file) as f:
        nakb_data = json.load(f)
        mapping = {}
        for key, val in nakb_data.items():
            if "standard_base" in val:
                try:
                    standard_base = val["standard_base"][0]
                except IndexError:
                    continue
                # RNA only
                if standard_base in ("A", "C", "G", "U"):
                    mapping[key] = standard_base
    return mapping


modified_to_unmodified_nakb = _load_nakb_mapping()

# fmt: off
rna_atom_groups = {
    "A": {
        "all": phosphate + sugar + \
               ["N9", "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N6"]
    },
    "U": {
        "all": phosphate + sugar + \
               ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
    },
    "G": {
        "all": phosphate + sugar + \
               ["N9", "N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", ]
    },
    "C": {
        "all": phosphate + sugar + \
                ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
    },
    "N": {
        "all": phosphate + sugar 
    }
}
