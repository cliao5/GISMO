import argparse
import pandas as pd

# Set up argparse
parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, required=True, help = "path to where all directories with results live")
parser.add_argument('--out', type = str, default = "merged_result", help = "path to output file (default: merged_result)")
parser.add_argument('--traits', type = str, required = True, help = "path to traits list (with PhenotypeC column)")
parser.add_argument('--groups', type = str, required = True, help = "path to groups list (the order in which -h2 was calculated)")
args = parser.parse_args()

# Read in traits table and get traits
traits_table = pd.read_table(args.traits)
traits = list(traits_table.PhenotypeC)

# Read in groups list
groups_table = pd.read_table(args.groups, header=None)
groups = list(groups_table.iloc[:, 0])

final = pd.DataFrame()
# Iterate through traits and construct merged table
for trait in traits:
    print(trait)
    d = pd.read_table(f"{args.indir}/{trait}/{trait}.results", header = 0)
    d = d.loc[d["Category"].isin([f"L2_{i}" for i in range(1, len(groups) + 1)])]
    d['group'] = groups
    d['trait'] = trait
    final = pd.concat([final, d])

final.to_csv(args.out, sep = "\t", index = False)
