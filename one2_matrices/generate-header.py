import os
rootdir = "human_hg38_reference/"

for f in os.listdir(rootdir):
    dir = f'{rootdir}{str(f)}/'
    if os.path.isdir(dir):
        species = [s for s in os.listdir(dir) if os.path.isdir(dir + s)]
        ss = '\t'.join(species)
        print(f"{f}\t{ss}")
