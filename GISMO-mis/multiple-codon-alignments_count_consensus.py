"""
Fasta files downloaded from: https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments/
Using: wget -r --no-parent -A '*fasta.gz' --no-check-certificate -e robots=off https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments/

Simple script to generate table counting number of synonymous/missense variants between each species and reference
Repeated for all 17439 variant fasta's


First, generate consensus reference sequence by codon (3 bp). Includes all species as unique votes (including those with duplicated names), as well as human reference (>REFERENCE)
Codon ties are permitted.

For synonymous and missense scoring:
If there is a gaps ("-") anywhere in either reference or query sequence, skip over it (don't count towards mis or syn)
In the case of multiple reference codons at a single position, none must be synonymous to be labelled as missense (otherwise, synonymous)
In the case of low confidence ("N"), consider all possible substitutions. Once again, will only be labelled as missense if no possible synonymous

Duplicated species names are averaged 
"""

# Set up
import argparse
import gzip
import time

## Codon mappings
codons = {'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
         'TTC': 'F', 'TTT': 'F',
         'TTA': 'L', 'TTG': 'L',
         'TAC': 'Y', 'TAT': 'Y',
         'TAA': '*', 'TAG': '*', 
         'TGC': 'C', 'TGT': 'C',
         'TGA': '*',
         'TGG': 'W',
         'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
         'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
         'CAC': 'H', 'CAT': 'H', 
         'CAA': 'Q', 'CAG': 'Q',
         'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
         'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 
         'ATG': 'M',
         'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
         'AAC': 'N', 'AAT': 'N',
         'AAA': 'K', 'AAG': 'K',
         'AGC': 'S', 'AGT': 'S',
         'AGA': 'R', 'AGG': 'R',
         'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
         'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
         'GAC': 'D', 'GAT': 'D', 
         'GAA': 'E', 'GAG': 'E',
         'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'
}

# Dictionary will store species:count, where count is a 17439-long list of ints (for each variant represented by each fasta)
# Each entry in order of processed variant Fasta
missense_counts = {}
synonymous_counts = {}
synonymous_counts_with_n = {}

# Also track all species and variants that we've come across
species_set = set()
genes_list = []
# Gene:Reference(s)
reference_dict = {}

"""
Function to generate consensus reference(s) given full list of lines from FASTA file
"""
def get_reference(lines):
    # Grab only the sequences
    seqs = [l for l in lines if ">" not in l]

    # Check that we have the right number
    assert(len(seqs) == len(lines)/2)
    # Check that sequences are the same and correct length
    l = len(seqs[0])
    assert(l%3 == 0)
    for s in seqs:
        assert(len(s) == l)

    # Iterate through codons (3's) to construct consensus sequence (each element is a list of most voted codon(s))
    sequences = []
    for i in range(l//3):
        counts = {}
        # Iterate through species sequences
        for s in seqs:
            codon = s[3*i:3*(i + 1)]
            if codon not in counts:
                counts[codon] = 1
            else:
                counts[codon] += 1

        # Find maximum 
        cs = [key for key, value in counts.items() if value == max(counts.values())]
        sequences.append(cs)


    assert(len(sequences) == l//3)
    for seq in sequences:
        assert(len(seq) > 0)

    return sequences

"""
Recursive function to consider all possibilities of replacing N's
Initial input before recursive calls is a list of a single codon (string) with possibly "N" characters to replace
"""
def replace_n(s):
    if any(['N' in codon for codon in s]):
        # Repeat for each first instance of 'N'
        new_codons = set()
        for codon in s:
            if 'N' not in codon:
                new_codons.add(codon)
                continue
            for ind in [i for i, nt in enumerate(codon) if nt == 'N']:
                new_codons.update([codon[:ind] + a + codon[ind+1:] for a in ['A', 'T', 'G', 'C']])
        return replace_n(list(new_codons))
    else:
        return s


"""
Query codon mappings dictionary (if possible) and return whether resulting AA match
"""
def is_synonymous(r: str, s: str):
    if (r in codons) and (s in codons):
        return codons[r] == codons[s]
    else:
        return False



"""
Function to compare reference to input sequence

NOTE: reference_seq should be a list of each possible codon at each position

Iterate through both by codon
Skip if equivalent codons OR if there are any deletions ("-")
Try all substitutions for "N" (conservatively) as synonymous if match

Return synonymous, missense, and synonymous with N counts
"""
def compare_to_reference(reference_seq: str, seq: str) -> int:
    # Check both sequences are the same length and multiple of 3
    assert(len(reference_seq) == len(seq)//3)
    assert(len(seq)%3 == 0)

    final_syn_count = 0
    final_mis_count = 0
    final_syn_with_n = 0
    # Iterate through sequence 3-characters at a time (each element of reference)
    for i in range(len(reference_seq)):
        r = reference_seq[i]
        s = seq[3*i:3*i+3]

        # For each possible codon in consensus reference
        syn_count = 0
        mis_count = 0
        syn_with_n = 0
        for c in r:
            # Nothing interesting if sequences are the same
            if c == s:
                continue

            # Skip over anything with a deletion as well
            if ("-" in c) or ("-" in s):
                continue
            
            # Account for N codons conservatively (if it can be synonymous, mark as synonymous)
            # ASSUME NO N's IN REFERENCE OTHER THAN NNN
            if (s == 'NNN' or c == 'NNN'):
                syn_count += 1
                syn_with_n += 1
                continue

            is_syn = False
            if ('N' in s):
                new_codons = replace_n([s])
                for new_codon in new_codons:
                    if (c == new_codon) or is_synonymous(c, new_codon):
                        syn_count += 1
                        syn_with_n += 1
                        is_syn = True
                        break
                if is_syn:
                    continue
                else:
                    mis_count += 1
                    continue

            # Automatically count for missense if no N's and either not in codon dict (excluding cases with '-')
            if is_synonymous(c, s):
                syn_count += 1
                continue
            else:
                mis_count += 1
                continue
    
        
        final_syn_count += (syn_count > 0)
        final_mis_count += (mis_count == len(r))
        final_syn_with_n += (syn_with_n > 0)
        


    return final_syn_count, final_mis_count, final_syn_with_n


"""
Function to count genes in a given fasta file. 
"""
def count_genes(fasta_file: str):
    # Read lines from fasta file
    with gzip.open(fasta_file, 'rt') as f:
        lines = f.read().splitlines()

    # Sanity check that variants are 2 lines at a time
    assert(len(lines)%2 == 0)

    # Go through the entire file (all lines) and generate reference (consensus majority)
    reference_seq = get_reference(lines)


    
    # Keep only ENST[#] for gene
    # From original, full ENST[#].[name].[#]
    gene = ""

    # Dictionaries will store species:count, where count is a list of ints (for each species sharing the same name represented in each fasta)
    # Keyed by unique species names
    species_ndup = {}
    species_syn_counts = {}
    species_syn_with_n_counts = {}
    species_mis_counts = {}
    species_set_var = set()
    # Iterate through all species now (INCLUDING reference)  
    for i in range(len(lines)//2):
         
        assert(lines[2*i][0] == ">")
        # Grab gene if not already grabbed, or check that it's consistent
        if "REFERENCE" not in lines[2*i]:
            g = lines[2*i].split()[1].split(".")[0] 
            if gene == "":
                gene = g
            else:
                # Check all have the same gene name prefix (up the the ENST[#])
                assert(gene == g)
        
        # Check all genes have the same length as reference for direct comparison
        assert(len(lines[2*i+1])//3 == len(reference_seq))
 
        # Grab sequence and compare to reference
        seq = lines[2*i+1]

        # Count regardless of whether greater than 50% indels (no more marking as NA)
        syn_count, mis_count, syn_with_n = compare_to_reference(reference_seq, seq)

        # Grab species and add to dictionaries if not already, or increment by 1 if duplicate
        if "REFERENCE" in lines[2*i]:
            species = "REFERENCE"
        else:
            species = lines[2*i].split()[0].split(">vs_")[1]

        if species not in species_ndup:
            species_ndup[species] = 1
            species_syn_counts[species] = [syn_count]
            species_syn_with_n_counts[species] = [syn_with_n]
            species_mis_counts[species] = [mis_count]
        else:
            species_ndup[species] += 1
            species_syn_counts[species].append(syn_count)
            species_syn_with_n_counts[species].append(syn_with_n)
            species_mis_counts[species].append(mis_count)

        species_set_var.add(species)



    # Now that we're done with all species for current gene, add to dictionary
    ## Go species by species
    for species in species_set_var:
        # If there are species with > 50% '-' (dropped, represented as None)
        if (None in species_syn_counts[species]) or (None in species_syn_with_n_counts[species]) or (None in species_mis_counts[species]):
            # Remove Nones
            species_syn_avg_complete = [s for s in species_syn_counts[species] if s != None] 
            species_syn_with_n_avg_complete = [s for s in species_syn_with_n_counts[species] if s != None] 
            species_mis_avg_complete = [s for s in species_mis_counts[species] if s != None]
            assert(len(species_syn_avg_complete) == len(species_syn_with_n_avg_complete) and len(species_syn_avg_complete) == len(species_mis_avg_complete))

            # If there was only one
            if len(species_syn_avg_complete) == 0:
                species_syn_avg = None
                species_syn_with_n_avg = None
                species_mis_avg = None
            # Otherwise, average as usual
            else:
                species_syn_avg = sum(species_syn_avg_complete) / len(species_syn_avg_complete)
                species_syn_with_n_avg = sum(species_syn_with_n_avg_complete) / len(species_syn_with_n_avg_complete)
                species_mis_avg = sum(species_mis_avg_complete) / len(species_mis_avg_complete)

        # Otherwise, no None entries
        else:
            ## Take average counts for each species(handling duplicates)
            species_syn_avg = sum(species_syn_counts[species]) / species_ndup[species]
            species_syn_with_n_avg = sum(species_syn_with_n_counts[species]) / species_ndup[species]
            species_mis_avg = sum(species_mis_counts[species]) / species_ndup[species]
        
        # Now that we have values, create entries in dictionary
        ## If we haven't seen species before 
        if (species not in synonymous_counts) or (species not in synonymous_counts_with_n) or (species not in missense_counts):
            assert((species not in synonymous_counts))
            assert((species not in synonymous_counts_with_n))
            assert((species not in missense_counts))
            synonymous_counts[species] = []
            synonymous_counts_with_n[species] = []
            missense_counts[species] = []
        
        # We need to pad with same number of genes that we've processed so far (so all species are the same length)
        ## Fill with Nones
        assert(len(synonymous_counts[species]) == len(synonymous_counts_with_n[species]) and len(synonymous_counts[species]) == len(missense_counts[species]))

        length = len(missense_counts[species])
        assert(length <= len(genes_list))
        for i in range(len(genes_list) - length):
            synonymous_counts[species].append(None)
            synonymous_counts_with_n[species].append(None)
            missense_counts[species].append(None)

        assert(len(synonymous_counts[species]) == len(synonymous_counts_with_n[species]) and len(synonymous_counts[species]) == len(missense_counts[species]))
        assert(len(synonymous_counts[species]) == len(genes_list))

        synonymous_counts[species].append(species_syn_avg)
        synonymous_counts_with_n[species].append(species_syn_with_n_avg)
        missense_counts[species].append(species_mis_avg)

        
    genes_list.append(gene)
    #assert(gene not in reference_dict)
    #reference_dict[gene] = reference_seq

    species_set.update(species_set_var)


def main(args):
    start = time.time()
    # Read in list of fasta files
    with open(args.fasta_list, "r") as f:
        fasta_files = f.read().splitlines()
    fasta_files = [args.fasta_prefix + f for f in fasta_files]
    
    # Process all fasta files
    for fasta_file in fasta_files:
        count_genes(fasta_file)
        print(f"{len(genes_list)}/{len(fasta_files)} variants processed ({time.time() - start} elapsed)")
    

    # Check for no duplicated genes
    assert(len(genes_list) == len(set(genes_list)))

    # Finish padding all species entries with Nones
    for species in species_set:
        length = len(missense_counts[species])
        assert(length <= len(genes_list))
        for i in range(len(genes_list) - length):
            synonymous_counts[species].append(None)
            synonymous_counts_with_n[species].append(None)
            missense_counts[species].append(None)
        assert(len(synonymous_counts[species]) == len(synonymous_counts_with_n[species]) and len(synonymous_counts[species]) == len(missense_counts[species]) and len(synonymous_counts[species]) == len(genes_list))
        assert(len(synonymous_counts[species]) == len(genes_list))

    species_list = list(species_set)
    species_list.sort()


    """
    reference_dict manipulation
    # Write consensus reference sequences in fasta format
    assert(len(fasta_files) == len(reference_list))
    with open(args.out+"_consensus-references", "w") as f:
        for i in range(len(fasta_files)):
            f.write(f'>{fasta_files[i]}\n{reference_list[i]}\n')
    """


    # Write synonymous-counts
    with open(args.out+"_synonymous-counts", "w") as f:
        line = "\t".join(genes_list) 
        f.write(f'species\t{line}\n') # Header
        for species in species_list:
            line = "\t".join([str(i) for i in synonymous_counts[species]])
            f.write(f'{species}\t{line}\n')

    # Write synonymous-counts-with-n
    with open(args.out+"_synonymous-counts-with-n", "w") as f:
        line = "\t".join([str(i) for i in genes_list])
        f.write(f'species\t{line}\n') # Header
        for species in species_list:
            line = "\t".join([str(i) for i in synonymous_counts_with_n[species]])
            f.write(f'{species}\t{line}\n')

    # Write missense-counts
    with open(args.out+"_missense-counts", "w") as f:
        line = "\t".join(genes_list)
        f.write(f'species\t{line}\n') # Header
        for species in species_list:
            line = "\t".join([str(i) for i in missense_counts[species]])
            f.write(f'{species}\t{line}\n')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta_list",
        help = "Path to list of fasta files (one per line)",
        type = str,
        required = True
    )
    parser.add_argument(
        "--fasta_prefix",
        help = "File prefix to fasta files",
        type = str,
        required = True
    )
    parser.add_argument(
        "--out",
        help="Output path prefix for final tables (csvs)",
        type=str,
        required = True
    )
    
    args = parser.parse_args()

    main(args)

