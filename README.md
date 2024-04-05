# GISMO

Gene Identity Score of Mammalian Orthologs

## GenerateScore
- Generating GISMO-mis and GISMO metrics. 
    #### inputs/
    - *consensus_v2_missense-counts.gz* and *consensus_v2_synonymous-counts.gz* output from scripts in **GISMO-mis**
    - *unmerged-species_combined_matrix_2023-07-26.tsv* output from scripts in **one2_matrices**


## GISMO-mis/
- Scripts to generate input files for GISMO-mis 
- Instead of comparing to human (">REFERENCE" in fasta), generate a consensus reference sequence by codon for each gene. Note that this permits multiple reference codons at a position (in the event of ties)
- If there is a deletion ("-") anywhere in either reference or query sequence for a gene, skip over it (not counted towards mis or syn totals)
- For any "N" bp, consider all possible substitutions. To be the most conservative, if it can be synonymous, mark is as synonymous. In other words, only classify if missense if all possibilities are missense.
- For synonymous and missense scoring in the case of multiple reference codons, again taking the most conservative route by only determining missense if all reference codons are only missense against queried codon.
    ### inputs/
    - All used *ENST[x].[gene].fasta.gz* files hosted at: https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments/
    ### outputs/
    - Primary outputs moved to: *GenerateScore/inputs/consensus_v2_missense-counts.gz* and *GenerateScore/inputs/consensus_v2_synonymous-counts.gz*

   


## one2_matrices/
- Generating GISMO one2 matrix (primary output: *GenerateScore/inputs/unmerged-species_combined_matrix_2023-07-26.tsv*)
- Input files hosted at: https://genome.senckenberg.de/download/TOGA/human_hg38_reference/



