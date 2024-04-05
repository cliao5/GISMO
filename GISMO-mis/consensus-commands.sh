# inputs/MultipleCodonAlignments/ fasta files hosted at: https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments/
# Outputs moved to *GenerateScore/inputs/consensus_v2_missense-counts.gz* and *GenerateScore/inputs/consensus_v2_synonymous-counts.gz*
python3 multiple-codon-alignments_count_consensus.py \
    --fasta_list inputs/multiple-codon-alignments_fasta-files \
    --fasta_prefix inputs/MultipleCodonAlignments/ \
    --out ../outputs/consensus_v2 \
    > ../outputs/consensus_v2_log