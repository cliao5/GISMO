#! /bin/bash

# First, generate cluster list and create gene-set/coord files
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial

# $1: input directory
# $2: input cluster list
# $3: input batch subfolder
# $4: output directory
# Example: 
# BUCKET=/home/rye/toga-human_hg38_reference/partitioned-heritability
# ./01_ldsc-run.sh "${BUCKET}"/inputs "${BUCKET}"/inputs/GISMO-mis/decile-list GISMO-mis/decile /home/rye/partitioned-heritability/outputs


# Activate environment
source activate ldsc

# Function to create annotation files and run ldsc
annot_ldsc(){
    # Make output dir for each chr
    mkdir -p chr"${6}"/; cd chr"${6}"/

    # Run make_annot.py to create annotation file
    python /home/rye/ldsc/make_annot.py \
        --gene-set-file "${1}"/"${3}"/"${5}"/gene-set-file \
        --gene-coord-file "${1}"/reference/ENSG_coord.txt \
        --windowsize 100000 \
        --bimfile "${1}"/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC."${6}".bim \
        --annot-file "${5}"."${6}".annot.gz

    # Run ldsc
    python /home/rye/ldsc/ldsc.py \
        --l2 \
        --bfile "${1}"/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC."${6}" \
        --ld-wind-cm 1 \
        --annot "${5}"."${6}".annot.gz \
        --thin-annot \
        --print-snps "${1}"/reference/list.txt \
        --out "${5}"."${6}"
        #--print-snps "${1}"/reference/hapmap3_snps/hm."${6}".snp
        

    # Switch back
    cd ../
}
export -f  annot_ldsc



inputDir=$1
listPath=$2
runName=$3
outputDir=$4

# Iterate through deciles
while IFS= read -r cluster; do
    start=$(date +%s)
    # Make output dir for each cluster
    mkdir -p "${outputDir}"/"${runName}"/"${cluster}"/; cd "${outputDir}"/"${runName}"/"${cluster}"/

    # Iterate parallel through chromosomes and call annot_ldsc
    parallel -j 15 annot_ldsc "${inputDir}" "${listPath}" "${runName}" "${outputDir}" "${cluster}" ::: {1..22} > "${cluster}"_annot_ldsc_log
    end=$(date +%s)
    echo "Elapsed Time: $(($end-$start)) seconds" >> "${cluster}"_annot_ldsc_log
done < "${listPath}"

conda deactivate