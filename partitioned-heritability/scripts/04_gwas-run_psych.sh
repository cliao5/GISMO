#! /bin/bash

# Create trait-list before and pass in as argument
# https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

# A bit messy and inelegant implementation... sorry
# $1: output (input for this part) directory
# $2: input cluster list
# $3: input batch subfolder
# $4: output directory
# $5: list of traits to run GWAS individually on
# $6: name for log
# Deciles: 
# BUCKET=/home/rye/toga-human_hg38_reference/partitioned-heritability
# ./04_gwas-run.sh /home/rye/partitioned-heritability/outputs "${BUCKET}"/inputs/GISMO-mis/decile-list GISMO-mis/decile /home/rye/partitioned-heritability/outputs "${BUCKET}"/inputs/sumstats/ukb31063_ldsc/traits-list_1 traits-list_1



# Activate environment 
source activate ldsc

# "${BUCKET}"/outputs
inputDir=$1
# "${BUCKET}"/inputs/GISMO-mis/decile-list
listPath=$2
# GISMO-mis/decile
runName=$3
# /home/rye/partitioned-heritability/outputs
outputDir=$4
# "${BUCKET}"/inputs/sumstats/ukb31063_ldsc/traits-list_1
traits=$5
# traits-list_1
name=$6



# Some hardcoding
# Assuming mounting gs://toga-human_hg38_reference to /home/rye/toga-human_hg38_reference/
baseline="/home/rye/toga-human_hg38_reference/partitioned-heritability/inputs/reference/1000G_Phase3_baselineLD_ldscores/baselineLD."
#baseline="/home/rye/toga-human_hg38_reference/partitioned-heritability/inputs/reference/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."
sumstatPrefix="/home/rye/toga-human_hg38_reference/partitioned-heritability/inputs/sumstats/sumstats_v2"
weights="/home/rye/toga-human_hg38_reference/partitioned-heritability/inputs/reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
frqfile="/home/rye/toga-human_hg38_reference/partitioned-heritability/inputs/reference/1000G_Phase3_frq/1000G.EUR.QC."



# Have to restructure ldsc file organization
while read -r cluster; do
    # Make temp dir for each cluster
    mkdir -p "${outputDir}"/"${runName}"/temp/"${cluster}"/; cd "${outputDir}"/"${runName}"/temp/"${cluster}"/

    # Iterate through chromosomes
    for chr in {1..22}; do
        # Make softlink for .l2.ldscore.gz and .l2.M_5_50 and .annot.gz
        ln -s "${inputDir}"/"${runName}"/"${cluster}"/chr"${chr}"/"${cluster}"."${chr}".l2.ldscore.gz .
        ln -s "${inputDir}"/"${runName}"/"${cluster}"/chr"${chr}"/"${cluster}"."${chr}".l2.M_5_50 .
        ln -s "${inputDir}"/"${runName}"/"${cluster}"/chr"${chr}"/"${cluster}"."${chr}".annot.gz .
    done 
done < "${listPath}"



# Create --ref-ld-chr input flag
d="${baseline}"
while read -r cluster; do
    d+=",${outputDir}"/"${runName}/temp/${cluster}/${cluster}."
done < "${listPath}"



# Function to run GWAS
# 1: sumstatPrefix
# 2: d
# 3: weights
# 4: frqfile
# 5: trait file name
gwas(){
    # Get trait name (cutting name from something like "1468_3.ldsc.imputed_v3.both_sexes.tsv.bgz")
    traitName="$(cut -d'.' -f1 <<< "${5}")"

    # Make output directory for trait
    echo "${traitName}"
    mkdir -p "${traitName}"; cd "${traitName}"

    # Run GWAS with trait
    python /home/rye/ldsc/ldsc.py \
        --h2 "${1}"/"${5}" \
        --print-coefficients \
        --ref-ld-chr "${2}" \
        --w-ld-chr "${3}" \
        --overlap-annot \
        --frqfile-chr "${4}" \
        --out "${traitName}"

    # Switch back
    cd ../

}
export -f  gwas



# Read traits
set -f
IFS=$'\n'
traitList=($(<"${traits}"))



# Parallel iterate through traitList
mkdir -p "${outputDir}"/"${runName}"/GWAS/; cd "${outputDir}"/"${runName}"/GWAS/
start=$(date +%s)

parallel -j 14 gwas "${sumstatPrefix}" "${d}" "${weights}" "${frqfile}" ::: "${traitList[@]}" > "${name}"_log

end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds" >> "${name}"_log



# Clean up
rm -r "${outputDir}"/"${runName}"/temp/
conda deactivate
