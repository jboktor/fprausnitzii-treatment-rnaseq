#!/bin/bash
#SBATCH --job-name=ParkQC
#SBATCH --time=12:00:00
#SBATCH --mem=16G

# Load conda
source ~/miniconda3/etc/profile.d/conda.sh

# Activate environment
conda activate MultiQCEnv

fastqc /groups/Ismagilov_Group/AW/Parkinsons/20250415_DataForAna/AllFastq/*.fastq.gz -o /groups/Ismagilov_Group/AW/Parkinsons/20250415_DataForAna/AllFastq/FastQC_Results/
multiqc /groups/Ismagilov_Group/AW/Parkinsons/20250415_DataForAna/AllFastq/FastQC_Results/* -o /groups/Ismagilov_Group/AW/Parkinsons/20250415_DataForAna/MultiQC/


#############################################################################################

#!/bin/bash
#SBATCH --time=24:00:00   # walltime
#SBATCH --nodes=2   # number of nodes
#SBATCH --mem-per-cpu=64G   # memory per CPU core
#SBATCH -J "LI-C1"   # job name CHANGED FOR EACH SAMPLE
#SBATCH --mail-user=awinnett@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Variables
sample="LI-C1" # CHANGED FOR EACH SAMPLE
fastqloc="/groups/Ismagilov_Group/AW/Parkinsons/AllFastq/"
outloc="/groups/Ismagilov_Group/AW/Parkinsons/"

mkdir -p "${outloc}/${sample}"

# STAR Alignment
STAR --genomeDir /groups/Ismagilov_Group/shared_resources/database_builds/STAR_gencode_GRCm39_genomic \
     --readFilesIn "${fastqloc}/${sample}_R1.fastq.gz" "${fastqloc}/${sample}_R2.fastq.gz" \
     --readFilesCommand zcat \
     --genomeLoad NoSharedMemory \
     --outFileNamePrefix "${outloc}/${batch}/${sample}/${sample}."

# Convert SAM to sorted BAM
samtools view -bS "${outloc}/${sample}/${sample}.Aligned.out.sam" | \
samtools sort -@ 8 -o "${outloc}/${sample}/${sample}.sorted.out.bam"

# Index the BAM file
samtools index "${outloc}/${sample}/${sample}.sorted.out.bam" "${outloc}/${batch}/${sample}/${sample}.sorted.bam.bai"

sbatch "${sample}_run_qualimap-featureCounts.sh"

#############################################################################################

#!/bin/bash
#SBATCH --time=6:00:00   # walltime
#SBATCH --nodes=2   # number of nodes
#SBATCH --mem-per-cpu=32G   # memory per CPU core
#SBATCH -J "LI-C1"   # job name CHANGED FOR EACH SAMPLE

# Specify sample and batch variables
sample="LI-C1" # CHANGED FOR EACH SAMPLE
outloc="/groups/Ismagilov_Group/AW/Parkinsons/"

# Qualimap RNA-seq analysis
qualimap rnaseq \
    -outdir "${outloc}/${sample}/" \
    -a proportional \
    -bam "${outloc}/${sample}/${sample}.sorted.out.bam" \
    -p non-strand-specific \
    -gtf "/groups/Ismagilov_Group/shared_resources/genomes/gencode_GRCm39_genomic/gencode.vM36.primary_assembly.annotation.gtf" \
    --java-mem-size=8G

# FeatureCounts
featureCounts \
    -p -O -J -T 8 \
    -g gene_id \
    -a "/groups/Ismagilov_Group/shared_resources/genomes/gencode_GRCm39_genomic/gencode.vM36.primary_assembly.annotation.gtf" \
    -o "${outloc}/${sample}/${sample}.counts.txt" \
    "${outloc}/${sample}/${sample}.sorted.out.bam"

#############################################################################################

#!/bin/bash

# Base directory where the SB# folders are located
BASE_DIR="/groups/Ismagilov_Group/AW/Parkinsons"  # Replace with the path to your base directory

# Directory where the quant.sf files will be copied
DEST_DIR="${BASE_DIR}/AllFeatureCounts"

# Create the destination directory if it doesn't exist
mkdir -p "${DEST_DIR}"

# Counter for the number of files copied
file_count=0

# Loop through the SB# directories
for avw_dir in *-*; do
    avw_name=$(basename "${avw_dir}")

    # Construct the path to the quant.sf file
    quant_file="${avw_dir}/${avw_name}.counts.txt"
    echo $quant_file

    # Check if the file exists
    if [[ -f "${quant_file}" ]]; then

        # Copy the file to the destination with the new name
        cp "${quant_file}" "${DEST_DIR}/${avw_name}.counts.txt"

        # Increment the file counter
        ((file_count++))
    else
        # error if file is not found
        echo "Error: counts.txt.summary file not found in ${avw_dir}"
    fi
done


#############################################################################################

#!/bin/bash

# File name
input_file="data.txt"
temp_file="temp.txt"

# Read each line from the file
while IFS=$'\t' read -r AVWName a b c d e f g h i; do
    echo $AVWName
    # Skip the header line
    if [ "$AVWName" == "AVWName" ]; then
        echo -e "AVWName\tR1FileSize\tR2FileSize\tSTAR_TotalReads\tSTAR_PercHuman\tSTAR_PercExonic\tSTAR_PercIntronic\tSTAR_PercIntergenic\tSTAR_PercOverlapExon\tFeatureCounts_Sum" > "$temp_file"
        continue
    fi
    counts_file="../$AVWName/${AVWName}.counts.txt"
    file_r1="${AVWName}_R1.fastq.gz"
    file_r2="${AVWName}_R2.fastq.gz"
    echo $file_r1 $file_r2
    if [ -f "$counts_file" ]; then
        sum_counts=$(awk '{sum += $NF} END {print sum}' "$counts_file")
    else
        sum_counts="N/A"
    fi

    new_r1_size=$( [ -f "$file_r1" ] && du -sh "$file_r1" | cut -f1 || echo "N/A" )
    new_r2_size=$( [ -f "$file_r2" ] && du -sh "$file_r2" | cut -f1 || echo "N/A" )
    filename=.Log.final.out
    STARfile=../$AVWName/$AVWName$filename
    num_reads=$(grep "Number of input reads" ../$AVWName/$AVWName$filename | awk '{print $NF}' || echo "N/A")
    STARmapping_rate=$(grep "Uniquely mapped reads %" ../$AVWName/$AVWName$filename | awk '{print $NF}' || echo "N/A")
    QUALIMAPfile=../$AVWName/rnaseq_qc_results.txt
    exonic=$(grep "exonic" "$QUALIMAPfile" | awk -F'[()]' '{print $2}' | tr -d '%' || echo "N/A")
    intronic=$(grep "intronic" "$QUALIMAPfile" | awk -F'[()]' '{print $2}' | tr -d '%' || echo "N/A")
    intergenic=$(grep "intergenic" "$QUALIMAPfile" | awk -F'[()]' '{print $2}' | tr -d '%' || echo "N/A")
    overlapexon=$(grep "overlapping exon" "$QUALIMAPfile" | awk -F'[()]' '{print $2}' | tr -d '%' || echo "N/A")

    echo -e "$AVWName\t$new_r1_size\t$new_r2_size\t$num_reads\t$STARmapping_rate\t$exonic\t$intronic\t$intergenic\t$overlapexon\t$sum_counts" >> "$temp_file"


done < "$input_file"



# Replace the original file with the updated data

mv "$temp_file" "$input_file"
