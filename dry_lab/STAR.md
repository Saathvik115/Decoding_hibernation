# STAR

## Purpose

**STAR** (Spliced Transcripts Alignment to a Reference) is a widely used RNA-seq aligner designed for ultrafast and accurate alignment of high-throughput sequencing reads to reference genomes. STAR is particularly well-suited for aligning reads that span exon-exon junctions, making it ideal for transcriptome studies, gene expression analysis, and splice junction detection.

STAR's alignment process involves two distinct phases:

1. **Genome Index Generation**: Creating a searchable index from a reference genome and annotation files
2. **Read Alignment**: Mapping sequencing reads to the indexed genome, producing alignment files for downstream analysis

The tool employs a **seed-and-extend** approach combined with **suffix arrays** to achieve rapid alignment speeds while maintaining high accuracy. STAR can detect novel splice junctions and handle reads with multiple exons, making it superior to traditional DNA aligners for RNA-seq data.

After quality trimming with tools like Trim Galore, the next critical step in RNA-seq analysis is aligning the cleaned reads to a reference genome. STAR produces alignment files in SAM/BAM format that can be used for downstream analyses such as:

- Gene expression quantification
- Differential expression analysis
- Splice junction discovery
- Variant calling from RNA-seq data
- Fusion gene detection

The bears in this project all have corresponding trimmed FASTQ files from the previous quality control step. These files will be aligned to the brown bear (*Ursus arctos*) reference genome using STAR's two-step workflow.
## Methods

All steps in this alignment pipeline utilize **Slurm scripts** for execution on high-performance computing environments. As described in previous workflows, Slurm scripts enable efficient resource management and job scheduling on HPC clusters (see SRA&FASTQ.md to get started with SLurm scripts).

STAR's alignment workflow consists of two essential steps that must be performed sequentially:

### Genome Index Generation

The first step involves creating a STAR genome index, which is a compressed, searchable representation of the reference genome. This process:

1. **Processes the reference genome FASTA file** to create suffix arrays for rapid sequence searching
2. **Incorporates annotation information** from GTF/GFF files to identify known splice junctions (in this project, a GTF file was used)
3. **Generates multiple index files** that enable fast alignment and junction detection
4. **Creates a directory structure** containing all necessary files for the alignment step

The genome index only needs to be generated once per reference genome and can be reused for multiple alignment jobs. This step is computationally intensive but significantly speeds up subsequent alignment processes. I've linked the reference genome FASTA file and GTF file for accessibility.

### Read Alignment

The second step performs the actual alignment of sequencing reads to the indexed genome:

1. **Loads the pre-generated genome index** into memory
2. **Processes paired-end reads** from trimmed FASTQ files
3. **Performs seed-and-extend alignment** with splice-aware mapping
4. **Detects splice junctions** both from annotation and novel discoveries
5. **Outputs alignment results** in SAM/BAM format with comprehensive alignment statistics

STAR's **two-pass alignment mode** can be enabled to improve alignment accuracy by first identifying splice junctions in an initial pass, then using these junctions to guide a second, more accurate alignment pass.

## Implementation

On my HPC (Hummingbird), I used two separate Slurm scripts to complete the STAR alignment workflow. The first script generates the genome index, while the second performs the actual read alignment.

### Genome Index Generation Script

```bash
#!/bin/bash
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/065/955/GCF_023065955.2_UrsArc2.0/ contains the reference genome for Ursus arctos

#slurm script to run STAR genome index creation
#SBATCH --job-name=rungenomeidx		# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=svoora@ucsc.edu   	# Where to send mail
#SBATCH --time=0-05:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                 		# Run a single task
#SBATCH --cpus-per-task=12                	# Use 12 threads for STAR
#SBATCH --output=/hb/groups/sip_eeb_01/saat/scripts/logs/rungenomeidx_%j.out    # Standard output and error log
#SBATCH --error=/hb/groups/sip_eeb_01/saat/scripts/logs/rungenomeidx_%j.err     # Standard output and error log
#SBATCH --mem=32G                    # Allocate memory for the job.

# moves to working directory
cd /hb/groups/sip_eeb_01/saat || exit 1

#downloading any tools/modules
module load star/2.7.10b   

#creates directory to store STAR data
mkdir -p data/star_genome_index

# create STAR genome index
STAR --runMode genomeGenerate --runThreadN 12 --genomeDir data/star_genome_index --genomeFastaFiles data/genome/GCF_023065955.2_UrsArc2.0_genomic.fna --sjdbGTFfile data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf --sjdbOverhang 100
```

This script requires substantial computational resources for genome index creation. The **32GB memory allocation** is necessary because STAR loads the entire genome into memory during index generation. **12 CPU cores** are utilized to parallelize the indexing process, significantly reducing runtime for large genomes like the brown bear reference.

The script loads `star/2.7.10b`, which provides access to the STAR aligner and its dependencies. Version 2.7.10b includes important improvements for splice junction detection and memory efficiency.

The script creates `data/star_genome_index/` to store the generated index files. This directory will contain multiple files including suffix arrays, splice junction databases, and genome parameters.

**STAR Command Parameters:**
- `--runMode genomeGenerate`: Specifies that this is a genome index generation run, not an alignment run
- `--runThreadN 12`: Uses 12 threads for parallel processing, matching the CPU allocation in the Slurm header
- `--genomeDir data/star_genome_index`: Specifies the output directory for the generated index files
- `--genomeFastaFiles data/genome/GCF_023065955.2_UrsArc2.0_genomic.fna`: Points to the reference genome FASTA file for *Ursus arctos*
- `--sjdbGTFfile data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf`: Provides the gene annotation file to incorporate known splice junctions into the index
- `--sjdbOverhang 100`: Sets the overhang length for splice junctions. This should be (read length - 1); for 101bp reads, 100 is optimal

### Read Alignment Script

```bash
#!/bin/bash

#slurm script to run STAR alignment on trimmed fastq data
#SBATCH --job-name=runSTAR		# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=svoora@ucsc.edu   	# Where to send mail
#SBATCH --time=0-10:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                 		# Run a single task
#SBATCH --cpus-per-task=12                	# Use 12 threads for STAR
#SBATCH --output=/hb/groups/sip_eeb_01/saat/scripts/logs/runstar_%j.out    # Standard output and error log
#SBATCH --error=/hb/groups/sip_eeb_01/saat/scripts/logs/runstar_%j.err     # Standard output and error log
#SBATCH --mem=32G                    # Allocate memory for the job.
#SBATCH --array=1-11					# array job

# moves to working directory
cd /hb/groups/sip_eeb_01/saat || exit 1

#downloading any tools/modules
module load star/2.7.10b   

#creates directory to store STAR data
mkdir -p analysis/2_star

#foldr with fastq files
fastq="data/fastq_names.tsv"

#call line in file we're processing
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$fastq")
read1=$(echo ${LINE} | awk '{ print $1; }')
read2=$(echo ${LINE} | awk '{ print $2; }')

#running star alignment
STAR --runThreadN 12 --genomeDir data/star_genome_index --sjdbGTFfile data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf --sjdbOverhang 100 --outFileNamePrefix analysis/2_star/${read1}_ --outFilterMultimapNmax 1 --readFilesIn analysis/1_trim/${read1}_val_1.fq analysis/1_trim/${read2}_val_2.fq --twopassMode Basic --outSAMtype BAM SortedByCoordinate
```

This script uses **array job functionality** to process all 11 samples efficiently. The **10-hour time limit** accommodates the potentially longer runtime required for alignment compared to index generation. **32GB memory** is allocated to handle the genome index loaded into memory during alignment.

Similar to previous workflows, the script uses `fastq_names.tsv` to manage paired-end FASTQ files. The `sed` and `awk` commands extract the appropriate R1 and R2 file names for each array job iteration.

The script creates `analysis/2_star/` to store alignment results, maintaining the organized pipeline structure where each analysis step has its own numbered directory.

**STAR Alignment Parameters:**
- `--runThreadN 12`: Uses 12 threads for parallel alignment processing
- `--genomeDir data/star_genome_index`: Points to the previously generated genome index directory
- `--sjdbGTFfile data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf`: Provides annotation for splice junction detection during alignment
- `--sjdbOverhang 100`: Must match the value used during index generation for consistency
- `--outFileNamePrefix analysis/2_star/${read1}_`: Sets the output file prefix using the sample name, creating unique files for each sample
- `--outFilterMultimapNmax 1`: Allows only uniquely mapping reads, filtering out multi-mapping reads that could introduce ambiguity
- `--readFilesIn analysis/1_trim/${read1}_val_1.fq analysis/1_trim/${read2}_val_2.fq`: Specifies the input trimmed FASTQ files from the previous Trim Galore step
- `--twopassMode Basic`: Enables two-pass alignment for improved accuracy by first identifying splice junctions, then realigning with this information
- `--outSAMtype BAM SortedByCoordinate`: Outputs directly to sorted BAM format, eliminating the need for separate SAM-to-BAM conversion and sorting steps

This comprehensive STAR workflow produces high-quality alignments ready for downstream RNA-seq analyses such as gene expression quantification, differential expression analysis, and splice junction characterization. The sorted BAM files generated by this pipeline can be directly used with tools like featureCounts or HTSeq (this project will utilize the former).
