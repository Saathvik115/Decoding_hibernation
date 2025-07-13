# Trim Galore

## Purpose
**Trim Galore** is a wrapper tool around two essential bioinformatics programs: **Cutadapt** and **FastQC**. It automates the process of quality and adapter trimming for FASTQ files in high-throughput sequencing data analysis. This tool is particularly valuable in bioinformatics pipelines as it ensures that downstream analyses, such as read alignment and variant calling, are performed on high-quality, clean sequencing data.

After obtaining FASTQ files from tools like `fasterq-dump` (as described in the SRA workflow), the next critical step is quality control and trimming. Raw sequencing reads often contain several issues that can compromise downstream analysis:

1. **Adapter contamination** - Sequencing adapters may remain attached to reads when the DNA fragment is shorter than the read length
2. **Low-quality bases** - Sequencing quality typically degrades toward the 3' end of reads  
3. **Artificial sequences** - Technical artifacts from library preparation
4. **Biased positions** - Certain positions may show systematic biases, especially in bisulfite sequencing

Trim Galore addresses these issues through a comprehensive process that combines quality trimming, adapter removal, and quality control assessment. The tool can automatically detect common adapter types (Illumina Universal, Nextera, and small RNA adapters) and handles both single-end and paired-end sequencing data. For paired-end data, it ensures that both read pairs maintain synchronization after trimming, which is crucial for downstream alignment and analysis.

The bears in this project all have corresponding FASTQ files generated from the SRA conversion process described earlier. The next step involves trimming these FASTQ files to remove low-quality bases and adapter sequences before proceeding with alignment and variant analysis.

## Methods
All steps in this quality trimming pipeline utilize a **Slurm script** for execution on high-performance computing environments. As described in the SRA workflow, Slurm scripts enable efficient resource management and job scheduling on HPC clusters like Hummingbird at UCSC.

In this step of the pipeline, we will be loading the *trimgalore* module via our SLURM script. Trim Galore combines several key functionalities that make it an ideal choice for automated quality control:

### Quality Trimming Algorithm
The tool employs a **sliding window approach** similar to the algorithm used by BWA:
1. Subtracts the quality threshold from all base quality scores
2. Calculates cumulative sums from each position to the end of the read
3. Identifies the position where the cumulative sum is minimal
4. Trims the read at that position

### Adapter Detection and Removal
Trim Galore uses **Cutadapt** for adapter trimming, which offers several advantages:
- Handles partial adapter matches with configurable stringency
- Tolerates sequencing errors within adapter sequences
- Processes both 5' and 3' adapter contamination
- Supports automatic detection of common adapter types

### Quality Control Integration
The tool automatically generates **FastQC reports** for trimmed files, providing immediate feedback on the effectiveness of the trimming process. This integration eliminates the need for separate quality control steps and ensures consistency in the analysis pipeline.

## Implementation
On my HPC (Hummingbird), I used the following Slurm script to run Trim Galore for quality and adapter trimming:

```bash
#!/bin/bash

#SBATCH --job-name=runtrimgalore            # Job name
#SBATCH --partition=128x24                  # Partition name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=svoora@ucsc.edu         # Where to send mail
#SBATCH --time=0-05:00:00                   # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=1                   # Use 1 thread for trimgalore
#SBATCH --output=scripts/logs/trimgalorerun_%j.out    # Standard output and error log
#SBATCH --error=scripts/logs/trimgalorerun_%j.err     # Standard output and error log
#SBATCH --mem=8G                            # Allocate memory for the job.
#SBATCH --array=1-11                        # array job

# moves to working directory
cd /hb/groups/sip_eeb_01/saat || exit 1

#downloading any tools/modules
module load trimgalore

#creates directory to store trimgalore data
mkdir -p analysis/1_trim/fastqc

#folder with fastq files
fastq="data/fastq_names.tsv"

#call line in file we're processing
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$fastq")
read1=$(echo ${LINE} | awk '{ print $1; }')
read2=$(echo ${LINE} | awk '{ print $2; }')

#using trimgalore to trim the fastq files
#using paired-end reads, quality score of 24, fastqc, stringency of 5, illumina adapter trimming, and minimum length of 50
#outputting to 1_trim folder
#clipping 12 bases from the start of each read
#fastqc output will be in 1_trim/fastqc folder
#using --paired option for paired-end reads
trim_galore --paired -q 24 --fastqc --fastqc_args "--noextract --nogroup --outdir analysis/1_trim/fastqc" --stringency 5 --illumina --length 50 -o analysis/1_trim --clip_R1 12 --clip_R2 12 data/fastq_files/${read1}.fastq data/fastq_files/${read2}.fastq
```

My Slurm header includes several important parameters for this trimming workflow. **Memory allocation (8G)** is specified because while Trim Galore typically has low memory requirements, 8GB ensures sufficient resources for processing large FASTQ files. Similar to the SRA workflow, I used an **array job** to process all 11 samples efficiently rather than creating separate scripts for each sample.

I loaded the `trimgalore` module, which provides access to both Trim Galore and its dependencies (Cutadapt and FastQC) in the HPC environment. The script then creates a well-organized directory structure with `analysis/1_trim/` as the main output directory for trimmed files and `analysis/1_trim/fastqc/` as a dedicated subdirectory for FastQC reports.

Similar to the SRA workflow, I used a TSV file (`fastq_names.tsv`) containing the names of paired FASTQ files. The script extracts the appropriate R1 and R2 file names for each array job using `sed` and `awk`. I then initialized variables `read1` and `read2` to store the first and second columns of the line, which contain the paired FASTQ file names.

**Trim Galore Paramaters** The main `trim_galore` command includes several critical parameters:

- `--paired`: Indicates paired-end reads, ensuring both R1 and R2 files are processed together and remain synchronized
- `-q 24`: Sets the quality score threshold. A Phred score of 24 corresponds to 99.6% base call accuracy, providing stringent quality filtering
- `--fastqc`: Automatically runs FastQC on trimmed files to generate quality reports
- `--stringency 5`: Controls the minimum overlap between read and adapter required for trimming
- `--illumina`: Specifically targets Illumina TruSeq adapters, the most common in modern sequencing
- `--length 50`: Sets minimum read length after trimming; shorter reads are discarded
- `--clip_R1 12 --clip_R2 12`: Removes 12 bases from the 5' end of both reads to eliminate potentially biased positions

The `--fastqc_args` parameter passes additional arguments to FastQC, including `--noextract` to keep output compressed, `--nogroup` to prevent base grouping in plots, and `--outdir` to specify the FastQC output directory.

This workflow efficiently processes paired-end FASTQ files, removing low-quality bases and adapters while generating comprehensive quality reports. The trimmed files are then ready for downstream analysis steps such as read alignment, variant calling, or differential expression analysis. After trimming, one should examine the FastQC reports to verify that quality scores have improved, adapter contamination has been removed, and overall data quality meets your requirements for downstream analysis.

I analyzed the FastQC reports my Trim Galore run gave me by using the **`multiqc`** package. Rather than loading each of the individual HTML files for FastQC that Trim Galore outputs per SRA number, the `multiqc` package enables the user to compile those FastQC reports in one place. In Hummingbird, I was able to load the module by typing `module load multiqc` in my console. Then, by changing into the directory with all my FastQC reports, I ran the tool by typing `multiqc .` in my console.
