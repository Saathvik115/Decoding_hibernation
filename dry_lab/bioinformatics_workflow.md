# Bioinformatics Processing and Quantification
To process raw short-read RNA-seq data from adipose tissue, the folloing bioinformatics processing and quantification steps were used.

**Note**: the analyses covered below were run on the *UCSC HPC (Hummingbird)*. In order to learn more about the specific paramters used and how to generalize them for implementation, consult the appropriate markdown files for that particular processing step.

## Pipeline Workflow

### Step 1: Quality Control and Trimming with TrimGalore

**Purpose**: Remove low-quality bases, adapter sequences, and artificial sequences from raw FASTQ files to ensure high-quality input for downstream analysis.

**Tool**: TrimGalore v0.4.2 (wrapper around Cutadapt and FastQC)

**Key Parameters**:
- Quality threshold: Phred score ≥ 24 (99.6% accuracy)
- Adapter detection: Illumina TruSeq adapters
- Minimum read length: 50 bases after trimming
- Hard clipping: 12 bases from 5' end of both reads
- Paired-end mode with synchronization

**Slurm Script Implementation**:
```bash
#using trimgalore to trim the fastq files
trim_galore --paired -q 24 --fastqc --fastqc_args "--noextract --nogroup --outdir analysis/1_trim/fastqc" --stringency 5 --illumina --length 50 -o analysis/1_trim --clip_R1 12 --clip_R2 12 data/fastq_files/${read1}.fastq data/fastq_files/${read2}.fastq
```

**Output**: Quality-trimmed FASTQ files and FastQC reports stored in `analysis/1_trim/`

---

### Step 2: Read Alignment with STAR

**Purpose**: Map trimmed sequencing reads to the brown bear reference genome, enabling detection of splice junctions and quantification of gene expression.

**Tool**: STAR v2.7.10b (Spliced Transcripts Alignment to a Reference)

**Reference Genome**: Brown bear (*Ursus arctos*) assembly GCF_023065955.2_UrsArc2.0

**Key Parameters**:
- Multi-threading: 12 cores
- Splice junction overhang: 100 bases
- Two-pass alignment mode for improved splice junction detection
- Unique mapping only (removes multi-mapping reads)
- Sorted BAM output

#### Step 2A: Genome Index Generation

**Slurm Script Implementation**:
```bash
# create STAR genome index
STAR --runMode genomeGenerate --runThreadN 12 --genomeDir data/star_genome_index --genomeFastaFiles data/genome/GCF_023065955.2_UrsArc2.0_genomic.fna --sjdbGTFfile data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf --sjdbOverhang 100
```

#### Step 2B: Read Alignment

**Slurm Script Implementation**:
```bash
#running star alignment
STAR --runThreadN 12 --genomeDir data/star_genome_index --sjdbGTFfile data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf --sjdbOverhang 100 --outFileNamePrefix analysis/2_star/${read1}_ --outFilterMultimapNmax 1 --readFilesIn analysis/1_trim/${read1}_val_1.fq analysis/1_trim/${read2}_val_2.fq --twopassMode Basic --outSAMtype BAM SortedByCoordinate
```

**Output**: Sorted BAM files containing aligned reads stored in `analysis/2_star/`

---

### Step 3: Gene-Level Read Quantification with featureCounts

**Purpose**: Count the number of reads mapping to each gene to generate a gene expression matrix suitable for downstream differential expression analysis.

**Tool**: featureCounts from Subread v1.6.3

**Key Parameters**:
- Paired-end mode with proper pair counting
- Multi-threading: 8 cores
- Feature type: exon-level counting
- Gene-level aggregation using gene_id attribute
- Annotation: GTF file matching the reference genome

**Environment Setup**:
```bash
# Load miniconda3
module load miniconda3

# Create conda environment with subread package
conda create -n featureCounts_env -c bioconda subread

# Activate environment
conda activate featureCounts_env
```

**Slurm Script Implementation**:
```bash
#running featureCounts
featureCounts -p -T 8 -t exon -g gene_id -a data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf -o analysis/3_featurecounts/bear_adipose_rawCounts.txt analysis/2_star/*_Aligned.sortedByCoord.out.bam
```

**Output**: Gene-level read count matrix (`bear_adipose_rawCounts.txt`) ready for differential expression analysis

---

## Pipeline Summary

This three-step pipeline transforms raw RNA-seq data into a quantified gene expression matrix:

1. **Quality Control** → Clean, high-quality reads
2. **Alignment** → Genomic positions of reads 
3. **Quantification** → Gene expression levels

The workflow is designed for scalability using Slurm array jobs, processing 11 brown bear samples in parallel across high-performance computing resources. Each step generates comprehensive output files and logs for quality assessment and troubleshooting. Below is a figure that visualizes the processing workflow used in this project.

<img width="930" height="517" alt="workflowimage" src="https://github.com/user-attachments/assets/a80957a2-5497-429a-ae66-218e08c8ae2f" />

## File Structure

Obviously, directories can have any structure so long as are able to store the results of each of these steps. I've attached an image below depicting my file structure:

<img width="194" height="413" alt="structureimage" src="https://github.com/user-attachments/assets/6a3058fa-e5b1-4729-9344-a37fab968f5a" />
