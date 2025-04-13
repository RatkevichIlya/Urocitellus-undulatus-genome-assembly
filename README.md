# Urocitellus-undulatus-genome-assembly

This repository includes scripts used in the assembly of the *Urocitellus undulatus* genome.

## Requirements

- Linux-based HPC cluster with Slurm
- 512G RAM
- 4 TB of disk space
- Conda (Tested with Miniconda3)
- FastQC - quality control of raw illumina reads
- MultiQC - visualization of FastQC and BUSCO results
- SeqKit - assembly of basic statistics for raw data and genomes
- FastP - filtering and trimming short illumina reads
- NanoQC - quality control of nanopore raw reads
- Porechop - filtering and trimming of nanopore reads
- MaSuRCA - genome assembly and analysis toolkit (installed via Conda)
- Flye - genome assembler used within MAECI-based pipeline
- Canu - genome assembler used in the MAECI-based pipeline
- wtdbg2 - genome assembler used in the MAECI-based pipeline
- NextPolish - polishing the assembly with high quality Illumina reads
- BUSCO - evaluation of genome assembly completeness
- BUSCO glires_odb10 dataset
- Quast - assessment the quality of genome assembly
- RepeatModeler - identification of tandem repeats in the genome and creation of a de novo library of them
- RepeatMasker - tandem repeat masking
- BRAKER - fully automated method for accurate gene structure annotation
- BlobToolKit – interactive quality assessment of genome assemblies
## 1. Installation


### 1.1. Create, activate and install tools in Conda environment

```bash
TOOLS=(
    fastqc
    seqkit
    fastp
    nanoqc
    porechop
    masurca
    flye
    canu
    wtdbg
    busco
    quast
    multiqc
    nextpolish
    repeatmasker
    repeatmodeler
)

# Loop through each tool, create conda environment, and install the package
for TOOL in "${TOOLS[@]}"; do
    echo "Creating conda environment for $TOOL..."
    conda create -n "$TOOL" -y
    conda activate "$TOOL"
    echo "Installing $TOOL from bioconda..."
    conda install -c bioconda "$TOOL" -y
    conda deactivate
    echo "$TOOL installation completed."
done

conda create -n blobtoolkit -y
conda activate blobtoolkit
echo "Installing Blobtoolkit from HCC..."
conda install hcc::blobtoolkit -y
conda deactivate
echo "Blobtoolkit installation completed."

echo "All tools installed successfully."
```
## 2. Evaluation and filtration of raw data

### 2.1. Estimation of baseline parameters of sequencing data using SeqKit

```bash
conda activate seqkit
seqkit stats <raw-data-files> > output.txt
```

### 2.2. Illumina data evaluation and filtration

Evaluation of illumina raw data with FastQC
```bash
conda activate fastqc
fastqc --outdir <raw-output-dir> <fastq-file1> <fastq-file2>
```
FastQC results visualisation using MultiQC
```bash
conda activate multiqc 
multiqc <raw-output-dir>
```

Filtration of illumina raw data
```bash
conda activate fastp
fastp --in1 <fastq-file1> \
--in2 <fastq-file2> \
--out1 <trimm-fastq-file1> \
--out2 <trimm-fastq-file2> \
-l 70 --cut_front --cut_tail --dedup \
--trim_tail1 1 \
-h fastp_report.html
```

Evaluation of illumina filtered data
```bash
conda activate fastqc
fastqc --outdir <trimm-output-dir> <trimm-fastq-file1> <trimm-fastq-file2>

conda activate multiqc 
multiqc <trimm-output-dir>
```

### 2.3. Nanopore data evaluation and filtration

Evaluation of nanopore raw data
```bash
conda activate nanoqc
nanoqc <nanopore-file>
```

Filtration of nanopore raw data
```bash
conda activate porechop
porechop -i <nanopore-file> -o <trimm-nanopore-file>
```

Evaluation of nanoopore filtered data
```bash
conda activate nanoqc
nanoqc <trimm-nanopore-file>
```
### 2.4 Statistics collection of cleansed data
```bash
conda activate seqkit
seqkit stats <filtered-data-files> > output.txt
```

## 3. Genome assembly and evaluation

### 3.1. MaSuRCA assembly with default parameters (based on Celera assembler)

Assembly
```bash
conda activate masurca
masurca -t 32 -i <fastq-file1>,<fastq-file2> -r <trimm-nanopore-file>
```

Assembly evaluation with BUSCO
```bash
conda activate busco
busco -i <assembly> -m genome -o ./output -c 64 -l glires_odb10 --offline
conda activate multiqc
mulitqc output/
```

Obtaining assembly statistics with SeqKit and Quast
```bash
conda activate seqkit
seqkit stats <assembly.fasta> > output.txt
conda activate quast
quast --threads 32 <assembly.fasta>
```

Improving assembly quality by polishing Illumina reads with NextPolish

1. Prepare sgs_fofn
```bash
ls reads1_R1.fq reads1_R2.fq > sgs.fofn
```

2. Create run.cfg
```bash
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
deltmp = yes
rerun = 3
parallel_jobs = 8
multithread_jobs = 16
genome = input.genome.fasta
genome_size = auto
workdir = ./01_rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa
```
3. Run
```bash
conda activate nextpolish
nextPolish run.cfg
```

4. SeqKit
```bash
conda activate seqkit
seqkit stats <assembly> > output.txt
```


5. BUSCO
```bash
conda activate busco
busco -i <assembly> -m genome -o ./busco_out -c 64 -l glires_odb10 --offline
multiqc busco_out/
```

### 3.2. MaSuRCA assembly using Flye builder and manually customized config file

#### Assembly

Config file data

```bash
DATA
PE= pe 500 50 /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq
NANOPORE=/FULL_PATH/nanopore.fa
END

PARAMETERS
EXTEND_JUMP_READS=0
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 1
USE_GRID=0
GRID_ENGINE=SLURM
GRID_QUEUE=all.q
GRID_BATCH_SIZE=500000000
LHE_COVERAGE=30
MEGA_READS_ONE_PASS=0
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.15
CLOSE_GAPS=1
NUM_THREADS = 64
JF_SIZE = 56000000000
SOAP_ASSEMBLY=0
FLYE_ASSEMBLY=1
END
```
Run the MaSuRCA script to generate a shell script from the config file, then run the resulting script.

```bash
conda activate masurca
masurca config.txt
./assembly.sh
```

Further steps to improve and evaluate the resulting assembly correspond to the steps given for running 
MaSuRCA with default parameters (see section 3.1.).

### 3.3. MAECI-based assembly pipeline

``` bash
#!/bin/bash -i

# Creating a directory for self-correction
mkdir -p self-correction

# Copying the assemblies
cp /FULL_PATH/wtdbg2/dbg.raw.fa self-correction/wtdbg2.fasta
cp /FULL_PATH/flye/assembly self-correction/flye.fasta
cp /FULL_PATH/canu/canu_out/NP01.contigs.fasta self-correction/canu.fasta

# Selecting the longest assembly by file size
ls -S self-correction/*.fasta | head -n 1 > self-correction/best_assembly.fasta

# Activating the environment and creating an index
conda activate masurca
minimap2 -d self-correction/reference.mmi $(cat self-correction/best_assembly.fasta)

# Alignment of Nanopore reads
minimap2 -a -x map-ont $(cat self-correction/best_assembly.fasta) \
    /FULL_PATH/<trimm-nanopore-file> -t 32 > self-correction/temp.sam

# First round of correction with Racon
conda activate racon
racon -t 32 /FULL_PATH/<trimm-nanopore-file> \
    self-correction/temp.sam $(cat self-correction/best_assembly.fasta) > self-correction/racon_1.fasta

# Second round of correction 
conda activate masurca
minimap2 -a -x map-ont self-correction/racon_1.fasta \
    /FULL_PATH/<trimm-nanopore-file> -t 32 > self-correction/temp.sam

conda activate racon
racon -t 32 /FULL_PATH/<trimm-nanopore-file> \
    self-correction/temp.sam self-correction/racon_1.fasta > self-correction/racon_2.fasta

# Third round of correction
conda activate masurca
minimap2 -a -x map-ont self-correction/racon_2.fasta \
    /FULL_PATH/<trimm-nanopore-file> -t 32 > self-correction/temp.sam

conda activate racon
racon -t 32 /FULL_PATH/<trimm-nanopore-file> \
    self-correction/temp.sam self-correction/racon_2.fasta > self-correction/racon_3.fasta

# Temporary files cleaning
rm -rf self-correction/temp.sam self-correction/racon_1.fasta \
       self-correction/racon_2.fasta self-correction/reference.mmi
```


Further steps to evaluate the resulting assembly correspond to the steps given for running 
MaSuRCA with default parameters (see section 3.1.).

## 4. 	Visualization assemblies statistics
Visualization of the obtained statistics for each assembly to identify the most complete and continuous assembly.

Self-written python scripts busqo\_visualization.py and quast\_visualization.py were used 
to visualize BUSCO and quast analysis results, respectively. 

The self-described python script busco_piechart.py was used to visualize the BUSCO scores 
of the selected genome with the best statistics.

All data and paths to them for self-written python scripts were written directly in the script.

Blobtoolkit tool was used to visualize the baseline statistics with snail-plot.

```bash
conda activate blobtoolkit
blobtools create --fasta <assembly.fasta> BlobDir
blobtools view --view snail BlobDir
```

## 5. Annotation

Further steps include annotation of the most complete genome

### The section is being finalized

1. *De novo* repeats identification with RepeatModeler utilizaton

2. Repeats soft masking with the usage of RepeatMasker

3. 
