# Urocitellus-undulatus-genome-assembly

This repository includes scripts used in the assembly of the Urocitellus undulatus genome.

## Requirements

- Linux-based HPC cluster with Slurm
- 512G RAM
- 3 TB of disk space
- Conda (tested with Miniconda3)
- FastQC (Quality control of raw illumina reads)
- MultiQC (Visualization of FastQC and BUSCO results)
- SeqKit (Assembly of basic statistics for raw data and genomes)
- FastP (Filtering and trimming short illumina reads)
- NanoQC (Quality control of nanopore raw reads)
- Porechop (Filtering and trimming of nanopore reads)
- MaSuRCA (installed via Conda)
- Flye (genome assembler used within MAECI-based pipeline)
- Canu (genome assembler used in the MAECI-based pipeline)
- wtdbg2 (genome assembler used in the MAECI-based pipeline)
- BUSCO (Evaluation of genome assembly completeness)
- BUSCO glires_odb10 dataset
- Quast (Assessment the quality of genome assembly)

## Installation

```bash
### Create, activate and install tools in Conda environment

#!/bin/bash

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
)

### Loop through each tool, create conda environment, and install the package
for TOOL in "${TOOLS[@]}"; do
    echo "Creating conda environment for $TOOL..."
    conda create -n "$TOOL" -y
    conda activate "$TOOL"
    echo "Installing $TOOL from bioconda..."
    conda install -c bioconda "$TOOL" -y
    conda deactivate
    echo "$TOOL installation completed."
done

echo "All tools installed successfully."


## Evaluation and filtration of raw data

### Estimation of baseline parameters of sequencing data using SeqKit

conda activate seqkit
seqkit stats <raw-data-files> > output.txt

## Illumina data evaluation and filtration

### Evaluation of illumina raw data with FastQC
conda activate fastqc
fastqc --outdir <raw-output-dir> <fastq-file1> <fastq-file2>

### FastQC results visualisation using MultiQC
conda activate multiqc 
multiqc <raw-output-dir>

### Filtration of illumina raw data
conda activate fastp
fastp --in1 <fastq-file1> \
--in2 <fastq-file2> \
--out1 <trimm-fastq-file1> \
--out2 <trimm-fastq-file2> \
-l 70 --cut_front --cut_tail --dedup \
--trim_tail1 1 \
-h fastp_report.html

### Evaluation of illumina filtered data
conda activate fastqc
fastqc --outdir <trimm-output-dir> <trimm-fastq-file1> <trimm-fastq-file2>

conda activate multiqc 
multiqc <trimm-output-dir>

## Nanopore data evaluation and filtration

### Evaluation of nanopore raw data
conda activate nanoqc
nanoqc <nanopore-file>

### Filtration of nanopore raw data
conda activate porechop
porechop -i <nanopore-file> -o <trimm-nanopore-file>

### Evaluation of nanoopore filtered data
conda activate nanoqc
nanoqc <trimm-nanopore-file>



## Genome assembly and evaluation

### MaSuRCA assembly with default parameters

** Assembly**
conda activate masurca
masurca -t 32 -i <fastq-file1>,<fastq-file2> -r <trimm-nanopore-file>

** Assembly evaluation with BUSCO**

Conda activate busco
busco -i <assembly> -m genome -o ./output -c 64 -l glires_odb10 --offline

** Obtaining assembly statistics with SeqKit**
seqkit stats <assembly> > output.txt

### MaSuRCA assembly using Flye builder and manually customized config file

** Assembly**
'config file data'
: '
DATA
PE= pe 500 50 /mnt/tank/scratch/iratkevich/squirrel/raw_data/L001_R1_001.fastq.gz /mnt/tank/scratch/iratkevich/squirrel/raw_data/L001_R2_001.fastq.gz
NANOPORE=/mnt/tank/scratch/iratkevich/squirrel/raw_data/nanopore.fastq.gz
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
'

conda activate masurca
masurca config.txt
./assembly.sh

** Assembly evaluation with BUSCO**

Conda activate busco
busco -i <assembly> -m genome -o ./output -c 64 -l glires_odb10 --offline

** Obtaining assembly statistics with SeqKit**
seqkit stats <assembly> > output.txt

### MAECI-based assembly pipeline

#!/bin/bash -i

** Creating a directory for self-correction**
mkdir -p self-correction

** Copying the assemblies**
cp /mnt/tank/scratch/iratkevich/squirrel/wtdbg2/dbg.raw.fa self-correction/wtdbg2.fasta
cp /mnt/tank/scratch/iratkevich/squirrel/flye/assembly self-correction/flye.fasta
cp /mnt/tank/scratch/iratkevich/squirrel/canu/canu_out/NP01.contigs.fasta self-correction/canu.fasta

** Selecting the longest assembly by file size **
ls -S self-correction/*.fasta | head -n 1 > self-correction/best_assembly.fasta

** Activating the environment and creating an index **
conda activate masurca
minimap2 -d self-correction/reference.mmi $(cat self-correction/best_assembly.fasta)

** Alignment of Nanopore reads **
minimap2 -a -x map-ont $(cat self-correction/best_assembly.fasta) \
    /mnt/tank/scratch/iratkevich/squirrel/trimmed/trimm_nanopore.fastq.gz -t 32 > self-correction/temp.sam

** First round of correction with Racon **
conda activate racon
racon -t 32 /mnt/tank/scratch/iratkevich/squirrel/trimmed/trimm_nanopore.fastq.gz \
    self-correction/temp.sam $(cat self-correction/best_assembly.fasta) > self-correction/racon_1.fasta

** Second round of correction **
conda activate masurca
minimap2 -a -x map-ont self-correction/racon_1.fasta \
    /mnt/tank/scratch/iratkevich/squirrel/trimmed/trimm_nanopore.fastq.gz -t 32 > self-correction/temp.sam

conda activate racon
racon -t 32 /mnt/tank/scratch/iratkevich/squirrel/trimmed/trimm_nanopore.fastq.gz \
    self-correction/temp.sam self-correction/racon_1.fasta > self-correction/racon_2.fasta

** Third round of correction **
conda activate masurca
minimap2 -a -x map-ont self-correction/racon_2.fasta \
    /mnt/tank/scratch/iratkevich/squirrel/trimmed/trimm_nanopore.fastq.gz -t 32 > self-correction/temp.sam

conda activate racon
racon -t 32 /mnt/tank/scratch/iratkevich/squirrel/trimmed/trimm_nanopore.fastq.gz \
    self-correction/temp.sam self-correction/racon_2.fasta > self-correction/racon_3.fasta

** Temporary files cleaning **
rm -rf self-correction/temp.sam self-correction/racon_1.fasta \
       self-correction/racon_2.fasta self-correction/reference.mmi

** Scaffolding, polishing, BUSCO, SeqKit etc. **
