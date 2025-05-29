# Urocitellus-undulatus-genome-assembly

This repository includes scripts used in the assembly of the *Urocitellus undulatus* genome.
The workflow presented is not a linked pipeline.
Variable and path names have sometimes been replaced to simplify command presentation.

## Requirements
System requirements
- Linux-based HPC cluster
- 512G RAM
- 4 TB of free disk space
- Conda (Tested with Miniconda3)
Raw data preprocessing
- FastQC - quality control of raw illumina reads
- MultiQC - visualization of FastQC and BUSCO results
- SeqKit - assembly of basic statistics for raw data and genomes
- FastP - filtering and trimming short illumina reads
- NanoQC - quality control of nanopore raw reads
- Porechop - filtering and trimming of nanopore reads
Genome assemblies
- MaSuRCA - genome assembly and analysis toolkit (installed via Conda)
- Flye - genome assembler used within MAECI-based pipeline
- Canu - genome assembler used in the MAECI-based pipeline
- wtdbg2 - genome assembler used in the MAECI-based pipeline
- NextPolish - polishing the assembly with high quality Illumina reads
Assemblies evaluation
- BUSCO - evaluation of genome assembly completeness
- BUSCO mammalia_odb12 dataset
- Quast - assessment the quality of genome assembly
- BlobToolKit – interactive quality assessment of genome assemblies
Genome annotation
- RepeatModeler - identification of tandem repeats in the genome and creation of a de novo library of them
- RepeatMasker - tandem repeat masking
- GeMoMa - a homology-based gene prediction program
- minimap2 - sequence alignment program
- ragtag - reference-based scaffolds orientation
- SYRI - predicting genomic differences between related genomes
Mitogenome assembly
- GetOrganelle - mitogenome assembly
Phylogenetic analysis
- MAFFT - sequence alignment program
- Gblocks - alignments filtering
- AMAS - alignments concatenation 
- IQ-TREE - phylogenetic tree reconstruction

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
    gemoma
    minimap2
    ragtag
    syri
    getorganelle
    mafft
    gblocks
    amas
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

# BlobToolKit installation is separated because of another channel
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

Self-written python scripts busco\_visualization.py and quast\_visualization.py were used 
to visualize BUSCO and quast analysis results, respectively. 

The self-described python script busco_piechart.py was used to visualize the BUSCO scores 
of the selected genome with the best statistics.

All data and paths to them for self-written python scripts were written directly in the script.

Blobtoolkit tool was utilized to visualize the baseline statistics with snail-plot.
```bash
conda activate blobtoolkit
blobtools create --fasta <assembly.fasta> BlobDir
blobtools view --view snail BlobDir
```

## 5. Annotation

Further steps include annotation of the most complete genome

### 1. *De novo* repeats identification with RepeatModeler utilizaton

Activating the conda environment
```bash
conda activate repeatmodeler
```

Create a new RepeatModeler BLAST database from the input genome.
```bash
BuildDatabase -name undulatus -engine ncbi genome.fasta
```

De novo repeats identification using RepeatModeler
```bash
RepeatModeler -pa 16 -engine ncbi -database ID_Genus-species
```

Adding species identificator as prefix to each fasta header line
```bash
cat consensi.fa.classified | seqkit fx2tab | awk '{ print "urocUnd1_"$0 }' | seqkit tab2fx > consensi.fa.prefix.fa.classified
```

Separating our repeats library into known and unknown elements
```bash
cat consensi.fa.prefix.fa.classified | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > consensi.fa.prefix.fa.classified.known
cat consensi.fa.prefix.fa.classified | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > consensi.fa.prefix.fa.classified.unknown
```

Quantify number of known and unknown classified elements
```bash
grep -c ">" consensi.fa.prefix.fa.classified.known
grep -c ">" consensi.fa.prefix.fa.classified.unknown
```

### 2. Repeats annotation and soft masking with the usage of RepeatMasker

conda activate repeatmasking

```bash
#Variables
GENOME=genome.fasta
KNOWN=consensi.fa.prefix.fa.classified.known
UNKNOWN=consensi.fa.prefix.fa.classified.unknown

# Output directories
mkdir -p logs 01_simple_out 02_tetrapoda_out 03_known_out 04_unknown_out 05_full_out

# Simple repeats
RepeatMasker -pa $THREADS -a -e ncbi -dir 01_simple_out -noint -xsmall $GENOME

# Tetrapoda repeats
RepeatMasker -pa $THREADS -a -e ncbi -dir 02_tetrapoda_out -nolow \
-species tetrapoda 01_simple_out/$(basename "$GENOME").masked

# Known repeats
RepeatMasker -pa $THREADS -a -e ncbi -dir 03_known_out -nolow \
-lib $KNOWN 02_tetrapoda_out/$(basename "$GENOME").masked.masked

# Unknown repeats
RepeatMasker -pa $THREADS -a -e ncbi -dir 04_unknown_out -nolow \
-lib $UNKNOWN 03_known_out/$(basename "$GENOME").masked.masked.masked

# Concatenation .out files
cat 01_simple_out/$(basename "$GENOME").out \
    <(tail -n +4 02_tetrapoda_out/$(basename "$GENOME").masked.out) \
    <(tail -n +4 03_known_out/$(basename "$GENOME").masked.masked.out) \
    <(tail -n +4 04_unknown_out/$(basename "$GENOME").masked.masked.masked.out) \
    > 05_full_out/genome.full_mask.out

# GFF3
rmOutToGFF3.pl 05_full_out/genome.full_mask.out > 05_full_out/genome.full_mask.gff3

# Genome soft masking
bedtools maskfasta -soft -fi "$GENOME" -bed 05_full_out/genome.full_mask.gff3 \
-fo 05_full_out/genome.full_masked.fasta
```

The repeats_vis.py script was used to visualize the results of the annotation of repeats.

### 3. Homology-based genome annotation using GeMoMa

The -w parameter was used to prioritize annotations of reference genomes of evolutionarily closer species for annotation.

```bash
mkdir -p gemoma_results

conda activate gemoma

GeMoMa -Xmx500G GeMoMaPipeline \
    threads=32 \
    AnnotationFinalizer.r=NO \
    p=false \
    o=true \
    t=genome.full_masked.fasta \
    outdir=gemoma_results/ \
    s=own w=3 i=u_parryii a=genomes/u_parryii/GCF_003426925.1/genomic.gff g=genomes/u_parryii/GCF_003426925.1/GCF_003426925.1_ASM342692v1_genomic.fna \
    s=own w=2 i=ictidomys a=genomes/Ictidomys_trid/GCF_016881025.1/genomic.gff g=genomes/Ictidomys_trid/GCF_016881025.1/GCF_016881025.1_HiC_Itri_2_genomic.fna \
    s=own w=2 i=mmm a=genomes/marmarmar/GCF_001458135.2/genomic.gff g=genomes/marmarmar/GCF_001458135.2/GCF_001458135.2_marMar_genomic.fna \
    s=own w=2 i=marm_flav a=genomes/marm_flav/GCF_047511675.1/genomic.gff g=genomes/marm_flav/GCF_047511675.1/GCF_047511675.1_mMarFla1.hap1_genomic.fna \
    s=own w=2 i=marm_monax a=genomes/Marmota_monax/GCF_021218885.2/genomic.gff g=genomes/Marmota_monax/GCF_021218885.2/GCF_021218885.2_Marmota_monax_Labrador192_F-V1.1_genomic.fna \
    s=own w=1 i=norw_rat a=genomes/norway_rat/GCF_036323735.1/genomic.gff g=genomes/norway_rat/GCF_036323735.1/GCF_036323735.1_GRCr8_genomic.fna \
    s=own w=1 i=black_rat a=genomes/black_rat/GCF_011064425.1/genomic.gff g=genomes/black_rat/GCF_011064425.1/GCF_011064425.1_Rrattus_CSIRO_v1_genomic.fna
```

The annotation\_vis.py and annotation\_vis_comp.py script were used to visualize the results of gene annotation to visualize 
the statistics of the obtained genomic annotation and to compare it with the reference genome annotation, respectively.

## 6. Verification of the representativeness of the genome assembly

To analyze structural similarity, we first need to create an extracted_chroms.fasta file that 
will contain only the assembled chromosomes of the reference genome.
To do this, create a file chrom_ids.txt containing only the names of the collected chromosomes. 
Then use seqkit to extract the records.

```bash
conda activate seqkit
seqkit grep -f chrom_ids.txt genome.fasta > extracted_chroms.fasta
```

Reference genome-based orientation of scaffolds. 
This step greatly accelerates the subsequent assembly of pseudochromosomes.
Gaps (stretches of "N" characters) are placed between adjacent query sequences to indicate the presence of unknown sequence.
```bash
conda activate ragtag
ragtag.py scaffold -t 32 extracted_chroms.fasta genome.full_masked.fasta
```

Reference-based assembly of pseudochromosomes for our genome. 
This is necessary for subsequent structural similarity analysis.
```bash
conda activate minimap2
minimap2 -t 128 -ax asm10 --eqx extracted_chroms.fasta ragtag_output/ragtag.scaffold.fasta > out.sam
conda activate syri
chroder -o chroder -F S out.sam extracted_chroms.fasta ragtag_output/ragtag.scaffold.fasta
```

After obtaining pseudochromosomes, the alignment must be done again. 
Then we run syri to identify structural differences and visualize the result using plotsr.
```bash
conda activate minimap2
minimap2 -t 32 -ax asm20 --eqx extracted_chroms.fasta chroder.qry.fasta > syri.sam
conda activate syri
syri -c syri.sam -r extracted_chroms.fasta -q chroder.qry.fasta -k -F S --nc 32
# The genomes.txt file contains the name of the assemblies being compared. 
# The order.txt file contains the names of chromosomes in the order in which they will be displayed by plotsr. 
plotsr --sr syri.out --genomes genomes.txt -o plotsr.png --chrord order.txt
```

Visualizations of the distribution of the number of scaffolds aligned to the reference chromosomes and 
the length of alignments attributable to each chromosome were obtained using the scripts 
map\_count\_visualization.py and map\_length\_visualization.py, respectively.

## 7. Mitogenome assembly

```bash
conda activate getorganelle
get_organelle_from_reads.py -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz \
  -o mitogenome_output -F animal_mt \
  -s u_undulatus.fasta,u_richardsonii.fasta,u_parryii1.fasta,u_parryii2.fasta \
  -t 32
```

## 8. Phylogenetics analysis

Establish the tools needed for phylogenetic analysis into one environment:

```bash
conda create -n phylo_env mafft iqtree amas gblocks biopython -c bioconda -c conda-forge
```

Phylogenetic analysis
```bash
conda activate phylo_env
for gene_set in sets_names.txt; do
    mafft --auto "$gene_set" > "${gene_set%.fasta}aligned.fasta"
done

for alignment in *_aligned.fasta; do
  Gblocks "$alignment" -t=d -b5=h -e=.gb
done

AMAS.py concat -i *_aligned.fasta.gb -u fasta -y nexus -d dna -t concatenated.fasta

iqtree2 -s concatenated.fasta -m MFP -B 1000 -alrt 1000 -T AUTO --prefix 15samples -o Castor_canadensis_NC_015108_1,Castor_fiber_NC_015072_1
```


## References

fastp – Chen S. (2023). fastp: an ultra-fast all-in-one FASTQ preprocessor. iMeta, 2(1), e107. https://doi.org/10.1002/imt2.107;

NanoPack – De Coster W., D’Hert S., Schultz D.T., Cruts M., Van Broeckhoven C. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666–2669. https://doi.org/10.1093/bioinformatics/bty149;

SeqKit 2.0 – Shen W., Le S., Li Y., Hu F. (2024). SeqKit 2.0: a powerful and ultrafast toolkit for FASTA/Q file manipulation. iMeta, 3(1), e191. https://doi.org/10.1002/imt2.191;
dep – slw287r. (2024). dep: Data Evaluation Pipeline. GitHub repository. https://github.com/slw287r/dep;

MaSuRCA – Zimin A.V., Puiu D., Luo M.-C., Zhu T., Koren S., Marçais G., Yorke J.A., Dvořák J., Salzberg S.L. (2017). Hybrid assembly of the large and highly repetitive genome of Aegilops tauschii, a progenitor of bread wheat, with the MaSuRCA mega-reads algorithm. Genome Research, 27(5), 787–792. https://doi.org/10.1101/gr.213405.116;

CABOG (в составе MaSuRCA) – Miller J.R., Delcher A.L., Koren S., Venter E., Walenz B.P., Brownley A., Johnson J., Li K., Mobarry C., Sutton G. (2008). Aggressive assembly of pyrosequencing reads with mates. Bioinformatics, 24(24), 2818–2824. https://doi.org/10.1093/bioinformatics/btn548;

wtdbg2 – Ruan J., Li H. (2019). Fast and accurate long-read assembly with wtdbg2. Nature Methods, 17(2), 155–158. https://doi.org/10.1038/s41592-019-0669-3;

Pilon – Walker B.J., Abeel T., Shea T., Priest M., Abouelliel A., Sakthikumar S., Cuomo C.A., Zeng Q., Wortman J., Young S.K., Earl A.M. (2014). Pilon: An integrated tool for comprehensive microbial variant detection and genome assembly improvement. PLoS ONE, 9(11), e112963. https://doi.org/10.1371/journal.pone.0112963;

QUAST-LG – Mikheenko A., Prjibelski A., Saveliev V., Antipov D., Gurevich A. (2018). Versatile genome assembly evaluation with QUAST-LG. Bioinformatics, 34(13), i142–i150. https://doi.org/10.1093/bioinformatics/bty266;

BUSCO – Manni M., Berkeley M.R., Seppey M., Simão F.A., Zdobnov E.M. (2021). BUSCO Update: Novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular Biology and Evolution, 38(10), 4647–4654. https://doi.org/10.1093/molbev/msab199;

BlobToolKit – Challis R., Richards E., Rajan J., Cochrane G., Blaxter M. (2020). BlobToolKit – interactive quality assessment of genome assemblies. G3: Genes|Genomes|Genetics, 10(4), 1361–1374. https://doi.org/10.1534/g3.119.400908;

RepeatModeler – Dfam-consortium. (2024). RepeatModeler. GitHub repository. https://github.com/Dfam-consortium/RepeatModeler;

RepeatMasker – Dfam-consortium. (2017). RepeatMasker. GitHub repository. https://github.com/Dfam-consortium/RepeatMasker;

RepBase – Bao W., Kojima K.K., Kohany O. (2015). Repbase Update, a database of repetitive elements in eukaryotic genomes. Mobile DNA, 6, 11. https://doi.org/10.1186/s13100-015-0041-9;

minimap2 – Li H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37(23), 4572–4574. https://doi.org/10.1093/bioinformatics/btab705;

SyRI – Goel M., Sun H., Jiao W.-B., Schneeberger K. (2019). SyRI: finding genomic rearrangements and local sequence differences from whole-genome assemblies. Genome Biology, 20, 277. https://doi.org/10.1186/s13059-019-1911-0;

RagTag – Alonge M., Lebeigle L., Kirsche M., Aganezov S., Wang X., Lippman Z.B., Schatz M.C. (2022). Automated assembly scaffolding elevates a new tomato system for high-throughput genome editing. Genome Biology, 23, 258. https://doi.org/10.1186/s13059-022-02823-7;

plotsr – Goel M., Schneeberger K. (2022). plotsr: visualizing structural similarities and rearrangements between multiple genomes. Bioinformatics, 38(4), 1192–1196. https://doi.org/10.1093/bioinformatics/btac196;

GetOrganelle – Jin J.-J., Yu W.-B., Yang J.-B., Song Y., dePamphilis C.W., Yi T.-S., Li D.-Z. (2020). GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. Genome Biology, 21, 241. https://doi.org/10.1186/s13059-020-02154-5;

ElasticBLAST – Camacho C., Hatcher E.L., Grigoriev I.V., etc. (2023). ElasticBLAST, a cloud-based BLAST tool. BMC Bioinformatics, 24, 95. https://doi.org/10.1186/s12859-023-05245-9;

MAFFT-DASH – Rozewicki J., Li S., Amada K.M., Standley D.M., Katoh K. (2019). MAFFT-DASH: integrated protein sequence and structure alignment. Nucleic Acids Research, 47(W1), W5–W10. https://doi.org/10.1093/nar/gkz342;

Gblocks – Talavera G., Castresana J. (2007). Improvement of phylogenies after removing divergent and ambiguously aligned blocks from protein sequence alignments. Systematic Biology, 56(4), 564–577. https://doi.org/10.1080/10635150701472164;

AMAS – Borowiec M.L. (2016). AMAS: a fast tool for alignment manipulation and computing of summary statistics. PeerJ, 4, e1660. https://doi.org/10.7717/peerj.1660.
