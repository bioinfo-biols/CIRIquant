## CIRquant ##

CIRIquant is a comprehensive analysis pipeline for circRNA detection and quantification in RNA-Seq data

### Author ###

Authors: Jinyang Zhang(zhangjinyang@biols.ac.cn), Fangqing Zhao(zhfq@biols.ac.cn)

Maintainer: Jinyang Zhang

### Release Notes ###

- Version 0.1: init commit

### License ###

The code is released under the MIT License. See the `LICENSE` file for more detail.

### Prerequisites ###

```
Softwares:
    bwa
    hisat2
    stringtie
    samtools>=1.3.1

Python packages:
    PyYAML
    argparse
    pysam
    numpy
    scipy
    scikit-learn
```

**NOTE**:

Samtools version should be higher than `1.3.1`, as older version of samtools may use deprecated parameters in `sort` and other commands

### 1. Installation ###

Use the setup.py for CIRIquant installation (clean install under virutalenv is highly recommended).

```bash
# create and activate virtual env
virtualenv venv
source ./venv/bin/activate

# Install CIRIquant and its requirement automatically
tar zxvf CIRIquant.tar.gz
cd CIRIquant
python setup.py install

# Manual installation of required pacakges is also supported
pip install -r requirements.txt
```

### 2. Running CIRIquant For circRNA quantifcation ###

```
Usage:
  CIRIquant [options] --config <config> -1 <m1> -2 <m2>

  <config>          Config file
  <m1>              Input mate1 reads (for paired-end data)
  <m2>              Input mate2 reads (for paired-end data)


Options (defaults in parentheses):

  -v                Run in verbose mode
  -o, -out          Output directory (default: current directory)
  -e, --log         Specific log file (default: sample_prefix.log)
  -p, --prefix      Output sample prefix (default: input sample name)
  -t, --threads     Number of CPU threads to use (defualt: 4)
  
  --bed             User provided Back-Spliced Junction Site in BED format
  --circ            User provided circRNA prediction results
  --tool            User provided tool name for circRNA prediction

  --RNaseR          CIRIquant output file of RNase R data (required for RNase R correction)
  --bam             Specific hisat2 alignment bam file against reference genome
  --no-gene         Skipe StringTie estimation of gene abundance
```

A YAML-formated config file is needed for CIRIquant to find software and reference needed. A valid example of config file is demonstrated below.

```YAML
// Example of config file
name: hg19
tools:
  bwa: /home/zhangjy/bin/bwa/
  hisat2: /home/zhangjy/bin/hisat2
  stringtie: /home/zhangjy/bin/stringtie
  samtools: /home/zhangjy/bin/samtools

reference:
  fasta: /home/zhangjy/Data/database/hg19.fa
  gtf: /home/zhangjy/Data/database/gencode.v19.annotation.gtf
  bwa_index: /home/zhangjy/Data/database/hg19/_BWAtmp/hg19
  hisat_index: /home/zhangjy/Data/database/hg19/_HISATtmp/hg19
```

Key | Description
----|-------------
name| the name of config file
bwa | the path of `bwa`
hisat2 | the path of `hisat2`
stringtie | the path of `stringite`
samtools | the path of `samtools`, samtools version below 1.3.1 is not supported
fasta | reference genome fasta, a fai index by `samtools faidx` is also needed under the same directory
gtf | annotation file of reference genome
bwa_index | index of reference genome using `bwa index -w bwtsw`
hisat_index | index of reference genome using `hisat2-build`

For quantification of user-provided circRNAs, a list of junction sites in bed format is required, for example:

```
chr1    10000   10099   chr1:10000|10099    .   +
chr1    31000   31200   chr1:31000|31200    .   -
```

**NOTE**: --circ and --tool options can only parse CIRI2 results temporarily, the support for other algorithms is under active development 
### 3. Output files ###

The main output of CIRIquant is a GTF file, that contains detailed information of BSJ and FSJ reads of circRNAs and annotation of circRNA back-spliced regions in the attribute columns

Description of each columns's value

| column | name | description |
|--------|------|-------------|
| 1 | chrom | chromosome / contig name |
| 2 | source | CIRIquant |
| 3 | type | circRNA
| 4 | start | 5' back-spliced junction site |
| 5 | end | 3' back-spliced junction site |
| 6 | score | CPM of circRNAs (#BSJ / #Mapped reads)
| 7 | strand | strand information
| 8 | . | .
| 9 | attributes | attributes seperated by semicolon |

The attributes containing several pre-defined keys and values:

| key | description|
| --- | -----------|
| circ_id | name of circRNA |
| circ_type | circRNA types: exon / intron / intergenic |
| bsj | number of bsj reads |
| fsj | number of fsj reads |
| junc_ratio | circular to linear ratio: 2 * bsj / ( 2 * bsj + fsj)
| rnaser_bsj | number of bsj reads in RNase R data (only when --RNaseR is specificed)
| rnaser_fsj | number of fsj reads in RNase R data (only when --RNaseR is specificed)
| gene_id | ensemble id of host gene |
| gene_name | HGNC symbol of host gene |
| gene_type | type of host gene in gtf file |

### 4. Example Usage ###

Test data set can be retrived under `test_data/quant` folder, you can replace the path of required software in the `chr1.yml` with your own version

```
cd test_data/quant
CIRIquant2 --config chr1.yml \
    -1 test_1.fq.gz \
    -2 test_2.fq.gz \
    -o ./test_output \
    -p test \
    -t 4
```

The output file `test.gtf` should be located under `test_data/quant/test_output/`


### 5. Generate RNase R effect corrected BSJ information ###

In order to remove effect for RNase R treatment, two steps of programs are needed

1. Run CIRIquant with RNase R treated sample
2. Use output gtf file in Step1 and run CIRIquant with `--RNaseR` option using output gtf in previous step

The output is in the same format as normal run, however the header line is appended with additional information of RNase R treatment

### 6. Run differential expression analysis for circRNAs

#### Study without biological replicate ####

For sample without replicate, the differential expression & differential splicing analysis is performed using `CIRI_DE`

```
Usage:
  CIRI_DE [options] -n <control> -c <case> -o <out>

  <control>         CIRIquant result of control sample
  <case>            CIRIquant result of study sample
  <out>             Output file

Options (defaults in parentheses):

  -p                p value threshold for DE and DS score calculation (default: 0.05)
  -t                numer of threads (default: 4)  
```

The output format `CIRI_DE` is in the format below:

| column | name | description |
|--------|------|-------------|
| 1 | circRNA_ID | circRNA identifier |
| 2 | Case_BSJ | number of BSJ reads in case |
| 3 | Case_FSJ | number of FSJ reads in case |
| 4 | Case_Ratio | junction ratio in case |
| 5 | Ctrl_BSJ | number of BSJ reads in case |
| 6 | Ctrl_FSJ | number of FSJ reads in case |
| 7 | Ctrl_Ratio | junction ratio  in case |
| 8 | DE_score | differential expression score |
| 9 | DS_score | differential splicing score |

#### Study with biological replicates ####

**NOTE**: Gene abundance are used for normalization analysis, so all samples should run without `--no-gene` option provided

For study with biological replicates, a customed analysis pipeline of edgeR is recommended, the `prep_CIRIquant.py` 
is provided to generate matrix of circRNA expression level / junction ratio 

A tab-seperated list of sample IDS and path to CIRIquant output GTF file should be specificed respectively

```
CONTROL1    <PATH_TO_CIRIquant_GTF>
CONTROL2    <PATH_TO_CIRIquant_GTF>
CONTROL3    <PATH_TO_CIRIquant_GTF>
CASE1    <PATH_TO_CIRIquant_GTF>
CASE2    <PATH_TO_CIRIquant_GTF>
CASE3    <PATH_TO_CIRIquant_GTF>
```

Then, run `prep_CIRIquant` to summarize the circRNA expression profile in all samples

```
Usage:
  prep_CIRIquant [options]

  -i                the file of sample list
  --bsj             where to output the circRNA expression matrix
  --ratio           where to output the circRNA junction ratio matrix
```

These count matrices (CSV files) can then be imported into R for use by DESeq2 and edgeR (using the DESeqDataSetFromMatrix and DGEList functions, respectively).

For example, we provide the demo scripts for differential expression using `edgeR`

```R
library('edgeR')

# Load raw data
lib_data <- read.csv("./lib_data.csv", header = T, row.names = 1)
lib_data <- lib_data[order(lib_data$Patient),]
gene_data <- read.csv("./gene_counts.csv", header = T, row.names = 1)
circ_data <- read.csv("./circ_bsj.csv", header = T, row.names = 1)

# TMM normalization using gene read counts
gene_DGE <- DGEList(counts = gene_data,
                    group = lib_data[colnames(gene_data), "Group"])
gene_DGE <- calcNormFactors(gene_DGE)

# Normalizing circRNA expression matrix
circ_DGE <- DGEList(counts = circ_data,
                    group = gene_DGE$samples[, "group"],
                    lib.size = gene_DGE$samples[, "lib.size"],
                    norm.factors = gene_DGE$samples[, "norm.factors"])

# Build design matrix
design <- model.matrix(~factor(lib_data[colnames(circ_data), "Sample"])
                       + factor(lib_data[colnames(circ_data), "Group"]))

# Estimate dispersion
circ_DGE <- estimateDisp(circ_DGE, design)

# Calculagte P-value
circ_fit <- glmFit(circ_DGE, design)
circ_lrt <- glmLRT(circ_fit)
topTags(circ_lrt)

# Output DE genes
o <- order(circ_lrt$table$PValue)
de <- decideTestsDGE(circ_lrt)
df <- circ_lrt$table
df$DE <- decideTestsDGE(circ_lrt)
df$FDR <- p.adjust(df$PValue, method = "fdr")
df <- df[o,]
write.csv(df, file = "./circ_de.csv", quote = FALSE)
```