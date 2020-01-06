## CIRIquant ##

[![Build Status](https://travis-ci.com/Kevinzjy/CIRIquant.svg?branch=master)](https://travis-ci.com/Kevinzjy/CIRIquant)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/Kevinzjy/CIRIquant)
[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](https://github.com/Kevinzjy/CIRIquant/blob/master/LICENSE)
![GitHub All Releases](https://img.shields.io/github/downloads/Kevinzjy/CIRIquant/total)
![SourceForge](https://img.shields.io/sourceforge/dm/ciri/CIRIquant)

CIRIquant is a comprehensive analysis pipeline for circRNA detection and quantification in RNA-Seq data

### Author ###

Authors: Jinyang Zhang(zhangjinyang@biols.ac.cn), Fangqing Zhao(zhfq@biols.ac.cn)

Maintainer: Jinyang Zhang

### Release Notes ###

- Version 1.0: The first released version of CIRIquant

### License ###

The code is released under the MIT License. See the `LICENSE` file for more detail.

### Citing CIRIquant

- Zhang, J., Chen, S., Yang, J. et al. Accurate quantification of circular RNAs identifies extensive circular isoform switching events. Nat Commun 11, 90 (2020) [doi:10.1038/s41467-019-13840-9](https://doi.org/10.1038/s41467-019-13840-9)

### Prerequisites ###

```
Softwares:
    bwa
    hisat2
    stringtie
    samtools>=1.9

Python packages:
    PyYAML
    argparse
    pysam
    numpy
    scipy
    scikit-learn
```

**NOTES**:

**1. Only python2 is supported**  

**2. Samtools version should be higher than `1.9`, as older version of samtools may use 
deprecated parameters in `sort` and `index` commands**

### 1. Installation ###

**Please use the latest released version from [GitHub](https://github.com/Kevinzjy/CIRIquant/releases) or [SourceForge](https://sourceforge.net/projects/ciri/files/CIRIquant/)**

Use the setup.py for CIRIquant installation (clean install under virutalenv is highly recommended).

```bash
# create and activate virtual env
pip install virtualenv
virtualenv venv
source ./venv/bin/activate

# Install CIRIquant and its requirement automatically
tar zxvf CIRIquant.tar.gz
cd CIRIquant
python setup.py install

# Manual installation of required pacakges is also supported
pip install -r requirements.txt
```

The package should take approximately 40 seconds to install on a normal computer.

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
  --no-gene         Skip StringTie estimation of gene abundance
```

A YAML-formated config file is needed for CIRIquant to find software and reference needed.
A valid example of config file is demonstrated below.

```YAML
// Example of config file
name: hg19
tools:
  bwa: /home/zhangjy/bin/bwa
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
bwa_index | prefix of BWA index for reference genome
hisat_index | prefix of HISAT2 index for reference genome

For quantification of user-provided circRNAs, a list of junction sites in bed format is required, for example:

```
chr1    10000   10099   chr1:10000|10099    .   +
chr1    31000   31200   chr1:31000|31200    .   -
```

**NOTE**: 
- For now, --circ and --tool options can only parse CIRI2 results.
- Gene expression values are needed for normalization, do not use `--no-gene` if you need to run DE analysis afterwards. 

### 3. Output files ###

The main output of CIRIquant is a GTF file, that contains detailed information of 
BSJ and FSJ reads of circRNAs and annotation of circRNA back-spliced regions in the attribute columns

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

Test data set can be retrived from `test_data.tar.gz`, you can replace the path of required software in
 the `chr1.yml` with your own version

```
tar zxvf test_data.tar.gz
cd test_data/quant
CIRIquant -t 4 \
          -1 ./test_1.fq.gz \
          -2 ./test_2.fq.gz \
          --config ./chr1.yml \
          --no-gene \
          -o ./test \
          -p test
```

The output file `test.gtf` should be located under `test_data/quant/test`  

The demo dataset should take approximately 5 minutes on a personal computer. It has been tested on 
my PC with Intel i7-8700 processor and 16G of memory, running Ubuntu 18.04 LTS.

### 5. Generate RNase R effect corrected BSJ information ###

In order to remove effect for RNase R treatment, two steps of programs are needed

1. Run CIRIquant with RNase R treated sample
2. Use output gtf file in Step1 and run CIRIquant with `--RNaseR` option using output gtf in previous step

The output is in the same format as normal run, however the header line is appended with additional 
information of RNase R treatment

### 6. Run differential expression analysis for circRNAs

#### Study without biological replicate ####

For sample without replicate, the differential expression & differential splicing analysis is 
performed using `CIRI_DE`

```
Usage:
  CIRI_DE [options] -n <control> -c <case> -o <out>

  <control>         CIRIquant result of control sample
  <case>            CIRIquant result of treatment cases
  <out>             Output file

Options (defaults in parentheses):

  -p                p value threshold for DE and DS score calculation (default: 0.05)
  -t                numer of threads (default: 4)

Example usage:
  CIRI_DE -n control.gtf -c case.gtf -o CIRI_DE.csv
```

The output format `CIRI_DE` is in the format below:

| column | name | description |
|--------|------|-------------|
| 1 | circRNA_ID | circRNA identifier |
| 2 | Case_BSJ | number of BSJ reads in case |
| 3 | Case_FSJ | number of FSJ reads in case |
| 4 | Case_Ratio | junction ratio in case |
| 5 | Ctrl_BSJ | number of BSJ reads in control |
| 6 | Ctrl_FSJ | number of FSJ reads in control |
| 7 | Ctrl_Ratio | junction ratio  in control |
| 8 | DE_score | differential expression score |
| 9 | DS_score | differential splicing score |

#### Study with biological replicates ####

For study with biological replicates, a customed analysis pipeline of edgeR is recommended and 
we provide `prep_CIRIquant` to generate matrix of circRNA expression level / junction ratio and `CIRI_DE_replicate` 
for DE analysis

**Step1**: Prepare CIRIquant output files

One should provide a text file listing sample information and path to CIRIquant output GTF files

```
CONTROL1 ./c1/c1.gtf C 1
CONTROL2 ./c2/c2.gtf C 2
CONTROL3 ./c3/c3.gtf C 3
CASE1 ./t1/t1.gtf T 1
CASE2 ./t2/t2.gtf T 2
CASE3 ./t3/t3.gtf T 3
```

The first three columns is required by default. For paired samples, you could also add a column of subject name.

| column | description |
|--------|-------------|
| 1 | sample name |
| 2 | path to CIRIquant output gtf |
| 3 | group ("C" for control, "T" for treatment) |
| 4 | subject (optional, only for paired samples) |

Then, run `prep_CIRIquant` to summarize the circRNA expression profile in all samples

```
Usage:
  prep_CIRIquant [options]

  -i                the file of sample list
  --lib             where to output library information
  --circ            where to output circRNA annotation information
  --bsj             where to output the circRNA expression matrix
  --ratio           where to output the circRNA junction ratio matrix

Example:
  prep_CIRIquant -i sample.lst \
                 --lib library_info.csv \
                 --circ circRNA_info.csv \
                 --bsj circRNA_bsj.csv \
                 --ratio circRNA_ratio.csv
```

These count matrices (CSV files) can then be imported into R for use by DESeq2 and edgeR 
(using the DESeqDataSetFromMatrix and DGEList functions, respectively).

**Step2**: Prepare StringTie output

The output of StringTie should locate under `output_dir/gene/prefix_out.gtf`. You need to use 
[prepDE.py](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py) from stringTie to
generate the gene count matrix for normalization.

For example, one can provide a text file `sample_gene.lst` containing sample IDs and path to StringTie outputs:

```text
CONTROL1 ./c1/gene/c1_out.gtf
CONTROL2 ./c2/gene/c2_out.gtf
CONTROL3 ./c3/gene/c3_out.gtf
CASE1 ./t1/gene/t1_out.gtf
CASE2 ./t2/gene/t2_out.gtf
CASE3 ./t3/gene/t3_out.gtf
```

Then, run `prepDE.py -i sample_gene.lst` and use `gene_count_matrix.csv` generated under current working directory 
for further analysis.

**Step3**: Differential expression analysis

For differential analysis using `CIRI_DE_replicate`, you need to install a R environment and `edgeR` package from Bioconductor.

```bash
Usage:
  CIRI_DE_replicate [options]

  --lib             library information by CIRIquant
  --bsj             circRNA expression matrix
  --gene            gene expression matrix
  --out             output differential expression result

Example:
  CIRI_DE.R --lib  library_info.csv \
            --bsj  circRNA_bsj.csv \
            --gene gene_count_matrix.csv \
            --out  circRNA_de.csv
```

Please be noted that the output results is **unfiltered**, 
and you could apply a more stringent filter on expression values to get a more convincing result.