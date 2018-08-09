## CIRquant2 ##

CIRIquant2 is a comprehensive analysis pipeline for circRNA detection and quantification in RNA-Seq data

### Author ###

Authors: Jinyang Zhang(zhangjinyang@biols.ac.cn), Fangqing Zhao(zhfq@biols.ac.cn)

Maintainer: Jinyang Zhang

### Release Notes ###

- Version 1.1: add help
- Version 1.0: add implementation of pysam

### License ###

The code is released under the MIT License. See the `LICENSE` file for more detail.

### Prerequisites ###

```
Softwares:
    bwa=0.7.12
    hisat2=2.1.0
    stringtie=1.3.3b
    samtools=1.5

Python packages:
    simplejson
    argparse
    pysam
    numpy
    scipy
    scikit-learn

R packages:
    optparse
    edgeR

```

**NOTE**:

Samtools version should be higher than `1.3.1`, as older version of samtools may use deprecated parameters in `sort` and other commands

### 1. Installation ###

Use the setup.py for CIRIquant installation (clean install under virutalenv is highly recommended).

```bash
# create and activate virtual env
virtualenv venv
source ./venv/bin/activate

# Install CIRIquant
git clone https://kevinzjy.github.com/CIRIquant2
cd CIRIquant2
python setup.py install

# Manual installation of required pacakges is also supported
pip install -r requirements.txt
```

### 2. Running CIRIquant For circRNA quantifcation ###

```
Usage:
  CIRIquant [options] -1 <m1> -2 <m2> --config <config>

  <m1>              Input mate1 reads (for paired-end data)
  <m2>              Input mate2 reads (for paired-end data)
  <config>          Config file

Options (defaults in parentheses):

  -v                Run in verbose mode
  --no-gene         Skipe StringTie estimation of gene abundance
  -e, --log         Specific log file (default: sample_prefix.log)
  --circ            User provide Back-Spliced Junction sites
  --bam             Specific hisat2 alignment bam file against reference genome
  -o, -out          Output directory (default: current directory)
  -p, --prefix      Output sample prefix (default: input sample name)
  -t, --threads     Number of CPU threads to use (defualt: 4)
  --RNaseR          CIRIquant output file of RNase R data (required for RNase R correction)

```

A JSON-formated config file is needed for CIRIquant to find software and reference needed. A valid example of config file is demonstrated below.

```json
// Example of config file
{
    "bwa": "bwa",
    "hisat2": "hisat2",
    "stringtie": "stringtie",
    "samtools": "samtools",
    "genome": "/histor/zhao/zhangjy/database/gencode_hg19/hg19.fa",
    "gtf": "/histor/zhao/zhangjy/database/gencode_hg19/gencode.v19.annotation.gtf",
    "bwa_index": "/histor/zhao/zhangjy/database/gencode_hg19/_BWAtmp/hg19",
    "hisat_index": "/histor/zhao/zhangjy/database/gencode_hg19/_HISATtmp/hg19"
}
```

Key | Description
----|-------------
bwa | the command of `bwa`, (default: bwa)
hisat2 | command of `hisat2`, (default: hisat2)
stringtie | command of `stringite`, (default: stringtie)
samtools | command of `samtools`, (default: samtools), samtools version below 1.3.1 is not supported
genome | reference genome fasta, a fai index by `samtools faidx` is also needed under the same directory
gtf | annotation file of reference genome
bwa_index | index of reference genome using `bwa index -w bwtsw`
hisat_index | index of reference genome using `hisat2-build`

For quantification of user-provided circRNAs, a list of junction sites in bed format is required

```
chr1    10000   10099   chr1:10000|10099    .   +
chr1    31000   31200   chr1:31000|31200    .   -
```


### 3. Output format of CIRIquant ###

The output of CIRIquant is in standard gtf format, containing information of BSJ and FSJ reads of circRNAs and annotation of circRNA back-spliced regions in the attribute columns

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
| gene_id | ensemble id of host gene |
| gene_name | HGNC symbol of host gene |
| gene_type | type of host gene in gtf file |

### 4. Example Usage ###

Test data set can be retrived under `test_data` folder, you can replace the index and reference genome in the `hg19.config` with your own genome information

```
cd test_data
CIRIquant2 --config human.config \
    -1 test_1.fq \
    -2 test_2.fq \
    -o ./test_output \
    -p test \
    -t 4
```

The output file `test.gtf` should be located under `test_data/test_output/`


### 5. Generate RNase R effect corrected BSJ information ###

In order to remove effect for RNase R treatment, two steps of programs are needed

1. Run CIRIquant with RNase R treated sample
2. Use output in Step1 and run CIRIquant with `--RNaseR` option using output gtf in previous step

The output is in the same format as normal run, however the header line is appended with additional information

### 6. Run differential expression analysis for circRNAs

The Differential expression analysis is performed using `CIRI_DE.py`

```
Usage:
  python [options] CIRI_DE.py {-n <control> -c <case> | -i <sample>} -o <out>

  <control>         CIRIquant result of control sample
  <case>            CIRIquant result of study sample
  <sample>          List of Sample files (for study with replicates)
  <out>             Output file
```

#### Study without biological replicate ####

For sample without replicate, `-c` and `-s` should be used to specific case and control data

The output format of sample without replicate is in the format below:

| column | name | description |
|--------|------|-------------|
| 1 | chrom | chromosome / contig name |
| 2 | start | 5' back-spliced junction site |
| 3 | end | 3' back-spliced junction site |


#### Group of samples ####

For group of samples, a list of sample IDS and their respective path should be provided

**NOTE**: Gene abundance are used for normalization analysis, so all samples should run without `--no-gene` option provided

```
CONTROL1    N   <PATH_TO_CONTROL1.bed>
CONTROL2    N   <PATH_TO_CONTROL2.bed>
CONTROL3    N   <PATH_TO_CONTROL3.bed>
CASE1    T   <PATH_TO_CASE1.bed>
CASE2    T   <PATH_TO_CASE2.bed>
CASE3    T   <PATH_TO_CASE3.bed>
```

The output in the same format as edgeR output
| column | name | description |
|--------|------|-------------|
| 1 | chrom | chromosome / contig name |
| 2 | start | 5' back-spliced junction site |
| 3 | end | 3' back-spliced junction site |
