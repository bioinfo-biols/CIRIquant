## CIRquant2 ##

CIRIquant2 is a comprehensive analysis pipeline for circRNA detection and quantification in RNA-Seq data

### Author ###

Authors: Jinyang Zhang(zhangjinyang@biols.ac.cn), Fangqing Zhao(zhfq@biols.ac.cn)

Maintainer: Jinyang Zhang

### Release Notes ###

- Version 1.1: add help
- Version 1.0: add implementation of pysam

### License ###

The code is released under the MIT License. See the `LICENSE` for more detail.

### Prerequisites ###

```
Softwares:
    bwa=0.7.12
    hisat2=2.1.0
    stringtie=1.3.3b
    samtools=1.5

Python packages:
    numpy
    simplejson
    pysam
    argparse

R packages:
    optparse
    edgeR

```

**NOTE**:

Samtools version should be higher than `1.3.1`, as different version of samtools may use different parameters in `sort` and other commands

The required python packages are listed in `requirements.txt`, and will be automatically installed, you can use `pip` for manual installation

```bash
pip install -r requirements.txt
```

### 1. Installation ###

Use the setup.py for CIRIquant installation (clean install under virutalenv is recommended).

```bash
# create and activate virtual env
virtualenv venv
source ./venv/bin/activate

# Install CIRIquant
git clone https://kevinzjy.github.com/CIRIquant2
cd CIRIquant2
python setup.py install
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
  --RNaseR          Run RNase R correction

```

A JSON-formated config file is needed for CIRIquant to find reference genome and their hisat2 index etc. A valid format of config file is demonstrate below.

```json
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

For user-provided circRNA junction sites, a circRNA list in bed format is needed

```
chr1    10000   10099   chr1:10000-10099    0   +
chr1    31000   31200   chr1:31000-31200    0   -
```


### 3. Output format of CIRIquant ###

The output of CIRIquant is in bed format, containing information of BSJ and FSJ reads of circRNAs

| column | name | description |
|--------|------|-------------|
| 1 | chrom | chromosome / contig name |
| 2 | start | 5' back-spliced junction site |
| 3 | end | 3' back-spliced junction site |

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

The output result should locate under `test_data/test_output/test.CIRI2_info.list`


### 5. Generate RNase R effect corrected BSJ information ###

In order to remove effect for RNase R treatment, two steps of programs are needed

1. Run CIRIquant with RNase R treated sample
2. Use output in Step1 and run CIRIquant with `--RNaseR` option using Total RNA data

The output should be in the same bed format, and BSJ should be the corrected expression value for circRNAs

| column | name | description |
|--------|------|-------------|
| 1 | chrom | chromosome / contig name |
| 2 | start | 5' back-spliced junction site |
| 3 | end | 3' back-spliced junction site |

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
