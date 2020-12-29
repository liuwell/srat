# srat
A small RNA analysis tool

## Dependency 
Sevral external software were depended for srat:
+ cutadapt >= 2.10
+ bowtie

## Mandatory
+ biopython
+ numpy
+ pandas
+ matplotlib
+ seaborn



## Usage

### Show help message of srat.py

```bash
srat.py -h
```

```
usage: srat.py [-h] -i INPUT [-o OUTDIR] -l LIBRARY
               [-t {common,testis,sperm,EV}] [-p THREADS] [--no_merge]
               [--spikein] [-v]

A small RNA analysis and visualization tool

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        the input directory of raw data
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -l LIBRARY, --library LIBRARY
                        the reference sequence for mapping and annotation
  -t {common,testis,sperm,EV}, --tissue {common,testis,sperm,EV}
                        the sample tissue type
  -p THREADS, --threads THREADS
                        number of threads to launch (default: 1)
  --no_merge            not merge the expression files
  --spikein             consider spikein in the library
  -v, --version         show program's version number and exit
```

### Merge small RNA profiles

```bash
merge.py -h
```

```
usage: merge.py [-h] -i INPUT -prefix PREFIX [-barplot]
                [-t {common,testis,sperm,EV}]

Merge small RNA profiles

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        the input directory of processed data
  -prefix PREFIX        the prefix of merging files
  -barplot              perform barplot of merged files
  -t {common,testis,sperm,EV}, --tissue {common,testis,sperm,EV}
                        the sample tissue type

```

### Perform correlation plot

```bash
corr_heatmap.py -h
```

```
usage: corr_heatmap.py [-h] -i INPUT

The correlation and heatmap of expression data, such as miRNA expression, gene
expression

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        the input data

```

