# KMA-cgMLST

## Usage

### Preparing your cgMLST scheme
First, prepare a single multi-fasta file containing all alleles. The [pubmlst_client](https://github.com/Public-Health-Bioinformatics/pubmlst_client)
tool may be helpful for finding and downloading MLST schemes.

The fasta deflines of each allele should follow the format: `locus_allele`. eg:
```
>MYCO000001_1
CGATCGATGCTATACTAGG.....
>MYCO000001_2
CGATGCTTAGCGATCTACGT....
```

Index the fasta using `kma`:
```
kma index -i <your_scheme.fa> -o <your_scheme>
```

### Running the pipeline

```
nextflow run BCCDC-PHL/kma-cgmlst \
  --fastq_input </path/to/fastqs> \
  --scheme </path/to/cgmlst_scheme> \
  --outdir </path/to/output_dir> 
```
