

# Phasing

This contains the pipeline used to run the genotyping experiments for the HGSVC3 paper.
See the steps below to re-produce the analysis.

## How to replicate the results


### Step 1: Downloading input files

Download PanGenie 1kg genotypes for HGSVC3 data from:
```bat
TODO
```


Download CHM13 reference genome:

```bat
wget -P inputs/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
```

 Download GRCh38 reference genome:

```bat
wget -P inputs/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz
gunzip inputs/hg38.no_alt.fa.gz
```


Download all genetic maps for CHM13 `` chr*.t2t.gmap.resorted.gmap.gz" `` from: https://github.com/JosephLalli/phasing_T2T/tree/master/t2t_lifted_chrom_maps and put them to ``inputs/``


Download and prepare Illumina-based SNP and indel calls for 1kg samples:


```

### Step 2: run snakemake pipeline
Dependencies:
* snakemake
* singularity
* conda

```bat
snakemake -j <nr_cores> --use-conda
```

