

# Genotyping

This contains the pipeline used to run the genotyping experiments for the HGSVC3 paper.
See the steps below to re-produce the analysis.

## How to replicate the results

### Step 1: Downloading input files

Download Minigraph-Cactus VCF and GFA files:
```bat
wget -P inputs/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13.gfa.gz
wget -P inputs/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13.vcf.gz
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

Download external genotyped sets:

* NYGC Illumina-based SV calls:
```bat
wget -P inputs/ https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz
wget -p inputs/ https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz.tbi
```
* HPRC PanGenie genotypes (filtered set)
```bat
wget -P inputs/ https://zenodo.org/records/6797328/files/all-samples_bi_all.vcf.gz?download=1
wget -P inputs/ https://zenodo.org/records/6797328/files/bi_all_filters.tsv.gz?download=1

gunzip inputs/bi_all_filters.tsv.gz

wget -P inputs/ https://zenodo.org/records/6797328/files/select_ids.py?download=1

zcat all-samples_bi_all.vcf.gz | python3 inputs/select_ids.py inputs/bi_all_filters.tsv filtered | bgzip -c > inputs/all-samples_bi_filtered.vcf.gz

tabix -p vcf inputs/all-samples_bi_filtered.vcf.gz
```

### Step 2: run snakemake pipeline
Dependencies:
* snakemake
* singularity
* conda

```bat
snakemake -j <nr_cores> --use-conda
```

