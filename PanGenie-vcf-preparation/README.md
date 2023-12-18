# Preparation of PanGenie input VCF


## Pipeline

We used our pipeline at: https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-PAV (commit: d9ff614).   
Config file used: `` config.yaml `` (stored in this folder)
 * **VCF**: `` pav_variants_batch1_alt.vcf.gz`` , produced from PAV calls as described below
 * **reference genome**: downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz
 * **sample info**: `` sample-info.tsv``  (stored in this folder)

## PAV data

PAV VCFs were downloaded from:

https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_snv_snv_alt.vcf.gz
https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_snv_snv_alt.vcf.gz.tbi

https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_indel_insdel_alt.vcf.gz
https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_indel_insdel_alt.vcf.gz.tbi

https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_sv_insdel_alt.vcf.gz
https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_sv_insdel_alt.vcf.gz.tbi

and merged into a single VCF using the command:

``` bat
bcftools concat -a <files> | bgzip > pav_variants_batch1_alt.vcf.gz
```
