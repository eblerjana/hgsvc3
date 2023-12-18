# Benchmarking Experiments


## Pipeline

We used our benchmarking pipeline at: https://github.com/eblerjana/genotyping-pipelines/tree/main/benchmarking-pipeline (commit: 4ff6af2).      
Config file used: `` config.yaml``  (stored in this folder)
 * **Biallelic VCF**: `` pav_variants_batch1_alt.vcf.gz`` , produced from PAV calls as described below
 * **Mulitallelic VCF**: `` pangenome.vcf.gz`` , produced from PAV calls as described in folder `` ../PanGenie-vcf-preparation`` 
 * **reference genome**: downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz
 * **Pilot samples / reads**: `` genotyping-pilot-reads.tsv``  (stored in this folder). "1000GP high-coverage cohorts" reads downloaded from EBI/ENA and merged into one file per sample.
 * **PanGenie**: version v3.0.1 (commit: e981195)


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
