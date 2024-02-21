# Benchmarking Experiments


## Pipeline


We used our benchmarking pipeline at: https://github.com/eblerjana/genotyping-pipelines/tree/main/benchmarking-pipeline (commit: 48e4d2e).      

Config file used: `` config.yaml``  (stored in this folder)

### PAV-hgsvc3-batch1-hg38 (PAV, Batch1 samples)

 * **Biallelic VCF**: `` pav_variants_batch1_alt.vcf.gz`` , produced from PAV calls as described below
 * **Mulitallelic VCF**: `` pangenome.vcf.gz`` , produced from PAV calls as described in folder `` ../PanGenie-vcf-preparation`` 
 * **reference genome**: downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz
 * **Pilot samples / reads**: `` genotyping-pilot-reads.tsv``  (stored in this folder). "1000GP high-coverage cohorts" reads downloaded from EBI/ENA and merged into one file per sample.
 * **PanGenie**: version v3.0.2 (commit:  7c171f3)

### PAV-hgsvc3-freeze1-hg38 (PAV, Freeze1 samples)

 * **Biallelic VCF**: `` pav_variants_freeze1_alt.vcf.gz`` , produced from PAV calls as described below
 * **Mulitallelic VCF**: `` pangenome.vcf.gz`` , produced from PAV calls as described in folder `` ../PanGenie-vcf-preparation`` 
 * **reference genome**: downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz
 * **Pilot samples / reads**: `` genotyping-pilot-reads.tsv``  (stored in this folder). "1000GP high-coverage cohorts" reads downloaded from EBI/ENA and merged into one file per sample.
 * **PanGenie**: version v3.0.2 (commit: 7c171f3)

### MC-hgsvc3-freeze1-chm13 (MC, Freeze1 samples)


 * **Biallelic VCF**: `` MC-hgsvc3-freeze1-chm13_filtered_ids_biallelic.vcf.gz`` , produced from MC calls as described in folder `` ../MC-vcf-preparation ``
 * **Mulitallelic VCF**: `` MC-hgsvc3-freeze1-chm13_filtered_ids.vcf.gz`` , produced from MC calls as described in folder `` ../MC-vcf-preparation ``
 * **reference genome**: downloaded from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
 * **Pilot samples / reads**: `` genotyping-pilot-reads.tsv``  (stored in this folder). "1000GP high-coverage cohorts" reads downloaded from EBI/ENA and merged into one file per sample.
 * **PanGenie**: version v3.0.2 (commit: 7c171f3)


### MC-hgsvc3-freeze1-chm13 (MC, Freeze1 + HPRC samples)


 * **Biallelic VCF**: `` MC-hgsvc3-freeze1-hprc-chm13_filtered_ids_biallelic.vcf.gz`` , produced from MC calls as described in folder `` ../MC-vcf-preparation ``
 * **Mulitallelic VCF**: `` MC-hgsvc3-freeze1-hprc-chm13_filtered_ids.vcf.gz`` , produced from MC calls as described in folder `` ../MC-vcf-preparation ``
 * **reference genome**: downloaded from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
 * **Pilot samples / reads**: `` genotyping-pilot-reads.tsv``  (stored in this folder). "1000GP high-coverage cohorts" reads downloaded from EBI/ENA and merged into one file per sample.
 * **PanGenie**: version v3.0.2 (commit: 7c171f3)


## PAV data


### Batch1

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

### Freeze1

PAV VCFs downloaded from:

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20231211_Freeze1/HGSVC3_GRCh38_Freeze1/variants_freeze1_snv_snv_alt.vcf.gz
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20231211_Freeze1/HGSVC3_GRCh38_Freeze1/variants_freeze1_snv_snv_alt.vcf.gz.tbi

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20231211_Freeze1/HGSVC3_GRCh38_Freeze1/variants_freeze1_indel_insdel_alt.vcf.gz
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20231211_Freeze1/HGSVC3_GRCh38_Freeze1/variants_freeze1_indel_insdel_alt.vcf.gz.tbi

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20231211_Freeze1/HGSVC3_GRCh38_Freeze1/variants_freeze1_sv_insdel_alt.vcf.gz
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20231211_Freeze1/HGSVC3_GRCh38_Freeze1/variants_freeze1_sv_insdel_alt.vcf.gz.tbi

and merged into a single VCF using the command:

``` bat
bcftools concat -a <files> | bgzip > pav_variants_freeze1_alt.vcf.gz
```

