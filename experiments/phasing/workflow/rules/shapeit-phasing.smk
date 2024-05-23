
rule shapeit_extract_chromosome:
	"""
	Extract specific chromosome.
	"""
	input:
		GENOTYPES
	output:
		"{results}/vcf/genotypes_{chrom}.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/vcf/genotypes_{chrom}.log"
	shell:
		"""
		bcftools view -r {wildcards.chrom} -Oz -o {output} &> {log}
		tabix -p vcf {output}
		"""
	

rule shapeit_set_low_qual_to_missing:
	"""
	Set low quality genotypes to missing.
	"""
	input:
		 lambda wildcards: MERGED_GENOTYPES[wildcards.chrom] if MERGED_GENOTYPES else "{results}/vcf/genotypes_{chrom}.vcf.gz" 
	output:
		temp("{results}/vcf/low-qual-missing_{chrom}.vcf.gz")
	log:
		"{results}/vcf/low-qual-missing_{chrom}.log"
	conda:
		"../envs/shapeit.yaml"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 5
	shell:
		"""
		bcftools +setGT {input}  -- -t q -n ./. -i 'FMT/GQ<20' 2> {log} | bgzip > {output}
		tabix -p vcf {output}
		"""


rule shapeit_phase_common:
	"""
	Phase a chromosome using shapeit.
	"""
	input:
		vcf = "{results}/vcf/low-qual-missing_{chrom}.vcf.gz" if config['low_qual_to_missing'] else config['vcf'],
		fam = FAM,
		map = lambda wildcards: MAPS[wildcards.chrom]
	output:
		"{results}/shapeit/phased_shapeit_{chrom}.bcf"
	log:
		"{results}/shapeit/phased_shapeit_{chrom}.log"
	benchmark:
		"{results}/shapeit/phased_shapeit_{chrom}-benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	threads: 24
        resources:
		mem_total_mb = 100000,
		runtime_hrs = 13
	shell:
		"""
		SHAPEIT5_phase_common --input {input.vcf} --pedigree {input.fam} --region {wildcards.chrom} --map {input.map} --output {output} --thread {threads} &> {log}
		"""


rule shapeit_phase_common_scaffold:
	"""
	Phase a chromosome using shapeit.
	"""
	input:
		vcf = "{results}/vcf/low-qual-missing_{chrom}.vcf.gz" if config['low_qual_to_missing'] else config['vcf'],
		fam = FAM,
		map = lambda wildcards: MAPS[wildcards.chrom],
		scaffold = lambda wildcards: SCAFFOLDS[wildcards.chrom]
	output:
		"{results}/shapeit-scaffold/phased_shapeit-scaffold_{chrom}.bcf"
	log:
		"{results}/shapeit-scaffold/phased_shapeit-scaffold_{chrom}.log"
	benchmark:
		"{results}/shapeit-scaffold/phased_shapeit-scaffold_{chrom}-benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	threads: 24
        resources:
		mem_total_mb = 100000,
		runtime_hrs = 13
	shell:
		"""
		SHAPEIT5_phase_common --input {input.vcf} --scaffold {input.scaffold} --pedigree {input.fam} --region {wildcards.chrom} --map {input.map} --output {output} --thread {threads} &> {log}
		"""
