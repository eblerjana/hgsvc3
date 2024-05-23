rule concat_vcfs:
	"""
	Combine the phased per-chromosome VCFs into a single one.
	"""
	input:
		expand("{{results}}/{{version}}/phased_{{version}}_{chrom}.bcf", chrom = [c for c in MAPS.keys()])
	output:
		"{results}/{version}/phased_{version}.vcf.gz"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		version = "shapeit|shapeit-scaffold"
	benchmark:
		"{results}/{version}/phased_{version}_benchmark.txt"
	threads: 24
        resources:
		mem_total_mb = 100000,
		runtime_hrs = 5,
		runtime_min = 1
	log:
		"{results}/{version}/phased_{version}.log"
	shell:
		"""
		bcftools concat -o {output} -O z --threads {threads} {input} &> {log}
		tabix -p vcf {output}
		"""


rule extract_sample_phasing:
	"""
	Extract phasing of a single sample. This is done to reduce the
	time/memory of whatshap compare.
	"""
	input:
		"{results}/{version}/phased_{version}.vcf.gz"
	output:
		vcf = temp("{results}/{version}/phased_{version}_{sample}.vcf.gz"),
		tbi = temp("{results}/{version}/phased_{version}_{sample}.vcf.gz.tbi")
	conda:
		"../envs/shapeit.yaml"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 1
	params:
		chrom = ','.join([c for c in MAPS.keys() if (not 'X' in c) and (not 'Y' in c)])
	shell:
		"""
		bcftools view --samples {wildcards.sample} --regions {params.chrom} {input} | bcftools view --min-ac 1 | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule extract_sample_truthset:
	"""
	Extract ground truth phasing of a single sample. Done to reduce
	the time/memory of whatshap compare.
	"""
	input:
		lambda wildcards: TRUTHSETS[wildcards.truthset]["vcf"]
	output:
		vcf = temp("{results}/{truthset}_{sample}.vcf.gz"),
		tbi = temp("{results}/{truthset}_{sample}.vcf.gz.tbi")
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		truthsets = "|".join(TRUTHSETS.keys()),
		sample = "|".join(EVAL_SAMPLES)
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 2
	params:
		chrom = ','.join([c for c in MAPS.keys() if (not 'X' in c) and (not 'Y' in c)])
	shell:
		"""
		bcftools view --samples {wildcards.sample} --regions {params.chrom} {input} | bcftools view --min-ac 1 | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""
		

rule evaluate:
	"""
	Evaluate Phasing based on truth sets.
	"""
	input:
		truthsets = "{results}/{truthset}_{sample}.vcf.gz",
		truthsets_tbi = "{results}/{truthset}_{sample}.vcf.gz.tbi",
		phasing = "{results}/{version}/phased_{version}_{sample}.vcf.gz",
		phasing_tbi = "{results}/{version}/phased_{version}_{sample}.vcf.gz.tbi"
	output:
		pair = "{results}/evaluation/{version}/evaluation_{version}_{sample}_{truthset}_pair.tsv",
		multi = "{results}/evaluation/{version}/evaluation_{version}_{sample}_{truthset}_multi.tsv"
	conda:
		"../envs/whatshap.yaml"
	wildcard_constraints:
		truthset = "|".join(TRUTHSETS.keys()),
		sample = "|".join(EVAL_SAMPLES)
	resources:
		mem_total_mb = 100000,
		runtime_hrs = 3
	params:
		names = lambda wildcards: wildcards.version + "," + wildcards.truthset
	log:
		"{results}/evaluation/{version}/evaluation_{version}_{sample}_{truthset}.log"
	shell:
		"""
		whatshap compare {input.phasing} {input.truthsets} --sample {wildcards.sample} --tsv-pairwise {output.pair} --tsv-multiway {output.multi} --names {params.names} &> {log}
		"""


rule plot_phasing_results:
	"""
	Plot the switch error rates.
	"""
	input:
		lambda wildcards: expand("{{results}}/evaluation/{{version}}/evaluation_{{version}}_{sample}_{{truthset}}_pair.tsv", sample = TRUTHSETS[wildcards.truthset]["evaluation_samples"])
	output:
		"{results}/evaluation/{version}/evaluation_{version}_{truthset}.pdf"
	wildcard_constraints:
		truthset='|'.join([c for c in TRUTHSETS.keys()])
	conda:
		"../envs/plotting.yml"
	shell:
		"""
		python3 workflow/scripts/plot-phasing-results.py -tsvfiles {input} -truthsetname {wildcards.truthset} -outname {output}
		"""


rule multiway_evaluation:
	input:
		truthsets = expand("{{results}}/{truthset}_{{sample}}.vcf.gz", truthset = [k for k in config["truthsets"].keys()]),
		truthsets_tbi = expand("{{results}}/{truthset}_{{sample}}.vcf.gz.tbi", truthset = [k for k in config["truthsets"].keys()]),
		phasing = "{results}/{version}/phased_{version}_{sample}.vcf.gz",
		phasing_tbi = "{results}/{version}/phased_{version}_{sample}.vcf.gz.tbi"
	output:
		pair = "{results}/evaluation/{version}/evaluation_{version}_{sample}_multiway_pair.tsv",
		multi = "{results}/evaluation/{version}/evaluation_{version}_{sample}_multiway_multi.tsv"
	conda:
		"../envs/whatshap.yaml"
	wildcard_constraints:
		sample = "|".join(EVAL_SAMPLES)
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 3
	params:
		names = "shapeit-phasing," + ",".join([k for k in TRUTHSETS.keys()])
	log:
		"{results}/evaluation/evaluation_{version}_{sample}_multiway.log"
	shell:
		"""
		whatshap compare {input.phasing} {input.truthsets} --sample {wildcards.sample} --tsv-pairwise {output.pair} --tsv-multiway {output.multi} --names {params.names} &> {log}
		"""

rule plot_multiway:
	input:
		tsv=expand("{{results}}/evaluation/{{version}}/evaluation_{{version}}_{sample}_multiway_multi.tsv", sample=JOINT_EVAL_SAMPLES),
		trios=FAM
	output:
		"{results}/evaluation/{version}/evaluation_{version}_multiway.pdf"
	conda:
		"../envs/plotting.yml"
	shell:
		"""
		python3 workflow/scripts/plot-multiway.py -tsv {input.tsv} -trios {input.trios} -outname {output} 
		"""
