configfile: "config-full.yaml"

# awk '{print $2 "\t" $4 "\t" $3}' {input} | awk -F '\t' '{OFS="\t"; if ($2 == "0"){$2="NA"}; OFS="\t"; if ($3 == "0"){$3="NA"}; print}' | grep -v "#"

chromosomes = [c for c in config['maps'].keys()]
outname = config['outname']

autosomes = [c for c in chromosomes if c not in ['X', 'Y', 'chrX', 'chrY']]

samples = set([])
for callset in config["truthsets"]:
	for sample in config["truthsets"][callset]["evaluation_samples"]:
		samples.add(sample)

samples = list(samples)

joint_samples = None
for callset in config["truthsets"]:
	if not joint_samples:
		joint_samples = set(config["truthsets"][callset]["evaluation_samples"])
	else:
		joint_samples.intersection(set(config["truthsets"][callset]["evaluation_samples"]))

print(joint_samples)


def output_pair_plots(wildcards):
	output = []
	for t in config['truthsets'].keys():
		for sample in config['truthsets'][t]['evaluation_samples']:
			output.append("{results}/evaluation/evaluation_{sample}_{truthset}_pair.tsv".format(results=outname, sample=sample, truthset=t))
	return output


rule all:
	input:
		expand("{results}/shapeit-phasing.vcf.gz", results=outname),
#		expand("{results}/evaluation/evaluation_{sample}_{truthset}_pair.tsv", results=outname, sample=config["evaluation_samples"], truthset = config['truthsets'].keys()),
		output_pair_plots,
		expand("{results}/evaluation/evaluation_{truthset}.pdf", results=outname, truthset = config['truthsets'].keys()),
		expand("{results}/evaluation/evaluation_{sample}_multiway_multi.tsv", results=outname, sample = ["HG00418", "HG00419", "HG00420", "HG01256", "HG01257", "HG01258", "NA19127", "NA19128", "NA19129", "NA19818", "NA19819", "NA19828"]),
		expand("{results}/evaluation/evaluation_multiway.pdf", results=outname)


rule set_low_qual_to_missing:
	"""
	Set low quality genotypes to missing.
	"""
	input:
		config['vcf']
	output:
		temp("{results}/vcf/low-qual-missing.vcf.gz")
	conda:
		"shapeit.yaml"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 5
	shell:
		"""
		bcftools +setGT {input}  -- -t q -n ./. -i 'FMT/GQ<20' | bgzip > {output}
		tabix -p vcf {output}
		"""


rule shapeit_phase_common:
	"""
	Phase a chromosome using shapeit.
	"""
	input:
		vcf = "{results}/vcf/low-qual-missing.vcf.gz" if config['low_qual_to_missing'] else config['vcf'],
		fam = config['fam'],
		map = lambda wildcards: config['maps'][wildcards.chrom]
	output:
		"{results}/shapeit-{chrom}.bcf"
	log:
		"{results}/shapeit-{chrom}.log"
	benchmark:
		"{results}/shapeit-{chrom}-benchmark.txt"
	conda:
		"shapeit.yaml"
	threads: 24
        resources:
		mem_total_mb = 100000,
		runtime_hrs = 20
	shell:
		"""
		SHAPEIT5_phase_common --input {input.vcf} --pedigree {input.fam} --region {wildcards.chrom} --map {input.map} --output {output} --thread {threads} &> {log}
		"""

rule concat_vcfs:
	"""
	Combine the phased per-chromosome VCFs into a single one.
	"""
	input:
		expand("{{results}}/shapeit-{chrom}.bcf", chrom = chromosomes)
	output:
		"{results}/shapeit-phasing.vcf.gz"
	conda:
		"shapeit.yaml"
	benchmark:
		"{results}/shapeit-phasing-benchmark.txt"
	threads: 24
        resources:
		mem_total_mb = 100000,
		runtime_hrs = 5,
		runtime_min = 1
	log:
		"{results}/shapeit-phasing.log"
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
		"{results}/shapeit-phasing.vcf.gz"
	output:
		vcf = temp("{results}/shapeit-phasing_{sample}.vcf.gz"),
		tbi = temp("{results}/shapeit-phasing_{sample}.vcf.gz.tbi")
	conda:
		"shapeit.yaml"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 1
	params:
		chrom = ','.join(autosomes)
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
		lambda wildcards: config["truthsets"][wildcards.truthset]["vcf"]
	output:
		vcf = temp("{results}/{truthset}_{sample}.vcf.gz"),
		tbi = temp("{results}/{truthset}_{sample}.vcf.gz.tbi")
	conda:
		"shapeit.yaml"
	wildcard_constraints:
		truthsets = "|".join(config["truthsets"].keys()),
		sample = "|".join(samples)
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 2
	params:
		chrom = ','.join(autosomes)
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
		phasing = "{results}/shapeit-phasing_{sample}.vcf.gz",
		phasing_tbi = "{results}/shapeit-phasing_{sample}.vcf.gz.tbi"
	output:
		pair = "{results}/evaluation/evaluation_{sample}_{truthset}_pair.tsv",
		multi = "{results}/evaluation/evaluation_{sample}_{truthset}_multi.tsv"
	conda:
		"whatshap.yaml"
	wildcard_constraints:
		truthset = "|".join(config["truthsets"].keys()),
		sample = "|".join(samples)
	resources:
		mem_total_mb = 100000,
		runtime_hrs = 3
	params:
		names = lambda wildcards: "shapeit-phasing," + wildcards.truthset
	log:
		"{results}/evaluation/evaluation_{sample}_{truthset}.log"
	shell:
		"""
		whatshap compare {input.phasing} {input.truthsets} --sample {wildcards.sample} --tsv-pairwise {output.pair} --tsv-multiway {output.multi} --names {params.names} &> {log}
		"""


rule plot_phasing_results:
	"""
	Plot the switch error rates.
	"""
	input:
		lambda wildcards: expand("{{results}}/evaluation/evaluation_{sample}_{{truthset}}_pair.tsv", sample = config["truthsets"][wildcards.truthset]["evaluation_samples"])
	output:
		"{results}/evaluation/evaluation_{truthset}.pdf"
	wildcard_constraints:
		truthset='|'.join([c for c in config["truthsets"].keys()])
	conda:
		"plotting.yml"
	shell:
		"""
		python3 plot-phasing-results.py -tsvfiles {input} -truthsetname {wildcards.truthset} -outname {output}
		"""


rule multiway_evaluation:
	input:
		truthsets = expand("{{results}}/{truthset}_{{sample}}.vcf.gz", truthset = [k for k in config["truthsets"].keys()]),
		truthsets_tbi = expand("{{results}}/{truthset}_{{sample}}.vcf.gz.tbi", truthset = [k for k in config["truthsets"].keys()]),
		phasing = "{results}/shapeit-phasing_{sample}.vcf.gz",
		phasing_tbi = "{results}/shapeit-phasing_{sample}.vcf.gz.tbi"
	output:
		pair = "{results}/evaluation/evaluation_{sample}_multiway_pair.tsv",
		multi = "{results}/evaluation/evaluation_{sample}_multiway_multi.tsv"
	conda:
		"whatshap.yaml"
	wildcard_constraints:
		sample = "|".join(samples)
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 3
	params:
		names = "shapeit-phasing," + ",".join([k for k in config["truthsets"].keys()])
	log:
		"{results}/evaluation/evaluation_{sample}_multiway.log"
	shell:
		"""
		whatshap compare {input.phasing} {input.truthsets} --sample {wildcards.sample} --tsv-pairwise {output.pair} --tsv-multiway {output.multi} --names {params.names} &> {log}
		"""

rule plot_multiway:
	input:
		tsv=expand("{{results}}/evaluation/evaluation_{sample}_multiway_multi.tsv", sample=joint_samples),
		trios=config['fam']
	output:
		"{results}/evaluation/evaluation_multiway.pdf"
	conda:
		"plotting.yml"
	shell:
		"""
		python3 plot-multiway.py -tsv {input.tsv} -trios {input.trios} -outname {output} 
		"""
