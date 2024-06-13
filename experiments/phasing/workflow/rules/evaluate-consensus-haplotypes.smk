configfile: "config/config.yaml"

from random import shuffle

sample_to_reads = {}
for line in open(READS, 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	sample_to_reads[fields[1]] = fields[7]

chroms = []
for line in open(REFERENCE, 'r'):
	if line.startswith('>'):
		fields = line.strip().split()
		chrom = fields[0][1:]
		if not chrom in ["Y", "chrY", "M", "chrM"]:
			chroms.append(chrom)

chroms = {}
for callset in PHASED_VCFS:
	chroms[callset] = []
	for line in open(PHASED_VCFS[callset]["reference"]):
		if line.startswith('>'):
			fields = line.strip().split()
			chrom = fields[0][1:]
			if not chrom in ["Y", "chrY", "M", "chrM"] and (not "random" in chrom) and not ("chrUn" in chrom):
				chroms[callset].append(chrom)

samples = [s for s in QV_SAMPLES]
shuffle(samples)
permuted_samples = [samples[i] for i in range(1, len(samples))] + [samples[0]]
randomized_samples = {samples[i] : permuted_samples[i] for i in range(len(samples))}

kmer_size = 21





rule meryl_counting_reads:
	"""
	Counting kmers using meryl
	"""
	input:
		lambda wildcards: sample_to_reads[wildcards.sample]
	output:
		directory("{results}/read-counts-meryl/{sample}.meryl")
	threads:
		24
	resources:
		mem_total_mb = 32768,
		mem_total_gb = 32,
		time_hrs = 2
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/read-counts-meryl/{sample}.log"
	shell:
		"""
		meryl count k={kmer_size} memory={resources.mem_total_gb} threads={threads} {input} output {output} &> {log}
		"""


rule meryl_counting_reads_assemblies:
	"""
	Counting kmers using meryl
	"""
	input:
		"{results}/haplotypes/{callset}/{callset}_{sample}_hap{haplotype}.fasta.gz"
	output:
		directory("{results}/read-counts-assemblies/{callset}/{callset}_{sample}_hap{haplotype}.meryl")
	threads:
		24
	resources:
		mem_total_mb = 32768,
		mem_total_gb = 32,
		time_hrs = 2
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/read-counts-assemblies/{callset}/{callset}_{sample}_hap{haplotype}.log"
	shell:
		"""
		meryl count k={kmer_size} memory={resources.mem_total_gb} threads={threads} {input} output {output} &> {log}
		"""





rule merqury_evaluate_assemblies:
	"""
	Compute assembly statistics using merqury.
	"""
	input:
		assemblies=expand("{{results}}/haplotypes/{{callset}}/{{callset}}_{{sample}}_hap{haplotype}.fasta.gz", haplotype = [1,2]),
		counts = "{results}/read-counts-meryl/{sample}.meryl"
	output:
		qv = "{results}/evaluation/merqury/assigned/{callset}/{sample}/{callset}_{sample}.qv",
		completeness = "{results}/evaluation/merqury/assigned/{callset}/{sample}/{callset}_{sample}.completeness.stats"
	threads:
		24
	resources:
		mem_total_mb = 50000,
		time_hrs = 2
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/evaluation/merqury/assigned/{callset}/{sample}/{callset}_{sample}.log"
	params:
		out_prefix = "{callset}_{sample}",
		outdir = "{results}/evaluation/merqury/assigned/{callset}/{sample}/"
	shell:
		"""
		cd {params.outdir} && merqury.sh {RUN_DIR}{input.counts} {RUN_DIR}{input.assemblies} {params.out_prefix} &> {RUN_DIR}{log}
		"""



rule merqury_evaluate_assemblies_randomized:
	"""
	Compute assembly statistics using merqury with randomized samples.
	"""
	input:
		assemblies=expand("{{results}}/haplotypes/{{callset}}/{{callset}}_{{sample}}_hap{haplotype}.fasta.gz", haplotype = [1,2]),
		counts = lambda wildcards: "{results}/read-counts-meryl/{s}.meryl".format(results=wildcards.results, s = randomized_samples[wildcards.sample])
	output:
		qv = "{results}/evaluation/merqury/randomized/{callset}/{sample}/{callset}_{sample}.qv",
		completeness = "{results}/evaluation/merqury/randomized/{callset}/{sample}/{callset}_{sample}.completeness.stats"
	threads:
		24
	resources:
		mem_total_mb = 50000,
		time_hrs = 2
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/evaluation/merqury/randomized/{callset}/{sample}/{callset}_{sample}.log"
	params:
		out_prefix = "{callset}_{sample}",
		outdir = "{results}/evaluation/merqury/randomized/{callset}/{sample}/"
	shell:
		"""
		cd {params.outdir} && merqury.sh {RUN_DIR}{input.counts} {RUN_DIR}{input.assemblies} {params.out_prefix} &> {RUN_DIR}{log}
		"""


rule merqury_plot_result_with_assemblies:
	"""
	Plot merged + single QV values.
	"""
	input:
		computed_qvs = expand("{{results}}/evaluation/merqury/{mode}/{callset}/{sample}/{callset}_{sample}.qv", mode = ["assigned", "randomized"], sample = QV_SAMPLES, callset = [c for c in PHASED_VCFS.keys()]),
		given_qvs = ASSEMBLY_QVS
	output:
		"{results}/evaluation/qv-values_merqury.pdf"
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		ls {input.computed_qvs} | python3 workflow/scripts/plot-qv-merqury.py {output} -assembly {input.given_qvs}
		"""
