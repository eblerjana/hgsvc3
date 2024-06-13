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
permuted_samples = [samples[i] for i in range(1, len(samples))] + [samples[0]]
randomized_samples = {samples[i] : permuted_samples[i] for i in range(len(samples))}

kmer_size = 21


rule compute_consensus:
	"""
	Insert all variants into the reference genome to produce haplotypes.
	"""
	input:
		vcf = lambda wildcards: PHASED_VCFS[wildcards.callset]["vcf"],
		reference = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"]
	output:
		fasta = temp("{results}/haplotypes/{callset}/{callset}_{sample}_hap{haplotype}_tmp.fasta.gz")
	log:
		"{results}/haplotypes/{callset}/{callset}_{sample}_hap{haplotype}.log"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		haplotype = "1|2"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 4
	shell:
		"""
		bcftools consensus --sample {wildcards.sample} --haplotype {wildcards.haplotype} -e 'ALT~\"<.*>\"' -f {input.reference} {input.vcf} 2> {log} | bgzip > {output.fasta} 
		"""


rule extract_chromosomes:
	"""
	Exclude chrY, since it is not present
	in all samples.
	"""
	input:
		"{results}/haplotypes/{callset}/{callset}_{sample}_hap{haplotype}_tmp.fasta.gz"
	output:
		"{results}/haplotypes/{callset}/{callset}_{sample}_hap{haplotype}.fasta.gz"
	conda:
		"../envs/shapeit.yaml"
	params:
		regions = lambda wildcards: " ".join(chroms[wildcards.callset]),
		haplotype = "1|2"
	shell:
		"""
		samtools faidx {input} {params.regions} | bgzip > {output}
		"""

rule meryl_counting:
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



rule merqury_evaluate_assemblies:
	"""
	Compute assembly statistics using merqury.
	"""
	input:
		assembly = "{results}/haplotypes/{callset}/{callset}_{sample}_hap{haplotype}.fasta.gz",
		counts = "{results}/read-counts-meryl/{sample}.meryl"
	output:
		qv = "{results}/evaluation/merqury/assigned/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}.qv",
		completeness = "{results}/evaluation/merqury/assigned/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}.completeness.stats"
	threads:
		24
	resources:
		mem_total_mb = 32768,
		time_hrs = 2,
		mem_total_gb = 32
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/evaluation/merqury/assigned/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}.log"
	params:
		out_prefix = "{results}/evaluation/merqury/assigned/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}"
	shell:
		"""
		./workflow/scripts/compute-qv.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_total_gb} &> {log}
		./workflow/scripts/compute-completeness.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_total_gb} 
		"""



rule merqury_evaluate_assemblies_randomized:
	"""
	Compute assembly statistics using merqury.
	"""
	input:
		assembly = "{results}/haplotypes/{callset}/{callset}_{sample}_hap{haplotype}.fasta.gz",
		counts = lambda wildcards: "{results}/read-counts-meryl/{s}.meryl".format(results=wildcards.results, s = randomized_samples[wildcards.sample])
	output:
		qv = "{results}/evaluation/merqury/randomized/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}.qv",
		completeness = "{results}/evaluation/merqury/randomized/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}.completeness.stats"
	threads:
		24
	resources:
		mem_total_mb = 32768,
		mem_total_gb = 32,
		time_hrs = 2
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/evaluation/merqury/randomized/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}.log"
	params:
		out_prefix = "{results}/evaluation/merqury/randomized/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}"
	shell:
		"""
		./workflow/scripts/compute-qv.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_total_gb} &> {log}
		./workflow/scripts/compute-completeness.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_total_gb} 
		"""



rule merqury_plot_result_with_assemblies:
	"""
	Plot merged + single QV values.
	"""
	input:
		computed_qvs = expand("{{results}}/evaluation/merqury/{mode}/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}.qv", mode = ["assigned", "randomized"], sample = QV_SAMPLES, callset = [c for c in PHASED_VCFS.keys()], haplotype = ["1", "2"]),
		given_qvs = ASSEMBLY_QVS
	output:
		"{results}/evaluation/qv-values_merqury.pdf"
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		ls {input.computed_qvs} | python3 workflow/scripts/plot-qv-merqury.py {output} -assembly {input.given_qvs}
		"""


rule merqury_plot_completeness:
	"""
	Plot merged + single completeness values.
	"""
	input:
		computed_qvs = expand("{{results}}/evaluation/merqury/{mode}/{callset}/{sample}_hap{haplotype}/{callset}_{sample}_hap{haplotype}.completeness.stats", mode = ["assigned", "randomized"], sample = QV_SAMPLES, callset = [c for c in PHASED_VCFS.keys()], haplotype = ["1", "2"])
	output:
		"{results}/evaluation/completeness_merqury.pdf"
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		ls {input.computed_qvs} | python3 workflow/scripts/plot-completeness-merqury.py {output}
		"""
