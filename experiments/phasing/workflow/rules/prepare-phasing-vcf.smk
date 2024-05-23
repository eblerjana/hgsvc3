
rule prepare_convert_vcf_to_bi:
	"""
	Convert VCF to biallelic representation and restrict to specific chromosome.
	"""
	input:
		lambda wildcards: GENOTYPES if wildcards.set == "genotyping" else PHASED_SET[wildcards.chrom]
	output:
		"{results}/vcfs/{set}_{chrom}.vcf.gz"
	log:
		"{results}/vcfs/{set}_{chrom}.log"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		set = "phasing|genotyping"
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} | bcftools norm -m- -Oz -o {output} &> {log}
		tabix -p vcf {output}
		"""	


rule prepare_get_panel_samples:
	"""
	Get list of panel samples from which variants for genotyping were called.
	"""
	input: 
		PANEL_VCF
	output:
		"{results}/panel-samples.tsv"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/panel-samples.log"
	shell:
		"""
		bcftools query -l {input} 2> {log} 1> {output}
		"""


rule prepare_get_samples_list:
	"""
	Get list of all samples present in the phased VCF.
	"""
	input:
		lambda wildcards: "{results}/vcfs/phasing_{chrom}.vcf.gz".format(results=wildcards.results, chrom = list(PHASED_SET.keys())[0])
	output:
		"{results}/samples.tsv"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/samples.log"
	shell:
		"""
		bcftools query -l {input} 1> {output} 2> {log}
		"""


rule prepare_extract_genotyped_samples:
	"""
	Restrict genotyped VCF to samples present in phased VCF.
	Order the remaining samples in same way.
	"""
	input:
		vcf = "{results}/vcfs/genotyping_{chrom}.vcf.gz",
		samples = "{results}/samples.tsv"
	output:
		"{results}/genotyping-intersected_{chrom}.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/genotyping-intersected_{chrom}.log"
	resources:
	shell:
		"""
		bcftools view -S {input.samples} {input.vcf} -Oz -o {output} &> {log}
		tabix -p vcf {output}
		"""


rule prepare_extract_rare_snp_ids:
	"""
	In the phased VCF, extract all SNPs with AF=0 across the
	panel samples (= rare variants to be merged into the genotyped set)
	"""
	input: 
		vcf = "{results}/vcfs/phasing_{chrom}.vcf.gz",
		samples = "{results}/panel-samples.tsv"
	output: "{results}/rare-snps_{chrom}.tsv"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/rare-snps_{chrom}.log"
	shell:
		"""
		bcftools view -S {input.samples} {input.vcf} --force-samples 2> {log} | bcftools view --max-ac 0 -H | awk '(length($4) == 1) && (length($5) == 1)' | cut -f1,2,4,5  1> {output}
		"""


rule filter_phasing_vcfs:
	"""
	Output subset of phased and genotyped VCFs to be
	combined later.
	"""
	input:
		phasing = "{results}/vcfs/phasing_{chrom}.vcf.gz",
		genotyped = "{results}/genotyping-intersected_{chrom}.vcf.gz",
		rare_snps = "{results}/rare-snps_{chrom}.tsv"
	output:
		tmp_genotypes = "{results}/filtered_genotypes_{chrom}.vcf",
		tmp_phasing = "{results}/filtered_phasing_{chrom}.vcf",
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/filter_phasing_{chrom}.log"
	resources:
	shell:
		"""
		python3 workflow/scripts/prepare-phasing-vcfs.py {input.phasing} {input.genotyped} {input.rare_snps} {output.tmp_genotypes} {output.tmp_phasing} &> {log}
		"""


rule prepare_add_tags:
	"""
	Add AC/AN tags to VCFs
	"""
	input:
		"{results}/filtered_{set}_{chrom}.vcf"
	output:
		"{results}/filtered_{set}_{chrom}.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/filtered_{set}_{chrom}.log"
	wildcard_constraints:
		set = "phasing|genotypes"
	shell:
		"""
		bcftools +fill-tags {input} -Oz -o {output}  -- -t AN,AC,AF &> {log}
		tabix -p vcf {output}
		"""



rule prepare_combine_vcfs:
	"""
	Combine both VCFs to obtain final phasing input.
	"""
	input:
		genotypes = "{results}/filtered_genotypes_{chrom}.vcf.gz",
		phasing = "{results}/filtered_phasing_{chrom}.vcf.gz"
	output:
		"{results}/phasing-input_{chrom}.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/phasing-input_{chrom}.log"
	resources:
	shell:
		"""
		bcftools concat -a {input.genotypes} {input.phasing} -Oz -o {output} &> {log}
		tabix -p vcf {output}
		"""	

