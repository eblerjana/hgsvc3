
NYGC consensus haplotypes:

cp /gpfs/project/projects/medbioinf/users/ebler/hgsvc3/experiments/phasing-1kg/results-with-rare/consensus-haplotypes/haplotypes/NYGC-GRCh38/all/NYGC-GRCh38_qv_consensus-haplotypes.agc .  
cp /gpfs/project/projects/medbioinf/users/ebler/hgsvc3/experiments/phasing-1kg/results-with-rare/consensus-haplotypes/haplotypes/NYGC-GRCh38/all/NYGC-GRCh38_all_summary.tsv .  
agc append NYGC-GRCh38_qv_consensus-haplotypes.agc /gpfs/project/projects/medbioinf/users/ebler/hgsvc3/data/reference/chm13v2.0.fa hg002v1.1.fasta.gz > NYGC-GRCh38_qv_consensus-haplotypes_chm13_HG002.agc  
echo -e 'reference\tchm13v2.0\tchm13v2.0\tnan\tnan\nHG002\thg002v1.1.pat\thg002v1.1.mat\tnan\tnan' >> NYGC-GRCh38_all_summary.tsv  

HG002 Q100 assembly downloaded via globus from:

/human-pangenomics/T2T/HG002/assemblies/
hg002v1.1.mat.fasta.gz
hg002v1.1.pat.fasta.gz


HPRC2 assemblies downloaded via globus from:

/human-pangenomics/submissions/6807247E-4F71-45D8-AECE-9E5813BA1D9F--verkko-v2.2.1-release2_asms/HG01255/verkko-hi-c/
HG01255.assembly.haplotype1.fasta.gz
HG01255.assembly.haplotype2.fasta.gz

/human-pangenomics/submissions/6807247E-4F71-45D8-AECE-9E5813BA1D9F--verkko-v2.2.1-release2_asms/HG04157/verkko-hi-c/
HG04157.assembly.haplotype1.fasta.gz
HG04157.assembly.haplotype2.fasta.gz
