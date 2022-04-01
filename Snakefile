# Hufsah , 12.3.2021


import math
from collections import defaultdict

path_to_bams = "path to all bam files"
path_to_snps_1toX="path to folder containing vcf file with SNPs belonging to autosomes and chromosome X"
path_to_snps_Y= "path to folder containing vcf file with SNPs belonging to Y chromosome"
path_to_finals ="path to folder containing for each sample, the mosaicatcher output files reporting background cell states"
SAMPLE = (glob_wildcards("path to all bam files/{sample}/all/"))
S2 = tuple(i for i in SAMPLE[0])
SAMPLES = sorted(set(S2))

    
print("Detected {} samples:".format(len(SAMPLES)))


rule all:
    input:
        expand('output/{sample}/snp_strand_counts.txt', sample = SAMPLES),
        expand('per_sample_configs/{sample}_snp_ref_inv.txt', sample = SAMPLES),
         expand('per_sample_configs_clean/{sample}_snp_ref_inv.txt', sample = SAMPLES),
        'one_table_per_snp.txt'


rule snp_strand_counts_1toX:
            input:
                bam = expand("{path}/bam/{{sample}}/all/", path= path_to_bams),
                bed = "invs_1toX.bed",
                vcf = expand("{path}/concat_chrs_44samples_bi_reheader.vcf.gz", path= path_to_snps_1toX)
            output:
                strand_counts="output/{sample}/snp_strand_counts_1toX.txt",
                other_info= "other/{sample}/read_info_1toX.txt"

            shell:
                "python3 recurrence_first_v2.py -i {input.bam} -b {input.bed} -o {output.strand_counts} -v {input.vcf} -s {wildcards.sample} > {output.other_info} "
                

rule snp_strand_counts_Y:
            input:
                bam = expand("{path}/bam/{{sample}}/all/", path= path_to_bams),
                bed = "invs_Y.bed",
                vcf = expand("{path}/chrY_bi.vcf.gz", path= path_to_snps_Y)
            output:
                strand_counts="output/{sample}/snp_strand_counts_Y.txt",
                other_info= "other/{sample}/read_info_Y.txt"

            shell:
                "python3 recurrence_first_v2.py -i {input.bam} -b {input.bed} -o {output.strand_counts} -v {input.vcf} -s {wildcards.sample} > {output.other_info} "

rule snp_strand_counts:
            input:
                strand_counts_Y= "output/{sample}/snp_strand_counts_Y.txt",
                strand_counts_1toX = "output/{sample}/snp_strand_counts_1toX.txt",
            output:
                strand_counts= "output/{sample}/snp_strand_counts.txt"

            shell:
                "awk 'FNR>1 || NR==1' {input.strand_counts_1toX} {input.strand_counts_Y} > {output.strand_counts}"


               
rule add_cell_states:
            input:
            	states=expand("{path}/{{sample}}/100000_fixed_norm.selected_j0.1_s0.1/final.txt", path= path_to_finals),
            	snp_counts_file = "output/{sample}/snp_strand_counts.txt",

   	    output:
   	    	out = 'per_sample_configs/{sample}_snp_ref_inv.txt'
   	    shell:
   	    	"""
   	    	Rscript cell_states_from_final.R \
        	-f {input.states} \
        	-b {input.snp_counts_file} \
        	-o {output.out}\
        	"""
rule clean_configs:
            input:
            	configs='per_sample_configs/{sample}_snp_ref_inv.txt'

   	    output:
   	    	out = 'per_sample_configs_clean/{sample}_snp_ref_inv.txt'
   	    shell:
   	    	"""
   	    	Rscript clean_reads.R
        	"""
rule one_table:
            input:
            	configs=expand("per_sample_configs_clean/{sample}_snp_ref_inv.txt",  sample = SAMPLES),

   	    output:
   	    	out = 'one_table_per_snp.txt'
   	    shell:
   	    	"""
   	    	Rscript 1snp_1table_allsamples.R
        	"""




