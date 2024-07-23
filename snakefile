REP = ["A1"]
SAMPLE_SIZE = ["M-1000", "M-500"]
THETA = ["t-0.1", "t-0.01"]
CHR = []
for i in range(1,101):
  CHR.append(str(i))
NSNP = ["L-500","L-1000", "L-10000"]

# Helper function
def get_params(x):
  out = x.split("-")[1]
  return out

#################
## Target Rule ##
#################

rule all:
    input:
        expand("results/all.txt", rep = REP, sample_size=SAMPLE_SIZE, theta=THETA, chr=CHR, nsnp=NSNP)


##########################
### Simluate Genotypes ###
##########################


rule simulate_genotypes:
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos-{chr}.vcf",
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos-{chr}.pop"
    params:
        M = lambda wildcards: get_params(wildcards.sample_size),
        theta = lambda wildcards: get_params(wildcards.theta),
        output_prefix = "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos-{chr}"
    shell:
        """
        python -u code/Simulate_Genotypes/generate_genotypes.py \
	      --outpre {params.output_prefix} \
       	--Nanc 1000 \
	      --NA 1000 \
	      --NB 1000 \
	      --NC 1000 \
	      --chr {wildcards.chr} \
  	    -M {params.M} \
  	    -t {params.theta} \
        -s1 2000 \
        -s2 1000 \
        -L 1 \
        -u 1e-4  \
        -nsnp 600
        rm -f {params.output_prefix}_*
        """

############################
### Convert VCF to plink ###
############################

rule format_VCF:
    input:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos-{chr}.vcf"
    output:
        gz="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos-{chr}.tmp.vcf.gz"
    params:
        header = "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/header_{chr}.txt"
    shell:
        """
	      head -n6 {input} > {params.header}
	      cat {params.header} <(cat {input} | awk -v OFS="\t" 'NR>6 {{$3=$1"_"$2"_A_T";$4="A"; $5="T"; print ;}}') | bgzip > {output.gz}
			  rm {params.header}
			  rm {input}
			  """

rule concat_vcfs:
    input:
        expand("output/Simulate_Genotypes/{{rep}}/{{sample_size}}/{{theta}}/genos-{chr}.tmp.vcf.gz", chr=CHR)
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.ids.vcf.gz"
    params:
        tmp = "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/temp.vcf.gz"
    shell:
        """
        bcftools concat {input} -o {params.tmp} -O z
        bcftools annotate --rename-chrs code/Simulate_Genotypes/convert_chr.txt {params.tmp} -o {output} -O z
        rm {params.tmp}
        rm {input}
        """

rule convert_vcf_to_plink:
    input:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.ids.vcf.gz"
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.psam",
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pgen",
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pvar"
    params:
        outprefix = "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos"
    shell:
        """
        plink2 \
        --double-id \
        --make-pgen \
        --out {params.outprefix} \
        --vcf {input}
        rm {input}
        """

rule concat_pop_files:
    input:
        all = expand("output/Simulate_Genotypes/{{rep}}/{{sample_size}}/{{theta}}/genos-{chr}.pop", chr=CHR),
        one = "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos-1.pop"
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pop"
    shell:
        """
        mv {input.one} {output}
        rm -f {input.all}
        """

####################
### Select sites ###
####################

rule get_AF:
    input:
        psam="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.psam",
        pgen="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pgen",
        pvar="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pvar"
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.afreq"
    params:
        plink_prefix = "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos"
    shell:
        """
        plink2 \
	      --pfile {params.plink_prefix} \
	      --freq \
	      --out {params.plink_prefix}
        """

rule get_snp_list:
    input:
        freq="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.afreq"
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/snp_list.txt"
    params:
        num_snp = lambda wildcards: get_params(wildcards.nsnp)
    shell:
        """
        Rscript code/Simulate_Genotypes/get_snp_list.R {input.freq} {output} {params.num_snp}
        """


###########################
### Compute Fst and PCA ###
###########################

rule calculate_fst:
    input:
        psam="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.psam",
        pgen="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pgen",
        pvar="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pvar",
        pop="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pop",
        snps="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/snp_list.txt"
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/genos.fst.summary"
    params:
        plink_prefix="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos",
        out_prefix="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/genos"
    shell:
        """
	plink2 \
	      --pfile {params.plink_prefix} \
	      --extract {input.snps} \
	      --pheno {input.pop} \
	      --fst POP \
	      --out {params.out_prefix}
		    """

rule PCA:
    input:
        psam="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.psam",
        pgen="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pgen",
        pvar="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pvar",
        snps="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/snp_list.txt"
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/genos.eigenvec",
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/genos.eigenval"
    params:
        plink_prefix = "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos",
        out_prefix = "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/genos"
    shell:
        """
        plink2 \
	      --pfile {params.plink_prefix} \
	      --extract {input.snps} \
	      --out {params.out_prefix} \
		    --pca 10 approx
		    """

###################
### Compute b^2 ###
###################

rule compute_b2:
    input:
        evec="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/genos.eigenvec",
        pop="output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/genos.pop"
    output:
        "output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/b2.txt"
    shell:
        """
        Rscript code/Simulate_Genotypes/compute_B2.R {input.evec} {input.pop} {output}
        """

#######################
### Compile Results ###
#######################


rule compile_results:
    input:
        b2=expand("output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/b2.txt", rep = REP, sample_size=SAMPLE_SIZE, theta=THETA, nsnp=NSNP),
        fst=expand("output/Simulate_Genotypes/{rep}/{sample_size}/{theta}/{nsnp}/genos.fst.summary", rep = REP, sample_size=SAMPLE_SIZE, theta=THETA, nsnp=NSNP)
    output:
        "results/all.txt"
    shell:
        """
        Rscript code/summarize_results.R {output} {input.b2}
        """





