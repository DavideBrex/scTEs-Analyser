

import pandas as pd

#read fastq files
table_fastq=pd.read_csv("Fastq_files.tsv", delimiter = "\t", index_col=False)


#read config file
configfile: "config.yaml"
#set image
singularity: "/shares/CIBIO-Storage/GROUPS/sharedLC/Davide/containers/tes-analyser-cont.sif"


rule all:
    input:
       "results/alignments/ESCs.bam"


rule StarSOLO_10X:
    input: 
        fq1 = table_fastq["Fq1"],
        fq2 = table_fastq["Fq2"]
    output: 
        "results/alignments/ESCs.bam"
    threads: 20
    params:
        type_sc=config["scRNAseqType"],
        genome_index=config["ref"]["genome_index"],
        out_dir="results/alignments/",
        whitelist=config["whitelist"]
    shell:
        """
        STAR --genomeDir {params.genome_index} \
        --readFilesIn {input.fq2} {input.fq1} \
        --runThreadN {threads} \
        --soloType {params.type_sc} \
        --soloCBwhitelist {params.whitelist} \
        --outSAMattributes NH HI AS nM CR CY UR UY \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 100 \
        --winAnchorMultimapNmax 100 \
        --outMultimapperOrder Random \
        --runRNGseed 777 \
        --outSAMmultNmax 1 \
        --limitBAMsortRAM 56736503447 \
        --outFileNamePrefix {params.out_dir}

        mv {params.out_dir}Aligned.sortedByCoord.out.bam {output}
        """


#build index for scTEs
rule build_index_scTEs:
    params:
        organism=config["organism"]
    output:
        gene_annot="resources/annot_file_genes.gtf.gz",
        rmsk="resources/rmsk.txt.gz"
    shell:
        """
        /home/davide.bressan-1/tools/scTE/bin/scTE_build -g {params.organism}

        mv *.gtf.gz {output.gene_annot}
        mv rmsk.txt.gz {output.rmsk}
        """




rule run_scTEs:
    input: rule.StarSOLO_10X.output
    output: 
    params:
    shell:
        """
        scTE_build -g {organism}
        """ 





# rule StarSOLO_other:
#     input: 
#         fq1 = table_fastq["Fq1"],
#         fq2 = table_fastq["Fq2"]
#     output: 
#         "results/alignments/ESCs.bam"
#     threads: 20
#     params:
#         type_sc=config["scRNAseqType"],
#         genome_index=config["ref"]["genome_index"],
#         out_dir="results/alignments/"
#     shell:
#         """
#         STAR --genomeDir {params.genome_index} \
#         --readFilesIn {input.fq1} {input.fq2} \
#         --runThreadN {threads} \
#         --soloType {params.type_sc} \
#         --outSAMattributes NH HI AS nM CR CY UR UY \
#         --readFilesCommand zcat \
#         --outSAMtype BAM SortedByCoordinate \
#         --outFilterMultimapNmax 100 \
#         --winAnchorMultimapNmax 100 \
#         --out-MultimapperOrder Random \
#         --runRNGseed 777 \
#         --outSAMmultNmax 1 \
#         --outFileNamePrefix {params.out_dir}

#         mv {params.out_dir}Aligned.sortedByCoord.out.bam {output}
#         """