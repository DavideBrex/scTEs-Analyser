

import pandas as pd

#read fastq files
table_fastq=pd.read_csv("Fastq_files.tsv", delimiter = "\t", index_col=False)


#read config file
configfile: "config.yaml"
#set image
singularity: "/shares/CIBIO-Storage/GROUPS/sharedLC/Davide/containers/tes-analyser-cont.sif"


organism_g=config["organism"]

rule all:
    input:
       "results/ESCs.hdf5"
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
        --limitBAMsortRAM 200000000000 \
        --outFileNamePrefix {params.out_dir}

        mv {params.out_dir}Aligned.sortedByCoord.out.bam {output}
        """


#build index for scTEs
rule build_index_scTEs:
    params:
        organism=config["organism"]
    output:
        "resources/index.exclusive.idx"
    shell:
        """
        /home/davide.bressan-1/tools/scTE/bin/scTE_build -g {params.organism}

        mv {params.organism}.exclusive.idx resources/index.exclusive.idx 
        """



rule run_scTEs:
    input: 
        bam_file=rules.StarSOLO_10X.output,
        index=rules.build_index_scTEs.output
    output: 
        "results/ESCs.hdf5"
    params:
        out="results/ESCs.hdf5",
        organism=config["organism"]
    threads: 5
    shell:
        """
        /home/davide.bressan-1/tools/scTE/bin/scTE \
        -i {input.bam_file} \
        -o {params.out} \
        -x {input.index}  \
        -p {threads} \
        --hdf5 True -CB CR -UMI UR    
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