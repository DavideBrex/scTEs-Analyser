

import pandas as pd

#read fastq files
table_fastq=pd.read_csv("Fastq_files.tsv", delimiter = "\t", index_col=False)


#read config file
configfile: "config.yaml"


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
        out_dir="results/alignments/"
    shell:
        """
        STAR --genomeDir {params.genome_index} \
        --readFilesIn {input.fq2} {input.fq1} \
        --runThreadN {threads} \
        --soloType {params.type_sc} \
        --outSAMattributes NH HI AS nM CR CY UR UY \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 100 \
        --winAnchorMultimapNmax 100 \
        --out-MultimapperOrder Random \
        --runRNGseed 777 \
        --outSAMmultNmax 1 \
        --outFileNamePrefix {params.out_dir}

        mv {params.out_dir}/Aligned.sortedByCoord.out.bam {output}
        """