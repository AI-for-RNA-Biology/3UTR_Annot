rule quant_utr:
    input:
        L2=config["results_dir"]+'annotation/L2.gtf',
        L2_300=config["results_dir"]+'annotation/L2_300.gtf',
        bam=config["bam_dir"]+'{SRR_ID}.Aligned.sortedByCoord.out.bam'
    output:
        UTR=config["results_dir"]+'counts_3UTR/{SRR_ID}.txt'
    params:
        output_dir_utr=config["results_dir"]+'counts_3UTR/'
    conda:
         "../../envs/subread.yml"
    shell:
        """
            mkdir -p {params.output_dir_utr}
            featureCounts -a {input.L2_300} -o {output.UTR} -g transcript_id -f -O -M -s 2 -p --verbose -T 4 {input.bam}
        """

rule quant_genes:
    input:
        gtf=config["ref_gtf"],
        bam=config["bam_dir"]+'{SRR_ID}.Aligned.sortedByCoord.out.bam'
    output:
        genes=config["results_dir"]+'counts_genes/{SRR_ID}.txt'
    params:
        output_dir_genes=config["results_dir"]+"counts_genes/"
    conda:
         "../../envs/subread.yml"
    shell:
        """
            mkdir -p {params.output_dir_genes}
            featureCounts -a {input.gtf} -o {output.genes} -g gene_id -O -M -s 2 -p --verbose -T 4 {input.bam}

        """
