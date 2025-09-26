#Create STAR Index
rule star_index:
    input:
        fasta=config["ref_fasta"],
        gtf=config["ref_gtf"]
    output:
        index_dir=directory(config["genome_dir"]+"Index/")
    conda:
         "../../envs/alignment.yml"
    threads:
        32
    shell:
        """
            STAR --runThreadN {threads} --runMode genomeGenerate \
            --genomeDir {output.index_dir} --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf}
        """
 
#Run Alignment
rule star_aln:
    input:
        read1=config["fastq_dir"]+"{SRR_ID}_1.fastq.gz",
        read2=config["fastq_dir"]+"{SRR_ID}_2.fastq.gz",
        index_dir=config["genome_dir"]+"Index/",
        gtf=config["ref_gtf"]
    log:
        config["results_dir"]+'BAM/{SRR_ID}.log'
    output:
        bam=config["results_dir"]+'BAM/{SRR_ID}.Aligned.sortedByCoord.out.bam'
    params:
        prefix=config["results_dir"]+"BAM/{SRR_ID}."
    threads:
        5
    conda:
        "../../envs/alignment.yml"
    shell:
        """
		    STAR --runThreadN {threads} --genomeDir {input.index_dir} \
		    --readFilesIn {input.read1} {input.read2} \
		    --readFilesCommand zcat --sjdbGTFfile {input.gtf} \
		    --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate \
		    --outReadsUnmapped Within --outSAMattributes Standard > {log}

            samtools index {output}
        """    
