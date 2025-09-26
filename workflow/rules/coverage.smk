#Create coverage file per sample for each strand (+ and -)
#Keep positions with at least 7 reads
rule stranded_coverage:
    input:
        bam=config["bam_dir"]+'{SRR_ID}.Aligned.sortedByCoord.out.bam'
    output:
        pos_cov=config["results_dir"]+'Cov/{SRR_ID}.pos.coverage.gz',
        neg_cov=config["results_dir"]+'Cov/{SRR_ID}.neg.coverage.gz'
    params:
        out_dir=config["results_dir"]+'Cov/',
        N_reads=7,
        Protocol='F2R1'
    conda:
         "../../envs/alignment.yml"
    shell:
        """
            #1 Seperate Mates (0x0040 = Mate 1 of pairs, 0x0080 = Mate2 of pairs)
            samtools view -b -f 0x0040 -h -o {params.out_dir}{wildcards.SRR_ID}.m1.bam {input.bam}  
            samtools view -b -f 0x0080 -h -o {params.out_dir}{wildcards.SRR_ID}.m2.bam {input.bam} 

            #Seperate strands
            #Mate 1, Watson (+) strand
            bedtools genomecov -ibam {params.out_dir}{wildcards.SRR_ID}.m1.bam -dz -split -strand + > {params.out_dir}{wildcards.SRR_ID}.m1.pos.coverage

            #Mate 1, Crick (-) strand
            bedtools genomecov -ibam {params.out_dir}{wildcards.SRR_ID}.m1.bam -dz -split -strand - > {params.out_dir}{wildcards.SRR_ID}.m1.neg.coverage

            #Mate 2, Watson (+) strand
            bedtools genomecov -ibam {params.out_dir}{wildcards.SRR_ID}.m2.bam -dz -split -strand + > {params.out_dir}{wildcards.SRR_ID}.m2.pos.coverage

            #Mate 2, Crick (-) strand
            bedtools genomecov -ibam {params.out_dir}{wildcards.SRR_ID}.m2.bam -dz -split -strand - > {params.out_dir}{wildcards.SRR_ID}.m2.neg.coverage

            #Cleanup
            rm {params.out_dir}{wildcards.SRR_ID}.m1.bam
            rm {params.out_dir}{wildcards.SRR_ID}.m2.bam

            #Merge Mate1 and Mate2, based on protocol
            PROTOCOL={params.Protocol}
            if [ $PROTOCOL == 'F2R1' ];then
                #Remove positions not covered by at least N_reads
                cat {params.out_dir}{wildcards.SRR_ID}.m1.neg.coverage {params.out_dir}{wildcards.SRR_ID}.m2.pos.coverage \
                | awk '{{if($3 > {params.N_reads}) print($0)}}' \
                > {params.out_dir}{wildcards.SRR_ID}.pos.coverage
            
                cat {params.out_dir}{wildcards.SRR_ID}.m1.pos.coverage {params.out_dir}{wildcards.SRR_ID}.m2.neg.coverage \
                | awk '{{if($3 > {params.N_reads}) print($0)}}' \
                > {params.out_dir}{wildcards.SRR_ID}.neg.coverage
                
            elif [ $PROTOCOL == 'F1R2' ];then 
                #Remove positions not covered by at least N_reads            
                cat {params.out_dir}{wildcards.SRR_ID}.m1.pos.coverage {params.out_dir}{wildcards.SRR_ID}.m2.neg.coverage \
                | awk '{{if($3 > {params.N_reads}) print($0)}}' \
                > {params.out_dir}{wildcards.SRR_ID}.pos.coverage

                cat {params.out_dir}{wildcards.SRR_ID}.m1.neg.coverage {params.out_dir}{wildcards.SRR_ID}.m2.pos.coverage \
                | awk '{{if($3 > {params.N_reads}) print($0)}}' \
                > {params.out_dir}{wildcards.SRR_ID}.neg.coverage
            fi            

            pigz {params.out_dir}{wildcards.SRR_ID}.pos.coverage
            pigz {params.out_dir}{wildcards.SRR_ID}.neg.coverage

            #Cleanup
            rm {params.out_dir}{wildcards.SRR_ID}.m2.neg.coverage
            rm {params.out_dir}{wildcards.SRR_ID}.m2.pos.coverage
            rm {params.out_dir}{wildcards.SRR_ID}.m1.neg.coverage
            rm {params.out_dir}{wildcards.SRR_ID}.m1.pos.coverage

        """
