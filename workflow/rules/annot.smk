rule merge_and_filter_utr:
    input:
        bam=config["bam_dir"]+'{SRR_ID}.Aligned.sortedByCoord.out.bam',
        pos_cov=config["results_dir"]+'Cov/{SRR_ID}.pos.coverage.gz',
        neg_cov=config["results_dir"]+'Cov/{SRR_ID}.neg.coverage.gz',
        annot=config["ref_annot"],
        rpm=config["rpm"]
    output:
        filter9=config["results_dir"]+'GrossSegments/GFF/{SRR_ID}.9.gff',
        filter2=config["results_dir"]+'GrossSegments/GFF/{SRR_ID}.2.gff',
        filter6=config["results_dir"]+'GrossSegments/GFF/{SRR_ID}.6.gff'
    params:
        pos_cov_dir=config["results_dir"]+'Cov/{SRR_ID}/pos/',
        neg_cov_dir=config["results_dir"]+'Cov/{SRR_ID}/neg/',
        out_gff_dir=config["results_dir"]+'GrossSegments/GFF/'
    group:
        "annot_1"
    conda:
         "../../envs/annot.yml"
    shell:
        """ 
            WD=$(pwd)

            #Make main directory
            mkdir -p {config[results_dir]}"Cov/"{wildcards.SRR_ID}

            #Positive strand coverage by chr (autosomes and X/Y)
            mkdir -p {params.pos_cov_dir}
            cd {params.pos_cov_dir}
            zcat {input.pos_cov} > {params.pos_cov_dir}/input.bed 
            awk '$1 ~ /^[0-9]+$/ || $1 ~ /^chr[0-9]+$/ || $1 == "X" || $1 == "Y" || $1 == "chrX" || $1 == "chrY" {{close(f);f=$1; print $2 > f".bed"}}' {params.pos_cov_dir}/input.bed 
            rm {params.pos_cov_dir}/input.bed 

            #Negative strand coverage by chr (autosomes and X/Y)
            mkdir -p {params.neg_cov_dir}
            cd {params.neg_cov_dir}
            zcat {input.neg_cov} > {params.neg_cov_dir}/input.bed 
            awk '$1 ~ /^[0-9]+$/ || $1 ~ /^chr[0-9]+$/ || $1 == "X" || $1 == "Y" || $1 == "chrX" || $1 == "chrY" {{close(f);f=$1; print $2 > f".bed"}}' {params.neg_cov_dir}/input.bed 
            rm {params.neg_cov_dir}/input.bed 

            cd $WD
            #Identify segments longer than 100 nt
            ls {params.pos_cov_dir}/*.bed | while read FILE
            do
                out_file=${{FILE%.bed}}.txt
                python2 ./scripts/findRegions.py -i ${{FILE}} -s 100 -t 80 -c 0 -o ${{out_file}}
            done

            ls {params.neg_cov_dir}/*.bed | while read FILE
            do
                out_file=${{FILE%.bed}}.txt
                python2 ./scripts/findRegions.py -i ${{FILE}} -s 100 -t 80 -c 0 -o ${{out_file}}
            done

            #Merge expressed regions gapped by low-mappable using script mergeNeighbors.R and filter 
            #Positive strand
            #Loop through chromosomes
            zcat {input.pos_cov} | awk '$1 ~ /^[0-9]+$/ || $1 ~ /^chr[0-9]+$/ || $1 == "X" || $1 == "Y" || $1 == "chrX" || $1 == "chrY" {{print $1}}' | sort | uniq | while read CHROM
            do
                Rscript ./scripts/merge_and_filter.R {params.pos_cov_dir} ${{CHROM}} {input.bam} + TRUE {input.annot} {input.rpm} TRUE
            done

            #Negative strand
            #Loop through chromosomes
            zcat {input.neg_cov} | awk '$1 ~ /^[0-9]+$/ || $1 ~ /^chr[0-9]+$/ || $1 == "X" || $1 == "Y" || $1 == "chrX" || $1 == "chrY" {{print $1}}' | sort | uniq | while read CHROM
            do
                Rscript ./scripts/merge_and_filter.R {params.neg_cov_dir} ${{CHROM}} {input.bam} - TRUE {input.annot} {input.rpm} TRUE
            done
            
            #Merge all GFF (across chromosomes)
            mkdir -p {config[results_dir]}"GrossSegments/"
            mkdir -p {params.out_gff_dir}

            cat {params.pos_cov_dir}/merge_and*/groiv1.gff | grep -v '#' > {params.out_gff_dir}{wildcards.SRR_ID}".2.gff"
            cat {params.neg_cov_dir}/merge_and*/groiv1.gff | grep -v '#' >> {params.out_gff_dir}{wildcards.SRR_ID}".2.gff"
            
            cat {params.pos_cov_dir}/merge_and*/gF2.gff | grep -v '#' > {params.out_gff_dir}{wildcards.SRR_ID}".6.gff"
            cat {params.neg_cov_dir}/merge_and*/gF2.gff | grep -v '#' >> {params.out_gff_dir}{wildcards.SRR_ID}".6.gff"

            cat {params.pos_cov_dir}/merge_and*/gF7.gff | grep -v '#' > {params.out_gff_dir}{wildcards.SRR_ID}".9.gff"
            cat {params.neg_cov_dir}/merge_and*/gF7.gff | grep -v '#' >> {params.out_gff_dir}{wildcards.SRR_ID}".9.gff"
        """

rule extend_utr:
    input:
        filter9=config["results_dir"]+'GrossSegments/GFF/{SRR_ID}.9.gff',
        filter2=config["results_dir"]+'GrossSegments/GFF/{SRR_ID}.2.gff',
        annot=config["ref_annot"],
    output:
        extended_gff=config["results_dir"]+'GrossSegments/Extended/{SRR_ID}.extended.gff'
    params:
        out_dir=config["results_dir"]+'GrossSegments/Extended/',
        LIM=300,
        MAXDIST=16243
    group:
        "annot_1"
    conda:
         "../../envs/annot.yml"
    shell:
        """
            mkdir -p {params.out_dir}
            Rscript ./scripts/id_far_extensions_general.R {input.filter2} {params.out_dir} {wildcards.SRR_ID} {params.LIM} {input.filter9} {params.MAXDIST} {input.annot}

            cp {params.out_dir}"new_extensions_"{wildcards.SRR_ID}"/gF9_"{wildcards.SRR_ID}.gff {params.out_dir}{wildcards.SRR_ID}.extended.gff
        """

rule merge_annotations:
    input:
        filter6=expand(config["results_dir"]+'GrossSegments/GFF/{SRR_ID}.6.gff',SRR_ID=SRR_ID),
        extended_gff=expand(config["results_dir"]+'GrossSegments/Extended/{SRR_ID}.extended.gff',SRR_ID=SRR_ID),
        annot=config["ref_annot"]
    output:
        gross_fragments=config["results_dir"]+"annotation/GrossFragments.gtf",
        gross_fragments_bed=config["results_dir"]+"annotation/GrossFragments.bed",
        gross_fragments_rds=config["results_dir"]+"annotation/GrossFragments.rds"
    params:
        out_dir=config["results_dir"]+"annotation/"
    conda:
         "../../envs/annot.yml"
    shell:
        """
            mkdir -p {params.out_dir}

            Rscript ./scripts/mergeAnnotations.R {output.gross_fragments} {input.annot} {input.filter6} {input.extended_gff}
        """
        


    

