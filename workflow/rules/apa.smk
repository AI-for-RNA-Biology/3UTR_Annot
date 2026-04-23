rule extract_cov_gross_segment:
    input:
        gross_fragments=config["results_dir"]+"annotation/GrossFragments.gtf",
        bam=config["bam_dir"]+'{SRR_ID}.Aligned.sortedByCoord.out.bam'
    params:
        WITHCOR='TRUE',
        PAIRS='TRUE',
        FRAC=0.05,
        PROTOCOL='reverse'
    output:
        gross_seg_cov=directory(config["results_dir"]+"utrCov/{SRR_ID}/")
    conda:
         "../../envs/annot.yml"
    shell:
        """
            cat {input.gross_fragments} | awk '$1 ~ /^[0-9]+$/ || $1 ~ /^chr[0-9]+$/ || $1 == "X" || $1 == "Y" || $1 == "chrX" || $1 == "chrY" {{print $1}}' | sort | uniq | while read CHROM
            do
                Rscript ./scripts/extractCovGrossSegments.R {input.bam} {output} ${{CHROM}} {wildcards.SRR_ID} {input.gross_fragments} {params.WITHCOR} {params.PAIRS} {params.FRAC} {params.PROTOCOL}
            done
        """

rule identify_APA_from_cov:
    input:
        gross_frag_cov=config["results_dir"]+"utrCov/{SRR_ID}/",
        gross_fragments=config["results_dir"]+"annotation/GrossFragments.gtf",
        gross_fragments_rds=config["results_dir"]+"annotation/GrossFragments.rds",
        annot=config["ref_annot"]
    params:
        METHOD='SEG',
        limD=50,
        winD=151,
        noSeg=10,
        isLog='TRUE',
        minCov=8,
        n_cores=12
    output:
        out_dir=directory(config["results_dir"]+'utrid/APA/{SRR_ID}/')
    conda:
         "../../envs/apa.yml"
    shell:
        """
            cat {input.gross_fragments} | awk '$1 ~ /^[0-9]+$/ || $1 ~ /^chr[0-9]+$/ || $1 == "X" || $1 == "Y" || $1 == "chrX" || $1 == "chrY" {{print $1}}' | sort | uniq | while read CHROM
            do
                Rscript ./scripts/identify_pA_sites.R {output.out_dir} ${{CHROM}} {input.gross_fragments_rds} {params.METHOD} {input.annot} {params.limD} {params.winD} {params.noSeg} {params.isLog} {params.minCov} {input.gross_frag_cov} &

                # At most as number of CPU cores
                [ $( jobs | wc -l ) -ge {params.n_cores} ] && wait

            done
            
        """

rule filter_UTR:
    input:
        apa_dir=expand(config["results_dir"]+'utrid/APA/{SRR_ID}',SRR_ID=SRR_ID),
        myPAS=config["results_dir"]+"annotation/myPAS.gtf",
        annot=config["ref_annot"],
        rpm=config["rpm"]
    params:
        out_dir=config["results_dir"]+"annotation/",
        polyA_DB=config['polyA_DB'],
        PolyA_Site=config['polyA_Site'],
        PolyA_Models=config['polyA_Models']
    output:
        L2=config["results_dir"]+'annotation/L2.gtf',
        L2_300=config["results_dir"]+'annotation/L2_300.gtf'
    conda:
         "../../envs/annot.yml"
    shell:
        """
            Rscript ./scripts/FilterUTR.R {params.out_dir} {input.myPAS} {params.polyA_DB} {params.PolyA_Site} {params.PolyA_Models} {input.annot} {input.rpm} {input.apa_dir}
        """

rule filter_UTR_gene:
    input:
        apa_dir=expand(config["results_dir"]+'utrid/APA/{SRR_ID}',SRR_ID=SRR_ID),
        myPAS=config["results_dir"]+"annotation/myPAS.gtf",
        annot=config["ref_annot"],
        rpm=config["rpm"]
    params:
        out_dir=config["results_dir"]+"annotation/",
        polyA_DB=config['polyA_DB'],
        PolyA_Site=config['polyA_Site'],
        PolyA_Models=config['polyA_Models'],
        utr_refseq=config["utr_refseq"],
        e75p_RefSeq=config["e75p_RefSeq"],
        assembly_report=config["assembly_report"]
    output:
        Focus = config["results_dir"]+'annotation/hg38_Denovo_Ensembl_RefSeq_Filt_Focus_300.gtf'
    conda:
         "../../envs/annot.yml"
    shell:
        """
            Rscript ./scripts/FilterUTRGeneLevel.R {params.out_dir} {input.myPAS} {params.polyA_DB} {params.PolyA_Site} {params.PolyA_Models} {input.annot} {params.utr_refseq} {params.e75p_RefSeq} {params.assembly_report} {input.rpm} {input.apa_dir}
        """

