rule generate_fasta:
    input:
        gross_fragments_bed=config["results_dir"]+"annotation/GrossFragments.bed",
        fasta=config["ref_fasta"],
    output:
        gross_fragments_fasta=config["results_dir"]+"annotation/GrossFragments.fasta"
    conda:
        "../../envs/bedtools.yml"
    shell:
        """
            #FASTA sequence of gross fragments
            bedtools getfasta -fi {input.fasta} -bed {input.gross_fragments_bed} -s -name -fo {output.gross_fragments_fasta}
        """

rule identify_pas:
    input:
        gross_fragments=config["results_dir"]+"annotation/GrossFragments.gtf",
        gross_fragments_fasta=config["results_dir"]+"annotation/GrossFragments.fasta",
        fasta=config["ref_fasta"],
        pas_motif=config["PAS_motif"]
    output:
        myPAS=config["results_dir"]+"annotation/myPAS.gtf"
    params:
        pas_fragments_dir=config["results_dir"]+"utrid/GrossSegments/"
    conda:
         "../../envs/annot.yml"
    shell:
        """            
            mkdir -p {params.pas_fragments_dir}
            #Find PAS elements in 3'UTR of the fragments
            echo -e  "1\n 2\n 3\n 4\n 5\n 6\n 7\n 8\n 9\n 10\n 11\n 12"| while read IX
            do
                Rscript ./scripts/findPASelement.R {params.pas_fragments_dir} ${{IX}} {input.gross_fragments_fasta} {input.gross_fragments} {input.pas_motif}
            done

            Rscript ./scripts/mergeGTF.R {params.pas_fragments_dir} {output.myPAS}
        """
