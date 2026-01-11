
rule all:
    input:
        'patient/patient_gene_scores.txt'
    
rule clin_to_phen:
    input:
        txt = "patient/hcm.txt"
    output:
        hpo_terms = "patient/patient_hpo_terms.txt"
    conda: "envs/hpo_gene_mapping.yaml"
    shell:
        """
            clinphen {input.txt} |awk 'NR > 1 {{print $1}}'> {output.hpo_terms}
        """
rule build_ic:
    input:
        obo = 'data/hp.obo',
        pheno = 'data/phenotype.hpoa'
    output:
        ics = 'data/hpo_ics.txt'
    conda: "envs/hpo_gene_mapping.yaml"
    shell:
        '''
            ./scripts/make_ics.py --obo {input.obo} \
            --phenotype {input.pheno} \
            --output {output.ics}
        '''

rule map_phenotype:
    input:
        obo = 'data/hp.obo',
        ics = 'data/hpo_ics.txt',
        gene_hpo = 'data/genes_to_phenotype.txt',
        patient_hpo = 'patient/patient_hpo_terms.txt'
    output:
        mapping = 'patient/patient_gene_scores.txt'
    conda: "envs/hpo_gene_mapping.yaml"
    shell:
        '''
            ./scripts/hpo_mapping.py \
            --obo {input.obo} \
            --ic {input.ics} \
            --genes {input.gene_hpo} \
            --patient {input.patient_hpo} \
            --output {output.mapping}
        '''
