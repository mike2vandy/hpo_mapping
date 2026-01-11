# hpo_mapping

This code mostly replicates the core of Exomiser to rank genes in order from highest to lowest chance of being associated with a rare Mendelian genetic disease. 

HPO information is in `data/`

Tool versions and packages used for snakemake are listed in `envs/hpo_gene_mapping.yaml`

Clinical paragraphs as the starting point are in `patients/`, along with extracted HPO terms, and ranked gene lists.

Scripts for generating IC scores and ranking genes are in scripts. 
