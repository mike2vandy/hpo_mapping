#! /usr/bin/env python

import argparse, math
from collections import defaultdict
from pronto import Ontology

def build_child_dict(obo):
    all_children_dict = {}  
    for term in obo.terms():
        descendants = term.subclasses() 
        descendant_ids = {d.id for d in descendants}        
        descendant_ids.add(term.id)
        all_children_dict[term.id] = descendant_ids

    return all_children_dict

def load_hpo_to_omims(phenotype_file):
    omims = set()
    hpos = defaultdict(set)

    with open(phenotype_file) as f:
        for line in f:
            if line.startswith('#') or line.startswith('database'):
                continue
            fields = line.strip().split('\t')
            omim_id, hpo_id = fields[0], fields[3]
            hpos[hpo_id].add(omim_id)
            omims.add(omim_id)
    return hpos, omims

def compute_ic(hpos, omims, children_dict):

    total_diseases = len(omims)
    for hpo_id, children in children_dict.items():
        disease_set = set()
        for child in children:
            disease_set |= hpos.get(child, set())
        raw = (len(disease_set) + 1) / (total_diseases + 1)
        ic = -math.log(raw)
        yield hpo_id, ic

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute HPO information content using disease annotations"
    )
    parser.add_argument(
        "--obo",
        required=True,
        help="Path to hp.obo file"
    )
    parser.add_argument(
        "--phenotype",
        required=True,
        help="phenotype.hpoa / HPO annotation file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="output file"
    )
    return parser.parse_args()

def main():
    args = parse_args()
    obo = Ontology(args.obo)
    children_dict = build_child_dict(obo)
    hpos, omims = load_hpo_to_omims(args.phenotype)
    with open(args.output, "w") as f:
        for hpo_id, ic in compute_ic(hpos, omims, children_dict):
            print(hpo_id, ic, sep="\t", file=f)

if __name__ == "__main__":
    main()