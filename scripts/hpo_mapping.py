#! /usr/bin/env python

import argparse
from pronto import Ontology

ancestor_cache = {}

def get_ancestors(term, obo):
    if term not in ancestor_cache:
        ancestor_cache[term] = {p.id for p in obo[term].superclasses(with_self=False)}
    return ancestor_cache[term]

def mica(p_term, d_term, obo, ics):
    a1 = {p_term} | get_ancestors(p_term, obo)
    a2 = {d_term} | get_ancestors(d_term, obo)
    common_ancestors = a1 & a2
    if not common_ancestors:
        return None
    return max(common_ancestors, key=lambda t: ics.get(t, 0.0))

def term_similarity(t1, t2, obo, ics):
    m = mica(t1, t2, obo, ics)
    if m is None:
        return 0
    return ics.get(m, 0)

def set_similarity(query_terms, ref_terms, obo, ics):
    scores = []
    for q in query_terms:
        term_scores = []
        for r in ref_terms:
            score = term_similarity(q, r, obo, ics)
            term_scores.append(score)
        best = max(term_scores)
        scores.append(best)
    return sum(scores) / len(scores)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Phenotype-based gene prioritization using HPO similarity"
    )

    parser.add_argument(
        "--obo",
        required=True,
        help="Path to hp.obo file"
    )

    parser.add_argument(
        "--ic",
        required=True,
        help="Information content file (HPO_ID <tab> IC)"
    )

    parser.add_argument(
        "--genes",
        required=True,
        dest="g2p",
        help="genes_to_phenotype.txt file"
    )

    parser.add_argument(
        "--patient",
        required=True,
        help="Patient HPO list (one HPO per line)"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output file path"
    )

    return parser.parse_args()

#laod ICs
def load_ic(ic_file):
    ics = {}
    with open(ic_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            hpo_id, ic_value = line.split('\t')[:2]
            ics[hpo_id] = float(ic_value)
    return ics

#load gene to phenotype
def load_gene_hpos(g2p_file):
    gene_hpos = {}
    with open(g2p_file) as f:
        next(f)
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            gene, hpo = fields[1], fields[2]
            gene_hpos.setdefault(gene, set()).add(hpo)
    return gene_hpos

#load patient hpos
def load_patient_hpos(patient_file):
    p_hpos = set()
    with open(patient_file) as f:
        next(f)
        for line in f:
            hpo = line.strip().split("\t")[0]
            if hpo.startswith("HP:"):
                p_hpos.add(hpo)
    return p_hpos


def main():
    args = parse_args()

    # Load ontology
    obo = Ontology(args.obo)

    # Load inputs
    ics = load_ic(args.ic)
    gene_hpos = load_gene_hpos(args.g2p)
    p_hpos = load_patient_hpos(args.patient)

    gene_scores = {}
    for gene, gene_hpo in gene_hpos.items():
        sim = set_similarity(p_hpos, gene_hpo, obo, ics)
        gene_scores[gene] = sim
    
    # Sort genes by similarity score in descending order
    with open(args.output, 'w') as out_file:
        sorted_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
        print("gene", 'sim_score', sep="\t", file=out_file)
        for gene, sim in sorted_genes:
            print(gene, round(sim, 4), sep="\t", file=out_file)

if __name__ == "__main__":
    main()