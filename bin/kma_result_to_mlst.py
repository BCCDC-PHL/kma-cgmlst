#!/usr/bin/env python

from __future__ import print_function

import argparse
import csv
import json
import os
import sys

from pprint import pprint

def parse_res_file(res_file_path, locus_allele_delimiter):
    
    res_fieldnames = [
        'template',
        'score',
        'expected',
        'template_length',
        'template_identity',
        'template_coverage',
        'query_identity',
        'query_coverage',
        'depth',
        'q_value',
        'p_value',
    ]
    
    with open(res_file_path, 'r') as f:
        loci = {}
        reader = csv.DictReader(f, fieldnames=res_fieldnames, dialect="excel-tab")
        next(reader) #skip header
        for row in reader:
            locus, allele = map(str.strip, row['template'].rsplit(locus_allele_delimiter, 1))
            if locus in loci:
                loci[locus][allele] = {
                    'locus_id': locus,
                    'allele_id': allele,
                    'score': int(row['score'].strip()),
                    'expected': int(row['expected'].strip()),
                    'template_length': int(row['template_length'].strip()),
                    'template_identity': float(row['template_identity'].strip()),
                    'template_coverage': float(row['template_coverage'].strip()),
                    'query_identity': float(row['query_identity'].strip()),
                    'query_coverage': float(row['query_coverage'].strip()),
                    'depth': float(row['depth'].strip()),
                    'q_value': float(row['q_value'].strip()),
                    'p_value': float(row['p_value'].strip()),
                }
            else:
                loci[locus] = {}
                loci[locus][allele] = {
                    'locus_id': locus,
                    'allele_id': allele,
                    'score': int(row['score'].strip()),
                    'expected': int(row['expected'].strip()),
                    'template_length': int(row['template_length'].strip()),
                    'template_identity': float(row['template_identity'].strip()),
                    'template_coverage': float(row['template_coverage'].strip()),
                    'query_identity': float(row['query_identity'].strip()),
                    'query_coverage': float(row['query_coverage'].strip()),
                    'depth': float(row['depth'].strip()),
                    'q_value': float(row['q_value'].strip()),
                    'p_value': float(row['p_value'].strip()),
                }

        return loci

def main(args):

    loci = parse_res_file(args.res, args.locus_allele_delimiter)
    print("\t".join([
        "locus_id",
        "allele_id",
        "template_identity",
        "template_coverage",
        "depth",
        ]))

    for locus, alleles in loci.items():
        best_allele = sorted(alleles.values(),
                             key=lambda x: x['score'], reverse=True)[0]['allele_id']

        print("\t".join([
            alleles[best_allele]['locus_id'],
            alleles[best_allele]['allele_id'],
            str(alleles[best_allele]['template_identity']),
            str(alleles[best_allele]['template_coverage']),
            str(alleles[best_allele]['depth']),
        ]))
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--res", dest="res", help="KMA result overview file")
    parser.add_argument("--locus_allele_delimiter", help="Delimiter separating locus id from allele id", default='_')
    args = parser.parse_args()
    main(args)
