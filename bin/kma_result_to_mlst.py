#!/usr/bin/env python

import argparse
import csv
import collections
import json
import os
import sys

from pprint import pprint

def parse_alleles_file(alleles_file_path, locus_allele_delimiter):
    """
    Take the '.name' file from kma database and return all locus names
    """
    loci = set()
    with open(alleles_file_path, 'r') as f:
        for line in f:
            locus = line.strip().split(locus_allele_delimiter)[0]
            loci.add(locus)
    loci = sorted(list(loci))
    return loci

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
        reader = csv.DictReader(f, fieldnames=res_fieldnames)
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

    all_loci = parse_alleles_file(args.alleles, args.locus_allele_delimiter)

    loci = parse_res_file(args.res, args.locus_allele_delimiter)
    print(",".join([
        "locus_id",
        "allele_id",
        "template_identity",
        "template_coverage",
        "depth",
        ]))

    mlst_output = collections.OrderedDict()
    mlst_output['sample_id'] = args.sample_id

    for locus in all_loci:
        if locus in loci:
            alleles = loci[locus]
            best_allele = sorted(alleles.values(),
                                 key=lambda x: x['score'], reverse=True)[0]['allele_id']

            print(",".join([
                alleles[best_allele]['locus_id'],
                alleles[best_allele]['allele_id'],
                str(alleles[best_allele]['template_identity']),
                str(alleles[best_allele]['template_coverage']),
                str(alleles[best_allele]['depth']),
            ]))

            if alleles[best_allele]['template_identity'] >= args.min_identity and alleles[best_allele]['template_coverage'] >= args.min_coverage:
                mlst_output[alleles[best_allele]['locus_id']] = alleles[best_allele]['allele_id']
            else:
                mlst_output[alleles[best_allele]['locus_id']] = '-'
        else:
                mlst_output[locus] = '-'


    with open(args.output, 'w') as f:
        fieldnames = ['sample_id'] + all_loci
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            dialect='unix',
            quoting=csv.QUOTE_MINIMAL,
            extrasaction='ignore',
            lineterminator='\n'
        )
        writer.writeheader()
        writer.writerow(mlst_output)
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("res", help="KMA result overview file")
    parser.add_argument("--locus-allele-delimiter", help="Delimiter separating locus id from allele id", default='_')
    parser.add_argument("-s", "--sample-id", help="Sample ID", default='unknown')
    parser.add_argument("-a", "--alleles", help="List of all alleles")
    parser.add_argument("-i", "--min-identity", type=float, default=100.0, help="Minimum identity to consider an allele match")
    parser.add_argument("-c", "--min-coverage", type=float, default=100.0, help="Minimum coverage to consider an allele match")
    parser.add_argument("-o", "--output", help="Output")
    args = parser.parse_args()
    main(args)
