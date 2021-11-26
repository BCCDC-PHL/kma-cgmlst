#!/usr/bin/env python3

import argparse
import csv

def main(args):
    alleles = None
    sample_id = None
    with open(args.cgmlst, 'r') as f:
        next(f) # skip header
        line = f.readline()
        line_split = line.strip().split(args.delimiter)
        sample_id = line_split[0]
        alleles = line_split[1:]

    num_alleles = len(alleles)
    num_called_alleles = len([x for x in alleles if x != '-'])
    output_fields = [
        'sample_id',
        'called_alleles',
        'total_loci',
        'percent_called',
    ]
    print(','.join(output_fields))
    print(','.join([
        sample_id,
        str(num_called_alleles),
        str(num_alleles),
        str(round(float(num_called_alleles) / float(num_alleles) * 100, 3)),
        ]))
        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cgmlst')
    parser.add_argument('-d', '--delimiter', default=',')
    args = parser.parse_args()
    main(args)
