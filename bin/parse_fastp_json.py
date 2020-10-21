#!/usr/bin/env python

import argparse
import csv
import json
import os
import sys


def parse_fastp_json(path_to_fastp_json):
    with open(path_to_fastp_json) as f:
        complete_fastp_report = json.load(f)
    return complete_fastp_report

def flatten_and_select_fields_from_summary(parsed_fastp):
    summary = {}
    # TODO: re-write this as a loop that iterates through parsed_fastp['summary']
    summary['total_reads_before_filtering'] = parsed_fastp['summary']['before_filtering']['total_reads']
    summary['total_bases_before_filtering'] = parsed_fastp['summary']['before_filtering']['total_bases']
    summary['q20_bases_before_filtering'] = parsed_fastp['summary']['before_filtering']['q20_bases']
    summary['q30_bases_before_filtering'] = parsed_fastp['summary']['before_filtering']['q30_bases']
    summary['q20_rate_before_filtering'] = parsed_fastp['summary']['before_filtering']['q20_rate']
    summary['q30_rate_before_filtering'] = parsed_fastp['summary']['before_filtering']['q30_rate']
    summary['total_reads_after_filtering'] = parsed_fastp['summary']['after_filtering']['total_reads']
    summary['total_bases_after_filtering'] = parsed_fastp['summary']['after_filtering']['total_bases']
    summary['q20_bases_after_filtering'] = parsed_fastp['summary']['after_filtering']['q20_bases']
    summary['q30_bases_after_filtering'] = parsed_fastp['summary']['after_filtering']['q30_bases']
    summary['q20_rate_after_filtering'] = parsed_fastp['summary']['after_filtering']['q20_rate']
    summary['q30_rate_after_filtering'] = parsed_fastp['summary']['after_filtering']['q30_rate']

    return summary
    

def main(args):

    fastp = parse_fastp_json(args.fastp_json)
    output_data = flatten_and_select_fields_from_summary(fastp)
    output_data['sample_id'] = args.sample_id
    
    output_fields = [
        "sample_id",
        "total_reads_before_filtering",
        "total_bases_before_filtering",
        "q20_bases_before_filtering",
        "q30_bases_before_filtering",
        "q20_rate_before_filtering",
        "q30_rate_before_filtering",
        "total_reads_after_filtering",
        "total_bases_after_filtering",
        "q20_bases_after_filtering",
        "q30_bases_after_filtering",
        "q20_rate_after_filtering",
        "q30_rate_after_filtering",
    ]
    
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, dialect='excel-tab')
    writer.writeheader()
    writer.writerow(output_data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", help="fastp json output")
    parser.add_argument("--fastp_json", help="fastp json output")
    args = parser.parse_args()
    main(args)
