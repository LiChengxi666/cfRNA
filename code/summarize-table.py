#!/usr/bin/env python
import argparse
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("summarize count tables")

def main():
    parser = argparse.ArgumentParser(description='Summarize counts tables into a single matrix')
    parser.add_argument('--indir', '-i', type=str, required=True, help='directory contains input tables')
    parser.add_argument('--formatter', type=str,default='{}',help="How to map sample id to file name of input table")
    parser.add_argument('--sample-ids', '-s', type=str, required=True, help='sample ids to include in the output matrix')
    parser.add_argument('--value-field', '-vf', type=int, required=True,help='index of field/column (start from 0) in the table to summarize as feature value')
    parser.add_argument('--row-field', '-rf', type=int, required=True, help='index of field/column (start from 0) in the table to summarize as feature name')
    parser.add_argument('--row-name', '-rn', type=str, default="feature", help='index name in the output matrix')
    parser.add_argument('--comment', '-c', default="#", type=str, help='line starts with this will be skipped')
    parser.add_argument('--first-line', '-f', type=int, default=0, help='only consider records after this line (start from first line (line 0) by default). all lines will be considered by default ')
    parser.add_argument('--output','-o', type=str, required=True, help='ouput matrix')
    parser.add_argument('--fillna',action="store_true", help='if specificied, fill na values with zero')
    args = parser.parse_args()
    records = []

    sample_ids = open(args.sample_ids).read().strip().split("\n")
    logger.info("Check inputs ...")
    not_presented = []
    for sample_id in sample_ids:
        path = os.path.join(args.indir, args.formatter.format(sample_id))
        if not os.path.exists(path):
            not_presented.append(path)
    if len(not_presented):
        logger.warning(f"{','.join(not_presented)} does not exists in input directory, exit.")
        sys.exit(1)
    logger.info("Load tables ...")
    for sample_id in tqdm(sample_ids):
        path = os.path.join(args.indir,args.formatter.format(sample_id))
        with open(path) as f:
            i = 0
            if args.first_line > 0:
                for xx in range(args.first_line):
                    i += 1
                    _ = next(f)
            for line in f:
                i += 1
                line = line.strip()
                if line.startswith(args.comment):
                    continue
                fields = line.split("\t")
                assert args.row_field < len(fields) and args.value_field < len(fields), f"in line {i}, only {len(fields)} available"
                feature_id, value = fields[args.row_field],fields[args.value_field]
                records.append((feature_id, sample_id, value))
    logger.info("Pivoting the table ...")
    table = pd.DataFrame.from_records(records)
    table.columns = [args.row_name,"sample_id","values"]
    matrix = table.pivot(index=args.row_name,columns="sample_id",values="values")
    if args.fillna is not None:
        matrix = matrix.fillna(0)
    matrix = matrix.astype(int)
    logger.info("Saving feature matrix ...")
    matrix.to_csv(args.output,sep="\t") 
    logger.info("All done .")


if __name__ == "__main__":
    main()

