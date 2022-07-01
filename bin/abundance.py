#!/usr/bin/env python

######
INFO = "anbundance.py"
__version__ = "0.1"
######

"""
Title:          abundance.py
Author:         J.P.M. Coolen
Date:           30-06-2022 (dd-mm-yyyy)
Description:    calculate abundance
"""

import pandas as pd
import argparse
import os
import re

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mapstat", type=str, required=True,
                        help=",mapstat output of KMA"),
    parser.add_argument("--res", type=str, required=True,
                        help=",.res output of KMA"),
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

def read_res(res):
    df = pd.read_csv(res, sep='\t', index_col=0, encoding='latin1')

    return(df)

def calc_abundance(input, res):
    with open(input) as mapfile:
        fragments_line =mapfile.readlines()[3]
    total_frags = re.split(r'(\t|\n)' ,fragments_line)[2]
    df_stats = pd.read_csv(input, sep='\t', index_col=0, header = 6, encoding='latin1')
    res['fragmentCount'] = df_stats['fragmentCount']
    res['abundance'] = df_stats['fragmentCount'] / int(total_frags)

    return(res)

if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()

    res = read_res(args.res)
    res = calc_abundance(args.mapstat, res)

    res.to_csv(args.res, sep="\t")

    print("Finished")