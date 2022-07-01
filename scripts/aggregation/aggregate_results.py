#!/usr/bin/env python

import pandas as pd
import glob
import argparse
import os

def arg_parser():
    """Handles the argument in and output on the command line, returns the
    arguments given by the user"""
    argp = argparse.ArgumentParser(description="Aggregates all the annotation tables to a single overview")
    argp.add_argument("-i", "--inputDir", default=str,
                      help="""location containing all annot_table.txt files to parse""")
    argp.add_argument("-o", "--outputDir", default=str, nargs="?",
                      help="""output file of the program, default = stdout""")
    return argp.parse_args()

def main():

    args = arg_parser()

    ls = glob.glob(f"{args.inputDir}/*.res")
    print(ls)

    frames = []

    # loop and open in panda
    for f in ls:
        print(f)
        sampleName = os.path.basename(f)
        print(sampleName)
        try:
            df = pd.read_csv(f, sep='\t', engine='python', comment='##')
                # add samplename column
            df['sampleName'] = sampleName
            
        except pd.errors.EmptyDataError:
            continue
        frames.append(df)

    # concat the files
    df_concat = pd.concat(frames)
    # filter to get quality mutations
    df_concat = df_concat[df_concat["Depth"] >= 1]
    df_concat = df_concat[df_concat["fragmentCount"] >= 100]
    df_concat = df_concat[df_concat["Template_Coverage"] >= 50]
    
    # create aggregation
    df_final = df_concat.pivot(index="#Template",
                                columns="sampleName",
                                values="Depth")

    # store result
    df_final.to_csv(f"{args.outputDir}/AGGREGATION_results.txt", sep="\t")

if __name__ == "__main__":
    # load arguments set global workdir
    # args = parse_args()
    # fill_html(args)
    print("Start")
    main()
    print("Finished")
