#!/usr/bin/env python

import pandas as pd
import glob
import argparse
import os
import random
import string

def arg_parser():
    """Handles the argument in and output on the command line, returns the
    arguments given by the user"""
    argp = argparse.ArgumentParser(description="Aggregates all the annotation tables to a single overview")
    argp.add_argument("-i", "--inputDir", type=str,
                      help="""location containing all annot_table.txt files to parse""")
    argp.add_argument("--depth", type=int, default=1,
                      help=".=> number of depth required to report [default=1]")
    argp.add_argument("--fragmentCount", type=int, default=100,
                      help=".=> number of fragmentCount required to report [default=100]")
    argp.add_argument("--templateCoverage", type=int, default=50,
                      help=".=> number of templateCoverage required to report [default=50]")
    argp.add_argument("--abundance", type=float, default=0.01,
                      help=".=> number of abundance required to report after groupby [default=0.01]")
    argp.add_argument("-o", "--outputDir", type=str,
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

    # filter to get quality hits
    df_concat = df_concat[df_concat["Depth"] >= args.depth]
    df_concat = df_concat[df_concat["fragmentCount"] >= args.fragmentCount]
    df_concat = df_concat[df_concat["Template_Coverage"] >= args.templateCoverage]

    tax = df_concat["#Template"].str.split(";", n=6, expand=True)
    species = tax[6].str.split(" ", n=2, expand=True)

    # change None splits to ""
    species = species.replace([None], "")

    df_concat["#Template"] = tax[5] + ";" + species[0] + " " + species[1]

    # extract columns #Template, sampleName, abundance
    df_concat = df_concat[["#Template", "sampleName", "abundance"]]

    # sum the same identification on abundance
    df_concat = df_concat.groupby(["#Template", "sampleName"]).sum().reset_index()

    print(df_concat)

    df_concat = df_concat[df_concat["abundance"] >= args.abundance]

    # rename header
    df_concat = df_concat.rename(columns={"#Template": "X"})

    # create aggregation
    df_final = df_concat.pivot(index="X",
                                columns="sampleName",
                                values="abundance")

    # store result
    df_final.to_csv(f"{args.outputDir}/AGGREGATION_results_"
                    f"D{args.depth}_"
                    f"F{args.fragmentCount}_"
                    f"C{args.templateCoverage}_"
                    f"A{args.abundance}.txt", sep="\t")

if __name__ == "__main__":
    # load arguments set global workdir
    # args = parse_args()
    # fill_html(args)
    print("Start")
    main()
    print("Finished")
