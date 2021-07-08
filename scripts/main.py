
######
INFO = "Process fastq.gz files and process primerUMI for RC-PCR"
__version__ = "0.0.1"
######

from pysam import FastxFile
from collections import Counter
import pandas as pd
import numpy as np
import argparse
import os
import sys

## tested with these settings of fastp
# fastp -i 18S-03_S99_L001_R1_001.fastq.gz -o 18S-03_S99_L001_R1_001.fastq_fastp.gz \
#                                 -I 18S-03_S99_L001_R2_001.fastq.gz -O 18S-03_S99_L001_R2_001.fastq_fastp.gz \
#                                 --umi --umi_loc=per_read --umi_len=24 --umi_prefix=UMI --html 18S-03_S99.html \
#                                 --json 18S-03_S99.json --length_required 100

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="R1 or R2 fastq.gz file processed by fastp with UMI"),
    parser.add_argument("-p", "--primers", type=str, required=True,
                        help="primer file in fasta format"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

#read fastq files
def read_fastq_primerUMI_to_list(fn):

    UMI_list = []

    # read all the fastq reads and extract info
    with FastxFile(fn) as fh:
        for read in fh:
            UMI_list.append(collect_primer_as_UMI(read.name))

    return(UMI_list)

#collect UMI data
def collect_primer_as_UMI(header):
    UMI = header.split("UMI_")
    return(UMI[1])

#measure UMI length from UMI_list
def retrieve_UMI_length(UMI_list):
    return(len(UMI_list[0].split("_")[0]))

#make counttable
def make_counttable(UMI_list):
    data = Counter(UMI_list)
    df = pd.DataFrame.from_dict(data, orient='index', columns=['Count'])
    df = df.sort_values('Count', ascending=False)
    return(df)

#read_primer_file]
def read_primer_file(primers):
    primerseq_list = []
    primername_list =[]

    # read all the fastq reads and extract info
    with FastxFile(primers) as fh:
        for read in fh:
            primerseq_list.append(read.sequence)
            primername_list.append(read.name)

    return(primerseq_list, primername_list)

#match_primers_with_primerUMI
def match_primers(primerseqs, primernames, counttable, UMI_length):

    forward_name = []
    reverse_name = []

    # loop over counttable to check primers
    for (columnName, columnData) in counttable.iteritems():
        if columnName == 'Forward':
            for primerUMI in columnData.values:
                for num, primerseq in enumerate(primerseqs, start=0):
                    if primerUMI.startswith(primerseq[0:UMI_length]):
                        present = 'Correct'
                        break
                    else:
                        present = 'Incorrect'

                # assess if there is a match, store primer name
                if present == 'Correct':
                    forward_name.append(primernames[num])
                else:
                    forward_name.append('no_match')

        if columnName == 'Reverse':
            for primerUMI in columnData.values:
                for num, primerseq in enumerate(primerseqs, start=0):
                    if primerUMI.startswith(primerseq[0:UMI_length]):
                        present = 'Correct'
                        break
                    else:
                        present = 'Incorrect'

                # assess if there is a match, store primer name
                if present == 'Correct':
                    reverse_name.append(primernames[num])
                else:
                    reverse_name.append('no_match')

    # create dataframe of match
    fw_df = pd.DataFrame({'fw_name': forward_name})
    rv_df = pd.DataFrame({'rv_name': reverse_name})

    return(fw_df,
           rv_df)

def calculate_correct_portion(counttable):

    correct = []
    incorrect = []

    # loop and count correct pairs
    for index, row in counttable.iterrows():
        fw_split = row['fw_name'].split('_')
        rv_split = row['rv_name'].split('_')

        if fw_split[0] == 'F' and rv_split[0] == 'R' and fw_split[1] == rv_split[1]:
            correct.append(row['percentage'])
        else:
            incorrect.append(row['percentage'])

    return((sum(correct)))

def extract_samplename(fullpath):

    # example of file name format
    # 18S-03_S99_L001_R1_001.fastq_fastp.gz
    # take first part before _
    # 18S-03
    basename = os.path.basename(fullpath)
    split = basename.split("_")
    samplename = split[0]

    return(samplename)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print("RC-PCR depth")
    args = parse_args()

    print(args)

    # obtain samplename for args.input which is full path to input fastq.gz
    samplename = extract_samplename(args.input)

    #read primers to list
    primerseq_list, primername_list = read_primer_file(args.primers)

    # read primerUMI headers in reads of fastq file
    UMI_list = read_fastq_primerUMI_to_list(args.input)

    UMIlen = retrieve_UMI_length(UMI_list)

    counttable = make_counttable(UMI_list)
    counttable.reset_index(level=0, inplace=True)

    # split the primersUMI to forward and reverse
    counttable = pd.concat([counttable["index"].str.split("_", n=1, expand=True), counttable], axis=1)
    counttable.columns = ["Forward", "Reverse", "primersUMI", "Count"]

    # retrieve matched fw and rv primers
    fw, rv = match_primers(primerseq_list, primername_list, counttable, UMIlen)

    # create new counttable
    counttable = pd.concat([counttable, fw, rv], axis=1)

    total = np.sum(counttable.loc[:, 'Count'].values)

    counttable['percentage'] = counttable.loc[:, 'Count'].values / total * 100

    perc_correct = calculate_correct_portion(counttable)

    print(f"UMI length: {UMIlen}")
    print(f"Percentage correct pairs: {perc_correct}")
    print(f"Number of reads: {sum(counttable['Count'])}")

    outputFile = os.path.join(args.outputDir, f"{samplename}_UMI_counttable.csv")
    print(f"File saved: {outputFile}")

    # save counttable to disk
    counttable.to_csv(outputFile)
    print("Finished")