#!/usr/bin/env python

######
INFO = "Convert results to PDF report"
__version__ = "0.2"
######

"""
Title:          final_report.py
Author:         J.P.M. Coolen
Date:           10-03-2022 (dd-mm-yyyy)
Description:    Convert results to PDF report
"""

import os
import argparse
import pandas as pd
import datetime

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sampleName", type=str, required=True,
                        help="name of sequence sample"),
    parser.add_argument("--annotation", type=str, required=False,
                        help="location to annotation file"),
    parser.add_argument("--params", type=str, required=False,
                        help="location to parameters.txt file"),
    parser.add_argument("--blast", type=str, required=False,
                        help="location to abricate blast result file"),
    parser.add_argument("--kma", type=str, required=False,
                        help="location to kma result file .res"),
    parser.add_argument("--shankey", type=str, required=False,
                        help="location to shankey.svg file"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

def fill_html(args):
    '''
        Code to fill in the placeholders in the html
        and generate a html and pdf

        :params JSON: JSON object containing all the results
        :params outputDir: directory to store the results

        :out pdf: pdf report of the results
        :out html: html report of the results
    '''

    import matplotlib
    matplotlib.use('Agg')
    from weasyprint import HTML
    from jinja2 import Environment, FileSystemLoader

    print('Start Filling')

    localdir = os.path.dirname(os.path.realpath(__file__))

    # create and render html file with tables
    env = Environment(loader=FileSystemLoader(localdir))
    template = env.get_template('report/final_report_template.html')

    # location of logo
    logo = os.path.join(localdir, "report/logo.png")
    logo = logo.replace(' ','%20')

    # date
    date = datetime.datetime.now().strftime('%Y%m%d%H%M%S')

    # load parameters.txt file
    params_df = pd.read_csv(args.params, sep='\t')

    # load blast .txt result file of abricate
    blast_df = pd.read_csv(args.blast, sep='\t')
    try:
        blast_df = blast_df.sort_values(by=['%COVERAGE'], ascending=False)
    except KeyError:
        print("no result")

    # load kma .res result file of kma
    kma_df = pd.read_csv(args.kma, sep='\t')
    kma_df = kma_df.sort_values(by=['q_value'], ascending=False)
    kma_df = kma_df.head(10)

    # fill html
    template_vars = {
        # pretty things
        "logo": logo,
        "version": __version__,

        # general info
        "sampleName": args.sampleName,
        "date": date,

        # classification results
        "kma": kma_df.to_html(index=False, header=True),

        # shankey plot of classification results
        "shankey": args.shankey,

        "blast": blast_df.to_html(index=False, header=True),

        # parameters
        "parameters": params_df.to_html(index=False, header=True),

    }

    # output pdf
    outfile = os.path.join(args.outputDir, '{}.pdf'.format(args.sampleName))

    # render html and write
    html_out = template.render(template_vars)
    with open(os.path.join(args.outputDir, '{}.html'.format(args.sampleName)), 'w') as html_file:
        html_file.write(html_out)

    # save html as pdf to disc
    HTML(string=html_out, base_url=__file__).write_pdf(outfile,
                                                       stylesheets=[os.path.join(localdir, 'report/style.css')])
    
if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    fill_html(args)
    print("Finished")
