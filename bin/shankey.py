#!/usr/bin/env python

######
INFO = "Shankey plot creation"
__version__ = "0.1"
######

"""
Title:          shankeyplot.py
Author:         J.P.M. Coolen
Date:           29-06-2022 (dd-mm-yyyy)
Description:    Creates skakey plot
"""

import pandas as pd
import plotly
import argparse
import os
import sys

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--res", type=str, required=True,
                        help="KMA res output"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

# extract the data from KMA.res to extract the lvls
def createdf(input):
    df = pd.read_csv(input, sep="\t")

    print(df["#Template"])

    lvls = df["#Template"].str.split(";", expand = True)

    df["lvl1"] = "Bacteria"
    df["lvl2"] = lvls[1]
    df["lvl3"] = lvls[2]
    df["lvl4"] = lvls[3]
    df["lvl5"] = lvls[4]
    df["lvl6"] = lvls[5]
    df["lvl7"] = lvls[6]

    return df

#source:https://medium.com/kenlok/how-to-create-sankey-diagrams-from-dataframes-in-python-e221c1b4d6b0
def genSankey(df, cat_cols=[], value_cols='', title='Sankey Diagram'):
    # maximum of 9 value cols -> 9 colors
    colorPalette = ['#4B8BBE', '#306998', '#FFE873', '#FFD43B', '#646464', '#218380', '#FBB13C' , '#D81159', '#73D2DE']
    labelList = []
    colorNumList = []
    for catCol in cat_cols:
        labelListTemp = list(set(df[catCol].values))
        colorNumList.append(len(labelListTemp))
        labelList = labelList + labelListTemp

    # remove duplicates from labelList
    labelList = list(dict.fromkeys(labelList))

    # define colors based on number of levels
    colorList = []
    for idx, colorNum in enumerate(colorNumList):
        colorList = colorList + [colorPalette[idx]] * colorNum

    # transform df into a source-target pair
    for i in range(len(cat_cols) - 1):
        if i == 0:
            sourceTargetDf = df[[cat_cols[i], cat_cols[i + 1], value_cols]]
            sourceTargetDf.columns = ['source', 'target', 'count']
        else:
            tempDf = df[[cat_cols[i], cat_cols[i + 1], value_cols]]
            tempDf.columns = ['source', 'target', 'count']
            sourceTargetDf = pd.concat([sourceTargetDf, tempDf])
        sourceTargetDf = sourceTargetDf.groupby(['source', 'target']).agg({'count': 'sum'}).reset_index()

    # add index for source-target pair
    sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: labelList.index(x))
    sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: labelList.index(x))

    # creating the sankey diagram
    data = dict(
        type='sankey',
        node=dict(
            pad=15,
            thickness=20,
            line=dict(
                color="black",
                width=0.5
            ),
            label=labelList,
            color=colorList
        ),
        link=dict(
            source=sourceTargetDf['sourceID'],
            target=sourceTargetDf['targetID'],
            value=sourceTargetDf['count']
        )
    )

    layout = dict(
        title=title,
        font=dict(
            size=10
        )
    )

    fig = dict(data=[data], layout=layout)
    return fig

if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()

    df = createdf(args.res)
    fig = genSankey(df, cat_cols=['lvl1', 'lvl2', 'lvl3', 'lvl4', 'lvl5', 'lvl6', 'lvl7'], value_cols='fragmentCount', title='16S')
    #plotly.offline.plot(fig, validate=False)

    plotly.io.write_image(fig, "shankey.png", scale=1, width=800, height=400)
    plotly.io.write_image(fig, "shankey.svg", scale=1, width=1200, height=600)

    print("Finished")
