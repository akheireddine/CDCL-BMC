#! /usr/bin/python3

import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

from cycler import cycler
from tabulate import tabulate
from constants import *
from utils     import *

def main():
    verbose = True

    parser = argparse.ArgumentParser(description='ranking EDACC CSV')
    parser.add_argument('filename_edacc_csv', metavar='<edacc_csv>',
                        help="File containing the benchmark")

    parser.add_argument('--interactive', dest='is_interactive',
                        action="store_true",
                        help='cactus plot interactive')

    parser.add_argument('--timeout', dest='timeout',
                        type=int,
                        help='timeout')

    parser.add_argument('-o', '--output', dest='output',
                        default="cactus.pdf",
                        help='output file of cactus plot')

    parser.add_argument('--title', dest='title',
                        default=None,
                        help='titleof cactus plot')

    parser.add_argument('-c', '--cumsum', dest='is_cumsum',
                        action="store_true",
                        help='output with cumulative time')

    parser.add_argument('-m', '--mezcal', dest='is_mezcal',
                        action="store_true",
                        help='output in mezcal form')


    args = parser.parse_args()

    filename_edacc_csv = args.filename_edacc_csv
    df = pd.read_csv(filename_edacc_csv)


    limit = None
    if args.timeout:
        limit = args.timeout

    cumsum = False
    if args.is_cumsum:
        cumsum = args.is_cumsum

    mezcal = False
    if args.is_mezcal:
        mezcal = args.is_mezcal


    cactus_plot(df, limit=limit, output=args.output, mezcal=mezcal,
                cumsum=cumsum, title=args.title)

    if args.is_interactive:
        plt.show()


def cactus_plot(df, output='cactus.pdf', mezcal=False, cumsum=False,
                limit=None, title=None, key=TIME_KEY):
    ax = None

    # Guess limit
    if limit == None:
        limit = guess_limit(df, key)

    my_color_list  = ['blue_1', 'red_1', 'green_1', 'orange_1', 'brown_2',
                      'pink_1', 'purple_1', 'silver']
    my_marker_list = ['x', '^', 'o', '+', 'v', '>', '<', '*']

    plt.rc('axes', prop_cycle=(cycler('color', my_color_list) +
                               cycler('marker', my_marker_list)))


    for solver in column_no_duplicate(df, SOLVER_KEY):
        is_solver = df[SOLVER_KEY] == solver

        is_sat     = df[RESULT_KEY] == SAT
        is_unsat   = df[RESULT_KEY] == UNSAT
        is_complete = (is_sat) | (is_unsat)

        full   = df[(is_solver) & (is_complete)][TIME_KEY].sort_values()

        if cumsum:
            full = full.cumsum().reset_index(drop=True)
        else:
            full = full.reset_index(drop=True)

        if mezcal:
            plt.plot(full, full.index, label=solver, markerfacecolor="None",
                     linewidth=.5, alpha=0.8)
            plt.xlabel("cumulative time (s)" if cumsum else "time (s)")
            plt.ylabel("#solved instances")

        else:
            plt.plot(full.index, full, label=solver, markerfacecolor="None",
                     linewidth=.5, alpha=0.8)
            plt.xlabel("#solved instances")
            plt.ylabel("cumulative time (s)" if cumsum else "time (s)")


    plt.grid(color='gray', linestyle='--', linewidth=1, alpha=0.6)
    plt.title(title)

    plt.legend(numpoints=1, markerfirst=False, loc="best")

    plt.savefig(output, transparent=False, bbox_inches='tight')
    print("Output in", output)

    plt.rcParams.update(plt.rcParamsDefault)

if __name__ == '__main__':
    main()
