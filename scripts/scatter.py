#! /usr/bin/python3

import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

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
                        help='scatter plot interactive')

    parser.add_argument('--timeout', dest='timeout',
                        type=int,
                        help='timeout')

    parser.add_argument('-o', '--output', dest='output',
                        default="scatter.pdf",
                        help='output file of scatter plot')

    parser.add_argument('--title', dest='title',
                        default=None,
                        help='titleof scatter plot')

    parser.add_argument('--key', dest='key',
                        help='scatter plot in this key (default wall time)')

    parser.add_argument('--logscale', dest='is_logscale',
                        action="store_true",
                        help='scatter plot in log scale')

    parser.add_argument('--s1', dest='s1', metavar='<s1>',
                        help='first solver')
    parser.add_argument('--s2', dest='s2', metavar='<s2>',
                        help='second solver')


    args = parser.parse_args()

    filename_edacc_csv = args.filename_edacc_csv
    df = pd.read_csv(filename_edacc_csv)

    key = TIME_KEY
    if args.key:
        key = args.key
    df[key] = pd.to_numeric(df[key], errors='coerce')

    if not args.s1 or not args.s2:
        print("Need solvers names option --s1 and --s2")
        exit()

    solver_one = args.s1
    solver_two = args.s2

    limit = None
    if args.timeout:
        limit = args.timeout

    scatter_plot(df, solver_one, solver_two, limit=limit,
                 output=args.output, title=args.title, key=key,
                 is_logscale=args.is_logscale)

    if args.is_interactive:
        plt.show()


def scatter_plot(df, s1, s2, output="scatter.pdf", limit=None, title=None,
                 key=TIME_KEY, is_logscale=False):
    df = df.copy()

    solvers = [s1, s2]
    keep_only_solvers(df, solvers)

    ax = None

    # Guess limit
    if limit == None:
        limit = guess_limit(df, key)

    is_fail = (df[RESULT_KEY] != SAT) & (df[RESULT_KEY] != UNSAT)
    df.loc[is_fail, key] = limit * 1.1

    table_sat   = []
    table_unsat = []

    for instance in column_no_duplicate(df, INSTANCE_KEY):
        is_instance = df[INSTANCE_KEY] == instance

        is_sat     = df[RESULT_KEY] == SAT
        is_unsat   = df[RESULT_KEY] == UNSAT

        is_s1      = df[SOLVER_KEY] == s1
        is_s2      = df[SOLVER_KEY] == s2

        val_s1 = df[(is_instance) & (is_s1)][key].item()
        val_s2 = df[(is_instance) & (is_s2)][key].item()

        if len(df[(is_instance) & (is_sat)].index) > 0:
            table_sat.append([val_s1, val_s2])
        elif len(df[(is_instance) & (is_unsat)].index) > 0:
            table_unsat.append([val_s1, val_s2])

    if len(table_unsat) > 0:
        df_unsat = pd.DataFrame(table_unsat, columns=[s1, s2])
        ax = df_unsat.plot.scatter(s1, s2, ax=ax,
        marker="+", color="red_1", linestyle="None", label="UNSAT",
                                   alpha=0.85)

    if len(table_sat) > 0:
        df_sat = pd.DataFrame(table_sat, columns=[s1, s2])
        ax = df_sat.plot.scatter(s1, s2, ax=ax,
        marker="x", color="blue_1", linestyle="None", label="SAT", alpha=0.85)


    plt.rcParams.update({'font.size': 13})

    plt.xlim(0, limit * 1.15)
    plt.ylim(0, limit * 1.15)

    if ax != None:
        ax.grid(color='gray', linestyle='--', linewidth=1, alpha=0.6)
        if not is_logscale:
            ax.set_aspect('equal')

    if  is_logscale:
        plt.yscale('symlog')
        plt.xscale('symlog')

    plt.plot([0, limit * 1.05], [0, limit * 1.05],
             color="black", linestyle="--", alpha=0.6)
    plt.plot([0, limit * 1.15], [limit * 1.05, limit * 1.05],
             linestyle="--", color="black", alpha=0.6)
    plt.plot([limit * 1.05, limit * 1.05], [limit * 1.15, 0],
             linestyle="--", color="black", alpha=0.6)
    plt.legend(numpoints=1, ncol=2, markerfirst=False, loc="center",
               bbox_to_anchor=(0.5, -0.20))


    if title != None:
        plt.title(title, loc='center')

    plt.savefig(output, transparent=False, bbox_inches='tight')
    plt.rcParams.update(plt.rcParamsDefault)

    print("Output in", output)

if __name__ == '__main__':
    main()
