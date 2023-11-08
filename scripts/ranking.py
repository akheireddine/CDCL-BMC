#! /usr/bin/python3
# Ranking for SAT benchmark

import csv
import pandas as pd
import numpy as np
import argparse

from tabulate import tabulate
from utils     import *
from constants import *


def main():
    verbose = True

    parser = argparse.ArgumentParser(description='ranking EDACC CSV')
    parser.add_argument('filename_edacc_csv', metavar='<edacc_csv>',
                     help="File containing the benchmark")

    parser.add_argument('--remove-memory', dest='rm_memout',
                        action="store_true",
                        help='remove memory exceeded instances')

    parser.add_argument('--remove-unknown', dest='rm_unknown',
                        action="store_true",
                        help='remove unknown instances')

    parser.add_argument('--remove-timeout', dest='rm_timeout',
                        action="store_true",
                        help='remove timeout instances')

    parser.add_argument('--add-vbs', dest='add_vbs',
                        action="store_true",
                        help='remove unknown instances')

    parser.add_argument('--solvers', nargs='?', dest='solvers',
                        help='list of wanted solvers')

    parser.add_argument('--query', nargs='?', dest='query',
                        help='filter df with query in form: '
                        'solver column operator value')

    parser.add_argument('--timeout', dest='timeout',
                        type=int, default=72000,
                        help='timeout')

    parser.add_argument('--filter-output', dest='filter_output',
                        help='filename of filtered output csv')

    parser.add_argument('--rename', nargs='?', dest='rename',
                        help='rename solver in format '
                        'solver=new_name')

    parser.add_argument('--ltl', nargs='?', dest='ltl',
                        help='list of wanted ltl: Recurrence,Persistence,Reactivity,...')
    parser.add_argument('--bound', nargs='?', dest='bound',
                        help='list of wanted bound: 20,40,60,80,100,200,...,4000')

    parser.add_argument('--model', nargs='?', dest='model',
                        help='list of wanted model: adding,krebs,...')
    args = parser.parse_args()

    filename_edacc_csv = args.filename_edacc_csv
    df = pd.read_csv(filename_edacc_csv)
    print(SOLVER_KEY)
    solvers = column_no_duplicate(df, SOLVER_KEY)
    if args.solvers != None:
        solvers = [s.strip() for s in args.solvers.split(',')]
    df = keep_only_solvers(df, solvers)

    if args.ltl != None:
        ltl = column_no_duplicate(df, LTL_KEY)
        ltl = [s.strip() for s in args.ltl.split(',')]
        df = keep_only_ltl(df, ltl)

    if args.bound != None:
        bound = [s.strip() for s in args.bound.split(',')]
        df = keep_only_bound(df, bound[0])

    if args.model != None:
        model = [s.strip() for s in args.model.split(',')]
        df = keep_only_model(df, model)

    if args.rm_memout:
        df = remove_memory(df)
    if args.rm_unknown:
        df = remove_unknown(df)
    if args.rm_timeout:
        df = remove_timeout(df)

    if args.rename:
        renames = [s.strip() for s in args.rename.split(',')]
        for r in renames:
            df = rename_solver(df, r)


    if args.query:
        queries = [s.strip() for s in args.query.split(',')]
        for query in queries:
            remove_with_query(df, query)

    if args.filter_output != None:
        out = df.to_csv(index=False, quotechar='"', quoting=csv.QUOTE_ALL)
        f = open(args.filter_output, "w+")
        print(out, file=f)
        f.close()
        if verbose:
            print("Output CSV written in", args.filter_output)


    df = add_parX(df, timeout=args.timeout)
    df = add_cti(df)

    if args.add_vbs:
        df = add_vbs_solver(df)

    ranking(df, verbose)

def add_parX(df, X=2, timeout=5000):
    is_not_sat     = df[RESULT_KEY] != SAT
    is_not_unsat   = df[RESULT_KEY] != UNSAT

    df[PAR2] = np.where((is_not_sat) & (is_not_unsat), timeout * X,
                          df[TIME_KEY])
    return df


def add_cti(df):
    num_solvers = len(column_no_duplicate(df, SOLVER_KEY))

    df[CTI] = df[TIME_KEY]

    for instance in column_no_duplicate(df, INSTANCE_KEY):
        is_instance = df[INSTANCE_KEY] == instance

        is_sat     = df[RESULT_KEY] == SAT
        is_unsat   = df[RESULT_KEY] == UNSAT

        table  = df[(is_instance) & ((is_sat) | (is_unsat))]
        if len(table.index) != num_solvers:
            df.loc[is_instance, CTI] = 0

    return df

def add_vbs_solver(df):
    output = []
    for instance in column_no_duplicate(df, INSTANCE_KEY):
        is_instance = df[INSTANCE_KEY] == instance

        is_sat     = df[RESULT_KEY] == SAT
        is_unsat   = df[RESULT_KEY] == UNSAT

        is_complete = (is_sat) | (is_unsat)

        if df[(is_instance) & (is_complete)][TIME_KEY].empty:
            idx = df[is_instance][TIME_KEY].idxmin()
        else:
            idx = df[(is_instance) & (is_complete)][TIME_KEY].idxmin()

        row = df.loc[idx].copy()
        row[SOLVER_KEY]   = "VBS"

        output += [row]

    return df.append(output)

def remove_with(df, column, at_leat_one=True, all_solvers=False):
    num_solvers = len(column_no_duplicate(df, SOLVER_KEY))
    for instance in column_no_duplicate(df, INSTANCE_KEY):
        is_instance = df[INSTANCE_KEY] == instance
        is_column  = df[RESULT_KEY] == column

        num_column = len(df[(is_instance) & (is_column)].index)
        if (all_solvers and num_column == num_solvers) or \
           (at_leat_one and num_column > 0):
            df.drop(df[is_instance].index, inplace=True)

    return df

def remove_memory(df):
    return remove_with(df, MEMORY, at_leat_one=True, all_solvers=False)

def remove_unknown(df):
    return remove_with(df, UNKNOWN, at_leat_one=True, all_solvers=False)

def remove_timeout(df):
    return remove_with(df, TIMEOUT, at_leat_one=False, all_solvers=True)

def remove_with_query(df, query_params):
    print("colums : ",df.columns,"\n query "+query_params+",  split :",query_params.split())
    solver, column, operator, value = query_params.split()
    if column == "ResultCode" :
        column = "Result Code"
    if column == "WallTime" :
        column = "Wall Time"
    # Some checks
    assert(solver and column and operator and value)
    assert(solver in column_no_duplicate(df, SOLVER_KEY))

    assert(column in df.columns)
    for c in operator: assert(c in "<>=")

    query = "(df[\"{col}\"] {op} {value}) | (df[\"{col}\"].isnull())".format(
        col=column, op=operator, value=value)

    for instance in column_no_duplicate(df, INSTANCE_KEY):
        is_instance = df[INSTANCE_KEY] == instance
        is_solver = df[SOLVER_KEY] == solver
        is_query = eval(query)

        if len(df[(is_instance) & (is_solver) & (is_query)].index) > 0:
            df.drop(df[is_instance].index, inplace=True)

def rename_solver(df, rename_params):
    solver, new_name = rename_params.split('=')
    assert(solver in column_no_duplicate(df, SOLVER_KEY))
    df[SOLVER_KEY] = df[SOLVER_KEY].replace(solver, new_name)
    return df

def ranking(df, verbose):
    num_solvers = len(column_no_duplicate(df, SOLVER_KEY))
    TOTAL = "TOTAL (" + str(int(len(df.index) / num_solvers)) + ")"
    keys = ["Solvers", "UNSAT", "SAT", TOTAL, "TIMEOUT", "MEMORY", "UNKNOWN"]

    if PAR2 in df.columns:
        keys += [PAR2]

    if CTI in df.columns:
        num = int(df[CTI].astype(bool).sum(axis=0) / num_solvers)
        name = CTI + " (" + str(num)  + ")"
        keys += [name]

    datas = pd.DataFrame(columns=keys)
    fmt='psql'

    ERR_WS     =  -11111111
    ERR_WS_MSG = 'wrong sat'

    ERR_WU     =  -2222222222
    ERR_WU_MSG = 'wrong unsat'

    ERR_AL     =  -33333333
    ERR_AL_MSG = 'wrong all'

    is_sat     = df[RESULT_KEY] == SAT
    is_unsat   = df[RESULT_KEY] == UNSAT
    is_timeout = df[RESULT_KEY] == TIMEOUT
    is_unknown = df[RESULT_KEY] == UNKNOWN
    is_memout  = df[RESULT_KEY] == MEMORY
    is_wrong   = df[RESULT_KEY] == WS


    # Compute inconsistent solvers
    # a solver answer UNSAT and the instance is SAT
    inconsitent_solvers = set()
    for instance in column_no_duplicate(df, INSTANCE_KEY):
        is_instance = df[INSTANCE_KEY] == instance

        sat  = df[(is_instance) & (is_sat)]
        unsat = df[(is_instance) & (is_unsat)]

        num_sat  = sat[RESULT_KEY].count()
        num_unsat = unsat[RESULT_KEY].count()

        if num_sat > 0 and num_unsat > 0:
            solvers = unsat[SOLVER_KEY].tolist()
            instance = unsat[EDACC_INSTANCE].tolist()
            inconsitent_solvers |=  set(solvers)
            if verbose:
                print(*solvers, sep=',', end='')
                print(' has wrong solution on UNSAT instance', *instance)

    # Compute number of success instance
    for solver in column_no_duplicate(df, SOLVER_KEY):
        print("Solver ",solver)
        is_solver = df[SOLVER_KEY] == solver

        unsat   = df[(is_solver) & (is_unsat)]  [RESULT_KEY].count()
        sat     = df[(is_solver) & (is_sat)]    [RESULT_KEY].count()
        total   = sat + unsat
        timeout = df[(is_solver) & (is_timeout)][RESULT_KEY].count()
        unknown = df[(is_solver) & (is_unknown)][RESULT_KEY].count()
        memout  = df[(is_solver) & (is_memout)] [RESULT_KEY].count()


        wrong   = df[(is_solver) & (is_wrong)][EDACC_INSTANCE]

        if verbose and  wrong.count() > 0:
                print(solver, "has wrong solution on SAT", *wrong.tolist())

        # Error handling
        if wrong.count() > 0 and solver in inconsitent_solvers:
            data =  [solver]
            for c in range(len(keys)-1):
                data += [ERR_AL]
        elif wrong.count() > 0:
            data =  [solver]
            for c in range(len(keys)-1):
                data += [ERR_WS]
        elif solver in inconsitent_solvers:
            data =  [solver]
            for c in range(len(keys)-1):
                data += [ERR_WU]
        else:
            data =  [solver, unsat, sat, total, timeout, memout, unknown]
            if PAR2 in df.columns:
                par2 = df[(is_solver)][PAR2].sum()
                data += [par2]

            if CTI in df.columns:
                cti = df[(is_solver)][CTI].sum()
                data += [cti]

        datas.loc[-1] = data
        datas.index = datas.index + 1

    datas = datas.sort_values(TOTAL, ascending=False)
    datas = datas.reset_index(drop=True)

    # Pretty print the output
    output = tabulate(datas, headers=keys, tablefmt=fmt, floatfmt='.0f')
    output = output.replace(str(ERR_WS), ERR_WS_MSG)
    output = output.replace(str(ERR_WU), ERR_WU_MSG)
    output = output.replace(str(ERR_AL), ERR_AL_MSG)

    print(output)

if __name__ == '__main__':
    main()
