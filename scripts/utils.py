from constants import *

import pandas as pd
import numpy as np

INSTANCE_KEY = EDACC_INSTANCE
SOLVER_KEY   = EDACC_SOLVER_CONFIGURATION
RESULT_KEY   = EDACC_RESULT_CODE
TIME_KEY     = EDACC_WALL_TIME
LTL_KEY      = EDACC_LTL_TYPE

def column_no_duplicate(df, column):
    assert(column in df.columns)
    return np.unique(df[column])


def guess_limit(df, key=TIME_KEY):
    limit = 0
    for instance in column_no_duplicate(df, INSTANCE_KEY):
        is_instance = df[INSTANCE_KEY] == instance

        is_sat     = df[RESULT_KEY] == SAT
        is_unsat   = df[RESULT_KEY] == UNSAT

        is_complete = (is_sat) | (is_unsat)

        if len(df[(is_instance) & (is_complete)].index) > 0:
            limit = max(limit, max(df[(is_instance) & (is_complete)][key]))
    return limit

def keep_only_solvers(df, solvers):
    for solver in solvers:
        assert(solver in column_no_duplicate(df, SOLVER_KEY))

    for solver in column_no_duplicate(df, SOLVER_KEY):
        is_solver = df[SOLVER_KEY] == solver
        if solver not in solvers:
            df.drop(df[is_solver].index, inplace=True)
    return df

def keep_only_ltl(df, ltls):
    for ltl in ltls:
        assert(ltl in column_no_duplicate(df, LTL_KEY))

    for ltl in column_no_duplicate(df, LTL_KEY):
        is_ltl = df[LTL_KEY] == ltl
        if ltl not in ltls:
            df.drop(df[is_ltl].index, inplace=True)
    return df

def keep_only_bound(df, bound):
    bound=str(bound)
    for inst in column_no_duplicate(df, INSTANCE_KEY):
        if bound not in inst:
            is_inst = ~df[INSTANCE_KEY].str.lower().str.contains(bound.lower())
            df.drop(df[is_inst].index, inplace=True)
    return df

def keep_only_model(df, models):
    for inst in column_no_duplicate(df, INSTANCE_KEY):
        keep=False
        for m in models:
           if m in inst:
                keep = True
                break
        if not(keep) :
            is_inst = df[INSTANCE_KEY] == inst
            df.drop(df[is_inst].index, inplace=True)
    return df

