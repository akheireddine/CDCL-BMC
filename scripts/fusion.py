#! /usr/bin/python3

import pandas as pd
import numpy as np

import csv
import numpy as np
import matplotlib.pyplot as plt
import argparse

from cycler import cycler
from tabulate import tabulate

import plotly.graph_objs as go
from plotly.offline import plot

from matplotlib.colors import ColorConverter


# Full stats on cnfs headers
YEAR               = 'year'
PATH               = 'path'
FILE               = 'file'
CATEGORY           = 'category'
MD5SUM             = 'md5sum'
IS_WELL_FORMED     = 'is_well_formed'
NUM_VARS           = 'num_vars'
NUM_CLAUSES        = 'num_clauses'
PERCENT_SYM_VARS   = 'percent_sym_vars'
IS_ONLY_INVOLUTION = 'is_only_involution'
IS_INVERTING       = 'is_inverting_perm'
NUM_ORBITS         = 'num_orbits'

NUM_ESBP = 'num_esbp'

# EDACC COLUMNS ID
EDACC_ID                     = 'ID'
EDACC_PRIORITY               = 'Priority'
EDACC_COMPUTE_QUEUE          = 'Compute Queue'
EDACC_COMPUTE_NODE           = 'Compute Node'
EDACC_COMPUTENODE_IP         = 'Compute Node IP'
EDACC_SOLVER                 = 'Solver'
EDACC_SOLVER_CONFIGURATION   = 'Solver Configuration'
EDACC_PARAMETERS             = 'Parameters'
EDACC_INSTANCE               = 'Instance'
EDACC_INSTANCE_MD5           = 'Instance MD5'
EDACC_RUN                    = 'Run'
EDACC_TIME                   = 'Time'
EDACC_WALL_TIME              = 'Wall Time'
EDACC_COST                   = 'Cost'
EDACC_SEED                   = 'Seed'
EDACC_STATUS                 = 'Status'
EDACC_RUN_TIME               = 'Run time'
EDACC_RESULT_CODE            = 'Result Code'
EDACC_CPU_TIME_LIMIT         = 'CPU Time Limit'
EDACC_WALL_CLOCK_TIME_LIMIT  = 'Wall Clock Time Limit'
EDACC_MEMORY_LIMIT           = 'Memory Limit'
EDACC_STACK_SIZE_LIMIT       = 'Stack Size Limit'
EDACC_SOLVER_OUTPUT          = 'Solver Output'
EDACC_LAUNCHER_OUTPUT        = 'Launcher Output'
EDACC_WATCHER_OTPUT          = 'Watcher Output'
EDACC_VERIFIER_OUTPUT        = 'Verifier Output'
EDACC_LTL_TYPE               = 'LTL'


# EDACC RESULT CODE
SAT     = 'SAT'
UNSAT   = 'UNSAT'
TIMEOUT = 'wall clock limit exceeded'
MEMORY  = 'memory limit exceeded'
UNKNOWN = 'unknown'
WS      = 'wrong solution'

# Custom
PAR2 = 'PAR2'
CTI  = 'CTI'



# Custom Colors
ColorConverter.cache = {}
ColorConverter.colors['UPMC_corporate_brown'] = (145/255, 120/255, 91/255)
ColorConverter.colors['beige'] = (96/255, 96/255, 86/255)
ColorConverter.colors['terminal_bg'] = (52/255, 56/255, 60/255)
ColorConverter.colors['UPMC_cool_gray'] = (97/255, 99/255, 101/255)
ColorConverter.colors['terminal_fg'] = (219/255, 219/255, 219/255)
ColorConverter.colors['white'] = (255/255, 255/255, 255/255)
ColorConverter.colors['black'] = (0/255, 0/255, 0/255)
ColorConverter.colors['turquoise'] = (26/255, 188/255, 156/255)
ColorConverter.colors['green_sea'] = (22/255, 160/255, 133/255)
ColorConverter.colors['emerald'] = (46/255, 204/255, 113/255)
ColorConverter.colors['nephritis'] = (39/255, 174/255, 96/255)
ColorConverter.colors['peter_river'] = (52/255, 152/255, 219/255)
ColorConverter.colors['belize_hole'] = (41/255, 128/255, 185/255)
ColorConverter.colors['amethyst'] = (155/255, 89/255, 182/255)
ColorConverter.colors['wisteria'] = (142/255, 68/255, 173/255)
ColorConverter.colors['wet_asphalt'] = (52/255, 73/255, 94/255)
ColorConverter.colors['midnight_blue'] = (44/255, 62/255, 80/255)
ColorConverter.colors['sun_flower'] = (241/255, 196/255, 15/255)
ColorConverter.colors['organge'] = (243/255, 156/255, 18/255)
ColorConverter.colors['carrot'] = (230/255, 126/255, 34/255)
ColorConverter.colors['pumpkin'] = (211/255, 84/255, 0/255)
ColorConverter.colors['alizarin'] = (231/255, 76/255, 60/255)
ColorConverter.colors['pomegranate'] = (192/255, 57/255, 43/255)
ColorConverter.colors['clouds'] = (236/255, 240/255, 241/255)
ColorConverter.colors['silver'] = (189/255, 195/255, 199/255)
ColorConverter.colors['concrete'] = (149/255, 165/255, 166/255)
ColorConverter.colors['asbestos'] = (127/255, 140/255, 141/255)
ColorConverter.colors['blue_1'] = (51/255, 102/255, 204/255)
ColorConverter.colors['red_1'] = (220/255, 57/255, 18/255)
ColorConverter.colors['yellow_1'] = (255/255, 153/255, 0/255)
ColorConverter.colors['green_1'] = (16/255, 150/255, 24/255)
ColorConverter.colors['purple_1'] = (153/255, 0/255, 153/255)
ColorConverter.colors['blue_2'] = (0/255, 153/255, 198/255)
ColorConverter.colors['pink_1'] = (221/255, 68/255, 119/255)
ColorConverter.colors['green_2'] = (102/255, 170/255, 0/255)
ColorConverter.colors['red_2'] = (184/255, 46/255, 46/255)
ColorConverter.colors['blue_3'] = (49/255, 99/255, 149/255)
ColorConverter.colors['purple_2'] = (153/255, 68/255, 153/255)
ColorConverter.colors['turquoise_1'] = (34/255, 170/255, 153/255)
ColorConverter.colors['green_olive_1'] = (170/255, 170/255, 17/255)
ColorConverter.colors['purple_3'] = (102/255, 51/255, 204/255)
ColorConverter.colors['orange_1'] = (230/255, 115/255, 0/255)
ColorConverter.colors['bordeau_1'] = (139/255, 7/255, 7/255)
ColorConverter.colors['purple_4'] = (101/255, 16/255, 103/255)
ColorConverter.colors['gray_green_1'] = (50/255, 146/255, 98/255)
ColorConverter.colors['blue_4'] = (85/255, 116/255, 166/255)
ColorConverter.colors['blue_5'] = (59/255, 62/255, 172/255)
ColorConverter.colors['brown_1'] = (183/255, 115/255, 34/255)
ColorConverter.colors['green_3'] = (22/255, 214/255, 32/255)
ColorConverter.colors['purple_5'] = (185/255, 19/255, 131/255)
ColorConverter.colors['pink_2'] = (244/255, 53/255, 158/255)
ColorConverter.colors['brown_2'] = (156/255, 89/255, 53/255)
ColorConverter.colors['green_olive_3'] = (169/255, 196/255, 19/255)
ColorConverter.colors['blue_6'] = (42/255, 119/255, 141/255)
ColorConverter.colors['green_4'] = (102/255, 141/255, 28/255)
ColorConverter.colors['green_olive_4'] = (190/255, 164/255, 19/255)
ColorConverter.colors['green_5'] = (12/255, 89/255, 34/255)
ColorConverter.colors['brown_3'] = (116/255, 52/255, 17/255)


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




def main_cactus(filename_edacc_csv, timeout, output="cactus.pdf", title=None):
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





def main_scatter_plotly():
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

    parser.add_argument('--ltl', nargs='?', dest='ltl',
                        help='list of wanted ltl: recurrence,persistence,reactivity,...')

    parser.add_argument('--bound', nargs='?', dest='bound',
                        help='list of wanted bound: 20,40,60,80,100,200,...,4000')

    parser.add_argument('-o', '--output', dest='output',
                        default=".",
                        help='output file of scatter plot')

    parser.add_argument('--title', dest='title',
                        default=None,
                        help='titleof scatter plot')

    parser.add_argument('--s1', dest='s1', metavar='<s1>',
                        help='first solver')
    parser.add_argument('--s2', dest='s2', metavar='<s2>',
                        help='second solver')

    parser.add_argument('--model', nargs='?', dest='model',
                        help='list of wanted model: adding,krebs,...')

    args = parser.parse_args()

    filename_edacc_csv = args.filename_edacc_csv
    df = pd.read_csv(filename_edacc_csv)
 

    if args.ltl != None:
        ltl = column_no_duplicate(df, LTL_KEY)
        ltl = [s.strip() for s in args.ltl.split(',')]
        print(ltl)
        df = keep_only_ltl(df, ltl)

    if args.bound != None:
        bound = [s.strip() for s in args.bound.split(',')]
        df = keep_only_bound(df, bound[0])

    if args.model != None:
        model = [s.strip() for s in args.model.split(',')]
        df = keep_only_model(df, model)

    if not args.s1 or not args.s2:
        print("Need solvers names option --s1 and --s2")
        exit()

    solver_one = args.s1
    solver_two = args.s2


    limit = None
    if args.timeout:
        limit = args.timeout

    scatter_plot(df, solver_one, solver_two, limit=limit,
                 output=args.output, title=args.title)



def scatter_plot(df, s1, s2, output=".", limit=None, title=None,
                 key=TIME_KEY):
    df = df.copy()

    solvers = [s1, s2]
    keep_only_solvers(df, solvers)

    ax = None

    # Guess limit
    if limit == None:
        limit = guess_limit(df)

    is_fail = (df[RESULT_KEY] != SAT) & (df[RESULT_KEY] != UNSAT)
    df.loc[is_fail, key] = limit * 1.1

    table_sat   = []
    table_unsat = []

    for instance in column_no_duplicate(df, INSTANCE_KEY):
        #print(instance)
        is_instance = df[INSTANCE_KEY] == instance
        is_unknown = False

        is_sat     = df[RESULT_KEY] == SAT
        is_unsat   = df[RESULT_KEY] == UNSAT
        is_unknown = df[RESULT_KEY] == UNKNOWN

        is_s1      = df[SOLVER_KEY] == s1
        is_s2      = df[SOLVER_KEY] == s2
        val_s1 = df[(is_instance) & (is_s1)][key].item()
        val_s2 = df[(is_instance) & (is_s2)][key].item()

        if len(df[(is_instance) & (is_unknown)].index) > 0 :
                print("...not plotted "+str(instance))
                pass
        elif len(df[(is_instance) & (is_sat)].index) > 0:
            table_sat.append([instance, val_s1, val_s2])
        elif len(df[(is_instance) & (is_unsat)].index) > 0:
            table_unsat.append([instance, val_s1, val_s2])

    df_unsat = pd.DataFrame(table_unsat, columns=[INSTANCE_KEY, s1, s2])
    df_sat = pd.DataFrame(table_sat, columns=[INSTANCE_KEY, s1, s2])

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df_unsat[s1])
        print(df_unsat[s2])
    fig_unsat = go.Scatter(
        x=df_unsat[s1],
        y=df_unsat[s2],
        text=df_unsat[EDACC_INSTANCE],
        mode='markers',
        name="UNSAT",
        marker=dict(symbol="cross", color="red", size=12),
    )

    fig_sat = go.Scatter(
        x=df_sat[s1],
        y=df_sat[s2],
        text=df_sat[EDACC_INSTANCE],
        mode='markers',
        name="SAT",
        marker=dict(symbol="x", color="blue", size=12),
    )
    

    linestyle = dict(color='rgba(0,0,0,0.6)', dash='dash')
    border_a = go.Scatter(
        x=[0, limit * 1.05],
        y=[0, limit * 1.05],
        showlegend = False,
        hoverinfo="none",
        mode='lines',
        line = linestyle
    )
    border_b = go.Scatter(
        x=[0, limit * 1.15],
        y=[limit * 1.05, limit * 1.05],
        showlegend = False,
        hoverinfo="none",
        mode='lines',
        line = linestyle
    )
    border_c = go.Scatter(
        x=[limit * 1.05, limit * 1.05],
        y=[limit * 1.15, 0],
        showlegend=False,
        hoverinfo="none",
        mode='lines',
        line = linestyle
    )

    data = [fig_unsat, fig_sat, border_a, border_b, border_c]

    layout =  go.Layout(
        autosize=False,
        height=800,
        width=800,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0.03)',#'rgba(77,210,255,0.03)',
        font=dict(
        #    family="Ariel",
            size=20,
        ),
        yaxis = dict(
            range=[0, limit*1.15],
            showgrid=True,
            zeroline=True,
            showline=True,
            gridcolor='rgba(0,0,0,0.1)',#'#bdbdbd',
            gridwidth=2,
            zerolinecolor='rgba(0,0,0,1)',
            zerolinewidth=3,
            linewidth=1,
 #           scaleanchor="x",
            scaleratio=1,
            title=s2,#r'$\Large\text{Kissat-MAB-H}_{\text{LP}-}LBD\leq3 \text{ (seconds)}$',#s2
        ),
        xaxis = dict(
            range=[0, limit*1.15],
            showgrid=True,
            zeroline=True,
            showline=True,
            gridcolor='rgba(0,0,0,0.1)',#'#bdbdbd',
            gridwidth=2,
            zerolinecolor='rgba(0,0,0,1)',#'#969696',
            zerolinewidth=3,
            linewidth=1,
#            scaleanchor="y",
            #titlefont=dict(size=30),
            scaleratio=10,
            title=s1,#r'$\Large\text{Kissat-MAB-Hordesat (seconds)}$',f
        )
    )
    fig = go.Figure(data=data, layout=layout)
    #fig.update_layout(plot_bgcolor='rgb(254,254,254)')
#    fig.show()
    fig.update_layout(
        autosize=False,
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=100,
            t=100,
            pad = 4),
        title=title,
        font=dict(size=18)
    )
    fig.write_image(output+"/"+s1+"_"+s2+".pdf")
    url = plot(fig, filename=output+"/"+s1+"_"+s2+".html")
    fig.write_image(output+".pdf")
    url = plot(fig, filename=output+".html")

#    print(url)


if __name__ == '__main__':
	pass
