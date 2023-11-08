#! /usr/env/python

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
