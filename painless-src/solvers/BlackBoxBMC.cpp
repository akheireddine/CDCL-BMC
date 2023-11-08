
#include <filesystem>

#include "BlackBoxBMC.h"
#include "../../mapleCOMSPS/mapleCOMSPS/core/Dimacs.h"
#include "../clauses/ClauseManager.h"
#include "../solvers/SolverFactory.h"
#include "../../mapleCOMSPS/mapleCOMSPS/utils/System.h"

#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/tl/parse.hh>
#include <spot/twa/bddprint.hh>
#include <spot/twa/formula2bdd.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/word.hh>
#include <spot/misc/common.hh>
#include <spot/misc/minato.hh>
#include <bddx.h>

using namespace MapleCOMSPS;

// Macros for minisat literal representation conversion
#define MINI_LIT(lit) lit > 0 ? mkLit(lit - 1, false) : mkLit((-lit) - 1, true)

#define INT_LIT(lit) sign(lit) ? -(var(lit) + 1) : (var(lit) + 1)

void makeMiniVec(std::vector<int> &cls, vec<Lit> &mcls)
{
    for (size_t i = 0; i < cls.size(); i++)
    {
        mcls.push(MINI_LIT(cls[i]));
    }
}

void copyMiniVecInt(std::vector<int> &cls, vec<int> &mcls)
{
    for (size_t i = 0; i < cls.size(); i++)
    {
        mcls.push(cls[i]);
    }
}

EnvBMC::EnvBMC(int id, bool inters = true,
               bool loopcond = true, bool edgecond = true,
               bool looprm = true) : SolverInterface(id, BMC),
                                     enable_intersection(inters),
                                     enable_loop_cond(loopcond),
                                     enable_edge_cond(edgecond),
                                     enable_loop_rm(looprm)
{
    seq = true;
    max_vars = 0;
    maxVarsPure = -1;
    minVarsPure = -1;
    nb_model_variable = 0;
    nb_property_variable = 0;
    nb_junction_variable = 0;
    k_and_fictif_state = 0;
    K = 0;
    aut = NULL;
    word_aut = NULL;
    numb_clauses.resize(4), size_clauses.resize(4), used_clauses_conflicts.resize(4);
    numb_call = 0, numb_clauses_falsified = 0, lbd_clauses = 0, numb_asserting_clauses = 0;
    time_init = 0, added_clauses = 0, time_blackbox = 0, time_sync_product = 0, remain_clauses = 0;
    prop_deci_rate_inf = Parameters::getIntParam("prop-rate-inf", 10);
    prop_deci_rate_sup = Parameters::getIntParam("prop-rate-sup", 100);
    depth_rate = Parameters::getIntParam("depth-rate", 100);
    lbd_limit = Parameters::getIntParam("lbd-limit-bmc", -1);
    solver = NULL;
    init_box_done = false;
    stop = false;
    reducer = NULL;

    new_trail = NULL;
    lockTrail = 1;
    oneTrailAtTime = Parameters::getBoolParam("oneTrail");
    noshare = Parameters::getBoolParam("no-share");
}

bool EnvBMC::ExportClause(vec<int> &cls, CONSTYPE const_type)
{
    //**************** USE REDUCER ******************//
    if (reducer && Parameters::getBoolParam("reducer"))
    {
        vec<Lit> strengthenedCls;
        SatResult res = reducer->solve_bmc(cls, strengthenedCls);
        if (res == UNSAT)
        {
            printf("c REDUNCER says UNSAT\n");
            return false;
        }
        ClauseExchange *ncls = ClauseManager::allocClause(strengthenedCls.size());
        std::unordered_set<int> levels;
        int undef = 0;
        int lvl;
        int sz = strengthenedCls.size();
        for (int i = 0; i < strengthenedCls.size(); i++)
        {
            int a = abs(INT_LIT(strengthenedCls[i]));
            ncls->lits[i] = INT_LIT(strengthenedCls[i]);
            lvl = curr_levels[abs((INT_LIT(strengthenedCls[i])))];
            if (lvl > -1)
                levels.insert(lvl);
            else
                undef++;
        }
        int lbd = levels.size() + undef;
        if (lbd_limit >= 0 && lbd > lbd_limit && const_type != EMPTY_INTERSECT)
        {
            ClauseManager::releaseClause(ncls);
            return true;
        }
        lbd_clauses += lbd;
        ncls->lbd = lbd;
        ncls->from = -1;
        ncls->size = sz;
        clausesToExport.addClause(ncls);
        added_clauses++;
    }
    else
    {
        ClauseExchange *ncls = ClauseManager::allocClause(cls.size());
        std::unordered_set<int> levels;
        // vec<Lit> cls_lit;
        int undef = 0;
        int lvl;

        for (int i = 0; i < cls.size(); i++)
        {
            // cls_lit.push(MINI_LIT(cls[i]));
            ncls->lits[i] = cls[i];
            lvl = curr_levels[abs(cls[i])];
            if (lvl > -1)
                levels.insert(lvl);
            else
                undef++;
        }
        int lbd = levels.size() + undef;
        if (noshare)
        {
            ClauseManager::releaseClause(ncls);
            return false;
        }
        if (lbd_limit >= 0 && lbd > lbd_limit && const_type != EMPTY_INTERSECT)
        {
            ClauseManager::releaseClause(ncls);
            return false;
        }
        if (Parameters::getBoolParam("fic-lbd"))
            lbd = std::min(cls.size(), 6); // cls.size();

        lbd_clauses += lbd;
        ncls->lbd = lbd;
        ncls->from = -1;
        ncls->size = cls.size();
        clausesToExport.addClause(ncls);
        added_clauses++;
    }
    //****************************************************

    return true;
}

bool EnvBMC::ImportClause()
{
    ClauseExchange *new_trail_ = NULL;

    if (trailToImport.getClause(&new_trail_) == false)
        return false;

    assert(new_trail_->from == -1);

    curr_trail.clear();
    // curr_levels.clear();
    reinitialize_trail_map();
    // Set their assignements if it exists
    for (int i = 0; i < new_trail_->lbd; i++)
    {
        int variable = new_trail_->lits[i];
        int abs_variable = abs(variable);
        assert(info_variables[abs_variable].in_property);
        trail_prop_spot[get_step(abs_variable)][get_name_variable(abs_variable)] = variable >= 0 ? l_True : l_False;
        curr_trail.emplace_back(variable);
        curr_levels[abs_variable] = new_trail_->lits[new_trail_->lbd + i];
    }
    return true;
}

bool EnvBMC::tryLockTrail(ClauseExchange *trail)
{
    int expected = 0;
    if (new_trail == NULL)
    {
        if (!lockTrail.compare_exchange_weak(expected, 2))
            return false;
        else
        {
            assert(lockTrail == 2);
            new_trail = trail;
            lockTrail = 0;
            return true;
        }
    }
    return false;
}

bool EnvBMC::ImportTrailAtATime()
{
    if (new_trail == NULL)
        return false;
    int expected = 0;
    while (!lockTrail.compare_exchange_weak(expected, 1))
        expected = 0;

    assert(lockTrail == 1);
    assert(new_trail->from == -1);

    curr_trail.clear();
    curr_levels.clear();
    reinitialize_trail_map();
    // Set their assignements if it exists
    for (int i = 0; i < new_trail->lbd; i++)
    {
        int variable = new_trail->lits[i];
        int abs_variable = abs(variable);
        assert(info_variables[abs_variable].in_property);
        trail_prop_spot[get_step(abs_variable)][get_name_variable(abs_variable)] = variable >= 0 ? l_True : l_False;
        curr_trail.emplace_back(variable);
        curr_levels[abs_variable] = new_trail->lits[new_trail->lbd + i];
    }
    return true;
}

bool EnvBMC::prepareClause(vec<int> &cls, vec<Lit> &out_cls, CONSTYPE const_type)
{
    out_cls.clear();
    std::unordered_set<int> levels;
    int undef = 0;
    int lvl;
    for (int i = 0; i < cls.size(); i++)
    {
        Lit l = MINI_LIT(cls[i]);
        out_cls.push(l),
            lvl = solver->getLevel(var(l));
        if (lvl > -1)
            levels.insert(lvl);
        else
            undef++;
    }
    int lbd = levels.size() + undef;
    if (lbd_limit >= 0 && lbd > lbd_limit && const_type != EMPTY_INTERSECT)
        return false;
    return true;
}

void EnvBMC::initialize_info_variables()
{
    maxVarsPure = -1;
    info_variables.resize(max_vars);
    k_and_fictif_state = pure_vars.size();
    K = k_and_fictif_state - 1;

    name_variable_id.resize(k_and_fictif_state);
    for (unsigned i = 0; i < pure_vars.size(); i++)
    {
        for (auto var_name : pure_vars[i])
        {
            int var = std::get<0>(var_name);
            std::string name = std::get<1>(var_name);
            nb_model_variable++;
            info_variables[var].num_step = i;
            info_variables[var].name = name;
            info_variables[var].in_property = false;
            name_variable_id[i].insert(std::make_pair(name, var));

            maxVarsPure = maxVarsPure < var ? var : maxVarsPure;
            minVarsPure = ((minVarsPure > var) || (minVarsPure == -1)) ? var : minVarsPure;
        }
    }
    for (int i = maxVarsPure + 1; i < max_vars; i++)
    {
        info_variables[i].num_step = -1;
        info_variables[i].num_partition = -1;
        info_variables[i].in_property = false;
    }
    nb_junction_variable += (max_vars - (maxVarsPure + 1));
}

void EnvBMC::init()
{
    loadFormula(Parameters::getFilename());
    initialize_info_variables();

    if (Parameters::getBoolParam("reducer"))
        reducer = (Reducer *)SolverFactory::createReducerSolver(SolverFactory::createMapleCOMSPSSolver(NULL));

    // Parse the LTL formula
    spot::parsed_formula pf = spot::parse_infix_psl(ltlspec);

    if (pf.format_errors(std::cerr))
    {
        std::cerr << "c fomula error" << std::endl;
        exit(1);
    }

    spot::translator trans;
    trans.set_type(spot::postprocessor::Buchi);
    aut = trans.run(pf.f);

    initialize_literals_map();
    init_box_done = true;
    std::cout << "c Number of formula states: " << aut->num_states() << "\n";
    std::cout << "c Number synchr states: " << aut->num_states() * K << "\n";
    std::cout << "c Number of property variables: " << nb_property_variable << "\n";
}

// Initialize literal's map
void EnvBMC::initialize_literals_map()
{
    std::string name;
    trail_prop_spot.clear();
    curr_levels.clear();
    trail_prop_spot.resize(K); // fictif state is not included
    // Set the literal at each step on the map
    for (unsigned step = 0; step < K; step++)
        for (spot::formula ap : aut->ap())
        {
            name = ap.get_child_of({}).ap_name();
            trail_prop_spot[step].insert(std::make_pair(name, l_Undef));
            int var = name_variable_id[step][name];
            curr_levels.insert(std::make_pair(var, -1));
            info_variables[var].in_property = true;
            nb_property_variable++;
            nb_model_variable--;
        }
    for (int i = 0; i < rename_loop.size(); i++)
    {
        int var = rename_loop[i];
        std::string name = "loop" + to_string(i);
        trail_prop_spot[i].insert(std::make_pair(name, l_Undef));
        curr_levels.insert(std::make_pair(var, -1));
        info_variables[var].in_property = true;
        info_variables[var].num_step = i;
        info_variables[var].name = name;
        nb_property_variable++;
        nb_model_variable--;
    }
}

void EnvBMC::reinitialize_trail_map()
{
    std::string name;
    // Set the literal at each step on the map
    for (unsigned step = 0; step < K; step++)
        for (auto [var_name, cnf_assign] : trail_prop_spot[step])
        {
            trail_prop_spot[step][var_name] = l_Undef;
            curr_levels[name_variable_id[step][var_name]] = -1;
        }
}

SatResult EnvBMC::solve(const vector<int> &cube)
{
    time_init = cpuTime();
    init();
    time_init = cpuTime() - time_init;
    printf("c Blackbox start solving...\n");
    if (!seq)
        explore_automata_each_loop_parallel();
    return NOTHING;
}

void EnvBMC::explore_automata_each_loop(vec<vec<Lit>> &out_refined)
{
    // double parsed_time = cpuTime();
    numb_call++;
    word_aut = NULL;
    int K_real = K;
    std::vector<bool> is_empty_inter(K, false);
    int nb_empty_inter = 0;
    std::vector<bdd> loop_constraints_aut(K, bddtrue);
    bdd loop_formula_all_aut = bddtrue;
    std::vector<std::vector<bdd>> edges_constraints_aut(K - 1);
    std::vector<bdd> edges_formula_all_aut(K - 1, bddtrue);
    vec<int> cls_result;
    double tmp_time;
    for (int i = 0; i < K; i++)
    {
        if (rename_loop[i] == 0)
        {
            continue;
        }
        if (solver->value(MINI_LIT(rename_loop[i])) == l_False)
            continue;

        create_word_automata(i);
        tmp_time = cpuTime();
        std::vector<std::vector<bdd>> constraints = loops_edges_constraints(i);
        time_sync_product += cpuTime() - tmp_time;
        if (constraints.empty() && enable_intersection)
        {
            is_empty_inter[i] = true;
            nb_empty_inter++;
        }
        if (!constraints.empty())
        {
            loop_constraints_aut[i] = constraints[0][0];
            if (loop_formula_all_aut == bddtrue) // init loop bdd
                loop_formula_all_aut = loop_constraints_aut[i];
            else if (loop_formula_all_aut != loop_constraints_aut[i]) // not equal at least once
                loop_formula_all_aut = bddfalse;
            if (i < K - 1)
            {
                edges_constraints_aut[i] = constraints[1];
                for (size_t j = 0; j < K - 1; j++)
                {
                    if (edges_formula_all_aut[j] == bddtrue) // init j-th edge bdd
                        edges_formula_all_aut[j] = edges_constraints_aut[i][j];
                    else if (edges_formula_all_aut[j] != edges_constraints_aut[i][j]) // not equal at least once on j-th edge
                        edges_formula_all_aut[j] = bddfalse;
                }
            }
        }
    }
    // last automaton, can create clauses if bdds are the same
    {
        if (nb_empty_inter == K_real && enable_intersection)
        {
            // std::cout << "c ALL INTERSECTION VIDE " << solver->curr_trail.size() << "\n";
            cls_result.clear();
            add_trail_to_cls(cls_result);
            out_refined.push();
            int id = out_refined.size() - 1;
            if (!prepareClause(cls_result, out_refined[id], EMPTY_INTERSECT))
                out_refined.pop();
            else
            {
                // //////////////////////// STATS ////////////////////
                numb_clauses[EMPTY_INTERSECT]++;
                size_clauses[EMPTY_INTERSECT] += out_refined[id].size();
                ///////////////////////////////////////////////////
            }
        }
        if (enable_loop_cond && loop_formula_all_aut != bddtrue && loop_formula_all_aut != bddfalse)
            createLoopEdgeClauses(loop_formula_all_aut, out_refined, K - 1, -1);
        if (enable_edge_cond)
        {
            for (size_t j = 0; j < K - 1; j++)
            {
                if (edges_formula_all_aut[j] != bddtrue)
                {
                    if (edges_formula_all_aut[j] != bddfalse)
                        createLoopEdgeClauses(edges_formula_all_aut[j], out_refined, j, -1);
                }
            }
        }
    }
    // try again if some bdd's transitions are not equivalent
    for (size_t i = 0; i < K; i++)
    {
        if (rename_loop[i] == 0)
            continue;
        if (solver->value(MINI_LIT(rename_loop[i])) == l_False)
            continue;

        cls_result.clear();
        if (nb_empty_inter < K_real && is_empty_inter[i] && enable_intersection)
        {
            create_trail_cls(cls_result, i);
            out_refined.push();
            int id = out_refined.size() - 1;
            if (!prepareClause(cls_result, out_refined[id], EMPTY_INTERSECT))
                out_refined.pop();
            else
            {
                // //////////////////////// STATS ////////////////////
                numb_clauses[EMPTY_INTERSECT]++;
                size_clauses[EMPTY_INTERSECT] += out_refined[id].size();
                ///////////////////////////////////////////////////
            }
        }
        else
        {
            bdd loop_constraint = loop_constraints_aut[i];
            std::vector<bdd> &edges_constraints = edges_constraints_aut[i];
            if (enable_loop_cond && loop_constraint != bddtrue && loop_constraint != bddfalse)
                if (loop_formula_all_aut == bddfalse) // bdd's are not equivalent between automata
                    createLoopEdgeClauses(loop_constraint, out_refined, K - 1, i);
            if (enable_edge_cond && i < K - 1)
            {
                if (edges_constraints.empty())
                    continue;
                for (unsigned j = 0; j < K - 1; j++)
                {
                    assert(edges_constraints[j] != bddfalse);
                    if (edges_constraints[j] != bddtrue)
                        if (edges_formula_all_aut[j] == bddfalse) // bdd's are not equivalent between automata
                            createLoopEdgeClauses(edges_constraints[j], out_refined, j, i);
                }
            }
        }
    }
    // time_blackbox += cpuTime() - parsed_time;
}

void EnvBMC::explore_automata_each_loop_parallel()
{
    int K_real = K;
    vec<vec<Lit>> out_refined;
    while (true)
    {
        lockTrail = 0;
        if (oneTrailAtTime)
            while (!ImportTrailAtATime())
                ;
        else
            while (!ImportClause())
                ;

        numb_call++;
        word_aut = NULL;
        std::vector<bool> is_empty_inter(K, false);
        int nb_empty_inter = 0;
        std::vector<bdd> loop_constraints_aut(K, bddtrue);
        bdd loop_formula_all_aut = bddtrue;
        std::vector<std::vector<bdd>> edges_constraints_aut(K - 1);
        std::vector<bdd> edges_formula_all_aut(K - 1, bddtrue);
        vec<int> cls_result;
        double tmp_time;
        bool compact = true;
        for (int i = 0; i < K; i++)
        {
            if (rename_loop[i] == 0)
            {
                compact = false;
                continue;
            }

            create_word_automata(i);
            std::vector<std::vector<bdd>> constraints = loops_edges_constraints(i);
            if (constraints.empty() && enable_intersection)
            {
                is_empty_inter[i] = true;
                nb_empty_inter++;
            }
            if (!constraints.empty())
            {
                loop_constraints_aut[i] = constraints[0][0];
                if (loop_formula_all_aut == bddtrue) // init loop bdd
                    loop_formula_all_aut = loop_constraints_aut[i];
                else if (loop_formula_all_aut != loop_constraints_aut[i]) // not equal at least once
                    loop_formula_all_aut = bddfalse;
                if (i < K - 1)
                {
                    edges_constraints_aut[i] = constraints[1];
                    if (compact)
                    {
                        for (size_t j = 0; j < K - 1; j++)
                        {
                            if (edges_formula_all_aut[j] == bddtrue) // init j-th edge bdd
                                edges_formula_all_aut[j] = edges_constraints_aut[i][j];
                            else if (edges_formula_all_aut[j] != edges_constraints_aut[i][j]) // not equal at least once on j-th edge
                                edges_formula_all_aut[j] = bddfalse;
                        }
                    }
                }
            }
        }
        // last automaton, can create clauses if bdds are the same

        if (nb_empty_inter == K_real && enable_intersection)
        {
            // std::cout << "c ALL INTERSECTION VIDE " << curr_trail.size() << "\n";
            cls_result.clear();
            add_trail_to_cls(cls_result);
            if (ExportClause(cls_result, EMPTY_INTERSECT))
            {
                // //////////////////////// STATS ////////////////////
                numb_clauses[EMPTY_INTERSECT]++;
                // size_clauses[EMPTY_INTERSECT] += out_refined[id].size();
                ///////////////////////////////////////////////////
            }
        }
        if (compact)
        {
            if (enable_loop_cond && loop_formula_all_aut != bddtrue && loop_formula_all_aut != bddfalse)
                createLoopEdgeClauses(loop_formula_all_aut, out_refined, K - 1, -1);
            if (enable_edge_cond)
            {
                for (size_t j = 0; j < K - 1; j++)
                    if (edges_formula_all_aut[j] != bddtrue)
                        if (edges_formula_all_aut[j] != bddfalse)
                        {
                            createLoopEdgeClauses(edges_formula_all_aut[j], out_refined, j, -1);
                        }
            }
        }
        // try again if some bdd's transitions are not equivalent
        for (size_t i = 0; i < K; i++)
        {
            if (rename_loop[i] == 0)
                continue;

            cls_result.clear();
            if (nb_empty_inter < K_real && is_empty_inter[i] && enable_intersection)
            {
                create_trail_cls(cls_result, i);
                if (ExportClause(cls_result, EMPTY_INTERSECT))
                {
                    // //////////////////////// STATS ////////////////////
                    numb_clauses[EMPTY_INTERSECT]++;
                    // size_clauses[EMPTY_INTERSECT] += out_refined[id].size();
                    ///////////////////////////////////////////////////
                }
            }
            else
            {
                bdd loop_constraint = loop_constraints_aut[i];
                std::vector<bdd> &edges_constraints = edges_constraints_aut[i];
                if (enable_loop_cond && loop_constraint != bddtrue && loop_constraint != bddfalse)
                    if (loop_formula_all_aut == bddfalse) // bdd's are not equivalent between automata
                        createLoopEdgeClauses(loop_constraint, out_refined, K - 1, i);
                if (enable_edge_cond && i < K - 1)
                {
                    if (edges_constraints.empty())
                        continue;
                    for (unsigned j = 0; j < K - 1; j++)
                    {
                        assert(edges_constraints[j] != bddfalse);
                        if (edges_constraints[j] != bddtrue)
                            if (!compact && edges_formula_all_aut[j] == bddfalse) // bdd's are not equivalent between automata
                                createLoopEdgeClauses(edges_constraints[j], out_refined, j, i);
                    }
                }
            }
        }
        ClauseManager::releaseClause(new_trail);
        new_trail = NULL;
    }
}

// Create automata of the new assignement
void EnvBMC::create_word_automata()
{
    spot::bdd_dict_ptr bdd_dict = aut->get_dict();

    word_aut = spot::make_twa_graph(bdd_dict);
    word_aut->copy_ap_of(aut);
    word_aut->new_states(K);

    for (unsigned step = 0; step < K; step++)
    {
        bdd cond = bddtrue;
        for (auto [var_name, cnf_assign] : trail_prop_spot[step])
        {
            if (cnf_assign == l_Undef || var_name == "loop" + to_string(step))
                continue;

            int var = bdd_dict->has_registered_proposition(spot::formula::ap(var_name), aut);
            if (var < 0)
                continue;

            if (cnf_assign == l_True)
                cond &= bdd_ithvar(var);
            else
                cond &= bdd_nithvar(var);
        }
        if (step != K - 1)
            word_aut->new_edge(step, step + 1, cond);
        else
            for (unsigned i = 0; i < K; i++)
            {
                if (rename_loop[i] > 0)
                    word_aut->new_edge(step, i, cond);
            }
    }
    // printf("\n\nMOT \n");
    // spot::print_hoa(std::cout, word_aut, "k");
}

// Create automata of the new assignement
void EnvBMC::create_word_automata(int loop)
{
    if (word_aut != NULL)
    {
        word_aut->edge_storage(num_edge_loop).dst = loop;
        return;
    }
    spot::bdd_dict_ptr bdd_dict = aut->get_dict();

    word_aut = spot::make_twa_graph(bdd_dict);
    word_aut->copy_ap_of(aut);
    word_aut->new_states(K);

    for (unsigned step = 0; step < K; step++)
    {
        bdd cond = bddtrue;
        for (auto [var_name, cnf_assign] : trail_prop_spot[step])
        {
            if (cnf_assign == l_Undef || var_name == "loop" + to_string(step))
                continue;
            int var = bdd_dict->has_registered_proposition(spot::formula::ap(var_name), aut);
            if (var < 0)
                continue;

            if (cnf_assign == l_True)
                cond &= bdd_ithvar(var);
            else
                cond &= bdd_nithvar(var);
        }
        if (step != K - 1)
            word_aut->new_edge(step, step + 1, cond);
        else
        {
            num_edge_loop = word_aut->new_edge(step, loop, cond);
            // std::cout << "\nc CONDITION "
            //   << " : " << spot::bdd_format_formula(aut->get_dict(), cond) << std::endl;
        }
    }
    // printf("\n\nMOT \n");
    // spot::print_hoa(std::cout, word_aut, "k");
}

std::vector<std::vector<bdd>> EnvBMC::loops_edges_constraints()
{

    spot::twa_graph_ptr p = spot::product(word_aut, aut);
    std::vector<std::vector<bdd>> constraints;
    spot::scc_info si(p);

    if (!si.is_useful_state(p->get_init_state_number()))
        return constraints;
    // std::vector<bool> usefull_loops(K, false);
    std::vector<bdd> loop_constraints(K, bddfalse);
    std::vector<bdd> edges_constraints(K - 1, bddfalse);
    constraints.resize(2);
    auto *ps = p->get_named_prop<spot::product_states>("product-states");
    int nb_scc = si.scc_count();
    unsigned step;

    for (int scc = 0; scc < nb_scc; scc++)
    {
        if (!si.is_useful_scc(scc))
            continue;
        for (unsigned s : si.states_of(scc))
        {
            step = (*ps)[s].first;
            if (step == K - 1)
            {
                for (auto &e : p->out(s))
                    if (si.is_accepting_scc(si.scc_of(e.dst)) && rename_loop[(*ps)[e.dst].first] > 0)
                    {
                        // usefull_loops[(*ps)[e.dst].first] = true;
                        loop_constraints[(*ps)[e.dst].first] |= e.cond;
                    }
            }
            else
            {
                for (auto &e : p->out(s))
                    if (si.is_useful_state(e.dst))
                        edges_constraints[step] |= e.cond;
            }
        }
    }
    constraints[0] = loop_constraints;
    constraints[1] = edges_constraints;
    return constraints;
}

std::vector<std::vector<bdd>> EnvBMC::loops_edges_constraints(int loop)
{
    spot::twa_graph_ptr p = spot::product(word_aut, aut);
    std::vector<std::vector<bdd>> constraints;
    spot::scc_info si(p);

    if (!si.is_useful_state(p->get_init_state_number()))
        return constraints;
    bdd loop_constraint = bddfalse;
    std::vector<bdd> edges_constraints(K - 1, bddfalse);
    constraints.resize(2, std::vector<bdd>());
    auto *ps = p->get_named_prop<spot::product_states>("product-states");
    int nb_scc = si.scc_count();
    unsigned step;

    for (int scc = 0; scc < nb_scc; scc++)
    {
        if (!si.is_useful_scc(scc))
            continue;
        for (unsigned s : si.states_of(scc))
        {
            step = (*ps)[s].first;
            if (step == K - 1)
            {
                for (auto &e : p->out(s))
                    if (si.is_accepting_scc(si.scc_of(e.dst))) // && ((*ps)[e.dst].first == loop))
                    {
                        // usefull_loops[(*ps)[e.dst].first] = true;
                        loop_constraint |= e.cond;
                    }
            }
            else
            {
                // if(step != loop) continue;
                for (auto &e : p->out(s))
                    if (si.is_useful_state(e.dst))
                        edges_constraints[step] |= e.cond;
            }
        }
    }

    constraints[0].emplace_back(loop_constraint);
    constraints[1] = edges_constraints;
    return constraints;
}

void EnvBMC::add_trail_to_cls(vec<int> &cls)
{
    if (seq)
        for (Lit l : solver->curr_trail)
            cls.push(INT_LIT(~l));
    else
        for (int v : curr_trail)
            cls.push(-v);
}

// ADD CLAUSE  (alpha -> -Li)  ---->    -alpha v -Li ----->  -alpha v -rename
void EnvBMC::constraints_to_remove_loop(vec<int> &cls_result, int loop)
{
    cls_result.clear();
    cls_result.push(-rename_loop[loop]);
    add_trail_to_cls(cls_result);
    //////////////////////// STATS ////////////////////
    numb_clauses[LOOP_EXCLUDED]++;
    size_clauses[LOOP_EXCLUDED] += cls_result.size();
    ///////////////////////////////////////////////////
}

// ADD CLAUSE  (L1 -> -L2 ^ -L3 ^ ...) ^ (L2 -> -L1 ^ -L3 ^ ...) ^ (L3 -> -L1 ^ -L2 ^ ...) ...
// ADD CLAUSE  (-L1 v (-L2 ^ -L3 ^ ...)) ^ (-L2 v (-L1 ^ -L3 ^ ...)) ^ (-L3 v (-L1 ^ -L2 ^ ...)) ...
// ADD CLAUSE  (-L1 v -L2) ^ (-L1 v -L3) ^ ...  ^ (-L2 v -L3) ^ ...
void EnvBMC::addLoopClauses()
{
    //     vec<Lit> cls;
    //     bool ok = true;
    //     for (int i = 0; i < rename_loop.size(); i++)
    //     {
    //         if (rename_loop[i] > 0)
    //         {
    //             for (int j = i + 1; j < rename_loop.size(); j++)
    //             {
    //                 if (rename_loop[j] > 0 && j != i)
    //                 {
    //                     cls.clear();
    //                     cls.push(MINI_LIT(-rename_loop[i]));
    //                     cls.push(MINI_LIT(-rename_loop[j]));
    //                     solver->addClause_(cls);
    //                 }
    //             }
    //         }
    //     }
}

// ADD CLAUSE  (alpha ^ Li -> constraints[i]) ---->  -alpha v -Li v constraints[i] ---> -alpha v -renamei v  constraints[i]
void EnvBMC::constraints_on_transition_loop(vec<int> &cls_result, std::vector<int> &constraints, int loop)
{
    cls_result.clear();
    // Lit rename = MINI_LIT(rename_loop[loop]);
    // add trail (alpha)
    add_trail_to_cls(cls_result);
    // add loop formula (rename variable)
    cls_result.push(-rename_loop[loop]);
    // for (int i = 0; i < rename_loop.size(); i++)
    //     if ((rename_loop[i] > 0) && (i != loop))
    //         cls_result[id].push(MINI_LIT(rename_loop[i]));
    // add bdd formula (constraint)
    copyMiniVecInt(constraints, cls_result);
    //////////////////////// STATS ////////////////////
    numb_clauses[LOOP_CONST]++;
    size_clauses[LOOP_CONST] += cls_result.size();
    //////////////////////////////////////////////////
}

// ADD CLAUSE (alpha -> constraints[i])  -----> -alpha v constraints[i]
void EnvBMC::constraints_on_edge(vec<int> &cls_result, std::vector<int> &constraints)
{
    cls_result.clear();
    // add trail (alpha)
    add_trail_to_cls(cls_result);
    // add bdd formula (constraint)
    copyMiniVecInt(constraints, cls_result);
    //////////////////////// STATS ////////////////////
    numb_clauses[EDGE_CONST]++;
    size_clauses[EDGE_CONST] += cls_result.size();
    //////////////////////////////////////////////////
}

// ADD CLAUSE (Li -> -alpha)   ----------->   (-Li v -alpha) ---------> -rename v -alpha --------> -rename v -alpha[0] v -alpha[1]...
void EnvBMC::create_trail_cls(vec<int> &cls_result, int loop)
{
    cls_result.clear();
    // Lit rename = MINI_LIT(rename_loop[loop]);
    // if (solver->value(rename) == l_True)
    // return;
    cls_result.push(-rename_loop[loop]);
    add_trail_to_cls(cls_result);
    //////////////////////// STATS ////////////////////
    numb_clauses[EMPTY_INTERSECT]++;
    size_clauses[EMPTY_INTERSECT] += cls_result.size();
    //////////////////////////////////////////////////
}

// ADD CLAUSE  (alpha ^ Li -> constraints[i]) ---->  -alpha v -Li v constraints[i] ---> -alpha v -renamei v  constraints[i]
void EnvBMC::createLoopEdgeClauses(bdd formula, vec<vec<Lit>> &out_cls, int curr_step, int loop)
{
    if (formula == bddtrue)
        return;

    spot::minato_isop isop(!formula);
    vec<int> cls_result;
    bdd cube, b;
    int id;

    // std::cout << "\nc EDGE "
    //   << " : " << spot::bdd_format_formula(aut->get_dict(), formula) << std::endl;

    while ((cube = isop.next()) != bddfalse)
    {
        bool already_sat = false, redundant = true;
        b = cube;
        cls_result.clear();

        if (b == bddfalse)
            return;
        // add bdd formula (constraint)
        while (b != bddtrue)
        {
            int var = bdd_var(b);
            const spot::bdd_dict::bdd_info &i = aut->get_dict()->bdd_map[var];
            int res_int = name_variable_id[curr_step][(i.f).get_child_of({}).ap_name()];
            bdd high = bdd_high(b);
            if (high == bddfalse)
                b = bdd_low(b);
            else
            {
                res_int = -res_int;
                // If bdd_low is not false, then b was not a conjunction.
                assert(bdd_low(b) == bddfalse);
                b = high;
            }
            assert(b != bddfalse);

            lbool value_sign = seq ? solver->value(MINI_LIT(res_int)) : trail_prop_spot[get_step(abs(res_int))][get_name_variable(abs(res_int))];
            if (value_sign == l_True)
            {
                already_sat = true;
                if (seq)
                    break;
            }
            // else if (value_sign == l_Undef)
            // redundant = false;
            cls_result.push(res_int);
        }
        // Unuseful information
        if (seq || Parameters::getBoolParam("sat-redund"))
            if (already_sat) // || redundant)
                continue;

        add_trail_to_cls(cls_result);
        if (loop != -1)
            cls_result.push(-rename_loop[loop]);

        if (seq)
        {
            out_cls.push();
            id = out_cls.size() - 1;
            if (!prepareClause(cls_result, out_cls[id], LOOP_CONST))
                out_cls.pop();
            else
            {
                //////////////////////// STATS ////////////////////
                if (curr_step == K - 1)
                {
                    numb_clauses[LOOP_CONST]++;
                    size_clauses[LOOP_CONST] += out_cls[id].size();
                }
                else
                {
                    numb_clauses[EDGE_CONST]++;
                    size_clauses[EDGE_CONST] += out_cls[id].size();
                }
            }
        }
        else
        {
            if (ExportClause(cls_result, LOOP_CONST))
            {
                if (curr_step == K - 1)
                {
                    numb_clauses[LOOP_CONST]++;
                    size_clauses[LOOP_CONST] += cls_result.size();
                }
                else
                {
                    numb_clauses[EDGE_CONST]++;
                    size_clauses[EDGE_CONST] += cls_result.size();
                }
            }
        }
        //////////////////////////////////////////////////
    }
}

bool EnvBMC::loadFormula(const char *filename)
{
    gzFile in = gzopen(filename, "rb");
    // parse_DIMACS(in, *solver);
    parse_DIMACS_BMC(*(this), in);

    gzclose(in);

    return true;
}

// Get the number of variables of the formula
int EnvBMC::getVariablesCount()
{
    return -1;
}

// Get a variable suitable for search splitting
int EnvBMC::getDivisionVariable()
{
    return -1;
}

// Set initial phase for a given variable
void EnvBMC::setPhase(const int var, const bool phase)
{
}

// Bump activity for a given variable
void EnvBMC::bumpVariableActivity(const int var, const int times)
{
}

// Interrupt the SAT solving, so it can be started again with new assumptions
void EnvBMC::setSolverInterrupt()
{
    stop = true;
}

void EnvBMC::unsetSolverInterrupt()
{
    stop = false;
}

// Diversify the solver
void EnvBMC::diversify(int id)
{
}

void EnvBMC::addClause(ClauseExchange *clause)
{
    clausesToAdd.addClause(clause);

    setSolverInterrupt();
}

void EnvBMC::addLearnedClause(ClauseExchange *clause)
{
    trailToImport.addClause(clause);
}

void EnvBMC::addClauses(const vector<ClauseExchange *> &clauses)
{
    clausesToAdd.addClauses(clauses);

    setSolverInterrupt();
}

void EnvBMC::addInitialClauses(const vector<ClauseExchange *> &clauses)
{
}

void EnvBMC::addLearnedClauses(const vector<ClauseExchange *> &clauses)
{
    for (size_t i = 0; i < clauses.size(); i++)
    {
        addLearnedClause(clauses[i]);
    }
}

void EnvBMC::getLearnedClauses(vector<ClauseExchange *> &clauses)
{
    clausesToExport.getClauses(clauses);
}

void EnvBMC::increaseClauseProduction()
{
}

void EnvBMC::decreaseClauseProduction()
{
}

SolvingStatistics
EnvBMC::getStatistics()
{
    SolvingStatistics stats;
    return stats;
}

vector<int>
EnvBMC::getModel()
{
    return {};
}

vector<int>
EnvBMC::getFinalAnalysis()
{
    return {};
}

vector<int>
EnvBMC::getSatAssumptions()
{
    return {};
}

void EnvBMC::printStats(double total_time)
{
    double mem_used = memUsedPeak();
    // printf("c restarts              : %" PRIu64 "\n", solver.starts);
    // printf("c conflicts             : %-12" PRIu64 "   (%.0f /sec)\n", solver.conflicts, solver.conflicts / cpu_time);
    // printf("c decisions             : %-12" PRIu64 "   (%4.2f %% random) (%.0f /sec)  (%.0f %% decision)\n", solver.decisions,
    // (float)solver.rnd_decisions * 100 / (float)solver.decisions,
    // solver.decisions / cpu_time,
    // (float)solver.decisions * 100 / (float)solver.nVars());
    // printf("c propagations          : %-12" PRIu64 "   (%.0f /sec)\n", solver.propagations, solver.propagations / cpu_time);
    // printf("c conflict literals     : %-12" PRIu64 "   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals) * 100 / (double)solver.max_literals);

    uint64_t totalGeneratedCls = getTotalClauses();
    printf("c calls                 : %-12" PRIu64 "\n", numb_call);
    printf("c total clauses         : %-12" PRIu64 "\n", totalGeneratedCls);
    printf("c total unit clauses    : %-12" PRIu64 "\n", unit_clauses);
    printf("c numb cls EMPTY_INTERS : %-12" PRIu64 "\n", numb_clauses[EMPTY_INTERSECT]);
    printf("c numb cls LOOP_REMOVED : %-12" PRIu64 "\n", numb_clauses[LOOP_EXCLUDED]);
    printf("c numb cls EDGE_CONST   : %-12" PRIu64 "\n", numb_clauses[EDGE_CONST]);
    printf("c numb cls LOOP_CONST   : %-12" PRIu64 "\n", numb_clauses[LOOP_CONST]);

    if (totalGeneratedCls > 0)
    {
        printf("c avg lbd               : %g\n", lbd_clauses * 1.0 / totalGeneratedCls);
        printf("c falsified cls         : %-12" PRIu64 "   (%4.2f %%)\n", numb_clauses_falsified, numb_clauses_falsified * 1.0 / totalGeneratedCls);
        printf("c asserting cls         : %-12" PRIu64 "   (%4.2f %%)\n", numb_asserting_clauses, numb_asserting_clauses * 1.0 / totalGeneratedCls);
        printf("c added cls             : %-12" PRIu64 "   (%4.2f %%)\n", added_clauses, added_clauses * 1.0 / totalGeneratedCls);
    }
    printf("c avg size EMPTY_INTERS : %g\n", (numb_clauses[EMPTY_INTERSECT] > 0) ? size_clauses[EMPTY_INTERSECT] * 1.0 / numb_clauses[EMPTY_INTERSECT] : -1);
    printf("c avg size LOOP_REMOVED : %g\n", (numb_clauses[LOOP_EXCLUDED] > 0) ? size_clauses[LOOP_EXCLUDED] * 1.0 / numb_clauses[LOOP_EXCLUDED] : -1);
    printf("c avg size EDGE_CONST   : %g\n", (numb_clauses[EDGE_CONST] > 0) ? size_clauses[EDGE_CONST] * 1.0 / numb_clauses[EDGE_CONST] : -1);
    printf("c avg size LOOP_CONST   : %g\n", (numb_clauses[LOOP_CONST] > 0) ? size_clauses[LOOP_CONST] * 1.0 / numb_clauses[LOOP_CONST] : -1);
    printf("c time blackbox         : %g s\n", time_blackbox);
    printf("c time init             : %g s\n", time_init);
    printf("c time synch product    : %g s\n", time_sync_product);
    if (mem_used != 0)
        printf("c Memory used           : %.2f MB\n", mem_used);
    printf("c Real time             : %g s\n", total_time - time_blackbox);
}