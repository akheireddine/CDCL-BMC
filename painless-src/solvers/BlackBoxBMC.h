#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <map>

#include <spot/twa/twagraph.hh>

#include "../solvers/SolverInterface.h"
#include "../solvers/Reducer.h"
#include "../clauses/ClauseBuffer.h"
#include "../utils/Threading.h"
#include "../utils/Parameters.h"
#include "../../mapleCOMSPS/mapleCOMSPS/core/SolverTypes.h"

typedef struct variableBelonging
{
    int num_step = -1;
    int num_partition = -1;
    bool in_property = false;
    std::string name="";
} VB;

namespace MapleCOMSPS
{
    class Lit;
    class Solver;
    template <class T>
    class vec;
    class lbool;
}

enum CONSTYPE
{
    EMPTY_INTERSECT = 0,
    LOOP_EXCLUDED = 1,
    EDGE_CONST = 2,
    LOOP_CONST = 3
};

class EnvBMC : public SolverInterface
{
public:
    EnvBMC(int id, bool inters = true,
           bool loopcond = true, bool edgecond = true,
           bool looprm = true);
    ~EnvBMC() {}

    //**************************  CNF < ---- > BMC ******************************************//

    void initialize_info_variables();
    inline bool is_property_variable(int v) { return info_variables[v].in_property; };
    inline int get_step(int v) { return info_variables[v].num_step; };
    inline bool is_ltl_property(int v) { return is_property_variable(v) && get_step(v) > -1; };
    inline bool is_model_variable(int v) { return get_step(v) > -1; };

    inline std::string get_name_variable(int v) { return info_variables[v].name; };

    /// Working strategies BMC parameters
    unsigned max_vars;
    int max_vars_property;
    int maxVarsPure;
    int minVarsPure;
    int nb_property_variable;
    int nb_model_variable;
    int nb_junction_variable;
    unsigned k_and_fictif_state, K;
    bool init_box_done;
    bool stop;
    bool oneTrailAtTime;

    std::vector<VB> info_variables;
    std::vector<std::map<std::string, int>> name_variable_id;
    std::vector<std::vector<std::tuple<int, std::string>>> pure_vars;

    //************************** CNF < ---- > AUT ******************************************//
    void init();
    SatResult solve(const vector<int> &cube);
    void create_word_automata();
    void create_word_automata(int loop);

    void explore_automata_each_loop(MapleCOMSPS::vec<MapleCOMSPS::vec<MapleCOMSPS::Lit>> &out_refined);
    void explore_automata_each_loop_parallel();
    void initialize_literals_map();
    void initialize_trail_prop();
    void reinitialize_trail_map();
    void constraints_to_remove_loop(MapleCOMSPS::vec<int> &cls_result, int loop);
    void constraints_on_transition_loop(MapleCOMSPS::vec<int> &cls_result, std::vector<int> &constraints, int loop);
    void constraints_on_edge(MapleCOMSPS::vec<int> &cls_result, std::vector<int> &constraints);
    void add_trail_to_cls(MapleCOMSPS::vec<int> &cls);

    std::vector<std::vector<bdd>> loops_edges_constraints();
    std::vector<std::vector<bdd>> loops_edges_constraints(int loop);
    inline void set_ltlspec(std::string l) { ltlspec = l; };
    int getTotalClauses() { return accumulate(numb_clauses.begin(), numb_clauses.end(), 0); };
    void addLoopClauses();
    void create_trail_cls(MapleCOMSPS::vec<int> &cls_result, int loop);

    //******************* EXPORT/IMPORT CLAUSES ***********************//
    /// Load formula from a given dimacs file, return false if failed.
    bool loadFormula(const char *filename);

    /// Get the number of variables of the current resolution.
    int getVariablesCount();

    /// Get a variable suitable for search splitting.
    int getDivisionVariable();

    /// Set initial phase for a given variable.
    void setPhase(const int var, const bool phase);

    /// Bump activity of a given variable.
    void bumpVariableActivity(const int var, const int times);

    /// Interrupt resolution, solving cannot continue until interrupt is unset.
    void setSolverInterrupt();

    /// Remove the SAT solving interrupt request.
    void unsetSolverInterrupt();

    /// Add a permanent clause to the formula.
    void addClause(ClauseExchange *clause);

    /// Add a list of permanent clauses to the formula.
    void addClauses(const vector<ClauseExchange *> &clauses);

    /// Add a list of initial clauses to the formula.
    void addInitialClauses(const vector<ClauseExchange *> &clauses);

    /// Add a learned clause to the formula.
    void addLearnedClause(ClauseExchange *clause);

    /// Add a list of learned clauses to the formula.
    void addLearnedClauses(const vector<ClauseExchange *> &clauses);

    /// Get a list of learned clauses.
    void getLearnedClauses(vector<ClauseExchange *> &clauses);

    void getTrailClauses(vector<ClauseExchange *> & clauses){ exit(3); };

    /// Request the solver to produce more clauses.
    void increaseClauseProduction();

    /// Request the solver to produce less clauses.
    void decreaseClauseProduction();

    /// Get solver statistics.
    SolvingStatistics getStatistics();

    /// Return the model in case of SAT result.
    vector<int> getModel();

    /// Native diversification.
    void diversify(int id);

    vector<int> getFinalAnalysis();

    vector<int> getSatAssumptions();

    bool ImportClause();
    bool ImportTrailAtATime();
    bool tryLockTrail(ClauseExchange *trail);

    bool ExportClause(MapleCOMSPS::vec<int> &cls, CONSTYPE const_type);
    bool prepareClause(MapleCOMSPS::vec<int> &cls, MapleCOMSPS::vec<MapleCOMSPS::Lit> &out_cls, CONSTYPE const_type);
    void createLoopEdgeClauses(bdd formula, MapleCOMSPS::vec<MapleCOMSPS::vec<MapleCOMSPS::Lit>> &out_cls, int curr_step, int loop);

    /// Buffer used to import clauses (units included).
    ClauseBuffer trailToImport;
    /// Buffer used to export clauses (units included).
    ClauseBuffer clausesToExport;

    /// Buffer used to add permanent clauses.
    ClauseBuffer clausesToAdd;

    Reducer *reducer;
    MapleCOMSPS::Solver *solver;
    //****************************************************************//

    atomic<int> lockTrail;
    ClauseExchange *new_trail;

    spot::twa_graph_ptr aut;
    spot::twa_graph_ptr word_aut;
    std::string ltlspec;

    std::vector<std::map<std::string, MapleCOMSPS::lbool>> trail_prop_spot;
    std::vector<int> curr_trail;
    std::map<int, int> curr_levels;
    int curr_lbd;

    int lbd_limit;

    std::vector<int> rename_loop;

    float prop_deci_rate_inf,prop_deci_rate_sup,depth_rate;
    int num_edge_loop;
    bool seq, noshare;

    //************************** STATS  ******************************************//
    std::vector<uint64_t> numb_clauses, size_clauses, used_clauses_conflicts;
    uint64_t numb_call, numb_clauses_falsified, lbd_clauses, numb_asserting_clauses, unit_clauses, added_clauses, remain_clauses;
    double time_init, time_blackbox, time_sync_product;
    bool enable_intersection, enable_loop_cond, enable_edge_cond, enable_loop_rm;

    void printStats(double total_time);
};