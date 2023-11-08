/****************************************************************************************[Dimacs.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef MapleCOMSPS_Dimacs_h
#define MapleCOMSPS_Dimacs_h

#include <stdio.h>
#include <unordered_set>
#include <tuple>
#include <string>

#include "../utils/ParseUtils.h"
#include "../core/SolverTypes.h"
#include "../core/Solver.h"



namespace MapleCOMSPS {

//=================================================================================================
// DIMACS Parser:

template<class B, class Solver>
static void readClause(B& in, Solver& S, vec<Lit>& lits) {
    int     parsed_lit, var;
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        while (var >= S.nVars()) S.newVar();
        lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
}

template<class B>
static void readClause_(B& in, vec<Lit>& lits, int& total_vars) {
    int     parsed_lit, var;
    using std::vector;
    std::vector<int> clause;
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        while (std::abs(parsed_lit) >= total_vars) total_vars++;
        clause.emplace_back(parsed_lit);
        var = abs(parsed_lit)-1;
        lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
}

// template<class B, class Solver, class EnvBMC>
// static void readClauseLoop(B& in, Solver& S, EnvBMC& env, int& total_vars) {
//     int     parsed_lit, var, id_cls;
//     int loop = parseInt(in);

//     while (env.clauses_loop_neg.size() <= loop)
//         env.clauses_loop_neg.push();

//     env.clauses_loop_neg[loop].push();
//     id_cls = env.clauses_loop_neg[loop].size() - 1;
//     while(*in != ':')
//         ++in;
//     ++in;
//     for (;;){
//         parsed_lit = parseInt(in);
//         if (parsed_lit == 0) break;
//         while (std::abs(parsed_lit) >= total_vars) total_vars++;
//         var = abs(parsed_lit)-1;
//         // while (var >= S.nVars()) S.newVar(true,false);
//         env.clauses_loop_neg[loop][id_cls].push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
//     }
// }

template<class B>
static bool operator_ltlspec(B& in, bool bis=false){
    if(bis){
        return (*in == '(' || *in == ')' || *in == '!' || 
                *in == '&' || *in == '|' || *in == '<' || 
                *in == '>' || *in == ' ' || *in == EOF);
    }
    return (*in == '(' || *in == ')' || *in == '!' || *in == 'G' || *in == 'F' || *in == 'U' || 
            *in == 'X' || *in == 'R' || *in == 'W' || *in == '&' || *in == '|' || *in == '<' || 
            *in == '>' || *in == '-' || *in == ' ' || *in == 'V');
}

template<class B, class EnvBMC>
static void readLTLSPEC(B& in, EnvBMC& env){
    std::string ltl="";
    for(;;){
        if (*in == EOF)
            break;
        if(operator_ltlspec(in)){
            ltl += *in;
            ++in;
        }
        else{
            ltl += "\"";
            while(!operator_ltlspec(in,true)){
                ltl += *in;
                ++in;
            }
            ltl += "\"";
        }
    }
    ltl.pop_back();ltl.pop_back();ltl.pop_back();
    env.set_ltlspec(ltl);
}


template<class B, class Solver>
static void parse_DIMACS_main(B& in, Solver& S) {
    vec<Lit> lits;
    int vars    = 0;
    int clauses = 0;
    int cnt     = 0;
    for (;;){
        skipWhitespace(in);
        if (*in == EOF) break;
        else if (*in == 'p'){
            if (eagerMatch(in, "p cnf")){
                vars    = parseInt(in);
                clauses = parseInt(in);
                // SATRACE'06 hack
                // if (clauses > 4000000)
                //     S.eliminate(true);
            }else{
                printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
        } else if (*in == 'c' || *in == 'p')
            skipLine(in);
        else{
            cnt++;
            readClause(in, S, lits);
            S.addClause_(lits); }
    }
    if (vars != S.nVars())
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of variables.\n");
    if (cnt  != clauses)
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of clauses.\n");
}

template<class B>
int parseVar(B& in)
{
    int parsed_val = -1, var = -1;
    for (;;)
    {
        skipWhitespace(in);
        if( *in < '0' || *in > '9' )
        {
            ++in;
            continue;
        }
        parsed_val = parseInt(in);
        if (parsed_val > 0)
        {
            var = parsed_val ;
            break;
        }
    }
    return var;
}


template<class B>
int parseVarLoop(B& in)
{
    int parsed_val = -1, var = -1;
    for (;;)
    {
        skipWhitespace(in);
        if( *in < '0' || *in > '9' )
        {
            ++in;
            continue;
        }
        parsed_val = parseInt(in);
        if (parsed_val >= 0)
        {
            var = parsed_val ;
            break;
        }
    }
    return var;
}

template<class B>
std::string parseNameVar(B& in)
{
    std::string name = "";
    for (;;)
    {
        skipWhitespace(in);
        if(eagerMatch(in,"Variable")){
            skipWhitespace(in);
            while(*in != '\n'){
                name += *in;
                ++in;
            }
            break;
        }
        else
        {
            ++in;
            continue;
        }
    }
    return name;
}

template<class B, class EnvBMC>
static void parse_NuSMV_DIMACS_main(B& in, EnvBMC& env)
{
    // gzFile in_gz = gzopen(filename, "rb");
    // StreamBuffer_ in(in_gz);
    vec<Lit> lits;
    int cnt        = 0;
    int total_vars = 0;
    int nb_vars    = 0;
    std::vector<std::tuple<int,std::string>> set_vars_time;
    env.max_vars = 0;
    // Get information about pure variables
    for(;;)
    {
        if ( *in == 'c' )
        {
            ++in; ++in;
            // c @@@@@@ Time
            if ( *in == '@' )
            {
                if( !set_vars_time.empty() )
                     env.pure_vars.emplace_back( set_vars_time );
                set_vars_time.clear();
            }
            // c CNF variable
            else if ( *in == 'C' )
            {
                int var_time = parseVar(in);
                std::string name = parseNameVar(in);
                set_vars_time.emplace_back( std::make_tuple(var_time,name) );
            }
            // c LOOP variable
            else if ( *in == 'L' )
            {
                env.rename_loop.emplace_back(parseVarLoop(in));
            }
            // c model
            else if ( *in == 'm' )
            {
                skipLine(in);
                break;
            }
            skipLine(in);
        }
        else
            ++in;
    }
    if( set_vars_time.empty() )
        throw std::runtime_error("dimacs: parse_NUMSV_DIMACS");
    env.pure_vars.emplace_back( set_vars_time );


    // Do normal parsing
    for (;;)
    {
        skipWhitespace(in);
        if (*in == EOF)
            break;
        else if (*in == 'p')
        {
            if (eagerMatch(in, "p cnf"))
            {
                env.max_vars    = parseInt(in);
                nb_vars = env.max_vars;
                env.max_vars_property = env.max_vars;
                parseInt(in);
            }
            else
                printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
        }
        else if (*in == 'c'){
            ++in;++in;
            // GET information about LTLSPEC
            if( *in == 'L'){
                ++in;
                if ( *in=='T' && eagerMatch(in, "TLSPEC "))
                    readLTLSPEC(in, env);
            }
            else
                skipLine(in);
        }
        else if (*in == 'p')
            skipLine(in);
        else
        {
            cnt++;
            readClause_(in,lits, total_vars);
        }
    }
    env.max_vars = total_vars;
}







// Inserts problem into solver.
//
template<class Solver>
static void parse_DIMACS(gzFile input_stream, Solver& S) {
    StreamBuffer in(input_stream);
    parse_DIMACS_main(in, S); }

template<class EnvBMC>
static void parse_DIMACS_BMC(EnvBMC& env, gzFile input_stream)
{
    StreamBuffer in(input_stream);

    parse_NuSMV_DIMACS_main(in, env);
}
//=================================================================================================

}

#endif
