// -----------------------------------------------------------------------------
// Copyright (C) 2017  Ludovic LE FRIOUX
//
// PaInleSS is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
// -----------------------------------------------------------------------------

#include "painless.h"

#include "utils/Logger.h"
#include "utils/Parameters.h"
#include "utils/System.h"
#include "utils/SatUtils.h"

#include "solvers/SolverFactory.h"

#include "clauses/ClauseManager.h"

#include "sharing/HordeSatSharing.h"
#include "sharing/SimpleSharing.h"
#include "sharing/StrengtheningSharing.h"
#include "sharing/Sharer.h"

#include "working/SequentialWorker.h"
#include "working/Portfolio.h"

#include <unistd.h>

using namespace std;

// -------------------------------------------
// Declaration of global variables
// -------------------------------------------
atomic<bool> globalEnding(false);

Sharer **sharers = NULL;

int nSharers = 0;

WorkingStrategy *working = NULL;

SatResult finalResult = UNKNOWN;

vector<int> finalModel;

// -------------------------------------------
// Main of the framework
// -------------------------------------------
int main(int argc, char **argv)
{
   Parameters::init(argc, argv);

   if (Parameters::getFilename() == NULL ||
       Parameters::getBoolParam("h"))
   {
      cout << "USAGE: " << argv[0] << " [options] input.cnf" << endl;
      cout << "Options:" << endl;
      cout << "\t-c=<INT>\t\t number of cpus, default is 24" << endl;
      cout << "\t-max-memory=<INT>\t memory limit in GB, default is 200" << endl;
      cout << "\t-t=<INT>\t\t timeout in seconds, default is no limit" << endl;
      cout << "\t-lbd-limit=<INT>\t LBD limit of exported clauses, default is 2" << endl;
      cout << "\t-shr-sleep=<INT>\t time in useconds a sharer sleep each "
              "round, default is 500000 (0.5s)"
           << endl;
      cout << "\t-shr-lit=<INT>\t\t number of literals shared per round, default is 1500" << endl;
      cout << "\t-v=<INT>\t\t verbosity level, default is 0" << endl;
      cout << "\t-no-blackbox \t\t disable cdcl-bmc calls." << endl;
      cout << "\t-prop-rate-inf=<INT>\t\t percentage of changes to call the blackbox (Low limit), default is 10." << endl;
      cout << "\t-prop-rate-sup=<INT>\t\t percentage of changes to call the blackbox (High limit), default is 10." << endl;
      cout << "\t-depth-rate=<INT>\t\t percentage of global decision to call the blackbox, default is 100." << endl;
      cout << "\t-reducer\t\t Enable reducer." << endl;
      cout << "\t-lbd-limit-bmc=<INT>\t LBD limit of exported clauses, default is -1 (unlimited)" << endl;

      cout << "\t-pin=<INT>\t\t Pin thread id>=0, to disable pining pin=-1, default is -1." << endl;

      cout << "\t-oneTrail\t\t Enable to treat only one trail at a time (no addition sharer)." << endl;
      cout << "\t-no-share\t\t Disable sharing of blackbox clauses to others solvers." << endl;
      cout << "\t-fic-lbd\t\t Creates fictif LBD of blackbox clauses = min(clause.size, 6)." << endl;
      cout << "\t-sat-redund\t\t Do not share satisfied clauses from blackbox constraints." << endl;

      return 0;
   }

   Parameters::printParams();

   int cpus = Parameters::getIntParam("c", 1);
   setVerbosityLevel(Parameters::getIntParam("v", 0));
   int pin = Parameters::getIntParam("pin", -1);

   // Create and init solvers
   vector<SolverInterface *> solvers;
   vector<SolverInterface *> solvers_VSIDS;
   vector<SolverInterface *> solvers_LRB;

   EnvBMC *bxbmc = NULL;
   if (!Parameters::getBoolParam("no-blackbox"))
   {
      bxbmc = SolverFactory::createBlackboxBMC();
      SolverFactory::createMapleCOMSPSSolvers(cpus, solvers, bxbmc);
      if (bxbmc->seq)
         solvers.push_back(bxbmc);
   }
   else
   {
      SolverFactory::createMapleCOMSPSSolvers(cpus, solvers);
   }
   int nSolvers = solvers.size();

   if(bxbmc && !bxbmc->seq){
      SolverFactory::nativeDiversification(solvers);
      for (int id = 0; id < nSolvers; id++)
      {
         if (id % 2)
         {
            solvers_LRB.push_back(solvers[id]);
         }
         else
         {
            solvers_VSIDS.push_back(solvers[id]);
         }
      }
      SolverFactory::sparseRandomDiversification(solvers_LRB);
      SolverFactory::sparseRandomDiversification(solvers_VSIDS);

      SolverFactory::sparseRandomDiversification(solvers);
   }
   // Init Sharing
   // 15 CDCL, 1 Reducer producers by Sharer
   vector<SolverInterface *> prod1;
   vector<SolverInterface *> prod2;
   vector<SolverInterface *> reducer1;
   vector<SolverInterface *> reducer2;
   vector<SolverInterface *> cons1;
   vector<SolverInterface *> cons2;
   vector<SolverInterface *> consCDCL;

   switch (Parameters::getIntParam("shr-strat", 0))
   {
   case 1:
      printf("c SimpleSharing without reducer\n");

      if (bxbmc)
      {
         nSharers = 2;
         prod1.insert(prod1.end(), solvers.begin(), solvers.end());
         prod1.insert(prod1.end(), bxbmc);
         cons1.insert(cons1.end(), solvers.begin(), solvers.end());

         prod2.insert(prod2.end(), solvers.begin(), solvers.end());
         cons2.insert(cons2.end(), bxbmc);

         solvers.push_back(bxbmc);
         nSolvers = solvers.size();
         sharers = new Sharer *[nSharers];

         sharers[0] = new Sharer(1, new SimpleSharing(), prod1, cons1); // share conflict clauses
         sharers[1] = new Sharer(2, new SimpleSharing(), prod2, cons2); // share trails
      }
      else
      {
         nSharers = 1;
         prod1.insert(prod1.end(), solvers.begin(), solvers.end());
         cons1.insert(cons1.end(), solvers.begin(), solvers.end());
         sharers = new Sharer *[nSharers];
         if (pin >= 0)
            sharers[0] = new Sharer(1, pin, new SimpleSharing(), prod1, cons1);
         else
            sharers[0] = new Sharer(1, new SimpleSharing(), prod1, cons1);
      }
      break;

   case 2:
      nSharers = 1;
      prod1.insert(prod1.end(), solvers.begin(), solvers.end());
      cons1.insert(cons1.end(), solvers.begin(), solvers.end());
      if (bxbmc && bxbmc->oneTrailAtTime)
      {
         std::cout << "c One trail evaluation: ON\n";
         solvers.push_back(bxbmc);
         nSolvers = solvers.size();
         prod1.insert(prod1.end(), bxbmc);
      }
      sharers = new Sharer *[nSharers];
      sharers[0] = new Sharer(1, new SimpleSharing(), prod1, cons1);
      break;

   case 5:
      if (cpus > 1)
      {
         solvers.push_back(SolverFactory::createReducerSolver(SolverFactory::createMapleCOMSPSSolver()));
         nSolvers = solvers.size();
      }
      nSharers = 1;
      prod1.insert(prod1.end(), solvers.begin(), solvers.end());
      cons1.insert(cons1.end(), solvers.begin(), solvers.end());
      sharers = new Sharer *[nSharers];
      sharers[0] = new Sharer(1, new SimpleSharing(), prod1, cons1);
      break;
   case 6:
      printf("c SimpleSharing without reducer for BMC\n");
      nSharers = 1;
      // prod1.emplace_back(solvers[0]);
      // prod1.emplace_back(solvers[1]);
      // cons1.emplace_back(solvers[0]);
      prod1.insert(prod1.end(), solvers.begin(), solvers.end());
      cons1.insert(cons1.end(), solvers.begin(), solvers.end());
      sharers = new Sharer *[nSharers];
      sharers[0] = new Sharer(1, new SimpleSharing(), prod1, cons1);
      break;

   default:
      break;
   }

   // Init working
   working = new Portfolio();
   if (pin >= 0)
   {
      for (size_t i = 0; i < nSolvers; i++)
      {
         working->addSlave(new SequentialWorker(solvers[i], i + 1 + pin));
      }
   }
   else
   {
      for (size_t i = 0; i < nSolvers; i++)
      {
         working->addSlave(new SequentialWorker(solvers[i]));
      }
   }

   // Init the management of clauses
   ClauseManager::initClauseManager();

   // Launch working
   vector<int> cube;
   working->solve(cube);

   // Wait until end or timeout
   int timeout = Parameters::getIntParam("t", -1);

   while (globalEnding == false)
   {
      sleep(1);

      if (timeout > 0 && getRelativeTime() >= timeout)
      {
         globalEnding = true;
         working->setInterrupt();
      }
   }

   // Delete sharers
   for (int id = 0; id < nSharers; id++)
   {
      sharers[id]->printStats();
      delete sharers[id];
   }
   delete sharers;

   // Print solver stats
   SolverFactory::printStats(solvers);

   // Delete working strategy
   delete working;

   // Delete shared clauses
   ClauseManager::joinClauseManager();

   // Print the result and the model if SAT
   double total_time = getRelativeTime();
   cout << "c Resolution time: " << total_time << " sec" << endl;
   if (bxbmc && bxbmc->init_box_done)
      bxbmc->printStats(total_time);
   if (finalResult == SAT)
   {
      cout << "s SATISFIABLE" << endl;

      if (Parameters::getBoolParam("no-model") == false)
      {
         printModel(finalModel);
      }
   }
   else if (finalResult == UNSAT)
   {
      cout << "s UNSATISFIABLE" << endl;
   }
   else
   {
      cout << "s UNKNOWN" << endl;
   }

   return 0;
}
