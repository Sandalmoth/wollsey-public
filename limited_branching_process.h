#ifndef __LIMITED_BRANCHING_PROCESS_H__
#define __LIMITED_BRANCHING_PROCESS_H__

#include <map>
#include <random>
#include <thread>

#include "branching_process.h"
#include "gene.h"
#include "medic.h"
#include "phenotypist.h"
#include "protocol.h"
#include "record.h"
#include "simulator.h"


class LimitedBranchingProcess : public BranchingProcess {
public:
  LimitedBranchingProcess(int max_iters, double mut_prob, Medic medic, Phenotypist phenotypist, int min_pop, int max_pop, double spring, int recordflags=0, int record_resolution=1000)
    : BranchingProcess(max_iters, mut_prob, medic, phenotypist, recordflags, record_resolution)
    , min_pop(min_pop)
    , max_pop(max_pop)
    , spring(spring) {
  }

  LimitedBranchingProcess(const LimitedBranchingProcess &bp)
    : BranchingProcess(bp.max_iters, bp.mutation_probability, bp.medic, bp.phenotypist, bp.recordflags, bp.record_resolution) {
    min_pop = bp.min_pop;
    max_pop = bp.max_pop;
    spring = bp.spring;
  }

  // general interface
  virtual std::thread run(Protocol protocol) override {
    return std::thread([=]{
        _run(protocol);
      });
  }

  virtual Record &get_record() override {
    return record;
  }

  virtual std::shared_ptr<Simulator> get_instance() override {
    auto bp =  std::make_shared<LimitedBranchingProcess>(max_iters
                                                         , mutation_probability
                                                         , medic
                                                         , phenotypist
                                                         , min_pop
                                                         , max_pop
                                                         , spring
                                                         , recordflags
                                                         , record_resolution
                                                         );
    for (auto &c: cells) {
      std::dynamic_pointer_cast<LimitedBranchingProcess>(bp)->add_cells(c.first, c.second);
    }
    return bp;
  }

protected:
  void _run(Protocol protocol);
  void _step();
  // double growth_probability(gene::dna d);

  // int max_iters;
  // double mutation_probability;
  // int recordflags;
  // int record_resolution;
  int min_pop;
  int max_pop;
  double spring;

  // int timestep = 0;
  // std::map<gene::dna, unsigned long long> cells;
  // std::map<gene::dna, int> newborn;
  // std::mt19937 rng;
  // Record record;
  // Medic medic;
  // Phenotypist phenotypist;
};


#endif
