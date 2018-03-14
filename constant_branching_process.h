#ifndef __CONSTANT_BRANCHING_PROCESS_H__
#define __CONSTANT_BRANCHING_PROCESS_H__

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


/*
  Variant of Branching Process that maintains population around a set number by manipulating the death rate.
*/


class ConstantBranchingProcess : public BranchingProcess {
public:
  ConstantBranchingProcess(int max_iters, double mut_prob, Medic medic, Phenotypist phenotypist, int population, double spring, int recordflags=0, int record_resolution=1000)
    : BranchingProcess(max_iters, mut_prob, medic, phenotypist, recordflags, record_resolution)
    , population(population)
    , spring(spring) {
  }

  ConstantBranchingProcess(const ConstantBranchingProcess &bp)
    : BranchingProcess(bp.max_iters, bp.mutation_probability, bp.medic, bp.phenotypist, bp.recordflags, bp.record_resolution) {
    population = bp.population;
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
    auto bp =  std::make_shared<ConstantBranchingProcess>(max_iters
                                                         , mutation_probability
                                                         , medic
                                                         , phenotypist
                                                         , population
                                                         , spring
                                                         , recordflags
                                                         , record_resolution
                                                         );
    for (auto &c: cells) {
      std::dynamic_pointer_cast<ConstantBranchingProcess>(bp)->add_cells(c.first, c.second);
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
  int population;
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
