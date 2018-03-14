#ifndef __BRANCHING_PROCESS_H__
#define __BRANCHING_PROCESS_H__

#include <map>
#include <random>
#include <thread>

#include "gene.h"
#include "medic.h"
#include "phenotypist.h"
#include "protocol.h"
#include "record.h"
#include "simulator.h"


class BranchingProcess : public Simulator {
public:
  BranchingProcess(int max_iters, double mut_prob, Medic medic, Phenotypist phenotypist, int recordflags=0, int record_resolution=1000)
    : max_iters(max_iters)
    , mutation_probability(mut_prob)
    , medic(medic)
    , phenotypist(phenotypist)
    , recordflags(recordflags)
    , record_resolution(record_resolution) {
    _init();
  }
  BranchingProcess(const BranchingProcess &bp) {
    max_iters = bp.max_iters;
    mutation_probability = bp.mutation_probability;
    medic = bp.medic;
    phenotypist = bp.phenotypist;
    recordflags = bp.recordflags;
    record_resolution = bp.record_resolution;

    _init();
  }

  int get_size();

  void add_cells(gene::dna d, int n);

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
    auto bp =  std::make_shared<BranchingProcess>(max_iters
                                                  , mutation_probability
                                                  , medic
                                                  , phenotypist
                                                  , recordflags
                                                  , record_resolution
                                                  );
    for (auto &c: cells) {
      std::dynamic_pointer_cast<BranchingProcess>(bp)->add_cells(c.first, c.second);
    }
    return bp;
  }

protected:
  void _run(Protocol protocol);
  void _step();
  void _init();
  double growth_probability(gene::dna d);

  void add_newborn(gene::dna d, int n);
  void merge_newborn();

  int max_iters;
  double mutation_probability;
  Medic medic;
  Phenotypist phenotypist;
  int recordflags;
  int record_resolution;

  int timestep = 0;
  std::map<gene::dna, unsigned long long> cells;
  std::map<gene::dna, int> newborn;
  std::mt19937 rng;
  Record record;
};


#endif
