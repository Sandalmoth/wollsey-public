#include "branching_process.h"


#include <algorithm>


using namespace std;


void BranchingProcess::_init() {
  random_device rd;
  rng.seed(rd());
  record.set_flags(recordflags);
  record.set_resolution(record_resolution);
}


void BranchingProcess::_run(Protocol protocol) {
  medic.reset();
  medic.set_protocol(protocol);
  timestep = 0;

  while (timestep < max_iters) {

    record.snapshot(cells, medic.get_drug_concentrations(timestep));

    _step();

    ++timestep;
    medic.advance();
  }

  record.snapshot(cells, medic.get_drug_concentrations(timestep));
  medic.drop_protocol();
}


double BranchingProcess::growth_probability(gene::dna d) {
  gene::id id(d);
  return phenotypist.get_reproduction_rate(id) * medic.get_fitness(id)
    + phenotypist.get_minimum_reproduction_rate(id);
}


int BranchingProcess::get_size() {
	return accumulate(cells.begin(), cells.end(), 0,
                    [](const int& left, const decltype(cells)::value_type &right) {
                      return left + right.second;
                    });
}


void BranchingProcess::_step(){
  // first, select some cells to kill
  // we need to avoid inducing drift when we kill first and reproduce second
  // so we scale the number of kills by how many we expect to reproduce
  double n_cells = get_size();
  double n_new_cells = accumulate(cells.begin(), cells.end(), 0.0
                                  , [this](const double &left, const decltype(cells)::value_type &right) {
                                    return left + right.second * growth_probability(right.first);
                                  });

  // kill cells
  for (auto &c: cells) {
    double death_rate = phenotypist.get_death_rate(gene::id(c.first));
    // scale it as killing before reproducing when they are supposed to be simultaneous
    // would induce drift
    death_rate *= n_cells / (n_cells + n_new_cells);
    c.second -= binomial_distribution<int>(c.second, death_rate)(rng);
  }

  // remove any genotype reduced to size 0 (or less (?))
  for (auto it = cells.begin(); it != cells.end(); /*no incrementation*/) {
    if (it->second <= 0) {
      it = cells.erase(it); // returns iterator to element that follows removed element
    } else {
      ++it;
    }
  }

  // now divide remaining cells
  for (auto &c: cells) {
    // find the number of cells (of this genotype) that divide
    int growth = binomial_distribution<int>(c.second, growth_probability(c.first))(rng);
    // if we actually get new cells, some of them might be mutants
    if (growth) {
      int mutants = binomial_distribution<int>(growth, mutation_probability)(rng);
      // if we actually get mutants, mutate some cells
      for (int i = 0; i < mutants; ++i) {
        gene::dna mutant = c.first;
        mutant.mutate_snv(rng);
        // cout << timestep << " Mutant " << gene::gene(mutant) << ((gene::id(mutant) == 0) ? " synonymous" : " missense") << endl;
        add_newborn(mutant, 1);
      }
      c.second += (growth - mutants); // add all non-mutant cells to cell group
    }
  }

  merge_newborn();
}


void BranchingProcess::add_cells(gene::dna d, int n) {
  auto it = cells.find(d);
  if (it == cells.end()) {
    cells.insert(make_pair(d, n));
  } else {
    it->second += n;
  }
}


void BranchingProcess::add_newborn(gene::dna d, int n) {
  auto it = newborn.find(d);
  if (it == newborn.end()) {
    newborn.insert(make_pair(d, n));
  } else {
    it->second += n;
  }
}


void BranchingProcess::merge_newborn() {
  for (auto &nb: newborn) {
    add_cells(nb.first, nb.second);
  }
  newborn.clear();
}
