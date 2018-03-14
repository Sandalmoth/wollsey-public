#include "constant_branching_process.h"


#include <algorithm>
#include <cmath>


using namespace std;


void ConstantBranchingProcess::_run(Protocol protocol) {
  medic.reset();
  medic.set_protocol(protocol);
  timestep = 0;

  while (timestep < max_iters) {

    record.snapshot(cells, medic.get_drug_concentrations(timestep));

    ConstantBranchingProcess::_step();

    ++timestep;
    medic.advance();
  }

  record.snapshot(cells, medic.get_drug_concentrations(timestep));
  medic.drop_protocol();
}


void ConstantBranchingProcess::_step() {
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
    // double min_death_rate = phenotypist.get_death_rate(gene::id(c.first));
    // as death rate is used to balance population, we can safely ignore the realistic death_rates
    double death_rate = n_new_cells / n_cells;
    int size_error = n_cells - population;
    double relative_size_error = static_cast<double>(size_error) / population;
    death_rate *= 1.0 + (relative_size_error) * spring;
    // cout << "size_error " << size_error << ' ' << population << endl;
    // cout << "relative_size_error " << relative_size_error << endl;
    // cout << "expected number of new cells " << n_new_cells << std::endl;
    // cout << "expected number of deaths " << death_rate * get_size() << std::endl;
    // cout << "n_cells " << n_cells << endl;
    // cout << "expected adjusted number of deaths " << death_rate * get_size() << std::endl;

    // scale it to compensate for killing before reproducing (when it should be simultaneous)
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
