#include "limited_branching_process.h"


#include <algorithm>
#include <cmath>


using namespace std;


void LimitedBranchingProcess::_run(Protocol protocol) {
  medic.reset();
  medic.set_protocol(protocol);
  timestep = 0;

  while (timestep < max_iters) {

    record.snapshot(cells, medic.get_drug_concentrations(timestep));

    LimitedBranchingProcess::_step();

    ++timestep;
    medic.advance();
  }

  record.snapshot(cells, medic.get_drug_concentrations(timestep));
  medic.drop_protocol();
}

// Deprecated in favour of altering dethrates alone for population control
/*
double LimitedBranchingProcess::growth_probability(gene::dna d) {
  double reproduction_rate =  phenotypist.get_reproduction_rate(gene::id(d)) * medic.get_fitness(gene::id(d));
  double death_rate = phenotypist.get_death_rate(gene::id(d));
  // We don't want to excessively exceed max_pop.
  // so, have reproduction_rate -> death_rate as size -> max_pop
  // What realistically limits cell growth is basically surface area
  // (look up Bertalanffy's equation)
  // Volume is basically proportional to size (number of cells)
  // so area ~ size^(2/3)
  // and the ratio size^(2/3) / size (= size^(-1/3)) should be what limits growth
  // This results in a Bertallanfy-type growth curve. (?)
  int sz = get_size();
  if (reproduction_rate > death_rate)
    reproduction_rate = reproduction_rate - (reproduction_rate - death_rate) * static_cast<double>(sz - min_pop*2) / static_cast<double>(max_pop - min_pop*2);
  return reproduction_rate;
}
*/


void LimitedBranchingProcess::_step() {
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
    double real_death_rate = phenotypist.get_death_rate(gene::id(c.first));
    double stable_death_rate = n_new_cells / n_cells;
    double death_rate = stable_death_rate;
    double relative_size_error = 0;
    if (n_cells > max_pop) {
      relative_size_error = static_cast<double>(n_cells - max_pop) / max_pop;
      death_rate = (real_death_rate > stable_death_rate) ? real_death_rate : stable_death_rate;
    } else if (n_cells < min_pop) relative_size_error = static_cast<double>(n_cells - min_pop) / min_pop;
    else death_rate = real_death_rate;

    death_rate *= 1.0 + relative_size_error * spring;

    // cout << n_cells << endl;
    // cout << "relative error " << relative_size_error << endl;
    // cout << "spring" << spring << endl;
    // cout << "real " << real_death_rate << endl;
    // cout << "stable " << stable_death_rate << endl;
    // cout << "final " << death_rate << endl;

    // As we don't want to go below our population minimum:
    // have death_rate -> 0 as size -> min_pop
    // Goes from 0 at size = min_pop to unmodified at size = 2*min_pop
    // int sz = get_size();
    // if (sz < min_pop) {
    //   death_rate = 0.0;
    // } else if (sz < min_pop*2 ){
    //   death_rate *= static_cast<double>(sz - min_pop) / static_cast<double>(min_pop);
    // }

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
