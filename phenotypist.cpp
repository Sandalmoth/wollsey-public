#include "phenotypist.h"

#include <algorithm>
#include <vector>
#include <utility>


using namespace std;


void Phenotypist::add_phenotype(gene::id id, Phenotype ph) {
  phenotypes[id] = ph;
}


double Phenotypist::get_reproduction_rate(gene::id id) {
  if (unified_reproduction_rate)
    return reproduction_rate;
  else {
    auto it = phenotypes.find(id);
    if (it == phenotypes.end()) {
      // genotype has unknown phenotype values. Interpolate to find an approximate.
      // first, gather up a list of all reproduction rates with how different they are
      vector<pair<double, double> > rates;
      gene::gene gene(id);
      for (auto &ph: phenotypes) {
        gene::gene g(ph.first);
        rates.push_back(make_pair(gene::difference(gene, g), ph.second.reproduction_rate));
      }
      double irate = 0.0; // So tempted to name it pirate
      switch (interpolation) {
      case InterpolationMethod::MAX:
        break;
      case InterpolationMethod::WORST:
        irate = min_element(rates.begin(), rates.end())->second;
        break;
      case InterpolationMethod::LINEAR:
        break;
      default:
        cerr << "Phenotype interpolation method selection error" << endl;
        break;
      }
      // Add to 'caching'
      // FIXME this will fail for non-uniform death rates
      // Best way is to break out interpolation of both death and reproduction rates into it's own function and run it when necessary from get_repr_rate or get_death_rate
      phenotypes[id] = Phenotype{irate, 0.0};
      return irate;
    } else
      return it->second.reproduction_rate;
  }
}


// TODO fix non-uniform functionality here too
double Phenotypist::get_death_rate(gene::id id) const {
  if (unified_death_rate)
    return death_rate;
  else
    return phenotypes.at(id).death_rate;
}


// TODO and here
double Phenotypist::get_minimum_reproduction_rate(gene::id id) {
  if (unified_death_rate)
    return death_rate;
  else
    return phenotypes.at(id).death_rate;
}


ostream &operator<<(ostream &out, const Phenotypist &ph) {
  if (ph.unified_death_rate)
    cout << "Unified death rate: " << ph.death_rate << endl;
  if (ph.unified_reproduction_rate)
    cout << "Unified reproduction rate: " << ph.reproduction_rate << endl;
  if (ph.unified_minimum_reproduction_rate)
    cout << "Unified minimum reproduction rate: " << ph.minimum_reproduction_rate << endl;
  if (ph.unified_death_rate && ph.unified_reproduction_rate && ph.unified_minimum_reproduction_rate)
    return out;
  cout << "  Id\trepr. r.\tdeath r." << endl;
  for (auto x: ph.phenotypes) {
    cout << "  " << x.first << '\t' << x.second.reproduction_rate << '\t' << x.second.death_rate << endl;
  }
  return out;
}
