#include "drug.h"

#include <algorithm>


using namespace std;


void Drug::set_wildtype_ic50(double ic50) {
  // As the wildtype (all 0 genotype) now has it's own special ic50
  // we should null out it's related ic50s to make sure they don't accidentally
  // show up as maximum values.
  // This way it becomes possible for a mutant to have a lower ic50
  // than the wildtype.
  for (auto &v: ic50s) {
    v[0] = 0.0;
  }
  set_ic50(0, ic50);
}


void Drug::set_ic50(gene::id id, double ic50) {
  cache.insert(make_pair(id, ic50));
}


double Drug::get_ic50(gene::id id) {
  // cout << id << endl;
  auto it = cache.find(id);
  double ic50 = 0;
  if (it == cache.end()) {
    // New genotype. We need to calculate ic50.
    // Double mutants and worse will simply retain the highest resistance of
    // their components.
    // TODO restructure to avoid inefficiencies in back calculation.
    auto ogene = gene::gene(id);
    auto gene = ogene.get_x();
    vector<double> gene_ic50s;
    for (size_t i = 0; i < gene.size(); ++i) {
      gene_ic50s.push_back(ic50s[i][gene[i]]);
    }
    ic50 = *max_element(gene_ic50s.begin(), gene_ic50s.end());

    cache.insert(make_pair(id, ic50));
  } else {
    ic50 = it->second;
  }
  return ic50;
}


ostream &operator<<(ostream &out, const Drug &drug) {
  try {
    cout << "wt " << drug.cache.at(0);
  } catch (exception &e) {
    // do nothing
  }
  cout << " snps | ";
  for (auto &x: drug.ic50s) {
    for (auto ic50: x) {
      cout << ic50 << ' ';
    }
    cout << "| ";
  }
  return out;
}
