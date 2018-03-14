#include "medic.h"

#include <cmath>
#include <numeric>


using namespace std;


template <typename Titerator>
double harmonic_mean(Titerator first, Titerator last) {
  int N = last - first;
  double rsum = accumulate(first, last, 0.0,
                           [](double a, double b) {
                             return a + 1.0/b;
                           });
  return N / rsum;
}


void Medic::add_drug(std::string name, Drug drug) {
  drugs.insert(make_pair(name, drug));
}


void Medic::set_protocol(Protocol p) {
  protocol = p;
}


void Medic::drop_protocol() {
  // Protocols can be quite big, so this is a convenient memory-saving function.
  protocol = Protocol();
}


double Medic::get_fitness(gene::id id) {
  auto drug_concentrations = protocol.get_drug_concentrations(timestep);
  // for (auto dc: drug_concentrations)
    // cout << dc.first << ' ' << dc.second << ' ';
  // cout << endl;
  vector<double> inhibition_coefficients;
  for (auto &dc: drug_concentrations) {
    if (dc.second == 0.0)
      continue;
    double ic50 = drugs.at(dc.first).get_ic50(id);
    // cout << -dc.second/ic50 << endl;
    // cout << dc.second << endl;
    // cout << ic50 << endl;
    // exit(0);
    inhibition_coefficients.push_back(exp2(-dc.second/ic50));
  }
  // for (auto dc: drug_concentrations)
  //   cout << dc.first << ' ' << dc.second << ' ';
  // cout << endl;
  // TODO Benchmark switch implementation against else-if implementation.
  // switch (inhibition_coefficients.size()) {
  // case 0:
  //   cout << "no drugs at " << timestep << endl;
  //   return 1.0;
  // case 1:
  //   cout << "1 drug with ic " << inhibition_coefficients[0] << " at " << timestep << endl;
  //   return inhibition_coefficients[0];
  // default:
  //   cout << inhibition_coefficients.size() << " drugs with ic " << harmonic_mean(inhibition_coefficients.begin(), inhibition_coefficients.end()) << " at " << timestep << endl;
  //   return harmonic_mean(inhibition_coefficients.begin(), inhibition_coefficients.end());
  // }
  if (inhibition_coefficients.size() == 0) {
    // no inhibitors
      // cout << "no drugs at " << timestep << endl;
    return 1.0;
  } else if (inhibition_coefficients.size() == 1) {
      // cout << "1 drug with ic " << inhibition_coefficients[0] << " at " << timestep << endl;
    return inhibition_coefficients[0];
  } else {
    double inhibition = harmonic_mean(inhibition_coefficients.begin(), inhibition_coefficients.end());
    // cout << inhibition_coefficients.size() << " drugs with ic " << inhibition << " at " << timestep << endl;
    // Inhibition is assumed to be proportional to fitness
    return inhibition;
  }
}

ostream &operator<<(ostream &out, const Medic &medic) {
  for (auto &d: medic.drugs) {
    cout << d.first << '\t' << d.second << '\n';
  }
  return out;
}
