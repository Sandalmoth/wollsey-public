#include "protocol.h"

#include <iostream>


using namespace std;


void Protocol::add_protocol(string name, vector<double> protocol) {
  protocols.insert(make_pair(name, protocol));
}


void Protocol::add_protocol(string name, double (*generator)(int), int length) {
  vector<double> protocol;
  for (int i = 0; i < length; ++i) {
    protocol.push_back(generator(i));
  }
  protocols.insert(make_pair(name, protocol));
}


vector<pair<string, double>> Protocol::get_drug_concentrations(int timestep) {
  vector<pair<string, double>> concentrations;
  for (auto &p: protocols) {
    // try {
    //   concentrations.push_back(make_pair(p.first, p.second.at(timestep)));
    // } catch (exception &e) {
    //   concentrations.push_back(make_pair(p.first, 0.0));
    // }
    if (timestep < static_cast<int>(p.second.size())) {
      concentrations.push_back(make_pair(p.first, p.second[timestep]));
    } else {
        concentrations.push_back(make_pair(p.first, 0.0));
    }
    //concentrations.push_back(make_pair(p.first, 0.0));
  }
  return concentrations;
}


// Approximate memory usage in kB
int Protocol::memory_usage() {
  int m = 0;
  for (auto &p: protocols) {
    m += p.second.size();
  }
  return sizeof(double) * m / 1000;
}


std::ostream &operator<<(std::ostream &out, const Protocol &p) {
  cout << ' ';
  for (size_t i = 0; i < (*p.protocols.begin()).second.size(); ++i) {
    cout << "( ";
    for (auto &x: p.protocols) {
      cout << x.second[i] << ' ';
    }
    cout << ") ";
  }
  return out;
}
