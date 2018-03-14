#ifndef __MEDIC_H__
#define __MEDIC_H__


#include <string>
#include <map>
#include <iostream>

#include "drug.h"
#include "protocol.h"


class Medic {
public:
  Medic() { }
  Medic(const Medic &medic) {
    drugs = medic.drugs;
  }

  void add_drug(std::string name, Drug drug);
  void set_protocol(Protocol p);
  void drop_protocol();
  double get_fitness(gene::id id);
  std::vector<std::pair<std::string, double>> get_drug_concentrations(int timestep) {
    return protocol.get_drug_concentrations(timestep);
  }


  void reset() {
    timestep = 0;
  }

  void advance() {
    ++timestep;
  }

  friend std::ostream &operator<<(std::ostream &out, const Medic &medic);

private:
  int timestep = 0;
  // Consider indexable datastructure to get rid of nlogn complexity.
  std::map<std::string, Drug> drugs;

  Protocol protocol;

};


#endif
