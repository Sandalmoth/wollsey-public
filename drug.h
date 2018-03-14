#ifndef __DRUG_H__
#define __DRUG_H__


#include <list>
#include <vector>
#include <map>
#include <iostream>

#include "gene.h"


class Drug {
public:
  Drug(std::vector<std::vector<double>> ic50s)
    : ic50s(ic50s) { }

  void set_ic50(gene::id id, double ic50);
  void set_wildtype_ic50(double ic50);
  double get_ic50(gene::id id);

  friend std::ostream &operator<<(std::ostream &out, const Drug &drug);

private:
  std::vector<std::vector<double>> ic50s;
  // Cache size unlimited, so we don't really care about order added or stuff like that
  std::map<gene::id, double> cache;

};


#endif
