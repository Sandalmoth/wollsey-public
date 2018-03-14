#ifndef __PHENOTYPIST_H__
#define __PHENOTYPIST_H__


#include <map>
#include <iostream>

#include "gene.h"


class Phenotypist {
public:
  Phenotypist(bool udr = false, double dr = 0.0, bool urr = false, double rr = 0.0, bool umrr = false, double mrr = 0.0)
    : unified_death_rate(udr)
    , death_rate(dr)
    , unified_reproduction_rate(urr)
    , reproduction_rate(rr)
    , unified_minimum_reproduction_rate(umrr)
    , minimum_reproduction_rate(mrr) { }

  struct Phenotype {
    double reproduction_rate;
    double death_rate;
    double minimum_reproduction_rate;
  };

  enum class InterpolationMethod {
    MAX,
    WORST,
    LINEAR
  };

  void add_phenotype(gene::id id, Phenotype ph);

  double get_reproduction_rate(gene::id id);
  double get_death_rate(gene::id id) const;
  double get_minimum_reproduction_rate(gene::id id);

  friend std::ostream &operator<<(std::ostream &out, const Phenotypist &ph);

private:
  bool unified_death_rate;
  double death_rate;
  bool unified_reproduction_rate;
  double reproduction_rate;
  bool unified_minimum_reproduction_rate;
  double minimum_reproduction_rate;
  InterpolationMethod interpolation = InterpolationMethod::WORST; // Do your worst!


  // TODO profile speed of using dna/gene/id for classification.
  // map works also as a cache (well, unlimited, so database?) for interpolated values
  std::map<gene::id, Phenotype> phenotypes;

};


#endif
