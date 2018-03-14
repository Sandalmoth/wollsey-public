#ifndef __RECORD_H__
#define __RECORD_H__


#include <map>
#include <vector>
#include <bitset>
#include <H5Cpp.h>

#include "gene.h"


enum Recordables {
  TIMELINE = 1,
  WT_HALF = 2,
  WT_NULL = 4
};


class Record {
public:
  void set_flags(int flags);
  void set_resolution(int _r) { resolution = _r; }

  void snapshot(const std::map<gene::id, unsigned long long> &cells, const std::vector<std::pair<std::string, double>> &drugs);
  void snapshot(const std::map<gene::dna, unsigned long long> &cells, const std::vector<std::pair<std::string, double>> &drugs);

  void store(H5::Group &gp_record);

  // History access functions
  int get_wt_half() { return wt_half; }

private:
  int resolution = 1000;
  int timestep = 0;

  // Recordables
  std::map<gene::id, std::vector<std::pair<int, int>>> timeline;
  std::map<std::string, std::vector<std::pair<int, double>>> protocol;
  int wt_half = -1;

  // Flags
  bool f_timeline = false;
  bool f_wt_half = false;
};


#endif
