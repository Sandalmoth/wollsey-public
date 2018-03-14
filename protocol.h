#ifndef __PROTOCOL_H__
#define __PROTOCOL_H__


#include <string>
#include <vector>
#include <utility>
#include <map>
#include <iostream>


class Protocol {
public:

  void add_protocol(std::string name, std::vector<double> protocol);
  void add_protocol(std::string name, double (*generator)(int), int length);
  std::vector<std::pair<std::string, double>> get_drug_concentrations(int timestep);

  int memory_usage();

  void set_name(std::string n) { name = n; }
  std::string get_name() { return name; }

  friend std::ostream &operator<<(std::ostream &out, const Protocol &p);

private:
  std::map<std::string, std::vector<double>> protocols;
  std::string name;
};


#endif
