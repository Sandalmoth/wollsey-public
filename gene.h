#ifndef __GENE_H__
#define __GENE_H__


#include <vector>
#include <cstdint>
#include <string>
#include <iostream>
#include <random>


namespace gene {

  // TODO maybe refactor such that comparator operators belong to a template.

  class dna;
  class gene;
  class id;

  enum class MutationMatrix {
    JC69,
    K80
  };

  void set_protein_options(std::vector<std::vector<char>> opt);
  void generate_translation_table();
  void set_mutation_matrix(MutationMatrix mm);

  class dna {
  public:
    dna(std::string s);
    operator gene() const;
    operator id() const;
    bool operator<(const dna &other) const { return x < other.x; }
    bool operator==(const dna &other) const { return x==other.x; }
    friend std::ostream &operator<<(std::ostream &out, const dna &d);
    void mutate_snv(std::mt19937 &rng);
    const std::vector<unsigned char> &get_x() const { return x; }
  private:
    std::vector<unsigned char> x;
  };

  class gene {
  public:
    gene(std::vector<unsigned char> x)
      : x(x) { }
    gene(std::vector<int> z);
    bool operator<(const gene &other) const { return x < other.x; }
    bool operator==(const gene &other) const { return x==other.x; }
    operator id() const;
    friend double difference(const gene &g1, const gene &g2);
    friend std::ostream &operator<<(std::ostream &out, const gene &g);
    const std::vector<unsigned char> &get_x() const { return x; }
  private:
    std::vector<unsigned char> x;
  };

  class id {
  public:
    id(uint64_t x)
      : x(x) { }
    operator gene() const;
    bool operator<(const id &other) const {return x < other.x; }
    bool operator==(const id &other) const {return x == other.x; }
    friend std::ostream &operator<<(std::ostream &out, const id &i);
    const uint64_t &get_x() const { return x; }
  private:
    uint64_t x;
  };

  bool is_allowed_option(int x, unsigned char codon);

  double difference(const gene &g1, const gene&g2);

  std::ostream &operator<<(std::ostream &out, const dna &d);
  std::ostream &operator<<(std::ostream &out, const gene &g);
  std::ostream &operator<<(std::ostream &out, const id &i);


} // end namespace gene


#endif
