#include "gene.h"

#include <algorithm>


using namespace std;


static vector<char> translation_table;
static vector<vector<char>> protein_options;
static vector<int> radix;
static vector<vector<double>> mutation_matrix;


void gene::set_protein_options(vector<vector<char>> opt) {
  protein_options = opt;
  for (auto x: protein_options) {
    radix.push_back(x.size());
  }
}


void gene::generate_translation_table() {
  // This hardcode is sort of annoying, but probably the easiest way
  // Simply use the codon as the index. should be fast
  // A = 00
  // C = 01
  // G = 10
  // T = 11

  translation_table.reserve(64);

  translation_table.push_back('K'); // AAA
  translation_table.push_back('N'); // AAC
  translation_table.push_back('K'); // AAG
  translation_table.push_back('N'); // AAT
  translation_table.push_back('T'); // ACA
  translation_table.push_back('T'); // ACC
  translation_table.push_back('T'); // ACG
  translation_table.push_back('T'); // ACT
  translation_table.push_back('R'); // AGA
  translation_table.push_back('S'); // AGC
  translation_table.push_back('R'); // AGG
  translation_table.push_back('S'); // AGT
  translation_table.push_back('I'); // ATA
  translation_table.push_back('I'); // ATC
  translation_table.push_back('M'); // ATG
  translation_table.push_back('I'); // ATT

  translation_table.push_back('Q'); // CAA
  translation_table.push_back('H'); // CAC
  translation_table.push_back('Q'); // CAG
  translation_table.push_back('H'); // CAT
  translation_table.push_back('P'); // CCA
  translation_table.push_back('P'); // CCC
  translation_table.push_back('P'); // CCG
  translation_table.push_back('P'); // CCT
  translation_table.push_back('R'); // CGA
  translation_table.push_back('R'); // CGC
  translation_table.push_back('R'); // CGG
  translation_table.push_back('R'); // CGT
  translation_table.push_back('L'); // CTA
  translation_table.push_back('L'); // CTC
  translation_table.push_back('L'); // CTG
  translation_table.push_back('L'); // CTT

  translation_table.push_back('E'); // GAA
  translation_table.push_back('D'); // GAC
  translation_table.push_back('E'); // GAG
  translation_table.push_back('D'); // GAT
  translation_table.push_back('A'); // GCA
  translation_table.push_back('A'); // GCC
  translation_table.push_back('A'); // GCG
  translation_table.push_back('A'); // GCT
  translation_table.push_back('G'); // GGA
  translation_table.push_back('G'); // GGC
  translation_table.push_back('G'); // GGG
  translation_table.push_back('G'); // GGT
  translation_table.push_back('V'); // GTA
  translation_table.push_back('V'); // GTC
  translation_table.push_back('V'); // GTG
  translation_table.push_back('V'); // GTT

  translation_table.push_back('*'); // TAA
  translation_table.push_back('Y'); // TAC
  translation_table.push_back('*'); // TAG
  translation_table.push_back('Y'); // TAT
  translation_table.push_back('S'); // TCA
  translation_table.push_back('S'); // TCC
  translation_table.push_back('S'); // TCG
  translation_table.push_back('S'); // TCT
  translation_table.push_back('*'); // TGA
  translation_table.push_back('C'); // TGC
  translation_table.push_back('W'); // TGG
  translation_table.push_back('C'); // TGT
  translation_table.push_back('L'); // TTA
  translation_table.push_back('F'); // TTC
  translation_table.push_back('L'); // TTG
  translation_table.push_back('F'); // TTT
}


void gene::set_mutation_matrix(MutationMatrix mm) {
  // Mutation matrix works as mutation_matrix[from][to]
  // Order of nucleotides same as before, i.e. A C G T
  switch (mm) {
  case MutationMatrix::JC69:
    mutation_matrix.push_back(vector<double>{0.0, 1.0, 1.0, 1.0});
    mutation_matrix.push_back(vector<double>{1.0, 0.0, 1.0, 1.0});
    mutation_matrix.push_back(vector<double>{1.0, 1.0, 0.0, 1.0});
    mutation_matrix.push_back(vector<double>{1.0, 1.0, 1.0, 0.0});
    break;
  case MutationMatrix::K80:
    // altered transistion/transversion ratio
    double tsr = 4.0;
    mutation_matrix.push_back(vector<double>{0.0, 1.0, tsr, 1.0});
    mutation_matrix.push_back(vector<double>{1.0, 0.0, 1.0, tsr});
    mutation_matrix.push_back(vector<double>{tsr, 1.0, 0.0, 1.0});
    mutation_matrix.push_back(vector<double>{1.0, tsr, 1.0, 0.0});
    break;
  }
}


gene::dna::dna(string s) {
  for (auto c: s) {
    switch (c) {
    case 'A':
      x.push_back(0);
      break;
    case 'C':
      x.push_back(1);
      break;
    case 'G':
      x.push_back(2);
      break;
    case 'T':
      x.push_back(3);
      break;
    default:
      cerr << "Unkown nucleotide: " << c;
    }
  }
}


gene::dna::operator gene() const {
  // generate list of codons
  std::vector<unsigned char> g;
  for (size_t i = 0; i < x.size(); i += 3) {
    unsigned char c = 0;
    for (int j = 0; j < 3; ++j) {
      c = c << 2;
      c += x[i+j];
    }
    g.push_back(c);
  }
  for (size_t i = 0; i < g.size(); ++i) {
    char aa = translation_table[g[i]];
    auto it = find(protein_options[i].begin(), protein_options[i].end(), aa);
    g[i] = (it - protein_options[i].begin());
  }
  return move(g);
}


gene::dna::operator id() const {
  return id(gene(*this));
}


void gene::dna::mutate_snv(mt19937 &rng) {
  // NOTE consider weighting position choice instead (might be faster?)
  auto posdist = uniform_int_distribution<int>(0, x.size() - 1);
  auto mut = x;
  int pos = 0;
  unsigned char newc;
  int snv;
  do {
    newc = 0;
    pos = posdist(rng);
    vector<unsigned char> n3(x.begin() + pos - pos%3, x.begin() + pos + 3 - pos%3);
    snv = discrete_distribution<int>(mutation_matrix[n3[pos%3]].begin(), mutation_matrix[n3[pos%3]].end())(rng);
    n3[pos%3] = static_cast<unsigned char>(snv);
    for (auto n: n3) {
      newc = newc << 2;
      newc += n;
    }
  } while (!is_allowed_option(pos / 3, newc));
  x[pos] = snv;
}


gene::gene::gene(std::vector<int> z) {
  for (auto y: z) {
    x.push_back(static_cast<unsigned char>(y));
  }
}

gene::gene::operator id() const {
  uint64_t i = 0;
  int k = 1;
  for (unsigned int j = 0; j < radix.size(); ++j) {
    i += x[j] * k;
    k *= radix[j];
  }
  return i;
}


bool gene::is_allowed_option(int x, unsigned char codon) {
  auto &v = protein_options[x];
  return find(v.begin(), v.end(), translation_table[codon]) != v.end();
}


gene::id::operator gene() const {
  uint64_t i = x;
  std::vector<unsigned char> g;
  for (unsigned int j = 0; j < radix.size(); ++j) {
    g.push_back(i % radix[j]);
    i /= radix[j];
  }
  return move(g);
}


// FIXME should be int?
double gene::difference(const gene &g1, const gene&g2) {
  double ndiff = 0.0;
  for (size_t i = 0; i < radix.size(); ++i) {
    if (g1.x[i] == g2.x[i]) {
      ndiff += 1.0;
    }
  }
  return ndiff;
}


ostream &gene::operator<<(ostream &out, const dna &d) {
  for (size_t i = 0; i < d.x.size() - 1; ++i) {
    switch (d.x[i]) {
    case 0:
      cout << "A";
      break;
    case 1:
      cout << "C";
      break;
    case 2:
      cout << "G";
      break;
    case 3:
      cout << "T";
      break;
    default:
      cerr << "Dna string contains non-nucleotides" << endl;
    }
    if (i % 3 == 2)
      cout << ' ';
  }
  switch (d.x.back()) {
  case 0:
    cout << "A";
    break;
  case 1:
    cout << "C";
    break;
  case 2:
    cout << "G";
    break;
  case 3:
    cout << "T";
    break;
  default:
    cerr << "Dna string contains non-nucleotides" << endl;
  }
  return out;
}


ostream &gene::operator<<(ostream &out, const gene &g) {
  for (size_t i = 0; i < g.x.size() - 1; ++i) {
    out << protein_options[i][g.x[i]] << ' ';
  }
  out << protein_options.back()[g.x.back()];
  return out;
}


ostream &gene::operator<<(ostream &out, const id &i) {
  out << i.x;
  return out;
}
