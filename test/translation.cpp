#include <iostream>
#include <cassert>
#include <chrono>
#include <random>

#include "../gene.h"


using namespace std;


int main() {

  // TODO add in assertions to make this into a proper test.

  gene::set_protein_options(vector<vector<char>>{
        vector<char>{'I', 'L', 'D', 'V'}
      , vector<char>{'Q', 'T', 'D', 'S', 'A'}
      , vector<char>{'I', 'C', 'F'}
    });
  gene::generate_translation_table();
  gene::set_mutation_matrix(gene::MutationMatrix::K80);

  gene::dna d("GATACCATT");
  gene::gene g = d;
  gene::id i = g;
  gene::gene g2 = i;
  cout << d << endl;
  cout << g << endl;
  cout << i << endl;
  cout << g2 << endl;

  cout << endl;

  gene::dna d2("ATTCAAATT");
  gene::gene g3 = d2;
  gene::id i2 = g3;
  gene::gene g4 = i2;
  cout << d2 << endl;
  cout << g3 << endl;
  cout << i2 << endl;
  cout << g4 << endl;

  cout << endl;

  random_device rd;
  mt19937 rng(rd());

  for (int i = 0; i < 10; ++i) {
    d.mutate_snv(rng);
    cout << d << '\t' << gene::gene(d) << '\t' << gene::id(d) << endl;
  }

  auto t1 = [val = chrono::system_clock::now()] { return chrono::system_clock::now() - val; };
  cout << chrono::duration_cast<chrono::milliseconds>(t1()).count() << " ms" << endl;;
}
