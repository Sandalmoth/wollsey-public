#include <iostream>
#include <cassert>
#include <chrono>

#include "../branching_process.h"
#include "H5Cpp.h"


using namespace std;


int main() {
  // Set up genetics
  gene::set_protein_options(vector<vector<char>>{
      vector<char>{'I', 'L', 'D', 'V'}
      , vector<char>{'Q', 'T', 'D', 'S', 'A'}
      , vector<char>{'I', 'C', 'F'}
    });
  gene::generate_translation_table();
  gene::set_mutation_matrix(gene::MutationMatrix::K80);

  // bp testing
  gene::dna d("GATACCATT");
  cout << "starting gene id " << gene::id(d) << endl;
  Medic m;
  Phenotypist ph (true, 0.09);
  ph.add_phenotype(gene::id(d), Phenotypist::Phenotype{0.1, 0.0});
  ph.add_phenotype(5, Phenotypist::Phenotype{0.2, 0.0});
  ph.add_phenotype(8, Phenotypist::Phenotype{0.3, 0.0});
  BranchingProcess bp(100, 0.1, m, ph, Recordables::TIMELINE | Recordables::WT_HALF, 1);
  bp.add_cells(d, 100);
  Protocol p;
  auto t1 = [val = chrono::system_clock::now()] { return chrono::system_clock::now() - val; };
  bp.run(p).join();
  cout << chrono::duration_cast<chrono::milliseconds>(t1()).count() << " ms" << endl;;
  Record &r = bp.get_record();
  try {
    H5::H5File outfile("out.hdf5", H5F_ACC_TRUNC);
    // medmo expo.py compatibility code
    H5::Group gp_record = outfile.createGroup("Record");
    H5::Group gp_generation = gp_record.createGroup("0");
    H5::Group gp_protocol = gp_generation.createGroup("0");
    H5::Group gp_this_record = gp_protocol.createGroup("0");
    r.store(gp_this_record);
  } catch (H5::Exception &e) {
    e.printError();
    return 1;
  }
}
