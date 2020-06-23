#include <iostream>
#include <tclap/CmdLine.h>
#include <H5Cpp.h>
#include <memory>
#include <thread>
#include <list>
#include <chrono>

#include "gene.h"
#include "medic.h"
#include "drug.h"
#include "protocol.h"
#include "extrah5.h"
#include "simulator.h"
#include "branching_process.h"
#include "limited_branching_process.h"
#include "constant_branching_process.h"
#include "phenotypist.h"


using namespace std;


const string VERSION = "1.01";


enum class Mode {
  SINGLE_RUN,
  GENETIC_ALGORITHM
};


enum class Sim {
  BRANCHING_PROCESS,
  LIMITED_BRANCHING_PROCESS,
  CONSTANT_BRANCHING_PROCESS
};


struct Arguments {
  string input_file;
  string run_file;
  string output_file;
  Mode mode;
  Sim sim;
  int verbosity;
  int max_threads;
  int record_resolution;
  int seed;
};


struct Parameters {
  int repeats;
};


int main(int argc, char **argv) {

  // ############################## //
  // ### ARGUMENT PARSING BLOCK ### //
  // ############################## //

  Arguments a;

  try {
    TCLAP::CmdLine cmd("General treatment simulator", ' ', VERSION);
    // input/output
		TCLAP::ValueArg<string> a_input_file("i", "input-file", "Path to input file", true, "", "*.hdf5", cmd);
		TCLAP::ValueArg<string> a_run_file("r", "run-file", "Path to run file", true, "", "*.hdf5", cmd);
		TCLAP::ValueArg<string> a_output_file("o", "output-file", "Path to output file", true, "", "*.hdf5", cmd); // TODO make non-obligatory?
    // Major modes
    vector<string> modestrings = {"single", "genetic"};
    TCLAP::ValuesConstraint<string> modes(modestrings);
		TCLAP::ValueArg<string> a_mode("m", "mode", "Program running mode", false, "single", &modes, cmd);
    vector<string> simstrings = {"branching_process"
                                 , "limited_branching_process"
                                 , "constant_branching_process"};
    TCLAP::ValuesConstraint<string> sims(simstrings);
		TCLAP::ValueArg<string> a_sim("e", "simulator", "Simulation engine", false, "branching_process", &sims, cmd);
    // Minor modes
    TCLAP::MultiSwitchArg a_verbosity("v", "verbose", "Increase verbosity of output", cmd);
    TCLAP::ValueArg<int> a_max_threads("j", "max-threads", "Max number of threads for simulations (+1 for main)", false, 1, "integer", cmd);
    // TODO add more recording options
    TCLAP::ValueArg<int> a_record_resolution("s", "recording-resolution", "How often the simulation state will be recorded (i.e. every x timesteps)", false, 100, "integer", cmd);
    TCLAP::ValueArg<int> a_seed("d", "seed", "RNG-seed", false, 2071, "integer", cmd);

    cmd.parse(argc, argv);

    a.input_file = a_input_file.getValue();
    a.run_file = a_run_file.getValue();
    a.output_file = a_output_file.getValue();
    a.verbosity = a_verbosity.getValue();
    a.max_threads = a_max_threads.getValue();
    a.record_resolution = a_record_resolution.getValue();
    a.seed = a_seed.getValue();

    string modestring = a_mode.getValue();
    if (modestring == "single") {
      a.mode = Mode::SINGLE_RUN;
    } else if (modestring == "genetic") {
      a.mode = Mode::GENETIC_ALGORITHM;
    } else {
      // shouldn't happen
      cerr << "Unknown mode " << modestring << endl;
      return 1;
    }
    string simstring = a_sim.getValue();
    if (simstring == "branching_process") {
      a.sim = Sim::BRANCHING_PROCESS;
    } else if (simstring == "limited_branching_process") {
      a.sim = Sim::LIMITED_BRANCHING_PROCESS;
    } else if (simstring == "constant_branching_process") {
      a.sim = Sim::CONSTANT_BRANCHING_PROCESS;
    } else {
      cerr << "Unknown simulator " << simstring << endl;
      return 1;
    }

  } catch (TCLAP::ArgException &e) {
    cerr << "TCLAP Error: " << e.error() << endl << "\targ: " << e.argId() << endl;
    return 1;
  }


  // ########################## //
  // ### PLAN PARSING BLOCK ### //
  // ########################## //

  Parameters p;
  vector<Protocol> protocols;

  try {
    H5::H5File runfile(a.run_file, H5F_ACC_RDONLY);

    // read parameters
    cout << "Reading parameters" << endl;
    H5::DataSet ds_parameters = runfile.openDataSet("Parameters");
    ds_parameters.openAttribute("repeats").read(H5::PredType::NATIVE_INT, &p.repeats);

    // read protocol specifications
    cout << "Reading protocols" << endl;
    H5::Group gp_protocols = runfile.openGroup("Protocol");
    int n_protocols = gp_protocols.getNumObjs();
    for (int i = 0; i < n_protocols; ++i) {
      string protocolname = gp_protocols.getObjnameByIdx(i);
      cout << "  Reading protocol " << protocolname << endl;
      H5::Group gp_protocol = gp_protocols.openGroup(protocolname);
      H5::DataSet ds_index = gp_protocol.openDataSet("Index");
      H5::DataSet ds_record = gp_protocol.openDataSet("Record");
      int length;
      gp_protocol.openAttribute("length").read(H5::PredType::NATIVE_INT, &length);
      int resolution;
      // TODO implement actually caring about resolution?
      gp_protocol.openAttribute("resolution").read(H5::PredType::NATIVE_INT, &resolution);

      Protocol tprot;
      tprot.set_name(protocolname);
      // Iterate over the drugs
      vector<string> drugnames = X5::read_vector<string>(ds_index, 20);
      for (size_t j = 0; j < drugnames.size(); ++j) {
        string d = drugnames[j];
        cout << "    Reading drug " << d << endl;
        hsize_t start[2] {j, 0};
        vector<double> trace = X5::read_column<double, 2>(ds_record, start, length, 1, H5::PredType::NATIVE_DOUBLE);
        tprot.add_protocol(d, trace);
      }
      protocols.push_back(move(tprot));
    }

  } catch (H5::Exception &e) {
    e.printErrorStack();
    return 1;
  }


  // ################################ //
  // ### INPUT FILE PARSING BLOCK ### //
  // ################################ //

  shared_ptr<Simulator> simulator;

  try {
    H5::H5File infile(a.input_file, H5F_ACC_RDONLY);

    // Read radix
    cout << "Reading radix" << endl;
    H5::DataSet ds_radix = infile.openDataSet("Gene/Radix");
    vector<int> radix = X5::read_vector<int>(ds_radix, H5::PredType::NATIVE_INT);
    // gene::set_radix(radix); // deprecated in favour of deriving from protein options
    // Read protein options (ragged matrix)
    cout << "Reading protein options" << endl;
    H5::DataSet ds_options = infile.openDataSet("Gene/Options");
    vector<vector<char>> protein_options;
    for (size_t i = 0; i < radix.size(); ++i) {
      hsize_t start[2] = {i, 0};
      protein_options.push_back(X5::read_column<char, 2>(ds_options, start, radix[i], 1, 20));
    }
    gene::set_protein_options(protein_options);
    gene::generate_translation_table();
    gene::set_mutation_matrix(gene::MutationMatrix::K80); // TODO add option

    Medic medic; // We need this now to have somewhere to put the drugs

    // Read drugs
    cout << "Reading drugs" << endl;
    H5::Group gp_drugs = infile.openGroup("Drugs");
    int n_drugs = gp_drugs.getNumObjs();
    for (int i = 0; i < n_drugs; ++i) {
      string drugname = gp_drugs.getObjnameByIdx(i);
      cout << "  Reading drug " << drugname << endl;
      H5::Group gp_drug = gp_drugs.openGroup(drugname);
      H5::DataSet ds_single = gp_drug.openDataSet("Single");
      // First read the single mutation (ragged) matrix
      vector<vector<double>> snv_ic50s;
      for (size_t j = 0; j < radix.size(); ++j) {
        hsize_t start[2] = {j, 0};
        snv_ic50s.push_back(
          X5::read_column<double, 2>(ds_single, start, radix[j], 1, H5::PredType::NATIVE_DOUBLE)
        );
      }
      Drug drug(snv_ic50s);
      // Now read unique ic50s
      vector<gene::gene> unique_genes;
      vector<double> unique_ic50s;
      H5::DataSet ds_unique_ic50s = gp_drug.openDataSet("Unique_IC50s");
      unique_ic50s = X5::read_vector<double>(ds_unique_ic50s, H5::PredType::NATIVE_DOUBLE);
      H5::DataSet ds_unique_genes = gp_drug.openDataSet("Unique_genes");
      for (size_t j = 0; j < unique_ic50s.size(); ++j) {
        hsize_t start[2] = {j, 0};
        hsize_t tdms[2];
        ds_unique_genes.getSpace().getSimpleExtentDims(tdms);
        unique_genes.push_back(gene::gene(X5::read_column<int, 2>(ds_unique_genes, start, radix.size(), 1, H5::PredType::NATIVE_INT)));
      }
      // Add to the drug
      for (size_t j = 0; j < unique_genes.size(); ++j) {
        gene::id id = gene::id(unique_genes[j]);
        if (id == 0) {
          drug.set_wildtype_ic50(unique_ic50s[j]);
        } else {
          drug.set_ic50(id, unique_ic50s[j]);
        }
      }
      medic.add_drug(drugname, drug);
    }
    if (a.verbosity) cout << medic << endl;

    // Read genotype properties
    cout << "Reading genotype properties" << endl;
    H5::Group gp_rates = infile.openGroup("Rates");
    // First, figure out whether to use unified deathrate/reproductionrate
    double unified_death_rate = -1.0;
    double unified_reproduction_rate = -1.0;
    double unified_minimum_reproduction_rate = -1.0;
    // if (gp_rates.attrExists("unified_reproduction_rate"))
    //   gp_rates.openAttribute("unified_reproduction_rate").read(H5::PredType::NATIVE_DOUBLE, &unified_reproduction_rate);
    // if (gp_rates.attrExists("unified_death_rate"))
    //   gp_rates.openAttribute("unified_death_rate").read(H5::PredType::NATIVE_DOUBLE, &unified_death_rate);
    // workaround since attrExists is unsupported in older versions of hdf5
    try {
        gp_rates.openAttribute("unified_reproduction_rate").read(H5::PredType::NATIVE_DOUBLE, &unified_reproduction_rate);
    } catch (...) { }
    try {
        gp_rates.openAttribute("unified_death_rate").read(H5::PredType::NATIVE_DOUBLE, &unified_death_rate);
    } catch (...) { }
    try {
      gp_rates.openAttribute("unified_minimum_reproduction_rate").read(H5::PredType::NATIVE_DOUBLE, &unified_minimum_reproduction_rate);
    } catch (...) { }

    Phenotypist phenotypist(unified_death_rate >= 0.0, unified_death_rate
                            , unified_reproduction_rate >= 0.0, unified_reproduction_rate
                            , unified_minimum_reproduction_rate >= 0.0, unified_minimum_reproduction_rate);
    vector<gene::gene> genes;
    vector<double> reproduction_rates;
    vector<double> death_rates;
    vector<double> minimum_reproduction_rates;
    // Now read individual rates (if relevant)
    if (unified_reproduction_rate < 0.0) {
      H5::DataSet ds_reproduction = gp_rates.openDataSet("Reproduction");
      reproduction_rates = X5::read_vector<double>(ds_reproduction, H5::PredType::NATIVE_DOUBLE);
      H5::DataSet ds_genes = gp_rates.openDataSet("Genes");
      for (size_t j = 0; j < reproduction_rates.size(); ++j) {
        hsize_t start[2] = {j, 0};
        genes.push_back(gene::gene(X5::read_column<int, 2>(ds_genes, start, radix.size(), 1, H5::PredType::NATIVE_INT)));
      }
    }
    for (size_t i = 0; i < genes.size(); ++i) {
      Phenotypist::Phenotype ph {i < reproduction_rates.size() ? reproduction_rates[i] : 0, i < death_rates.size() ? death_rates[i] : 0, i < minimum_reproduction_rates.size() ? minimum_reproduction_rates[i] : 0};
      phenotypist.add_phenotype(gene::id(genes[i]), ph);
    }
    if (a.verbosity) cout << phenotypist << endl;

    // Read starting cells
    cout << "Reading starting cells" << endl;
    H5::Group gp_cells = infile.openGroup("Cells");
    H5::DataSet ds_genomes = gp_cells.openDataSet("Genomes");
    H5::DataSet ds_counts = gp_cells.openDataSet("Counts");
    vector<string> cell_genome_strings = X5::read_vector<string>(ds_genomes, 3*radix.size() + 1);
    vector<int> cell_counts = X5::read_vector<int>(ds_counts, H5::PredType::NATIVE_INT);
    vector<gene::dna> cell_genomes;
    for (auto gs: cell_genome_strings) {
      cell_genomes.push_back(gene::dna(gs));
      cout << gs << endl;
      auto td1 = gene::dna(gs);
      auto tg1 = gene::gene(td1);
      auto ti1 = gene::id(td1);
      auto ti2 = gene::id(tg1);
      auto tg2 = gene::gene(ti1);
      auto tg3 = gene::gene(ti2);
      cout << td1 << endl;
      cout << tg1 << endl;
      cout << ti1 << endl;
      cout << ti1 << endl;
      cout << tg1 << endl;
      cout << tg1 << endl;
    }

    // Read simulation parameters
    cout << "Reading simulation parameters" << endl;
    H5::DataSet ds_parameters = infile.openDataSet("Parameters");
    int max_iters;
    double mutation_probability;
    int max_cells;
    int min_cells;
    int population_size;
    double spring_force;
    ds_parameters.openAttribute("max_iters").read(H5::PredType::NATIVE_INT, &max_iters);
    ds_parameters.openAttribute("mutation_probability").read(H5::PredType::NATIVE_DOUBLE, &mutation_probability);
    ds_parameters.openAttribute("min_cells").read(H5::PredType::NATIVE_INT, &min_cells);
    ds_parameters.openAttribute("max_cells").read(H5::PredType::NATIVE_INT, &max_cells);
    ds_parameters.openAttribute("population_size").read(H5::PredType::NATIVE_INT, &population_size);
    ds_parameters.openAttribute("spring_force").read(H5::PredType::NATIVE_DOUBLE, &spring_force);
    if (a.verbosity > 1) {
      cout << "max_iters\t" << max_iters << endl;
      cout << "mutation_probability\t" << mutation_probability << endl;
      cout << "max_cells\t" << max_cells << endl;
      cout << "min_cells\t" << min_cells << endl;
      cout << "population_size\t" << population_size << endl;
      cout << "spring_force\t" << spring_force << endl;
    }

    // Set up the simulator
    switch (a.sim) {
    case Sim::BRANCHING_PROCESS:
      simulator = make_shared<BranchingProcess>(max_iters
                                                , mutation_probability
                                                , medic
                                                , phenotypist
                                                , Recordables::TIMELINE | Recordables::WT_HALF
                                                , a.record_resolution
                                                );
      // And also someshing like
      // reinterpret_cast<BranchingProcess>(simulator)->add_cells
      // dynamic_pointer_cast<BranchingProcess>(simulator)->add_cells(gene::dna("ATGCTGGGGCAGTACGAGGACGAGGTGATTTTCATGATGTTCCTGCATTTT"), 100);  // T315I
      // dynamic_pointer_cast<BranchingProcess>(simulator)->add_cells(gene::dna("ATGCTGGGGCAGTACGAGGACGAGGTGACTTTCATGATGTTCCTGCATTTT"), 100); // WT
      for (size_t i = 0; i < cell_genomes.size(); ++i) {
        dynamic_pointer_cast<BranchingProcess>(simulator)->add_cells(cell_genomes[i], cell_counts[i]);
      }
      break;
    case Sim::LIMITED_BRANCHING_PROCESS:
      simulator = make_shared<LimitedBranchingProcess>(max_iters
                                                       , mutation_probability
                                                       , medic
                                                       , phenotypist
                                                       , min_cells
                                                       , max_cells
                                                       , spring_force
                                                       , Recordables::TIMELINE | Recordables::WT_HALF
                                                       , a.record_resolution
                                                       );
      // dynamic_pointer_cast<LimitedBranchingProcess>(simulator)->add_cells(gene::dna("ATGCTGGGGCAGTACGAGGACGAGGTGACTTTCATGATGTTCCTGCATTTT"), 100); // WT
      for (size_t i = 0; i < cell_genomes.size(); ++i) {
        dynamic_pointer_cast<LimitedBranchingProcess>(simulator)->add_cells(cell_genomes[i], cell_counts[i]);
      }
      break;
      break;
    case Sim::CONSTANT_BRANCHING_PROCESS:
      simulator = make_shared<ConstantBranchingProcess>(max_iters
                                                        , mutation_probability
                                                        , medic
                                                        , phenotypist
                                                        , population_size
                                                        , spring_force
                                                        , Recordables::TIMELINE | Recordables::WT_HALF
                                                        , a.record_resolution
                                                        );
      // dynamic_pointer_cast<ConstantBranchingProcess>(simulator)->add_cells(gene::dna("ATGCTGGGGCAGTACGAGGACGAGGTGACTTTCATGATGTTCCTGCATTTT"), 100); // WT
      for (size_t i = 0; i < cell_genomes.size(); ++i) {
        dynamic_pointer_cast<ConstantBranchingProcess>(simulator)->add_cells(cell_genomes[i], cell_counts[i]);
      }
      break;
      break;
    default:
      cerr << "Unknown (or unimplemented) simulator" << endl;
      return 1;
    }

  } catch (H5::Exception &e) {
    e.printErrorStack();
    return 1;
  }

  cout << "Loaded simulation system" << endl;


  // ################# //
  // ### RUN BLOCK ### //
  // ################# //

  cout << "\nRunning simulation(s)" << endl;
  auto total_timer = [start = chrono::system_clock::now()] { return chrono::system_clock::now() - start; };

  try {
    H5::H5File outfile(a.output_file, H5F_ACC_TRUNC);
    H5::Group gp_record = outfile.createGroup("Record");
    X5::make_attr(gp_record, "n_repeats", p.repeats, H5::PredType::NATIVE_INT);

    switch (a.mode) {
    case Mode::GENETIC_ALGORITHM:
      cerr << "Genetic algorithms unimplemented as of now" << endl;
      // TODO yea
      break;

    case Mode::SINGLE_RUN:
      { // scope the groups
        if (a.verbosity) cout << "Single run mode" << endl;
        H5::Group gp_this_generation(gp_record.createGroup("0"));
        // Run every protocol p.repeats times.
        for (size_t j = 0; j < protocols.size(); ++j) {
          if (a.verbosity) cout << "  Running protocol " << j+1 << "/" << protocols.size() << '\t' << protocols[j].get_name() << endl;
          list<thread> threads;
          vector<shared_ptr<Simulator>> sims;
          sims.reserve(p.repeats); // check if this is neccessary with ptrs (as simulators shouldnt move when vector expands)
          int i = 0;
          for (int k = 0; k < a.max_threads; ++k) {
            if (k >= p.repeats)
              break;
            sims.emplace_back(simulator->get_instance());
            threads.emplace_front(sims.back()->run(protocols[j]));
            ++i;
          }
          while (int(i) < p.repeats) {
            threads.back().join();
            threads.pop_back();
            sims.emplace_back(simulator->get_instance());
            threads.emplace_front(sims.back()->run(protocols[j]));
            if (a.verbosity > 1) {
              cout << "\r    Simulation repeat " << i+1 << "/" << p.repeats;
              cout.flush();
            }
            ++i;
          }
          // Join from oldest (back) to newest (front).
          for (auto rit = threads.rbegin(); rit != threads.rend(); ++rit) {
            rit->join();
          }
          threads.clear();
          // Consider merging writing to file inside loop to save memory.
          H5::Group gp_this_protocol(gp_this_generation.createGroup(to_string(j)));
          for (int w = 0; w < p.repeats; ++w) {
            H5::Group gp_this_record(gp_this_protocol.createGroup(to_string(w)));
            sims[w]->get_record().store(gp_this_record);
          }
          if (a.verbosity > 1) cout << endl;
        }
      }
      // TODO turn into a general function?
      break;

    default:
      cerr << "Unknown or unimplemented program mode" << endl;
    }

    // Copy input file into output file
    H5::H5File infile(a.input_file, H5F_ACC_RDONLY);
    X5::rcopy(infile, outfile);
    // Also copy runfile
    // verbose as rcopy cant be used as it wont merge attributes of the Parameters dataset
    H5::H5File runfile(a.run_file, H5F_ACC_RDONLY);
    X5::copy(runfile, outfile, "Protocol");
    H5::DataSet ds_opar = outfile.openDataSet("Parameters");
    H5::DataSet ds_rpar = runfile.openDataSet("Parameters");
    X5::copy_attrs(ds_rpar, ds_opar);

  } catch (H5::Exception &e) {
    e.printErrorStack();
    return 1;
  }

  cout << "Total simulation time: " << chrono::duration_cast<chrono::seconds>(total_timer()).count() << " s" << endl;

}
