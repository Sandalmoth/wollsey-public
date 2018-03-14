#include "record.h"


using namespace std;


const hsize_t INDEX_CHUNKSIZE[1] {16};
const hsize_t RECORD_CHUNKSIZE[2] {16, 256};
const int CELL_INDEX_FILL[1] {0};
const int CELL_RECORD_FILL[1] {0};
const string DRUG_INDEX_FILL[1] {"UNKNOWN"};
const double DRUG_RECORD_FILL[1] {0.0};
const int DEFLATE_LEVEL = 4;


void Record::set_flags(int flags) {
  f_timeline = flags & Recordables::TIMELINE;
  f_wt_half = flags & Recordables::WT_HALF;
}


void Record::snapshot(const std::map<gene::dna, unsigned long long> &cells,
                      const std::vector<std::pair<std::string, double>> &drugs) {
  if (f_timeline) {
    if (timestep % resolution == 0) {
      std::map<gene::id, unsigned long long> cell_ids;
      for (auto &c: cells) {
        gene::id id = gene::id(gene::gene(c.first));
        auto it = cell_ids.find(id);
        if (it == cell_ids.end()) {
          cell_ids.insert(make_pair(id, c.second));
        } else {
          it->second += c.second;
        }
      }
      snapshot(cell_ids, drugs);
      --timestep; // prevent double-incrementing timestep as it is incremented in the other snapshot function
    }
  }

  if (f_wt_half and wt_half < 0) {
    // This could be a one-time calculation to save time
    int n_cells = 0;
    int n_wildtype_cells = 0;
    for (auto &c: cells) {
      gene::id id = gene::id(c.first);
      n_cells += c.second;
      if (id == 0)
        n_wildtype_cells += c.second;
    }
    if (n_wildtype_cells*2 < n_cells) {
      wt_half = timestep;
    }
    // cout << timestep << ' ' << n_cells << ' ' << n_wildtype_cells << endl;
    // auto it = cells.begin();
    // while (it != cells.end()) {
    //   gene::id id = gene::id(it->first);
    //   cout << id << endl;
    //   if (id == 0)
    //     break;
    //   ++it;
    // }
    // auto it = cells.find(0); // 0 is wildtype
    // if (it != cells.end()) {
    //   cout << timestep << ' ' << n_cells << ' ' << it->first << ' ' << it->second << endl;
    //   if (it->second < static_cast<unsigned long long>(n_cells/2)) {
    //     wt_half = timestep;
    //   }
  }
  ++timestep;
}


void Record::snapshot(const std::map<gene::id, unsigned long long> &cells,
                      const std::vector<std::pair<std::string, double>> &drugs) {
  if (f_timeline) {
    if (timestep % resolution == 0) {
      // record cell situation
      for (auto &c: cells) {
        auto it = timeline.find(c.first);
        if (it == timeline.end()) {
          timeline.insert(make_pair(c.first, vector<pair<int, int>>{make_pair(timestep, c.second)}));
        } else {
          it->second.push_back(make_pair(timestep, c.second));
        }
      }
      // record drug situation
      for (auto &d: drugs) {
        auto it = protocol.find(d.first);
        if (it == protocol.end()) {
          protocol.insert(make_pair(d.first, vector<pair<int, double>>{make_pair(timestep, d.second)}));
        } else {
          it->second.push_back(make_pair(timestep, d.second));
        }
      }
    }
  }

  if (f_wt_half) {
    // This could be a one-time calculation to save time
    int n_cells = 0;
    for (auto &c: cells) {
      n_cells += c.second;
    }
    auto it = cells.find(0); // 0 is wildtype
    if (it != cells.end()) {
      if (it->second < static_cast<unsigned long long>(n_cells/2) and wt_half < 0) {
        wt_half = timestep;
      }
    }
  }
  ++timestep;
}


void Record::store(H5::Group &gp) {
  if (f_timeline) {
    H5::Group gp_cells = gp.createGroup("Cells");
    H5::Group gp_drugs = gp.createGroup("Drugs");
    // ### CELLS ### //
    // TODO consider merging setup and writing
    // Set up the data sets
    {
      H5::DataSet ds_index;
      {
        hsize_t dims[1] {timeline.size()};
        H5::DataSpace sp(1, dims);
        H5::DSetCreatPropList pl;
        pl.setFillValue(H5::PredType::NATIVE_INT, CELL_INDEX_FILL);
        pl.setDeflate(DEFLATE_LEVEL);
        hsize_t chunk_dims[1] {INDEX_CHUNKSIZE[0]};
        if (INDEX_CHUNKSIZE[0] > dims[0])
          chunk_dims[0] = dims[0];
        pl.setChunk(1, chunk_dims);
        ds_index = gp_cells.createDataSet("Index", H5::PredType::NATIVE_INT, sp, pl);
      }
      H5::DataSet ds_record;
      {
        hsize_t dims[2] {timeline.size(), static_cast<hsize_t>(timestep / resolution) + 1};
        H5::DataSpace sp(2, dims);
        H5::DSetCreatPropList pl;
        pl.setFillTime(H5D_FILL_TIME_ALLOC);
        pl.setFillValue(H5::PredType::NATIVE_INT, CELL_RECORD_FILL);
        pl.setDeflate(DEFLATE_LEVEL);
        hsize_t chunk_dims[2] {RECORD_CHUNKSIZE[0], RECORD_CHUNKSIZE[1]};
        if (RECORD_CHUNKSIZE[0] > dims[0])
          chunk_dims[0] = dims[0];
        if (RECORD_CHUNKSIZE[1] > dims[1])
          chunk_dims[1] = dims[1];
        pl.setChunk(2, chunk_dims);
        ds_record = gp_cells.createDataSet("Record", H5::PredType::NATIVE_INT, sp, pl);
      }
      // Write data
      H5::DataSpace sp_index = ds_index.getSpace();
      H5::DataSpace sp_record = ds_record.getSpace();
      size_t row = 0;
      for (auto &c: timeline) {
        // Write id to index
        {
          hsize_t start[1] {row};
          hsize_t count[1] {1};
          sp_index.selectHyperslab(H5S_SELECT_SET, count, start);
          H5::DataSpace sp_mem(1, count);
          ds_index.write(&c.first, H5::PredType::NATIVE_INT, sp_mem, sp_index);
        }
        // Parse record vector into normal timeline vector
        // Memory should be fine since we're doing one at a time.
        // Beware that simulations longer than ~1e9 steps could mean trouble
        // as RAM usage could go into several gigabytes
        // (assuming they're run with resolution 1).
        {
          vector<int> trace;
          int i = 0;
          for (auto t: c.second) {
            while (i < t.first) {
              trace.push_back(0);
              i += resolution;
            }
            trace.push_back(t.second);
            i += resolution;
          }
          hsize_t start[2] {row, 0};
          hsize_t count[2] {1, trace.size()};
          sp_record.selectHyperslab(H5S_SELECT_SET, count, start);
          hsize_t dims_mem[2] {1, trace.size()};
          H5::DataSpace sp_mem(2, dims_mem);
          ds_record.write(trace.data(), H5::PredType::NATIVE_INT, sp_mem, sp_record);
        }
        ++row;
      }
    }
    // ### DRUGS ### //
    // We can run without drugs (but not without cells) so this section needs a condition
    if (protocol.size() > 0) {
      H5::DataSet ds_index;
      H5::StrType stringtype(0, 20);
      {
        hsize_t dims[1] {protocol.size()};
        H5::DataSpace sp(1, dims);
        H5::DSetCreatPropList pl;
        pl.setFillValue(stringtype, DRUG_INDEX_FILL);
        pl.setDeflate(DEFLATE_LEVEL);
        hsize_t chunk_dims[1] {INDEX_CHUNKSIZE[0]};
        if (INDEX_CHUNKSIZE[0] > dims[0])
          chunk_dims[0] = dims[0];
        pl.setChunk(1, chunk_dims);
        ds_index = gp_drugs.createDataSet("Index", stringtype, sp, pl);
      }
      H5::DataSet ds_record;
      {
        hsize_t dims[2] {protocol.size(), static_cast<hsize_t>(timestep / resolution) + 1};
        H5::DataSpace sp(2, dims);
        H5::DSetCreatPropList pl;
        pl.setFillTime(H5D_FILL_TIME_ALLOC);
        pl.setFillValue(H5::PredType::NATIVE_DOUBLE, DRUG_RECORD_FILL);
        pl.setDeflate(DEFLATE_LEVEL);
        hsize_t chunk_dims[2] {RECORD_CHUNKSIZE[0], RECORD_CHUNKSIZE[1]};
        if (RECORD_CHUNKSIZE[0] > dims[0])
          chunk_dims[0] = dims[0];
        if (RECORD_CHUNKSIZE[1] > dims[1])
          chunk_dims[1] = dims[1];
        pl.setChunk(2, chunk_dims);
        ds_record = gp_drugs.createDataSet("Record", H5::PredType::NATIVE_DOUBLE, sp, pl);
      }
      H5::DataSpace sp_index = ds_index.getSpace();
      H5::DataSpace sp_record = ds_record.getSpace();
      size_t row = 0;
      for (auto &d: protocol) {
        {
          hsize_t start[1] {row};
          hsize_t count[1] {1};
          sp_index.selectHyperslab(H5S_SELECT_SET, count, start);
          H5::DataSpace sp_mem(1, count);
          ds_index.write(d.first.data(), stringtype, sp_mem, sp_index);
        }
        {
          vector<double> trace;
          int i = 0;
          for (auto t: d.second) {
            while (i < t.first) {
              trace.push_back(0);
              i += resolution;
            }
            trace.push_back(t.second);
            i += resolution;
          }
          hsize_t start[2] {row, 0};
          hsize_t count[2] {1, trace.size()};
          sp_record.selectHyperslab(H5S_SELECT_SET, count, start);
          hsize_t dims_mem[2] {1, trace.size()};
          H5::DataSpace sp_mem(2, dims_mem);
          ds_record.write(trace.data(), H5::PredType::NATIVE_DOUBLE, sp_mem, sp_record);
        }
        ++row;
      }
    }
  } // endif (f_timeline)

  if (f_wt_half) {
    hsize_t dims[1] {1};
    H5::DataSpace sp(1, dims);
    H5::Attribute at = gp.createAttribute("wt_half", H5::PredType::NATIVE_INT, sp);
    at.write(H5::PredType::NATIVE_INT, &wt_half);
  }
}
