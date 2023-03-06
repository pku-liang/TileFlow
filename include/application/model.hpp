
#pragma once

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/bitset.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include "mapping/parser.hpp"
#include "mapping/arch-properties.hpp"
#include "mapping/constraints.hpp"
#include "compound-config/compound-config.hpp"
#include "model/sparse-optimization-parser.hpp"

//--------------------------------------------//
//                Application                 //
//--------------------------------------------//

class Application
{
 public:
  std::string name_;

  struct Stats
  {
    double energy;
    double cycles;
  };

 protected:
  // Critical state.
  std::vector<problem::Workload> workloads_;
  model::Engine::Specs arch_specs_;

  // Many of the following submodules are dynamic objects because
  // we can only instantiate them after certain config files have
  // been parsed.

  // The mapping.
  Mapping* mapping_;

  // Abstract representation of the architecture.
  ArchProperties* arch_props_;

  // Constraints.
  mapping::Constraints* constraints_;  
  
  // Application flags/config.
  bool verbose_ = false;
  bool auto_bypass_on_failure_ = false;
  std::string out_prefix_;

  // Sparse optimization
  sparse::SparseOptimizationInfo* sparse_optimizations_;

 private:

  // Serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version = 0);

 public:

  Application(config::CompoundConfig* config,
              std::string output_dir = ".",
              std::string name = "timeloop-model");

  // This class does not support being copied
  Application(const Application&) = delete;
  Application& operator=(const Application&) = delete;

  ~Application();

  // Run the evaluation.
  Stats Run();
};

