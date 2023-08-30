#include <iostream>
#include <csignal>
#include <cstring>

#include "application/model.hpp"
#include "compound-config/compound-config.hpp"
#include "util/args.hpp"

#include "tileflow/problem/problem.hpp"
#include "tileflow/mapping/mapping.hpp"
#include "tileflow/loop-analysis/nest-analysis.hpp"
#include "tileflow/model/topology.hpp"
#include "tileflow/mapper/checker.hpp"
#include "tileflow/mapper/mapper.hpp"

extern bool gTerminateEval;

//--------------------------------------------//
//                    MAIN                    //
//--------------------------------------------//

void show_energy(
  const model::TileFlow::Topology& topology,
  std::ostream& o = std::cout) {
  auto arith = topology.GetArithmeticLevel();
  o << "==========AccessEnergy===========" << std::endl;
  o << "metric, energy" << std::endl;
  o << "Arith::energy_per_op," 
    << arith->GetSpecs().op_energy_map.at("random_compute") << std::endl;
  for (unsigned i = 0; i < topology.NumStorageLevels(); i++){
    auto buffer = topology.GetStorageLevel(i);
    auto& specs = buffer->GetSpecs();
    o << "Buffer::" << buffer->Name() << "::energy_per_op::read," << specs.op_energy_map.at("random_read") << std::endl;
    o << "Buffer::" << buffer->Name() << "::energy_per_op::update," << specs.op_energy_map.at("random_update") << std::endl;
    o << "Buffer::" << buffer->Name() << "::energy_per_op::fill," << specs.op_energy_map.at("random_fill") << std::endl;
  }
  o << "========End AccessEnergy=========" << std::endl;
}

int main(int argc, char* argv[])
{
  assert(argc >= 2);

  std::vector<std::string> input_files;
  std::string output_dir = ".";
  bool success = ParseArgs(argc, argv, input_files, output_dir);
  if (!success)
  {
    std::cerr << "ERROR: error parsing command line." << std::endl;
    exit(1);
  }

  auto config = new config::CompoundConfig(input_files);

  auto root = config->getRoot();

  if (root.exists("macro")) 
    TileFlow::macros = root.lookup("macro");
  
  if (root.exists("verbose"))
    root.lookupValue("verbose", TileFlow::verbose_level);
  
  auto problem = root.lookup("problem");
  problem::TileFlow::Workloads workloads;

  config::CompoundConfigNode arch;

  if (root.exists("arch"))
  {
    arch = root.lookup("arch");
  }
  else if (root.exists("architecture"))
  {
    arch = root.lookup("architecture");
  }
  
  bool is_sparse_topology = root.exists("sparse_optimizations");

  model::Engine::Specs arch_specs_ = model::Engine::ParseSpecs(arch, is_sparse_topology);

  if (root.exists("ERT"))
  {
    std::cout << "Found Accelergy ERT (energy reference table), replacing internal energy model." << std::endl;
    auto ert = root.lookup("ERT");
    arch_specs_.topology.ParseAccelergyERT(ert);
    if (root.exists("ART")){ // Nellie: well, if the users have the version of Accelergy that generates ART
      auto art = root.lookup("ART");
      arch_specs_.topology.ParseAccelergyART(art);  
    }
  }

  std::cout << "Begin ParseWorkload..." << std::endl;
  problem::TileFlow::ParseWorkloads(problem, workloads);
  problem::Workload::SetCurrShape(&workloads.get_shape());

  if (TileFlow::verbose_level)
    workloads.Print();

  model::TileFlow::Topology topology;

  for (unsigned storage_level_id = 0; storage_level_id < arch_specs_.topology.NumStorageLevels();
   ++ storage_level_id){
    auto buffer = arch_specs_.topology.GetStorageLevel(storage_level_id);
    TILEFLOW_COND_WARNING(buffer->size.IsSpecified(), "No memory size specified at " << buffer->name.Get());
    if (verbose_level) {
      std::cout << buffer->name.Get() << ": ";
      std::cout << buffer->size.Get() << "words" << std::endl;
    }
   }

  std::cout << "Begin Spec..." << std::endl; 
  topology.Spec(arch_specs_.topology);
  if (verbose_level)
    show_energy(topology,std::cout);

  auto mapping = 
    mapping::TileFlow::ParseAndConstruct(root.lookup("mapping"), arch_specs_, workloads);
  
  if (TileFlow::verbose_level)
    mapping.Print();

  bool enable_mem_check_ = true;
  bool enable_spatial_check_ = true;
  bool enable_loopcount_check_ = true;
  if (root.exists("check")) {
    auto checknode = root.lookup("check");
    checknode.lookupValue("mem", enable_mem_check_);
    checknode.lookupValue("spatial", enable_spatial_check_);
    checknode.lookupValue("loopcount", enable_loopcount_check_);
  }

  TileFlow::Checker checker(workloads, mapping, topology
  , enable_mem_check_, enable_spatial_check_, enable_loopcount_check_);
  
  checker.check();

  if (verbose_level)
    checker.display();

  TileFlow::mapper::Objective obj = TileFlow::mapper::CYCLE;
  unsigned timeout = 600;
  unsigned topk = 1;
  std::string search_alg = "random";
  if (root.exists("tileflow-mapper")) {
    auto mapper = root.lookup("tileflow-mapper");
    std::string objective;
    if (mapper.lookupValue("objective", objective)){
      if (objective == "cycle") obj = TileFlow::mapper::CYCLE;
      else if (objective == "energy") obj = TileFlow::mapper::ENERGY;
    }
    mapper.lookupValue("timeout", timeout);
    mapper.lookupValue("alg", search_alg);
    mapper.lookupValue("topk", topk);
  }

  TileFlow::mapper::Mapper mapper(checker.get_constraints(), workloads, mapping, arch_specs_, topology, obj, timeout, search_alg, topk);

  auto result = mapper.search();
  assert(result);

  TILEFLOW_LOG("Verify result...");
  checker.check(result);
  TILEFLOW_LOG("Check passed!");
  
  mapper.report();

  if (root.exists("output")) {
    std::string filename;
    root.lookupValue("output", filename);
    mapper.dump(filename);
  }
  
  return 0;
}