#include <iostream>
#include <csignal>
#include <cstring>

#include "application/model.hpp"
#include "compound-config/compound-config.hpp"
#include "util/args.hpp"

#include "timeloopX/problem/problem.hpp"
#include "timeloopX/mapping/mapping.hpp"
#include "timeloopX/loop-analysis/nest-analysis.hpp"

extern bool gTerminateEval;

//--------------------------------------------//
//                    MAIN                    //
//--------------------------------------------//

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
  
  auto problem = root.lookup("problem");
  problem::TimeloopX::Workloads workloads;
  

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

  std::cout << "level names: ";
  for (auto name: arch_specs_.topology.LevelNames())
    std::cout << name << ",";
  std::cout << std::endl;

  problem::TimeloopX::ParseWorkloads(problem, workloads);

  
  auto mapping = mapping::TimeloopX::ParseAndConstruct(root.lookup("mapping"), arch_specs_, workloads);
  
  mapping.Print();
  
  workloads.Print();

  analysis::TimeloopX::NestAnalysis analysis(workloads, mapping);
  analysis.get_tilewise_workloads();
  analysis.Print();

  std::cout << "Parser check passed!" << std::endl;

  return 0;
}
