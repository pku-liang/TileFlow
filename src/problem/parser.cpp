#include "timeloopX/problem/problem.hpp"

namespace problem { 

namespace TimeloopX{


std::string ParseWorkload(config::CompoundConfigNode config, problem::Workload& workload) {
  std::string name;
  if (config.exists("op")) {
    auto op = config.lookup("op");
    workload.ParseShape(op);
    assert(op.lookupValue("name", name));
  }
  else {
    std::cerr << "ERROR: no op property passed for an operation." << std::endl;
    exit(1);
  }
  if (config.exists("instance")) {
    auto factorized_bound = config.lookup("instance");
    problem::ParseWorkloadInstance(factorized_bound, workload);  
  }
  else { 
    std::cerr << "ERROR: no instance passed for an op. Please make sure an instance is passed for each op." << std::endl;
    exit(1);
  }

  return name;
}

void ParseWorkloads(config::CompoundConfigNode config, Workloads& workloads) {
  std::cerr << config.isList() << "," << config.isArray() << std::endl;
  if (config.isList()) {
    for (int i = 0; i < config.getLength(); ++i) {
      std::shared_ptr<problem::Workload> p_workload(new Workload());
      std::string name = problem::TimeloopX::ParseWorkload(config[i], *p_workload);
      assert(workloads.add_workload(name, p_workload));
    }
  }
  else {
    std::shared_ptr<problem::Workload> p_workload(new Workload());
    std::string name = problem::TimeloopX::ParseWorkload(config, *p_workload); 
    assert(workloads.add_workload(name, p_workload));
  }
}

bool Workloads::add_workload(const std::string & name, std::shared_ptr<problem::Workload>& workload) {
  if (opname2id.count(name)) {
    TIMELOOPX_WARNING("Duplicate op named " << name << ". Drop all but the first.");
    return false;
  }
  opname2id[name] = workloads.size();
  workloads.push_back(std::move(workload));
  return true;
}


} // namespace TimeloopX 

} // namespace problem 