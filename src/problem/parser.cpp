#include <regex>

#include "tileflow/problem/problem.hpp"

namespace problem { 

namespace TileFlow{


std::vector<std::string> Split(
    const std::string& buffer) {
    std::vector<std::string> retval;
    std::regex re("([A-Za-z]+)[,]*[[:space:]]*", std::regex::extended);
    std::smatch sm;
    std::string str = std::string(buffer);
    str = str.substr(0, str.find("#"));

    while (std::regex_search(str, sm, re))
    {
        retval.push_back(sm[1]);
        str = sm.suffix().str();
    }

    return retval;
}

std::string ParseWorkload(config::CompoundConfigNode config, problem::TileFlow::Workload& workload) {
  std::string name;
  TILEFLOW_ASSERT(config.lookupValue("name", name), "No name specified for an op");
  workload.set_name(name);
  workload.ParseShape(config);
  std::string ins, out;
  TILEFLOW_ASSERT(config.lookupValue("ins", ins), "No ins specified for op " << name << ".");
  TILEFLOW_ASSERT(config.lookupValue("out", out), "No out specified for op " << name << ".");

  workload.set_io(Split(ins), Split(out));

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
  if (config.exists("io")){
    auto io = config.lookup("io");
    std::string ins, outs;
    TILEFLOW_ASSERT(io.lookupValue("ins", ins), "No ins property in problem::io");
    TILEFLOW_ASSERT(io.lookupValue("outs", outs), "No outs property in problem::io");
    workloads.set_io(Split(ins), Split(outs));
  }
  else {
    TILEFLOW_ERROR("No io found in problem.");  
  }
  
  if (config.exists("ops")) {
    auto ops = config.lookup("ops");
    for (int i = 0; i < ops.getLength(); ++i){
      std::shared_ptr<Workload> p_workload(new Workload(workloads));
      std::string name = problem::TileFlow::ParseWorkload(ops[i], *p_workload);
      assert(workloads.add_workload(name, p_workload));
    }
  }
  else {
    std::shared_ptr<Workload> p_workload(new Workload(workloads));
    std::string name = problem::TileFlow::ParseWorkload(config, *p_workload);
    assert(workloads.add_workload(name, p_workload));
  }
}

bool Workloads::add_workload(const std::string & name, std::shared_ptr<Workload>& workload) {
  if (workloads_.count(name)) {
    TILEFLOW_WARNING("Duplicate op named " << name << ". Drop all but the first.");
    return false;
  }

  workloads_[name] = std::move(workload);
  return true;
}

void Workloads::set_io(const std::vector<std::string>& ins, const std::vector<std::string>& outs) {
  ins_ = ins;
  outs_ = outs;
}

void Workload::set_io(const std::vector<std::string>& ins, const std::vector<std::string>& outs) {
  ins_ = ins;
  TILEFLOW_ASSERT(outs.size() == 1, "Currently, we only support op with single output");
  out_ = outs.front();
}

void Workloads::Print() {
  std::cout << "--------------Workloads------------" << std::endl;
  std::cout << "ins: ";
  for (auto t: ins_) {
    std::cout << t << ",";
  }
  std::cout << std::endl;
  std::cout << "outs: ";
  for (auto t: outs_) {
    std::cout << t << ",";
  }
  std::cout << std::endl;
  std::cout << "Tensors:" << std::endl;
  for (int i = 0; i < common_shape_.NumDataSpaces; ++i) {
    std::cout << "  " << common_shape_.DataSpaceIDToName[i];
    auto& proj = common_shape_.Projections[i];
    for (auto& expr: proj) {
      std::cout << "[";
      int tmp = 0;
      for (auto& term: expr) {
        if (term.first != -1)
          std::cout << common_shape_.CoefficientIDToName[term.first] << "*";
        std::cout << common_shape_.FactorizedDimensionIDToName[term.second];
        if (++tmp != expr.size())
          std::cout << "+";
      }
      std::cout << "]";
    }
    std::cout << std::endl;
  }
  for (auto& ptr: workloads_) {
    ptr.second->Print();
  }
  std::cout << "------------End Workloads----------" << std::endl;
}

void Workload::Print() {
  std::cout << "Op: " << name_;
  std::cout << "(";
  for (auto t: ins_) std::cout << t << ",";
  std::cout << ")->" << out_ << std::endl;
}

void Workload::apply_binding(const std::unordered_map<std::string, std::string>& binding) {
  if (binding_applied) 
    return;
  auto& common_shape_ = workloads_.common_shape_;
  
  for (auto& kv: shape_.FactorizedDimensionNameToID){
    TILEFLOW_ASSERT(binding.count(kv.first), kv.first << " of op " << name_ << " is not bond to any runtime iteration.");
    std::string iter = binding.at(kv.first);
    if (!common_shape_.FactorizedDimensionNameToID.count(iter)){
      common_shape_.FactorizedDimensionIDToName[common_shape_.NumFactorizedDimensions] = iter;
      common_shape_.FactorizedDimensionNameToID[iter] = common_shape_.NumFactorizedDimensions++;
    }
  }

  for (auto& kv: shape_.CoefficientNameToID) {
    std::string coeff_name = name_ + "::" + kv.first;
    common_shape_.CoefficientIDToName[common_shape_.NumCoefficients] = coeff_name;
    common_shape_.CoefficientNameToID[coeff_name] = common_shape_.NumCoefficients++;
  }
  
  for (auto& kv: shape_.DataSpaceIDToName) {
    std::string access_pattern_name = name_ + "::" + kv.second;
    auto & id = common_shape_.NumDataSpaces;
    common_shape_.DataSpaceOrder[id] = shape_.DataSpaceOrder[kv.first];
    common_shape_.IsReadWriteDataSpace[id] = shape_.IsReadWriteDataSpace[kv.first];
    common_shape_.DataSpaceIDToName[id] = access_pattern_name;
    common_shape_.DataSpaceNameToID[access_pattern_name] = id++;
  }

  for (auto& proj: shape_.Projections) {
    common_shape_.Projections.emplace_back();
    auto & new_proj = common_shape_.Projections.back();
    for (auto& expr: proj) {
      new_proj.emplace_back();
      auto & new_expr = new_proj.back();
      for (auto& term: expr) {
        Shape::CoefficientID new_coeff_id = term.first == shape_.NumCoefficients? -1:
        common_shape_.CoefficientNameToID[name_+"::"+shape_.CoefficientIDToName[term.first]];
        Shape::FactorizedDimensionID new_factorized_dim = common_shape_.FactorizedDimensionNameToID[binding.at(shape_.FactorizedDimensionIDToName[term.second])];
        new_expr.emplace_back(new_coeff_id, new_factorized_dim);
      }
    }
  }

  // provide a walk around; 
  common_shape_.DefaultCoefficients[-1] = 1;
  binding_applied = true;
}

} // namespace TileFlow 

} // namespace problem 