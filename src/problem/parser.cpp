#include <regex>

#include "timeloopX/problem/problem.hpp"

namespace problem { 

namespace TimeloopX{


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

std::string ParseWorkload(config::CompoundConfigNode config, problem::TimeloopX::Workload& workload) {
  std::string name;
  TIMELOOPX_ASSERT(config.lookupValue("name", name), "No name specified for an op");
  workload.set_name(name);
  workload.ParseShape(config);
  std::string ins, out;
  TIMELOOPX_ASSERT(config.lookupValue("ins", ins), "No ins specified for op " << name << ".");
  TIMELOOPX_ASSERT(config.lookupValue("out", out), "No out specified for op " << name << ".");

  workload.set_io(Split(ins), Split(out));

  workload.ParseShape(config);

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
    TIMELOOPX_ASSERT(io.lookupValue("ins", ins), "No ins property in problem::io");
    TIMELOOPX_ASSERT(io.lookupValue("outs", outs), "No outs property in problem::io");
    workloads.set_io(Split(ins), Split(outs));
  }
  else {
    TIMELOOPX_ERROR("No io found in problem.");  
  }
  
  if (config.exists("ops")) {
    auto ops = config.lookup("ops");
    for (int i = 0; i < ops.getLength(); ++i){
      std::shared_ptr<Workload> p_workload(new Workload());
      std::string name = problem::TimeloopX::ParseWorkload(ops[i], *p_workload);
      assert(workloads.add_workload(name, p_workload));
    }
  }
  else {
    std::shared_ptr<Workload> p_workload(new Workload());
    std::string name = problem::TimeloopX::ParseWorkload(config, *p_workload);
    assert(workloads.add_workload(name, p_workload));
  }
}

bool Workloads::add_workload(const std::string & name, std::shared_ptr<Workload>& workload) {
  if (opname2id_.count(name)) {
    TIMELOOPX_WARNING("Duplicate op named " << name << ". Drop all but the first.");
    return false;
  }
  opname2id_[name] = workloads_.size();
  

  for (auto& t: workload->ins_)
    if (!tensor2id_.count(t)) {
      tensor2id_[t] = tensor2id_.size();
      id2tensor_.push_back(t);
    }
  if (!tensor2id_.count(workload->out_)){
    tensor2id_[workload->out_] = tensor2id_.size();
    id2tensor_.push_back(workload->out_);
  }

  workloads_.push_back(std::move(workload));
  return true;
}

void Workloads::set_io(const std::vector<std::string>& ins, const std::vector<std::string>& outs) {
  ins_ = ins;
  outs_ = outs;
  for (auto& in: ins) {
    if (tensor2id_.count(in) == 0) {
      tensor2id_[in] = tensor2id_.size();
      id2tensor_.push_back(in);
    }
  }
  for (auto& out: outs) {
    if (tensor2id_.count(out) == 0) {
      tensor2id_[out] = tensor2id_.size();
      id2tensor_.push_back(out);
    }
  }
}

void Workload::set_io(const std::vector<std::string>& ins, const std::vector<std::string>& outs) {
  ins_ = ins;
  TIMELOOPX_ASSERT(outs.size() == 1, "Currently, we only support op with single output");
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
  for (auto& ptr: workloads_) {
    ptr->Print();
  }
  std::cout << "------------End Workloads----------" << std::endl;
}

void Workload::Print() {
  std::cout << "Op: " << name_;
  std::cout << "(";
  for (auto t: ins_) std::cout << t << ",";
  std::cout << ")->" << out_ << std::endl;
}

} // namespace TimeloopX 

} // namespace problem 