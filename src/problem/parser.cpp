#include <regex>

#include "tileflow/problem/problem.hpp"

using TileFlow::macros;
using TileFlow::verbose_level;

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

  Workload::Densities densities;
  std::string density_distribution;

  // shared pointer for parsed density distribution specs
  std::shared_ptr<DensityDistributionSpecs> density_distribution_specs;
  YAML::Node ynode;

  // 1) shared density specification for all dataspaces
  double common_avg_density;
  if (config.exists("commonDensity")){
    config::CompoundConfigNode density_config;
    if (! config.lookup("commonDensity").isMap()){
      config.lookupValue("commonDensity", common_avg_density);
      ynode["distribution"] = "fixed-structured";
      ynode["density"] = common_avg_density;
      density_config = config::CompoundConfigNode(nullptr, ynode, new config::CompoundConfig("dummy.yaml"));
    } else {
      density_config = config.lookup("commonDensity");
    }
    auto density_specs = DensityDistributionFactory::ParseSpecs(density_config);
    // assign all dataspaces the same density value
    for (unsigned i = 0; i < GetShape()->NumDataSpaces; i++){
      densities[i]= DensityDistributionFactory::Construct(density_specs);
      // make sure the density model is correctly set
      assert (densities[i] != NULL);
    }
  }

  // 2) density specifications for each dataspace
  else if (config.exists("densities"))
  {
    auto config_densities = config.lookup("densities");
    for (unsigned i = 0; i < GetShape()->NumDataSpaces; i++){
      double dataspace_avg_density;
      config::CompoundConfigNode density_config;
      std::string dataspace_name = GetShape()->DataSpaceIDToName.at(i);

      if (config_densities.exists(GetShape()->DataSpaceIDToName.at(i)))
      {
		config_densities.lookupValue(GetShape()->DataSpaceIDToName.at(i), dataspace_avg_density);
        
		// if the specific dataspace's density is specified
        if (!config_densities.lookup(GetShape()->DataSpaceIDToName.at(i)).isMap())
        {
          // single number for density is given, default to fixed density distribution
          assert(config_densities.lookupValue(GetShape()->DataSpaceIDToName.at(i), dataspace_avg_density));
          ynode["distribution"] = "fixed-structured";
          ynode["density"] = dataspace_avg_density;
          density_config = config::CompoundConfigNode(nullptr, ynode, new config::CompoundConfig("dummy.yaml"));
        } else
        {
          density_config = config_densities.lookup(GetShape()->DataSpaceIDToName.at(i));
        }
        auto density_specs = DensityDistributionFactory::ParseSpecs(density_config);
        densities[i] = DensityDistributionFactory::Construct(density_specs);
      }
      else
      {
        // no density specified, roll back to default
        ynode["distribution"] = "fixed-structured";
        ynode["density"] = 1.0;
        density_config = config::CompoundConfigNode(nullptr, ynode, new config::CompoundConfig("dummy.yaml"));
        auto density_specs = DensityDistributionFactory::ParseSpecs(density_config);
        densities[i]= DensityDistributionFactory::Construct(density_specs);
      }

      // make sure the density model is correctly set
      assert (densities[i] != NULL);
    }

    // 3) no density specification -> dense workload tensors
  } else {
    config::CompoundConfigNode density_config;
    for (unsigned i = 0; i < GetShape()->NumDataSpaces; i++){
      ynode["distribution"] = "fixed-structured";
      ynode["density"] = 1.0;
      density_config = config::CompoundConfigNode(nullptr, ynode, new config::CompoundConfig("dummy.yaml"));
      auto density_specs = DensityDistributionFactory::ParseSpecs(density_config);
      densities[i]= DensityDistributionFactory::Construct(density_specs);

      // make sure the density model is correctly set
      assert (densities[i] != NULL);
    }
  }
  workload.SetDensities(densities);

  // 3) set common shape
  // workload.set_common_shape();

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

  std::vector<std::string> dims;
  if (config.exists("dimensions")) {
    config.lookupArrayValue("dimensions", dims);
  }
  else {
    TILEFLOW_ERROR("no dimensions passed for the workloads.");
  }
  workloads.set_dims(dims);
  
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

  if (config.exists("instance")) {
    auto factorized_bounds = config.lookup("instance");
    for (auto dim: dims) {
      TILEFLOW_ASSERT(factorized_bounds.exists(dim), "no instance passed for axis " << dim << ".");
      std::string _tmp;
      factorized_bounds.lookupValue(dim, _tmp);
      int bound;
      if (macros.exists(_tmp)) {
        TILEFLOW_ASSERT(macros.lookupValue(_tmp, bound), 
        _tmp << " is not a specified MACRO");
      }
      else bound = std::stoi(_tmp);
      workloads.set_factorized_bound(dim, bound); 
    }
  }
  else { 
    std::cerr << "ERROR: no instance passed for an op. Please make sure an instance is passed for each op." << std::endl;
    exit(1);
  }

  if (config.exists("coefficient")) {
    workloads.set_coeffs(config);
  }
}

void Workloads::set_coeffs(const config::CompoundConfigNode& coeffs){
    for (unsigned i = 0; i < GetShape()->NumCoefficients; i++)
    {
      coefficients_[i] = GetShape()->DefaultCoefficients.at(i);
      coeffs.lookupValue(GetShape()->CoefficientIDToName.at(i), coefficients_[i]);
    }
}

bool Workloads::add_workload(const std::string & name, std::shared_ptr<Workload>& workload) {
  if (workloads_.count(name)) {
    TILEFLOW_WARNING("Duplicate op named " << name << ". Drop all but the first.");
    return false;
  }


  auto shape_ = workload->GetShape();
  for (auto& kv: shape_->FlattenedDimensionNameToID) {
    TILEFLOW_ASSERT(common_shape_.FlattenedDimensionNameToID.count(kv.first), 
      "Op:" << name << "'s dimension " << kv.first <<  " is not declared in global scope.");
  }
  for (auto& kv: shape_->DataSpaceIDToName) {
    std::string access_pattern_name = kv.second;
    if (common_shape_.DataSpaceNameToID.count(access_pattern_name)) {
      auto id = common_shape_.DataSpaceNameToID.at(access_pattern_name);
      // TODO: there should be a check on the access pattern; 
      TILEFLOW_ASSERT(common_shape_.DataSpaceOrder[id] == shape_->DataSpaceOrder.at(kv.first),
      "tensor " << access_pattern_name << " has mismatched order");
      continue;
    }
    auto& proj = shape_->Projections[kv.first];
    auto & id = common_shape_.NumDataSpaces;
    assert(id == common_shape_.Projections.size());
    common_shape_.Projections.emplace_back();
    auto & new_proj = common_shape_.Projections.back();
    for (auto& expr: proj) {
      new_proj.emplace_back();
      auto & new_expr = new_proj.back();
      for (auto& term: expr) {
        Shape::CoefficientID new_coeff_id = term.first == shape_->NumCoefficients? -1:
        common_shape_.CoefficientNameToID[shape_->CoefficientIDToName.at(term.first)];
        Shape::FactorizedDimensionID new_factorized_dim = 
          common_shape_.FactorizedDimensionNameToID[shape_->FactorizedDimensionIDToName.at(term.second)];
        new_expr.emplace_back(new_coeff_id, new_factorized_dim);
      }
    }
    common_shape_.DataSpaceOrder[id] = shape_->DataSpaceOrder.at(kv.first);
    common_shape_.IsReadWriteDataSpace[id] = shape_->IsReadWriteDataSpace.at(kv.first);
    common_shape_.DataSpaceIDToName[id] = access_pattern_name;
    common_shape_.DataSpaceNameToID[access_pattern_name] = id;
    densities_[id] = workload->GetDensity(kv.first);
    id++;
  }

  for (auto& kv: shape_->CoefficientNameToID) {
    std::string coeff_name = kv.first;
    common_shape_.CoefficientIDToName[common_shape_.NumCoefficients] = coeff_name;
    common_shape_.CoefficientNameToID[coeff_name] = common_shape_.NumCoefficients++;
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
  std::cout << "dimensions:";
  for (auto& kv: factorized_bounds_) 
    std::cout << "[" << common_shape_.FactorizedDimensionIDToName[kv.first] 
      << "," << kv.second << "]";
  std::cout << std::endl;
  // std::cout << "ins: ";
  // for (auto t: ins_) {
  //   std::cout << t << ",";
  // }
  // std::cout << std::endl;
  // std::cout << "outs: ";
  // for (auto t: outs_) {
  //   std::cout << t << ",";
  // }
  // std::cout << std::endl;
  std::cout << "Tensors:" << std::endl;
  for (int i = 0; i < (int)common_shape_.NumDataSpaces; ++i) {
    std::cout << "  " << common_shape_.DataSpaceIDToName[i];
    auto& proj = common_shape_.Projections[i];
    for (auto& expr: proj) {
      std::cout << "[";
      int tmp = 0;
      for (auto& term: expr) {
        if (term.first != (unsigned)-1)
          std::cout << common_shape_.CoefficientIDToName[term.first] << "*";
        std::cout << common_shape_.FactorizedDimensionIDToName[term.second];
        if (++tmp != (int)expr.size())
          std::cout << "+";
      }
      std::cout << "]";
    }
    std::cout << std::endl;
  }
  for (auto& ptr: workloads_) {
    ptr.second->Print();
  }
  std::cout << "UseFlattening: " << common_shape_.UsesFlattening << std::endl;
  std::cout << "Shape:" << std::endl;
  common_shape_.show();
  std::cout << "------------End Workloads----------" << std::endl;
}

void Workload::Print(std::ostream& o) {
  o << "Op: " << name_;
  o << "(";
  for (auto t: ins_) o << t << ",";
  o << ")->" << out_ << std::endl;
}

void Workloads::set_factorized_bound(const std::string& dim, int bound) {
  assert(common_shape_.FactorizedDimensionNameToID.count(dim));
  factorized_bounds_[common_shape_.FactorizedDimensionNameToID[dim]] = bound;
}

void Workloads::set_dims(const std::vector<std::string>& dims) {
  assert(common_shape_.NumFlattenedDimensions == 0);
  assert(common_shape_.NumFactorizedDimensions == 0);
  for (auto& dim: dims) {
    common_shape_.FactorizedDimensionNameToID[dim] = common_shape_.NumFactorizedDimensions;
    common_shape_.FactorizedDimensionIDToName[common_shape_.NumFactorizedDimensions] = dim;
    common_shape_.FlattenedDimensionIDToName[common_shape_.NumFlattenedDimensions] = dim;
    common_shape_.FlattenedDimensionNameToID[dim] = common_shape_.NumFlattenedDimensions;
    common_shape_.FlattenedToFactorized.push_back({common_shape_.NumFactorizedDimensions});
    common_shape_.FactorizedToFlattened[common_shape_.NumFlattenedDimensions] = common_shape_.NumFactorizedDimensions;
    common_shape_.FactorizedToFlattened[common_shape_.NumFactorizedDimensions] = common_shape_.NumFlattenedDimensions;
    common_shape_.FlattenedToFactorized.push_back({common_shape_.NumFactorizedDimensions});
    common_shape_.NumFlattenedDimensions ++;
    common_shape_.NumFactorizedDimensions ++;
  }
}

const problem::Workload& Workloads::get_workload() const {
  if (!workload_constructed_) {
    workload_.SetShape(common_shape_);
    workload_.SetCoefficients(coefficients_);
    workload_.SetFactorizedBounds(factorized_bounds_);
    workload_.SetDensities(densities_);
    workload_constructed_ = true;
  }
  return workload_;
}

} // namespace TileFlow 

} // namespace problem 

/**
 * 1. GetCoefficient
 * 2. vector_stride_ scale
 * 3. lifting loopnest computation;
*/