#include <regex>

#include "mapping/arch-properties.hpp"

#include "tileflow/mapping/mapping.hpp"

namespace mapping {

namespace TileFlow {

ArchProperties arch_props_;
const problem::TileFlow::Workloads* p_workloads_ = nullptr;

void tolower(std::string& str){
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) {return std::tolower(c);});
}

ScopeNode::ScopeNode(config::CompoundConfigNode config): Node(Node::Scope){
    std::string type_s = "sequential";
    config.lookupValue("type", type_s);
    tolower(type_s);
    if (type_s.find("seq") != std::string::npos) {
        type = Sequential;
    }
    else if (type_s.find("para") != std::string::npos) {
        type = Parallel;
    }
    else if (type_s.find("pipe") != std::string::npos) {
        type = Pipeline;
    }
    else {TILEFLOW_ERROR("ScopeNode type error. Should has type sequential/parallel");}

}   

TileNode::TileNode(config::CompoundConfigNode config): Node(Node::Tile) {
    std::string type_s = "temporal";
    config.lookupValue("type", type_s);
    tolower(type_s);    
    if (type_s.find("temp") != std::string::npos) {
        type_ = Temporal;
    }
    else if (type_s.find("spatial") != std::string::npos) {
        type_ = Spatial;
    }
    else {
        TILEFLOW_ERROR("Unknown Tile type" << type_s);
    }

    std::unordered_map<std::string, std::pair<int, int> > loop_bounds;
    std::string buffer;
    if (config.lookupValue("factors", buffer)) {
        loop_bounds = ParseFactors(buffer);
    }
    else {
        TILEFLOW_ERROR("No factors");
    }

    std::vector<std::string> iters;
    if (config.lookupValue("permutation", buffer)) {
        iters = ParsePermutations(buffer);
    }
    else {
        for (auto& kv: loop_bounds) iters.push_back(kv.first);
        TILEFLOW_WARNING("No permutation specified. Infer instead.");
    }
    
    TILEFLOW_ASSERT(iters.size() == loop_bounds.size(), "permutation & factor iter mismatch");

    ParseStorageLevel(config);

    unsigned split = iters.size();
    config.lookupValue("split", split);
    for (int i = (int)iters.size()-1; i >= 0; --i) {
        loopnests_.emplace_back();
        loop::TileFlow::Descriptor& loop = loopnests_.back();
        loop.name_ = iters[i];
        loop.start = 0;
        loop.end = loop_bounds[iters[i]].first;
        loop.residual_end = loop_bounds[iters[i]].second;
        loop.stride = 1;
        loop.spacetime_dimension = type_s == "spatial"? 
        ((unsigned)i < split? spacetime::Dimension::SpaceX : spacetime::Dimension::SpaceY) 
        : spacetime::Dimension::Time; 
    }
    
}

std::unordered_map<std::string, std::pair<int, int> > Node::ParseFactors(
    const std::string& buffer) {
    std::unordered_map<std::string, std::pair<int, int> > loop_bounds;
    std::regex re("([A-Za-z]+)[[:space:]]*[=]*[[:space:]]*([0-9]+)(,([0-9]+))?", std::regex::extended);
    std::smatch sm;
    std::string str = std::string(buffer);
    str = str.substr(0, str.find("#"));

    while (std::regex_search(str, sm, re))
    {
        std::string dimension_name = sm[1];

        int end = std::stoi(sm[2]);

        int residual_end = end;
        if (sm[4] != "")
        {
            residual_end = std::stoi(sm[4]);
        }

        loop_bounds[dimension_name] = {end, residual_end};

        str = sm.suffix().str();
    }

    return loop_bounds;
}

std::vector<std::string> Node::ParsePermutations(
    const std::string & buffer 
){
    std::vector<std::string> iters;
    
    std::istringstream iss(buffer);
    char token;
    while (iss >> token) {
        iters.push_back(std::string(1, token));
    }
    
    return iters;
}

void Node::ParseStorageLevel(config::CompoundConfigNode directive)
{
  auto num_storage_levels = arch_props_.StorageLevels();
    
  //
  // Find the target storage level. This can be specified as either a name or an ID.
  //
  std::string storage_level_name;
  unsigned storage_level_id;
    
  if (directive.lookupValue("target", storage_level_name))
  {
    // Find this name within the storage hierarchy in the arch specs.
    for (storage_level_id = 0; storage_level_id < num_storage_levels; storage_level_id++)
    {
      if (arch_props_.Specs().topology.GetStorageLevel(storage_level_id)->level_name == storage_level_name)
        break;
    }
    if (storage_level_id == num_storage_levels)
    {
      std::cerr << "ERROR: target storage level not found: " << storage_level_name << std::endl;
      exit(1);
    }
  }
  else
  {
    int id;
    assert(directive.lookupValue("target", id));
    assert(id >= 0  && id < int(num_storage_levels));
    storage_level_id = static_cast<unsigned>(id);
  }

  assert(storage_level_id < num_storage_levels);

  storage_level_name_ = storage_level_name;
  storage_level_ = storage_level_id;
}

OpNode::OpNode(config::CompoundConfigNode config): Node(Node::Op) {
    assert(config.lookupValue("name", name_));
    p_workload = p_workloads_->get_workload(name_);

    std::string str;
    assert(config.lookupValue("binding", str));
    std::regex re("([A-Za-z])[[:space:]]*[:][[:space:]]*([A-Za-z])", std::regex::extended);
    std::smatch sm;

    while (std::regex_search(str, sm, re))
    {
        std::string iter_name = sm[1];

        std::string dim_name = sm[2];

        binding_[dim_name] = iter_name;

        str = sm.suffix().str();
    }
    p_workload->apply_binding(binding_);
}


Node* RecursiveParse(config::CompoundConfigNode config) {
    std::string node_type; 

    if (!config.lookupValue("node-type", node_type)) {
        TILEFLOW_ERROR("No node-type is specified.");
    }
    tolower(node_type);

    Node * node = nullptr;
    if (node_type == "tile") {
        node = new TileNode(config);
    }
    else if (node_type == "op") {
        node = new OpNode(config);
    }
    else if (node_type == "scope") {
        node = new ScopeNode(config);
    }
    else {
        TILEFLOW_ERROR(node_type << " is not a valid type.");
    }
        
    assert(node != nullptr);
    
    if (config.exists("subtree")) {
        config = config.lookup("subtree");
        if (config.isList()){
            for (int i = 0; i < config.getLength(); i++){
                node->add_child(RecursiveParse(config[i]));
            }
        }
        else {
            node->add_child(RecursiveParse(config));
        }
    }
    else if (node_type != "op") {
        TILEFLOW_ERROR("Exist non-Op leaf node.");
    }

    return node;
}

void TileNode::display(std::string prefix, bool recursive) const{
    for (auto& loop: loopnests_) {
        std::cout << prefix;
        loop.Print(std::cout, true);
        std::cout << ", " << arch_props_.Specs().topology.GetStorageLevel(storage_level_)->level_name;
        std::cout << std::endl;
        prefix += "  ";
    }
    if (recursive)
        for (auto child: children_)
            child->display(prefix);
}

void ScopeNode::display(std::string prefix, bool recursive) const{
    std::cout << prefix << "Scope: ";
    if (type == Sequential) std::cout << "Sequential";
    else if (type == Parallel) std::cout << "Parallel";
    else if (type == Pipeline) std::cout << "Pipeline";
    if (recursive) {
        std::cout << "{" << std::endl;
        for (auto child: children_) 
            child->display(prefix + "  ");
        std::cout << prefix << "}" << std::endl;
    }
}

void OpNode::display(std::string prefix, bool) const {
    std::cout << prefix;
    p_workload->Print();
    std::cout << prefix << ", binding: "; 
    for (auto& bind: binding_) {
        std::cout << bind.first << ":" << bind.second;
        std::cout << ",";
    } 
    std::cout << std::endl;
}

Mapping ParseAndConstruct(config::CompoundConfigNode config,
                          model::Engine::Specs& arch_specs,
                          const problem::TileFlow::Workloads& workloads)
{
    arch_props_ = ArchProperties();
    arch_props_.Construct(arch_specs);
    p_workloads_ = &workloads;
    
    Mapping mapping;
    mapping.root = RecursiveParse(config);
    mapping.fanoutX_map = arch_props_.FanoutX();
    mapping.fanoutY_map = arch_props_.FanoutY();

    return mapping;
}

void Mapping::Print() {
    std::cout << "-----------------Mapping---------------" << std::endl;
    std::cout << "root: " << root << std::endl;
    root->display("");
    std::cout << "---------------------------------------" << std::endl;

}


} // namespace TileFlow 

} // namespace mapping 