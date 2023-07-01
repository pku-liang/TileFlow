#include <regex>

#include "mapping/arch-properties.hpp"

#include "tileflow/mapping/mapping.hpp"
#include "tileflow/mapper/mapper.hpp"

using TileFlow::macros;
using TileFlow::verbose_level;
using TileFlow::global_symbol_table_;

namespace mapping {

namespace TileFlow {

ArchProperties arch_props_;
const problem::TileFlow::Workloads* p_workloads_ = nullptr;

void tolower(std::string& str){
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) {return std::tolower(c);});
}

Node::Node(
    Node::type_t t, 
    config::CompoundConfigNode config):type_(t){
    name_ = type2name_.at(type_);

    if (config.exists("bypass"))
        config.lookupArrayValue("bypass", bypassed_);
    
    config.lookupValue("profile", profile_);
    std::string tag = "";
    config.lookupValue("tag", tag);
    if (tag.size())
        name_ += "::" + tag;
}

ScopeNode::ScopeNode(config::CompoundConfigNode config): Node(Node::Scope, config){
    std::string type_s = "sequential";
    config.lookupValue("type", type_s);
    tolower(type_s);
    if (type_s.find("seq") != std::string::npos) {
        type = Sequential;
        name_ += "::Sequential";
    }
    else if (type_s.find("para") != std::string::npos) {
        type = Parallel;
        name_ += "::Parallel";
    }
    else if (type_s.find("pipe") != std::string::npos) {
        type = Pipeline;
        name_ += "::Pipeline";
    }
    else if (type_s.find("shar") != std::string::npos) {
        type = Sharing;
        name_ += "::Sharing";
    }
    else {TILEFLOW_ERROR("ScopeNode type error. Should has type sequential/parallel");}
    

}   

TileNode::TileNode(config::CompoundConfigNode config): Node(Node::Tile, config) {
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
    
    TILEFLOW_ASSERT(iters.size() == loop_bounds.size(), "permutation " << buffer << " & factor iter mismatch");

    ParseStorageLevel(config);

    if (config.exists("multicast")) {
        config.lookupValue("multicast", multicast_enabled_);
    }

    unsigned split = iters.size();
    config.lookupValue("split", split);
    for (int i = (int)iters.size()-1; i >= 0; --i) {
        loopnests_.emplace_back();
        loop::TileFlow::Descriptor& loop = loopnests_.back();
        loop.name_ = iters[i];
        loop.dimension = problem::GetShape()->FactorizedDimensionNameToID.at(loop.name_);
        loop.start = 0;
        loop.end = loop_bounds[iters[i]].first;
        loop.residual_end = loop_bounds[iters[i]].second;
        loop.stride = 1;
        loop.spacetime_dimension = type_s == "spatial"? 
        ((unsigned)i < split? spacetime::Dimension::SpaceX : spacetime::Dimension::SpaceY) 
        : spacetime::Dimension::Time; 
    }
    
    name_ += type_ == Temporal? "::Temporal" : "::Spatial"; 
}


OpNode::OpNode(config::CompoundConfigNode config): Node(Node::Op, config) {
    assert(config.lookupValue("name", op_name_));
    p_workload = p_workloads_->get_workload(op_name_);
    name_ += "::" + op_name_;
}

std::unordered_map<std::string, std::pair<int, int> > Node::ParseFactors(
    const std::string& buffer) {
    std::unordered_map<std::string, std::pair<int, int> > loop_bounds;
    std::regex re("([A-Za-z]+)[[:space:]]*[=]*[[:space:]]*([0-9A-Za-z_?]+)(,([0-9]+))?", std::regex::extended);
    std::smatch sm;
    std::string str = std::string(buffer);
    str = str.substr(0, str.find("#"));

    while (std::regex_search(str, sm, re))
    {
        std::string dimension_name = sm[1];

        int end;
        if (macros.exists(sm[2])){
            macros.lookupValue(sm[2], end);
        }
        else {
            char* ptr = nullptr;
            end = std::strtol(sm[2].str().c_str(), &ptr, 10);
            if (ptr && *ptr) {
                end = global_symbol_table_.insert(sm[2]);
            }
        }

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
  name_ += "::" + storage_level_name_;
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

void TileNode::display(std::string prefix, bool recursive, const SymbolTable* symbol_table, std::ostream& o) const{
    display_active_tensors(prefix, o);
    if (symbol_table == nullptr) symbol_table = & global_symbol_table_;
    for (auto& loop: loopnests_) {
        o << prefix;
        o << "for " << loop.name_ << " in [" << loop.start << ":";
        if (loop.end < 0) {
            auto& entry = symbol_table->lookup(loop.end);
            o << entry.name_;
            if (entry.fixed_) o << "(" << entry.value_ << ")";
        }
        else o << loop.end; 
        o << ")";
        if (loop::IsSpatial(loop.spacetime_dimension))
        {
            if (loop::IsSpatialX(loop.spacetime_dimension))
                o << " (Spatial-X)";
            else
                o << " (Spatial-Y)";
        }
        o << ", " << arch_props_.Specs().topology.GetStorageLevel(storage_level_)->level_name;
        o << std::endl;
        prefix += "  ";
    }

    if (recursive)
        for (auto child: children_)
            child->display(prefix, recursive, symbol_table, o);
}

void ScopeNode::display(std::string prefix, bool recursive, const SymbolTable* symb_table, std::ostream& o) const{
    display_active_tensors(prefix, o);
    o << prefix << "Scope: ";
    if (type == Sequential) o << "Sequential";
    else if (type == Parallel) o << "Parallel";
    else if (type == Pipeline) o << "Pipeline";
    if (recursive) {
        o << "{" << std::endl;
        for (auto child: children_) 
            child->display(prefix + "  ", recursive, symb_table, o);
        o << prefix << "}" << std::endl;
    }
    else o << std::endl;
}

void OpNode::display(std::string prefix, bool, const SymbolTable*, std::ostream& o) const {
    display_active_tensors(prefix, o);
    o << prefix;
    p_workload->Print(o);
    o << std::endl;
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

void Mapping::Print() const{
    std::cout << "-----------------Mapping---------------" << std::endl;
    std::cout << "root: " << root << std::endl;
    root->display("");
    std::cout << "---------------------------------------" << std::endl;
}


} // namespace TileFlow 

} // namespace mapping 