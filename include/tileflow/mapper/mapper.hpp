#pragma once 

#include <vector> 

#include "tileflow/mapper/checker.hpp"
#include "tileflow/loop-analysis/nest-analysis.hpp"

using mapping::TileFlow::Node;
using mapping::TileFlow::OpNode;
using mapping::TileFlow::TileNode;
using mapping::TileFlow::ScopeNode;
using mapping::TileFlow::Visitor;

namespace TileFlow {

class Mapper {
    const std::vector<Constraint>& constraints_;
    const problem::TileFlow::Workloads& workloads_;
    const mapping::TileFlow::Mapping& mapping_;
    const model::Engine::Specs& arch_specs_;
    const model::Topology& topology_;
    bool enable_mem_check_ = true;
    bool enable_spatial_check_ = true;
    SymbolTable optimum_;

public:
    Mapper(const std::vector<Constraint>& constraints_,
        const problem::TileFlow::Workloads& workloads_, 
        const mapping::TileFlow::Mapping& mapping_, 
        const model::Engine::Specs& arch_specs_, 
        const model::Topology& topology_,
        bool enable_mem_check_ = true, 
        bool enable_spatial_check_ = true):
        constraints_(constraints_), workloads_(workloads_), 
        mapping_(mapping_), arch_specs_(arch_specs_), 
        topology_(topology_), enable_mem_check_(enable_mem_check_), 
        enable_spatial_check_(enable_spatial_check_){}
    
    SymbolTable search();

    void dump(const std::string& filename);

    void report();

};

} // namespace TileFlow 