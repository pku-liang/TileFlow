#pragma once 

#include "tileflow/loop-analysis/nest-analysis.hpp"

#include "model/topology.hpp"

namespace model {

namespace TileFlow {


class Topology: public model::Topology {
public: 
    void eval(
        const mapping::TileFlow::Mapping& mapping, 
        const analysis::TileFlow::NestAnalysis& analysis);    
    friend class StatCalculator;
};

class StatCalculator: public mapping::TileFlow::Visitor {
    void visitTile(const TileNode*) override;
    void visitScope(const ScopeNode*) override;
    void visitOp(const OpNode*) override;
    bool break_on_failure;
    std::stack<std::uint64_t> cycles_;
    double energy_;
    model::TileFlow::Topology& topology_;
    const mapping::TileFlow::Mapping& mapping_;
    const analysis::TileFlow::NestAnalysis& analysis_;
public: 
    StatCalculator(model::TileFlow::Topology& topology, 
    const mapping::TileFlow::Mapping& mapping,
    const analysis::TileFlow::NestAnalysis& analysis)
    : topology_(topology), mapping_(mapping), analysis_(analysis){}
    void run(const Node* root) override;
};

}

} // namespace model 