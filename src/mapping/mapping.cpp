#include "tileflow/mapping/mapping.hpp"

namespace mapping {
    
namespace TileFlow {

void Visitor::visitScope(const ScopeNode* node){
    for (auto child: node->children_) 
        child->accept(this);
}

void Visitor::visitOp(const OpNode* node){
    for (auto child: node->children_) 
        child->accept(this);
}

void Visitor::visitTile(const TileNode* node){
    for (auto child: node->children_) 
        child->accept(this);
}

void Visitor::run(const Node* root) {
    root->accept(this);
}

loop::Nest TileNode::constructLoopNest(
    const std::map<std::string, problem::Shape::FactorizedDimensionID> & name2id) const{
    std::uint64_t storage_level = 0;
    loop::Nest loop_nest;
    uint64_t num_subnests_added = 0;
    for (auto loop: loopnests_)
    {
        loop.dimension = name2id.at(loop.name_);
        // Ignore trivial factors
        // This reduces computation time by 1.5x on average.
        if (loop.start + loop.stride < loop.end){
            loop_nest.AddLoop(loop);
            num_subnests_added ++;
        }
    }
    if (num_subnests_added == 0) {
        loop_nest.AddLoop(0, 0, 1, 1, spacetime::Dimension::Time);
    }
    loop_nest.AddStorageTilingBoundary();
    return loop_nest;
}

} // namespace TileFlow

} // namespace mapping 