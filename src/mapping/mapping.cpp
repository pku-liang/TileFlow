#include "tileflow/mapping/mapping.hpp"

using TileFlow::global_symbol_table_;

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

loop::Nest TileNode::constructLoopNest(const SymbolTable* symbol_table_) const{
    loop::Nest loop_nest;
    uint64_t num_subnests_added = 0;
    for (auto loop: loopnests_)
    {
        // Ignore trivial factors
        // This reduces computation time by 1.5x on average.
        if (loop.end <= 0) {
            assert(symbol_table_);
            loop.end = symbol_table_->lookup(loop.end).value_;
        }
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

void Node::display_active_tensors(std::string prefix) const {
    bool isEmpty = active_tensors_.read_tensors.size()
        + active_tensors_.update_tensors.size()
        + active_tensors_.fill_tensors.size()
        + active_tensors_.wb_tensors.size();
    if (!isEmpty) return;
    std::cout << prefix;
    if (active_tensors_.read_tensors.size()) {
        std::cout << "read: ";
        for (auto id: active_tensors_.read_tensors) 
            std::cout << problem::GetShape()->DataSpaceIDToName.at(id) << " ";
    }
    if (active_tensors_.update_tensors.size()) {
        std::cout << "update: ";
        for (auto id: active_tensors_.update_tensors) 
            std::cout << problem::GetShape()->DataSpaceIDToName.at(id) << " ";
    }
    if (active_tensors_.fill_tensors.size()){
        std::cout << "fill: ";
        for (auto id: active_tensors_.fill_tensors) 
            std::cout << problem::GetShape()->DataSpaceIDToName.at(id) << " ";
    }
    if (active_tensors_.wb_tensors.size()) {
        std::cout << "write-back: ";
        for (auto id: active_tensors_.wb_tensors) 
            std::cout << problem::GetShape()->DataSpaceIDToName.at(id) << " ";
    }
    std::cout << std::endl;
}

} // namespace TileFlow

} // namespace mapping 