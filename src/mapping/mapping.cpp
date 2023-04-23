#include "tileflow/mapping/mapping.hpp"

using TileFlow::global_symbol_table_;

namespace mapping {
    
namespace TileFlow {

const std::unordered_map<Node::type_t, std::string> Node::type2name_ = {
    {Node::Tile, "Tile"},
    {Node::Op, "Op"},
    {Node::Scope, "Scope"}
};

void Node::add_child(const Node* child){
    if (type_ == Node::Scope) {
        
        unsigned storage_level;
        std::string storage_level_name = "Unknown";
        if (child->get_type() == Node::Tile) {
            // if (static_cast<const TileNode*>(child)->get_tile_type() == TileNode::Temporal){
            //     storage_level = child->get_storage_level() + 1;
            // }
            // else {
                storage_level = child->get_storage_level();
                storage_level_name = child->get_storage_name();
            // }
        }
        else if (child->get_type() == Node::Scope) {
            storage_level = child->get_storage_level();
            storage_level_name = child->get_storage_name();
        }
        else {
            TILEFLOW_ERROR("Scope Node should not have a op child");
        }
        assert(storage_level_ == unsigned(-1) || storage_level_ == storage_level);
        storage_level_ = storage_level;
        storage_level_name_ = storage_level_name;
    }
    assert(child != nullptr); 
    children_.push_back(child); 
    child->set_parent(this);
}

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
            loop.residual_end = loop.end = symbol_table_->lookup(loop.end).value_;
        }
        if (loop.start + loop.stride < loop.end){
            assert((type_==TileNode::Spatial && loop::IsSpatial(loop.spacetime_dimension))
            || (type_==TileNode::Temporal && !loop::IsSpatial(loop.spacetime_dimension)));
            loop_nest.AddLoop(loop);
            num_subnests_added ++;
        }
    }
    if (num_subnests_added == 0) {
        loop_nest.AddLoop(0, 0, 1, 1, type_ == TileNode::Spatial? spacetime::Dimension::SpaceX : spacetime::Dimension::Time);
    }
    loop_nest.AddStorageTilingBoundary();
    return loop_nest;
}

void Node::display_active_tensors(std::string prefix, std::ostream&o) const {
    bool isEmpty = active_tensors_.read_tensors.size()
        + active_tensors_.update_tensors.size()
        + active_tensors_.fill_tensors.size()
        + active_tensors_.wb_tensors.size();
    if (!isEmpty) return;
    o << prefix;
    if (active_tensors_.read_tensors.size()) {
        o << "read: ";
        for (auto id: active_tensors_.read_tensors) 
            o << problem::GetShape()->DataSpaceIDToName.at(id) << " ";
    }
    if (active_tensors_.update_tensors.size()) {
        o << "update: ";
        for (auto id: active_tensors_.update_tensors) 
            o << problem::GetShape()->DataSpaceIDToName.at(id) << " ";
    }
    if (active_tensors_.fill_tensors.size()){
        o << "fill: ";
        for (auto id: active_tensors_.fill_tensors) 
            o << problem::GetShape()->DataSpaceIDToName.at(id) << " ";
    }
    if (active_tensors_.wb_tensors.size()) {
        o << "write-back: ";
        for (auto id: active_tensors_.wb_tensors) 
            o << problem::GetShape()->DataSpaceIDToName.at(id) << " ";
    }
    o << std::endl;
}

} // namespace TileFlow

} // namespace mapping 