#include "timeloopX/mapping/mapping.hpp"

namespace mapping {
    
namespace TimeloopX {

void Visitor::visitScope(ScopeNode* node){
    for (auto child: node->children_) 
        child->accept(this);
}

void Visitor::visitOp(OpNode* node){
    for (auto child: node->children_) 
        child->accept(this);
}

void Visitor::visitTile(TileNode* node){
    for (auto child: node->children_) 
        child->accept(this);
}

} // namespace TimeloopX

} // namespace mapping 