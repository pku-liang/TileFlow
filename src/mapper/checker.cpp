#include <stack>

#include "tileflow/mapper/checker.hpp"

namespace TileFlow {

void ShapeConstraintParser::visitTile(const TileNode* node){
    for (auto loop: node->get_loops()) {
        TILEFLOW_ASSERT(
            problem::GetShape()->FactorizedDimensionIDToName.at(loop.dimension) == loop.name_,
            "mismatch between " <<  problem::GetShape()->FactorizedDimensionIDToName.at(loop.dimension) 
            << " and " << loop.name_);
        auto& expr = exprs_[loop.dimension];
        if (loop.end > 0) {
            expr.first *= loop.end;
        }
        else {
            expr.second.push_back(
                std::make_shared<VariableExpr>(loop.end));
        }
    }
    for (auto child:node->get_children()) child->accept(this);
    for (auto loop: node->get_loops()) {
        auto& expr = exprs_[loop.dimension];
        if (loop.end > 0) {
            expr.first /= loop.end;
        }
        else {
            expr.second.pop_back();
        }
    }
}

void ShapeConstraintParser::visitOp(const OpNode* node) {
    for(auto& kv: node->get_workload()->GetShape()->FactorizedDimensionNameToID){
        int dim = problem::GetShape()->FactorizedDimensionNameToID.at(kv.first);
        int scale = workload_.GetFactorizedBound(dim);
        auto & expr = exprs_[dim];
        assert(scale % expr.first == 0);
        if (!expr.second.size()) {
            TILEFLOW_ASSERT(scale == expr.first, "At"; node->display("", false); 
            std::cout <<  "mismatch for " << kv.first << ": " 
            << scale << "v.s." << expr.first);
        }
        else {
            constraints.push_back({std::make_shared<CondExpr>(
                std::make_shared<ProductExpr>(expr.second),
                std::make_shared<ParameterExpr>(scale / expr.first), 
                CondExpr::EQU), 
            " loopcount constraint for tiling of dim " + kv.first});
        }
    }
}

std::vector<Constraint> ShapeConstraintParser::parse(const Node* root) {
    exprs_.clear();
    exprs_.resize(problem::GetShape()->NumFlattenedDimensions);
    constraints.clear();
    for (auto& x: exprs_) x.first = 1; 
    root->accept(this);
    return constraints;
}

void SpatialScopeSwapper::visitScope(const ScopeNode* node){
    auto parent = node->get_parent();
    if (parent && parent->get_type() == Node::Tile) {
        auto parent_ = static_cast<const TileNode*>(parent);
        if (parent_->get_tile_type() == TileNode::Spatial) {
            TILEFLOW_ASSERT(node->get_scope_type() == ScopeNode::Sequential
            || node->get_scope_type() == ScopeNode::Sharing, 
            "A spatial tile's child cannot be a parallel/pipeline scope");
            
            auto node_ = const_cast<ScopeNode*>(node);
            
            std::vector<const Node*> new_nodes;
            for (auto child: node_->get_children()) {
                Node * new_node = new TileNode(*parent_);
                assert(new_node);
                new_nodes.push_back(new_node);
                new_node->reset_children();
                new_node->set_parent(node_);
                new_node->add_child(child);
            }

            node_->set_children(new_nodes);
            node_->set_parent(parent_->get_parent());

            if (parent_->get_parent()) {
                const_cast<Node*>(parent_->get_parent())->replace_child(parent_, node_);
            }

            // erase all child for safe deconstruct 
            const_cast<TileNode*>(parent_)->reset_children();
            delete parent_;
        }
    }
    for (auto child: node->get_children())
        child->accept(this);
}

void Checker::check(){
    swap_spatial_scope();
    sanity_check();
    get_active_tensors();
    get_shape_constraints();
    // this depends on get_active_tensors;
    get_memory_constraints();
}

void Checker::get_active_tensors(){
    std::vector<const OpNode *> opnodes = 
    mapping::TileFlow::CollectOpNode().collectOpNodes(mapping_.root);
    std::unordered_map<std::string, const Node *> tensor2producer;
    for (auto &t : workloads_.get_ins()){
        tensor2producer[t] = mapping_.root;
    }

    std::unordered_map<const Node *, std::vector<problem::Shape::DataSpaceID>> access_pattern;

    for (auto node : opnodes)
    {
        auto &ptr = node->get_workload();
        for (auto &t : ptr->get_ins())
        {
            TILEFLOW_ASSERT(tensor2producer.count(t), 
                "Op " << node->get_name() << "'s input " << t << " is unclear");
            auto producer = tensor2producer[t];
            problem::Shape::DataSpaceID tensor_id = workloads_.get_shape().DataSpaceNameToID.at(t);
            problem::Shape::DataSpaceID producer_id = producer->get_type() == Node::Op ? 
                tensor_id:problem::Shape::DataSpaceID(-1);
            if (verbose_level) {
                if (producer->get_type() == Node::Op) {
                    std::cout << node->get_name() << " consumes " <<
                         static_cast<const OpNode*>(producer)->get_name() << "'s " << t << std::endl;
                }
                else std::cout << node->get_name() << " consumes global input " << t << std::endl;
            }
            add_access_pattern(producer_id, producer, tensor_id, node);
        }
        tensor2producer[ptr->get_out()] = node;
    }

    for (auto &t : workloads_.get_outs())
    {
        TILEFLOW_ASSERT(tensor2producer.count(t) && tensor2producer[t]->get_type() == Node::Op, "Output " << t << " is unclear");
        int id = workloads_.get_shape().DataSpaceNameToID.at(t);
        add_access_pattern(id, tensor2producer[t], problem::Shape::DataSpaceID(-1), mapping_.root);
        if (verbose_level) {
            std::cout << "Output uses " << static_cast<const OpNode*>(tensor2producer[t])->get_name() << "'s " << t << std::endl;
        }
    }
}

void Checker::add_access_pattern(
            problem::Shape::DataSpaceID producer_id,
            const Node *producer,
            problem::Shape::DataSpaceID consumer_id,
            const Node *consumer)
{
    
    std::stack<const Node *> sp, sc;
    while (producer)
    {
        sp.push(producer);
        producer = producer->get_parent();
    }
    while (consumer)
    {
        sc.push(consumer);
        consumer = consumer->get_parent();
    }
    std::stack<const Node *> common;
    while (!sc.empty() && !sp.empty() && sc.top() == sp.top())
    {
        common.push(sc.top());
        sc.pop();
        sp.pop();
    }

    while (!common.empty()) {
        auto node = common.top();
        common.pop();
        sc.push(node); sp.push(node);
        if (node->get_type() == Node::type_t::Tile && 
            static_cast<const TileNode *>(node)->get_tile_type() == TileNode::Temporal)
            break;
    }

    if (producer_id != problem::Shape::DataSpaceID(-1)) {
        while (!sp.empty()) {
            auto node = sp.top();
            sp.pop();
            if (problem::GetShape()->IsReadWriteDataSpace.at(producer_id)) {
                node->get_active_tensors().read_tensors.push_back(producer_id);
                if (!sp.empty()) {
                    sp.top()->get_active_tensors().fill_tensors.push_back(producer_id);
                }
            }
            node->get_active_tensors().update_tensors.push_back(producer_id);
            if (!sp.empty()) {
                sp.top()->get_active_tensors().wb_tensors.push_back(producer_id);
            }
        }
    }

    if (consumer_id != problem::Shape::DataSpaceID(-1)) {
        while (!sc.empty()) {
            auto node = sc.top();
            sc.pop();
            node->get_active_tensors().read_tensors.push_back(consumer_id);
            if (!sc.empty()) {
                sc.top()->get_active_tensors().fill_tensors.push_back(consumer_id);
            }
        }
    }
}

void Checker::display() {
    std::cout << "==============Checker BEG================" << std::endl;
    mapping_.Print();
    std::cout << "constraints:" << std::endl;
    for (auto& cons: constraints) {
        std::cout << "\t";
        
        cons.expr->display(global_symbol_table_); 
        if (verbose_level>1)
            std::cout << "\t# " << cons.msg;
        std::cout << std::endl;}
    std::cout << "==============Checker END================" << std::endl;
}

void Checker::swap_spatial_scope() {
    SpatialScopeSwapper pass;
    pass.run(mapping_.root);
}

void Checker::sanity_check() {
    SanityChecker checker(topology_);
    checker.run(mapping_.root);
}

void Checker::get_shape_constraints() {
    ShapeConstraintParser parser(workloads_.get_workload());
    auto cons = parser.parse(mapping_.root);
    constraints.insert(constraints.end(), cons.begin(), cons.end());
}

void Checker::get_memory_constraints() {
    MemoryConstraintParser parser(workloads_.get_workload(), topology_);
    auto cons = parser.parse(mapping_.root);
    constraints.insert(constraints.end(), cons.begin(), cons.end());
}

std::shared_ptr<Expr> MemoryConstraintParser::cal_footprint(unsigned pv) {
    std::vector<std::shared_ptr<Expr>> prod;
    for (unsigned data_space_dim = 0; data_space_dim < problem::GetShape()->DataSpaceOrder.at(pv); 
    data_space_dim++)
    {
        std::vector<std::shared_ptr<Expr> > sum;
        for (auto& term : problem::GetShape()->Projections.at(pv).at(data_space_dim))
        {
            int coeff = 1;
            if (term.first != problem::GetShape()->NumCoefficients)
            {
                // If Coefficient is negative, flip high/low.
                coeff = workload_.GetCoefficient(term.first);
                if (coeff < 0) coeff = -coeff;
            }
            auto & factor = factors_[term.second];
            coeff *= factors_[term.second].first;
            std::vector<std::shared_ptr<Expr> > exprs;
            if (coeff != 1) {
                exprs.push_back(std::make_shared<ParameterExpr>(coeff));
            }
            if (factor.second.size() == 1) {
                exprs.push_back(std::make_shared<VariableExpr>(factor.second.front()));
            }
            else {
                for (auto id: factor.second)
                    exprs.push_back(std::make_shared<VariableExpr>(id));
            }
            if (exprs.size() == 1) sum.push_back(exprs.front());
            else sum.push_back(std::make_shared<ProductExpr>(exprs));
        }
        if (sum.size() == 1) prod.push_back(sum.front());
        else prod.push_back(std::make_shared<SumExpr>(sum));
    }
    if (prod.size() == 1) return prod.front();
    return std::make_shared<ProductExpr>(prod);
}

void MemoryConstraintParser::visitTile(
    const TileNode* node) {
    for (auto child: node->get_children())
        child->accept(this);
    for (auto& loop: node->get_loops()) {
        auto& factor = factors_[loop.dimension];
        if (loop.end > 0) 
            factor.first *= loop.end;
        else factor.second.push_back(loop.end);
    }

    if (node->get_tile_type() == TileNode::Temporal && 
    topology_.GetStorageLevel(node->get_storage_level())->GetSpecs().size.IsSpecified()) {
        std::set<int> active_tensors;
        active_tensors.insert(node->get_active_tensors().read_tensors.begin(), 
            node->get_active_tensors().read_tensors.end());
        active_tensors.insert(node->get_active_tensors().update_tensors.begin(), 
            node->get_active_tensors().update_tensors.end());
        std::vector<std::shared_ptr<Expr> > operands;
        for (auto pv: active_tensors)
            operands.push_back(cal_footprint(pv));
        
        std::shared_ptr<Expr> footprints = operands.front();
        if (operands.size() > 1) {
            footprints = std::make_shared<SumExpr>(operands); 
        }
        constraints_.push_back({std::make_shared<CondExpr>(
            footprints,
            std::make_shared<ParameterExpr>(
                (int)topology_.GetStorageLevel(node->get_storage_level())->GetSpecs().size.Get()),
            CondExpr::LEQ
        ), "Tile node memory constraint at " + node->get_storage_name()});
    }
}

void MemoryConstraintParser::visitOp(
    const OpNode* 
){
    factors_.clear();
    factors_.resize(problem::GetShape()->NumFlattenedDimensions);
    for (auto& factor: factors_) factor.first = 1;
}

void MemoryConstraintParser::visitScope(
    const ScopeNode* node 
) {
    for (auto child: node->get_children()){
        child->accept(this);
    }
    if (node->get_scope_type() == ScopeNode::Sharing
    && topology_.GetStorageLevel(node->get_storage_level())->GetSpecs().size.IsSpecified()) {
        std::set<int> active_tensors;
        active_tensors.insert(node->get_active_tensors().read_tensors.begin(), 
            node->get_active_tensors().read_tensors.end());
        active_tensors.insert(node->get_active_tensors().update_tensors.begin(), 
            node->get_active_tensors().update_tensors.end());
        std::vector<std::shared_ptr<Expr> > operands;
        for (auto pv: active_tensors)
            operands.push_back(cal_footprint(pv));
        
        std::shared_ptr<Expr> footprints = operands.front();
        if (operands.size() > 1) {
            footprints = std::make_shared<SumExpr>(operands); 
        }
        constraints_.push_back({std::make_shared<CondExpr>(
            footprints,
            std::make_shared<ParameterExpr>(
                (int)topology_.GetStorageLevel(node->get_storage_level())->GetSpecs().size.Get()),
            CondExpr::LEQ
        ), "Sharing Scope constraint at " + node->get_storage_name()});
    }
}

std::vector<Constraint> MemoryConstraintParser::parse(const Node* root){
    root->accept(this);
    return constraints_;
} 

void SanityChecker::run(const Node* root){
    storage_level_ = topology_.NumStorageLevels();
    root->accept(this);
}

void SanityChecker::visitTile(const TileNode* node){
    auto old_storage_level = storage_level_;
    TILEFLOW_ASSERT(node->get_tile_type() != TileNode::Spatial || 
        node->get_storage_level() == storage_level_, 
        "bad spatial tile level at " << node->get_storage_name() 
        << ", should be the same with upper temporal node.");
    TILEFLOW_ASSERT(node->get_tile_type() != TileNode::Temporal ||
        node->get_storage_level() ==  (storage_level_) || 
        node->get_storage_level() ==  (--storage_level_), 
        "skip or reverse storage level mapping at " << node->get_storage_name()
        << "," << node->get_storage_level() << ":" << storage_level_);
    auto& children = node->get_children();
    TILEFLOW_ASSERT(children.size() == 1, "Tile node should have only child.");
    auto child = children.front();
    child->accept(this);
    storage_level_ = old_storage_level;

    if (node->get_tile_type() == TileNode::Spatial){
        TILEFLOW_ASSERT(child->get_type() == Node::Tile, 
        "Spatial's child must be temporal tile.");
        TILEFLOW_ASSERT(static_cast<const TileNode*>(child)->get_tile_type() == TileNode::Temporal, 
        "Spatial's child must be temporal tile.");
    }

}

void SanityChecker::visitScope(const ScopeNode* node){
    auto& children = node->get_children();
    TILEFLOW_ASSERT(children.size() > 1, "Scope node should have more than one child.");
    for (auto child: children) child->accept(this);
}

void SanityChecker::visitOp(const OpNode* ) {
    TILEFLOW_ASSERT(storage_level_ == 0, 
    " missing temporal tiles for storage level under " << storage_level_);
} 

} // namespace TileFlow 