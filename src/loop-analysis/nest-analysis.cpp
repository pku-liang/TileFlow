#include <stack>

#include "tileflow/loop-analysis/nest-analysis.hpp"

namespace analysis {

namespace TileFlow {

    std::vector<const OpNode*> CollectOpNode::collectOpNodes(Node* root) {
        opnodes_.clear();
        root->accept(this);
        return std::move(opnodes_);
    }

    void CollectOpNode::visitOp(const OpNode* node) {
        opnodes_.push_back(node);
    }

    void NestAnalysis::get_tilewise_workloads(){
        std::vector<const OpNode*> opnodes = CollectOpNode().collectOpNodes(mapping_.root);
        std::cout << "Ops: " << std::endl;
        for (auto node: opnodes) {
            node->display("", false);
        }
        std::unordered_map<std::string, const Node*> tensor2producer;
        for (auto& t: workloads_.get_ins()) {
            tensor2producer[t] = mapping_.root;
        }

        std::unordered_map<const Node*, std::vector<problem::Shape::DataSpaceID> > access_pattern;

        for (auto node: opnodes) {
            auto & ptr = node->get_workload();
            for (auto& t: ptr->get_ins()){
                TILEFLOW_ASSERT(tensor2producer.count(t), "Op " << node->get_name() << "'s input " << t << " is unclear");
                auto producer = tensor2producer[t];
                problem::Shape::DataSpaceID producer_id = producer->get_type() == Node::Op? 
                    workloads_.get_shape().DataSpaceNameToID.at(static_cast<const OpNode*>(producer)->get_name() + "::" + t):
                    problem::Shape::DataSpaceID(-1);
                problem::Shape::DataSpaceID consumer_id = workloads_.get_shape().DataSpaceNameToID.at(node->get_name() + "::" + t);
                add_access_pattern(producer_id, producer, consumer_id, node, access_pattern);
            }
            tensor2producer[ptr->get_out()] = node;
        }

        for (auto& t: workloads_.get_outs()) {
            TILEFLOW_ASSERT(tensor2producer.count(t), "Output " << t << " is unclear");
            auto producer_op_name = static_cast<const OpNode*>(tensor2producer[t])->get_name();
            int id = workloads_.get_shape().DataSpaceNameToID.at(producer_op_name + "::" + t);
            add_access_pattern(id, tensor2producer[t], problem::Shape::DataSpaceID(-1), mapping_.root, access_pattern);
        }

        for (auto& kv: access_pattern) {
            problem::Shape shape(workloads_.get_shape());
            shape.NumDataSpaces = kv.second.size();
            for (auto iter = shape.DataSpaceIDToName.begin(); iter != shape.DataSpaceIDToName.end(); ){
                if (find(kv.second.begin(), kv.second.end(), iter->first) == kv.second.end()){
                    shape.DataSpaceNameToID.erase(iter->second);
                    iter = shape.DataSpaceIDToName.erase(iter);
                }
                else {iter ++;}
            }
            configs[kv.first].workload.SetShape(shape);
        }
    }

    void NestAnalysis::add_access_pattern (
        problem::Shape::DataSpaceID producer_id, 
        const Node* producer, 
        problem::Shape::DataSpaceID consumer_id,
        const Node* consumer, 
        std::unordered_map<const Node*, std::vector<problem::Shape::DataSpaceID> >& access_pattern) {
        std::stack<const Node*> sp, sc;
        while (producer) {
            sp.push(producer);
            producer = producer->get_parent();
        }
        while (consumer) {
            sc.push(consumer);
            consumer = consumer->get_parent();            
        }
        std::stack<const Node*> common;
        while(!sc.empty() && !sp.empty() && sc.top() == sp.top()) {
            common.push(sc.top());
            sc.pop(); sp.pop();
        }
        while (!common.empty()){
            auto node = common.top();
            sc.push(node); sp.push(node);
            common.pop();
            if (node->get_type() == Node::type_t::Tile && !static_cast<const TileNode*>(node)->is_spatial())
                break;
        }   
        
        while (!sc.empty()){
            auto node = sc.top();
            access_pattern[node].push_back(consumer_id);
            sc.pop();
        }
        
        while (!sp.empty()) {
            auto node = sp.top();
            access_pattern[node].push_back(producer_id);
            sp.pop();
        }
    }

    void NestAnalysis::Print() {
        std::cout << "-----------------Nest Analysis----------------" << std::endl;
        Displayer(*this).display();
        std::cout << "--------------END Nest Analysis---------------" << std::endl;
    }

    void Displayer::visitTile(const TileNode* node) {
        std::cout << prefix_ << "Tile:";
        auto& configs_ = analysis_.configs;
        if (configs_.count(node)) {
            auto shape = configs_.at(node).workload.GetShape();
            for (auto t: shape->DataSpaceNameToID) {
                std::cout << t.first << ",";
            }
            std::cout << std::endl;
        }
        node->display(prefix_, false);
        auto old_prefix = prefix_;
        for (int i = 0; i < (int)node->n_level(); ++i) 
            prefix_ += "   ";
        for (auto child: node->get_children()) 
            const_cast<Node*>(child)->accept(this);
        prefix_ = old_prefix;
    }
    void Displayer::visitScope(const ScopeNode*node) {
        node->display(prefix_, false);
        std::cout << "{" << std::endl;
        std::string old_prefix = prefix_;
        prefix_ += "   ";
        for (auto child: node->get_children()) 
            const_cast<Node*>(child)->accept(this);
        prefix_ = old_prefix;
        std::cout << prefix_ << "}" << std::endl;
    }
    void Displayer::visitOp(const OpNode* node) {
        node->display(prefix_, false);
    }

    void DatamovementCalculator::visitOp(const OpNode* node){}

    void DatamovementCalculator::visitTile(const TileNode* node){
        analysis::TileFlow::PerfectLoopnestAnalyzer analyzer(*this);
        auto& config = analysis_.configs[node];
        auto loopnest = node->constructLoopNest(analysis_.workloads_.get_shape().FactorizedDimensionNameToID);
        analyzer.init(&config.workload, &loopnest, analysis_.mapping_.fanoutX_map, analysis_.mapping_.fanoutY_map);
        // analyzer.calculateDataMovement();
    }

    void DatamovementCalculator::visitScope(const ScopeNode* node){
        std::vector<problem::OperationSpace> deltas;
        for (auto& child: node->get_children()) {
            child->accept(this);
            deltas.emplace_back(deltas_.top());
            deltas_.pop();
        }
        deltas_.push(combineDeltas(deltas, node->type));
    }

    problem::OperationSpace DatamovementCalculator::combineDeltas(
        const std::vector<problem::OperationSpace>& deltas, 
        ScopeNode::type_t type){
        if (type == ScopeNode::Parallel || type == ScopeNode::Pipeline){
            
        }
        else if (type == ScopeNode::Sequential) {

        }
    }

    problem::OperationSpace DatamovementCalculator::computeDelta(Node* node){
        node->accept(this);
        assert(!deltas_.empty());
        auto ret = deltas_.top();
        deltas_.pop();
        return ret;
    }

    void PerfectLoopnestAnalyzer::init(problem::Workload*wl, loop::Nest*nest, 
            std::map<unsigned, std::uint64_t> fanoutX_map,
            std::map<unsigned, std::uint64_t> fanoutY_map){
        Init(wl, nest, fanoutX_map, fanoutY_map);
        num_epochs_ = dm.screen_shot_.num_epochs_;
    }

    void PerfectLoopnestAnalyzer::calculateDataMovement() {
        if (nest_state_.size() != 0)
        {
            InitializeNestProperties();
            InitializeLiveState();
            DetectImperfectFactorization();

            // Recursive call starting from the last element of the list.
            dm.deltas_.push(ComputeDeltas(nest_state_.rbegin()));
            CollectWorkingSets();
        }

        // Done.
        working_sets_computed_ = true;
    }

} // namespace TileFlow 

} // namespace analysis 