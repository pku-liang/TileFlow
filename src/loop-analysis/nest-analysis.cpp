#include <stack>

#include "timeloopX/loop-analysis/nest-analysis.hpp"

namespace analysis {

namespace TimeloopX {

    std::vector<OpNode*> CollectOpNode::collectOpNodes(Node* root) {
        opnodes_.clear();
        root->accept(this);
        return std::move(opnodes_);
    }

    void CollectOpNode::visitOp(OpNode* node) {
        opnodes_.push_back(node);
    }

    void NestAnalysis::get_alive_tensors(){
        std::vector<OpNode*> opnodes = CollectOpNode().collectOpNodes(mapping_.root);
        std::cout << "Ops: " << std::endl;
        for (auto node: opnodes) {
            node->display("", false);
        }
        std::unordered_map<std::string, Node*> tensor2producer;
        for (auto& t: workloads_.get_ins()) {
            tensor2producer[t] = mapping_.root;
        }

        for (auto node: opnodes) {
            auto & ptr = node->get_workload();
            for (auto& t: ptr->get_ins()){
                TIMELOOPX_ASSERT(tensor2producer.count(t), "Op " << node->get_name() << "'s input " << t << " is unclear");
                parse_alive_tensors(t, tensor2producer[t], node);
            }
            tensor2producer[ptr->get_out()] = node;
        }

        for (auto& t: workloads_.get_outs()) {
            TIMELOOPX_ASSERT(tensor2producer.count(t), "Output " << t << " is unclear");
            parse_alive_tensors(t, tensor2producer[t], mapping_.root);
        }
    }

    void NestAnalysis::parse_alive_tensors (const std::string & tensor, const Node* producer, const Node* consumer) {
        int id = workloads_.tensor2id(tensor);
        std::stack<const Node*> sp, sc;
        while (producer) {
            sp.push(producer);
            producer = producer->get_parent();
        }
        while (consumer) {
            sc.push(consumer);
            consumer = consumer->get_parent();            
        }
        const Node* root = nullptr;
        while(sc.top() == sp.top()) {
            root = sc.top();
            sc.pop(); sp.pop();
        }
        assert(root);
        sc.push(root);
        while (!sp.empty()) {
            sc.push(sp.top());
            sp.pop();
        }

        while(!sc.empty()) {
            auto node = sc.top();
            sc.pop();
            configs[node].alive_tensors.set(id);
        }
    }

    void NestAnalysis::Print() {
        std::cout << "-----------------Nest Analysis----------------" << std::endl;
        Displayer(*this).display();
        std::cout << "--------------END Nest Analysis---------------" << std::endl;
    }

    void Displayer::visitTile(TileNode* node) {
        std::cout << prefix_ << "Tile:";
        auto& configs_ = analysis_.configs;
        if (configs_.count(node)) {
            auto& config = configs_.at(node);
            auto& tensors = analysis_.workloads_.get_tensors();
            if (config.alive_tensors.any()) {
                std::cout << "alive tensors:";
                for (int i = 0; i < (int)tensors.size(); ++i) {
                    if (config.alive_tensors.test(i)) {
                        std::cout << tensors[i] << ",";
                    }
                }
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
    void Displayer::visitScope(ScopeNode*node) {
        node->display(prefix_, false);
        std::cout << "{" << std::endl;
        std::string old_prefix = prefix_;
        prefix_ += "   ";
        for (auto child: node->get_children()) 
            const_cast<Node*>(child)->accept(this);
        prefix_ = old_prefix;
        std::cout << prefix_ << "}" << std::endl;
    }
    void Displayer::visitOp(OpNode* node) {
        node->display(prefix_, false);
    }

} // namespace TimeloopX 

} // namespace analysis 