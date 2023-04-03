
#include "tileflow/model/topology.hpp" 

namespace model {

namespace TileFlow {

    void Topology::eval(
        const mapping::TileFlow::Mapping& mapping, 
        const analysis::TileFlow::NestAnalysis& analysis){
        
        StatCalculator pass(*this, mapping, analysis);
        stats_.cycles = total_network_latency_;
        pass.run(mapping.root);

        std::cout << "Energy: " << stats_.energy << std::endl;
        std::cout << "Cycles: " << stats_.cycles << std::endl;
    }

    void StatCalculator::visitTile(const TileNode* node) {
        for (auto child: node->get_children()) 
            child->accept(this);
        if (!node->is_spatial()) {
            auto cycle = cycles_.top();
            cycles_.pop();
            auto& tile = analysis_.get_tile(node);
            auto storage_id = node->get_storage_level();
            auto storage_level = std::static_pointer_cast<model::BufferLevel>(topology_.GetStorageLevel(storage_id)->Clone());
            tiling::CompoundMask mask = {};
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; ++pv) 
                mask[pv] = true;
            storage_level->Evaluate(tile, mask, 
                0, 
                cycle, break_on_failure);
            storage_level->FinalizeBufferEnergy();
            cycles_.push(storage_level->Cycles());
            energy_ += storage_level->Energy();

            auto connection = topology_.connection_map_[storage_id];
            auto rf_net = connection.read_fill_network->Clone();
            rf_net->Evaluate(tile, break_on_failure);
            energy_ += rf_net->Energy(); 
            auto du_net = connection.drain_update_network->Clone();
            du_net->Evaluate(tile, break_on_failure);
            energy_ += du_net->Energy();

             std::cout << "Storage<" << storage_id << ">:" << std::endl 
                << *storage_level; 
            std::cout << "Connect<" << storage_id << ">:";
            std::cout << "rf_net: " << rf_net->Energy();
            std::cout << "du_net: " << du_net->Energy();
            std::cout << std::endl;
        }
    }

    void StatCalculator::visitScope(const ScopeNode* node) {

        auto type = node->get_scope_type();
        std::uint64_t cycle;
        if (type == ScopeNode::Sequential || type == ScopeNode::Sharing) {
            cycle = 0;
            for (auto child: node->get_children()) {
                child->accept(this);
                cycle += cycles_.top();
                cycles_.pop();
            }
        }
        else if (type == ScopeNode::Parallel || type == ScopeNode::Pipeline) {
            cycle = 0;
            for (auto child: node->get_children()) {
                child->accept(this);
                cycle = std::max(cycle, cycles_.top());
                cycles_.pop();
            }    
        }
        cycles_.push(cycle);
    }

    void StatCalculator::visitOp(const OpNode* node) {
        auto level = topology_.GetArithmeticLevel()->Clone();
        auto &tile = analysis_.get_tile(node);
        tiling::CompoundMask mask = {};
        for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; ++pv) 
            mask[pv] = true;
        level->Evaluate(tile, mask, 0, 
            tile.compute_info.accesses, break_on_failure);
        cycles_.push(level->Cycles());
        energy_ += level->Energy();
        std::cout << "Arithmetic::" << node->get_name() << ":" << level->Energy() << std::endl;
    }

    void StatCalculator::run(const Node* root) {
        break_on_failure = false;
        energy_ = 0.0;
        root->accept(this);
        assert(cycles_.size() == 1);
        topology_.stats_.energy = energy_;
        topology_.stats_.cycles += cycles_.top();
    }

} // namespace TileFlow

} // namespace model