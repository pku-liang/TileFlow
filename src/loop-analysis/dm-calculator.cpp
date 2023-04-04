#include <stack>

#include "tileflow/loop-analysis/nest-analysis.hpp"

namespace analysis
{
    namespace TileFlow
    {

        RetVal DatamovementCalculator::eval(const Node *root)
        {
            break_on_failure = false;
            MemoryState::set_workload(&workload_);
            InputParam input;
            input.num_epochs_ = 1;
            for (unsigned dim = 0;
                 dim < problem::GetShape()->NumFlattenedDimensions; dim++)
            {
                input.cur_transform_[dim] = 0;
            }
            input.curr_node_ = root;
            input_stack_.push(input);
            root->accept(this);
            assert(ret_stack_.size() == 1);
            return ret_stack_.top();
        }

        void DatamovementCalculator::visitOp(const OpNode *node)
        {
            auto input = input_stack_.top();
            assert(input.curr_node_ == node);
            input_stack_.pop();
            // Add a mask here.
            RetVal ret;
            auto &active_tensors = analysis_.configs.at(node).active_read_tensors;
            problem::OperationSpace point_set(MemoryState::workload_,
                                              input.cur_transform_, input.cur_transform_); // A single point
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
            {
                if (find(active_tensors.begin(), active_tensors.end(), pv) == active_tensors.end())
                    point_set.GetDataSpace(pv).Reset();
                else {
                    ret.access_stat_[pv](1,1).accesses = point_set.GetSize(pv);
                }
            }
            ret.last_working_set_.insert(0, point_set);
            // The compute level does not have data persistance.
            ret.deltas_ = ret.last_working_set_; // - input.init_working_set_;

            // std::cout << "--------Op---------" << std::endl;
            // std::cout << input << std::endl;
            // std::cout << ret << std::endl;
            // std::cout << "-------------------" << std::endl;
            tiling::CompoundTile tile;
            auto & compute_info = tile.compute_info;
            compute_info.accesses = input.num_epochs_;
            compute_info.replication_factor = analysis_.configs[node].replication_factor;
        
            for (unsigned op_id = 0; op_id < tiling::arithmeticOperationTypes.size();
                 op_id++)
            {
                auto op_name = tiling::arithmeticOperationTypes[op_id];
                compute_info.fine_grained_accesses[op_name] = 0;
            }
            compute_info.fine_grained_accesses["random_compute"] 
                = compute_info.accesses * compute_info.replication_factor;
            compute_info.avg_replication_factor = compute_info.replication_factor;

            auto &acc_compute_info = analysis_.configs[node].stats_[input.space_stamp_].compute_info_;
            acc_compute_info.accesses += compute_info.accesses;
            acc_compute_info.replication_factor = compute_info.replication_factor;

            tiling::CompoundMask mask = {};
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; ++pv) 
                mask[pv] = true;
            auto level = topology_.GetArithmeticLevel()->Clone();
            level->Evaluate(tile, mask, 0, tile.compute_info.accesses, break_on_failure);
            
            ret.cycle_ = level->Cycles();
            ret_stack_.push(ret);
        }

        void DatamovementCalculator::visitTile(const TileNode *node)
        {
            auto &input = input_stack_.top();
            assert(input.curr_node_ == node);
            analysis::TileFlow::PerfectLoopnestAnalyzer analyzer(*this, input, analysis_.configs[node]);
            auto &config = analysis_.configs[node];
            analyzer.init(&analysis_.common_workload_, &config.loop_nest);
            auto ret = analyzer.calculateDataMovement();
            if (input.num_epochs_ && (node->get_parent() == nullptr || 
                node->get_parent()->get_type() != Node::Tile))
                finalizeStat(node, ret);
            ret_stack_.push(ret);
        }

        void DatamovementCalculator::visitScope(const ScopeNode *node)
        {
            auto input = input_stack_.top();
            assert(input.curr_node_ == node);
            input_stack_.pop();
            auto type = node->get_scope_type();
            auto ret_ = RetVal();

            for (auto& child: node->get_children()) {
                input.curr_node_ = child;
                input_stack_.push(input);
                child->accept(this);
                auto ret = ret_stack_.top();
                ret_stack_.pop();
                if (type == ScopeNode::Pipeline || type == ScopeNode::Sequential)
                    input.init_working_set_ = ret.last_working_set_;
                else 
                    ret_.last_working_set_.Add(ret.last_working_set_);
                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces;
                ++pv) {
                    ret_.access_stat_[pv].Accumulate(ret.access_stat_[pv]);
                }
                if (type == ScopeNode::Sharing || type == ScopeNode::Sequential)
                    ret_.cycle_ += ret.cycle_;
                else ret_.cycle_ = std::max(ret.cycle_, ret_.cycle_);
            }

            if (type == ScopeNode::Pipeline || type == ScopeNode::Sequential)
                ret_.last_working_set_ = input.init_working_set_;
            ret_stack_.push(ret_);
        }
 
        RetVal DatamovementCalculator::computeDelta(const InputParam &input)
        {
            input_stack_.push(input);
            input.curr_node_->accept(this);
            auto ret = ret_stack_.top();
            ret_stack_.pop();
            std::cout << "======ComputeDelta========" << std::endl;
            std::cout << "------Node--------" << std::endl;
            input.curr_node_->display("", false);
            std::cout << "------Input-------" << std::endl;
            std::cout << input;
            std::cout << "------Output------" << std::endl;
            std::cout << ret;
            std::cout << "============END===========" << std::endl;

            return ret;
        }

        void DatamovementCalculator::finalizeStat(const Node* node, RetVal& ret) {
            auto storage_id = node->get_storage_level();
            auto storage_level = std::static_pointer_cast<model::BufferLevel>(
                topology_.GetStorageLevel(storage_id)->Clone());
            tiling::CompoundMask mask = {};
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; ++pv) 
                mask[pv] = true;
            storage_level->Evaluate(*ret.p_tile_, mask, 
                0, 
                ret.cycle_, break_on_failure);
            storage_level->FinalizeBufferEnergy();
            ret.cycle_ = storage_level->Cycles();

            auto connection = topology_.connection_map_.at(storage_id);
            auto rf_net = connection.read_fill_network->Clone();
            rf_net->Evaluate(*ret.p_tile_, break_on_failure);
            auto du_net = connection.drain_update_network->Clone();
            du_net->Evaluate(*ret.p_tile_, break_on_failure);
            assert(ret.p_tile_.unique());
            ret.p_tile_.reset();

             std::cout << "Storage<" << storage_id << ">:" << std::endl 
                << *storage_level; 
            std::cout << "Connect<" << storage_id << ">:";
            std::cout << "rf_net: " << rf_net->Energy();
            std::cout << "du_net: " << du_net->Energy();
            std::cout << std::endl;
        }

    } // namespace TileFlow
} // namespaec analysis