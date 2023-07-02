#include <stack>

#include "tileflow/loop-analysis/nest-analysis.hpp"

namespace analysis
{
    namespace TileFlow
    {

        void DataMovements::report(std::ostream& o) const {
            for(auto& kv: data_) {
                kv.second.report(o, kv.first);
            }
        }
    
        void DataMovement::report(std::ostream&o, const std::string& prefix) const {
            for (auto& kv: data_) {
                o << prefix << "::" << kv.first << ", " << kv.second;
                o << std::endl;
            }
        }

        RetVal DatamovementCalculator::eval(const Node *root)
        {
            break_on_failure = false;
            energy_ = 0.0;
            MemoryState::set_workload(&workload_);
            InputParam input;
            input.num_epochs_ = 1;
            assert(input.cur_transform_.Order() == problem::GetShape()->NumFlattenedDimensions);
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
            assert(input_stack_.size() == 1);
            auto input = input_stack_.top();
            assert(input.curr_node_ == node);
            input_stack_.pop();
            // Add a mask here.
            RetVal ret;
            if (input.num_epochs_ == 0) {
                ret_stack_.push(ret);
                return;
            }

            auto &fill_tensors = node->get_active_tensors().fill_tensors;
            auto &wb_tensors = node->get_active_tensors().wb_tensors;
            problem::OperationSpace point_set(MemoryState::workload_,
                                              input.cur_transform_, input.cur_transform_); // A single point
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
            {
                if (find(fill_tensors.begin(), fill_tensors.end(), pv) == fill_tensors.end() && 
                find(wb_tensors.begin(), wb_tensors.end(), pv) == wb_tensors.end())
                    point_set.GetDataSpace(pv).Reset();
                else {
                    ret.access_stat_[pv](1,1).accesses = point_set.GetSize(pv) * input.num_epochs_;
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
            energy_ += level->Energy();
            ret_stack_.push(ret);

            auto& data_movement_ = analysis_.data_movements_[level->Name()];
            data_movement_["Flops"] += level->Cycles() * compute_info.replication_factor; 
            data_movement_["Energy"] += level->Energy(); 
            
            if (verbose_level > 1) {
                std::cout << "========BEG Compute Stat=========" << std::endl;
                std::cout << input;
                std::cout << ret;
                level->Print(std::cout);
                std::cout << "========END Compute Stat=========" << std::endl;
            }
        }

        void DatamovementCalculator::visitTile(const TileNode *node)
        {
            assert(input_stack_.size() == 1);
            auto input = input_stack_.top();
            input_stack_.pop();
            assert(input.curr_node_ == node);
            analysis::TileFlow::PerfectLoopnestAnalyzer analyzer(*this, input, analysis_.configs[node]);
            auto &config = analysis_.configs[node];
            analyzer.init(&analysis_.common_workload_, &config.loop_nest);
            auto ret = analyzer.calculateDataMovement();
            // if current node is root or its parent is of another storage level 
            // we finalize the stat 
            if (input.num_epochs_ && (node->get_parent() == nullptr || 
                node->get_parent()->get_type() != Node::Tile || 
                static_cast<const TileNode*>(node->get_parent())->get_tile_type() != TileNode::Spatial)){
                finalizeStat(node->is_spatial()? 
                    node->get_children().front()->get_storage_level():node->get_storage_level(),
                     ret, node->is_profile());
            }
            ret_stack_.push(ret);
        }

        void DatamovementCalculator::visitScope(const ScopeNode *node)
        {
            auto input = input_stack_.top();
            assert(input.curr_node_ == node);
            input_stack_.pop();
            auto type = node->get_scope_type();
            auto ret_ = RetVal();

            std::set<unsigned> active_tensors;
            active_tensors.insert(node->get_active_tensors().fill_tensors.begin(), 
                                node->get_active_tensors().fill_tensors.end());
            active_tensors.insert(node->get_active_tensors().wb_tensors.begin(), 
                                node->get_active_tensors().wb_tensors.end());

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
                for (auto pv: active_tensors) {
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
            if (verbose_level > 1) {
                std::cout << "======ComputeDelta========" << std::endl;
                std::cout << "------Node--------" << std::endl;
                input.curr_node_->display("", false);
                std::cout << "------Input-------" << std::endl;
                std::cout << input;
                std::cout << "------Output------" << std::endl;
                std::cout << ret;
                std::cout << "============END===========" << std::endl;
            }

            return ret;
        }

        void DatamovementCalculator::finalizeStat(unsigned storage_id, RetVal& ret, bool profile) {
            auto storage_level = std::static_pointer_cast<model::BufferLevel>(
                topology_.GetStorageLevel(storage_id)->Clone());
            tiling::CompoundMask mask = {};
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; ++pv) 
                mask[pv] = true;
            storage_level->Evaluate(*ret.p_tile_, mask, 
                0, 
                ret.cycle_, break_on_failure);
            if (!(storage_level->Cycles() >= ret.cycle_) && verbose_level) {
                TILEFLOW_WARNING(storage_level->Name() << " runs " << storage_level->Cycles() << " > sublevel's " << ret.cycle_);
            }
            double slow_down = storage_level->Cycles() / (0.0 + ret.cycle_);
            ret.cycle_ = storage_level->Cycles();

            if (profile) {
                storage_level->FinalizeBufferEnergy();
                energy_ += storage_level->Energy();

                auto connection = topology_.connection_map_.at(storage_id);
                auto rf_net = connection.read_fill_network->Clone();
                rf_net->Evaluate(*ret.p_tile_, break_on_failure);
                auto du_net = connection.drain_update_network->Clone();
                du_net->Evaluate(*ret.p_tile_, break_on_failure);
                assert(ret.p_tile_.unique());
                energy_ += rf_net->Energy() + du_net->Energy();
                auto& data_movement = analysis_.data_movements_[storage_level->Name()];
                auto& stat = storage_level->GetStats();
                auto& specs = storage_level->GetSpecs();
                // double all_accesses = 0;
                data_movement["Energy"] += storage_level->Energy() + rf_net->Energy() + du_net->Energy();
                data_movement["Accesses"] += storage_level->Accesses();
                data_movement["SlowDown"] = std::max(data_movement["SlowDown"], slow_down);
                data_movement["CapUtil"] = std::max(data_movement["CapUtil"], storage_level->CapacityUtilization());
                data_movement["SpatialUtil"] = std::max(data_movement["SpatialUtil"], stat.utilized_instances.Max() / specs.instances.Get());
                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; ++pv){
                    std::string suffix = "::" + problem::GetShape()->DataSpaceIDToName.at(pv);
                    data_movement["Read" + suffix] += (double)stat.reads.at(pv);
                    data_movement["Update" + suffix] += (double)stat.updates.at(pv);
                    data_movement["Fill" + suffix] += (double)stat.fills.at(pv);
                    data_movement["Write" + suffix] += (double)stat.fills.at(pv) + (double)stat.updates.at(pv);
                    data_movement["Read"] += (double)stat.reads.at(pv);
                    data_movement["Update"] += (double)stat.updates.at(pv);
                    data_movement["Fill"] += (double)stat.fills.at(pv);
                    data_movement["Write"] += (double)stat.updates.at(pv) + (double)stat.fills.at(pv);
                    // all_accesses += (double)stat.reads.at(pv) + (double)stat.updates.at(pv) + (double)stat.fills.at(pv);
                }
                // auto param = (storage_level->Energy() + rf_net->Energy() + du_net->Energy()) / all_accesses;
                // std::cout << "Buffer::" << storage_level->Name() << ", " << param << std::endl;
            }

            if (verbose_level > 1) {
                std::cout << "============BEG finalizeStat==============" << std::endl; 
                std::cout << "storage_id:" << storage_id << std::endl;
                storage_level->Print(std::cout); 
                std::cout << ret;
                std::cout << "============END finalizeStat=============" << std::endl; 
            }
            ret.p_tile_.reset();
        }

    } // namespace TileFlow
} // namespaec analysis