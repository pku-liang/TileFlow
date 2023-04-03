#include <stack>

#include "tileflow/loop-analysis/nest-analysis.hpp"

namespace analysis
{
    namespace TileFlow
    {

        void DatamovementCalculator::run(const Node *root)
        {
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
        }

        void DatamovementCalculator::visitOp(const OpNode *node)
        {
            auto input = input_stack_.top();
            assert(input.curr_node_ == node);
            input_stack_.pop();
            auto &compute_info = analysis_.configs[node].stats_[input.space_stamp_].compute_info_;
            compute_info.accesses += input.num_epochs_;
            compute_info.replication_factor = analysis_.configs[node].replication_factor;
            // Add a mask here.
            RetVal ret;
            auto &active_tensors = analysis_.configs.at(node).active_read_tensors;
            problem::OperationSpace point_set(MemoryState::workload_,
                                              input.cur_transform_, input.cur_transform_); // A single point
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
            {
                if (find(active_tensors.begin(), active_tensors.end(), pv) == active_tensors.end())
                    point_set.GetDataSpace(pv).Reset();
            }
            ret.last_working_set_.insert(0, point_set);
            // The compute level does not have data persistance.
            ret.deltas_ = ret.last_working_set_; // - input.init_working_set_;

            // std::cout << "--------Op---------" << std::endl;
            // std::cout << input << std::endl;
            // std::cout << ret << std::endl;
            // std::cout << "-------------------" << std::endl;
            tile_generator_.visitOp(node, ret.tile_);
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
            tile_generator_.visitTile(node, ret.tile_);
            ret_stack_.push(ret);
        }

        void DatamovementCalculator::visitScope(const ScopeNode *node)
        {
            auto input = input_stack_.top();
            assert(input.curr_node_ == node);
            input_stack_.pop();
            auto type = node->get_scope_type();
            auto ret_ = RetVal();
            if (type == ScopeNode::Sequential)
            {
                for (auto &child : node->get_children())
                {
                    input.curr_node_ = child;
                    input_stack_.push(input);
                    child->accept(this);
                    auto ret = ret_stack_.top();
                    ret_stack_.pop();
                    ret_.deltas_.Add(ret.deltas_);
                    input.init_working_set_ = ret.last_working_set_;
                }
                ret_.last_working_set_ = input.init_working_set_;
            }
            else if (type == ScopeNode::Parallel)
            {
                for (auto &child : node->get_children())
                {
                    input.curr_node_ = child;
                    input_stack_.push(input);
                    child->accept(this);
                    auto ret = ret_stack_.top();
                    ret_stack_.pop();
                    ret_.deltas_.Add(ret.deltas_);
                    ret_.last_working_set_.Union(ret.last_working_set_);
                }
            }
            else if (type == ScopeNode::Sharing)
            {
                for (auto &child : node->get_children())
                {
                    input.curr_node_ = child;
                    input_stack_.push(input);
                    child->accept(this);
                    auto ret = ret_stack_.top();
                    ret_stack_.pop();
                    ret_.last_working_set_.Union(ret.last_working_set_);
                }
                ret_.deltas_ = ret_.last_working_set_.Substract(input.init_working_set_);
            }
            else if (type == ScopeNode::Pipeline)
            { // other schema's children shares the memory
                // 0, 1, ..., children.size() - 1
                /*
                0  1  2  3  4
                p1          0
                p2 p1       1
                p3 p2 p1    2
                   p3 p2 p1 3
                      p3 p2 4
                         p3 5
                [:i]
                1 p0 p1 p0 ... pk pk-1, ..., p1;
                2 pk-1, ..., p0;
                ...
                N+1-k pk-1 ... p0
                ...
                N-1 pk-1, pk-2
                N pk-1
                */
                for (auto &child : node->get_children())
                {
                    input.curr_node_ = child;
                    input_stack_.push(input);
                    child->accept(this);
                    auto ret = ret_stack_.top();
                    ret_stack_.pop();
                    ret_.deltas_.Union(ret.deltas_);
                    ret_.last_working_set_.Union(ret.last_working_set_);
                }
            }

            ret_stack_.push(ret_);
        }

        RetVal DatamovementCalculator::computeDelta(const InputParam &input)
        {
            input_stack_.push(input);
            input.curr_node_->accept(this);
            auto ret = ret_stack_.top();
            ret_stack_.pop();
            return ret;
        }

        void DatamovementCalculator::visitCollect(const CollectNode* node) {
            auto input = input_stack_.top();
            auto ret = ret_stack_.top(); 
            ret_stack_.pop();
            
        }

        void TileStatGenerator::visitTile(const TileNode* node, tiling::CompoundTile& tile) {
            if (!node->is_spatial()) {
                auto & infos = tile.data_movement_info;
                const auto & config = analysis_.configs[node];
                for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; pv++){
                    auto& info = infos[pv];
                    info.subnest = config.loop_nest.loops;
                    for (auto& kv: config.stats_) {
                        if ((int)kv.second.access_stat_.size() <= pv) {
                            std::cout << "ERROR: size " << kv.second.access_stat_.size()
                                << ", pv: " << pv << std::endl; 
                        }
                        info.access_stats.Accumulate(kv.second.access_stat_[pv]);
                        // Is this correct?
                        info.size += kv.second.max_size_[pv];
                        info.link_transfers += kv.second.link_transfer_[pv];
                    }
                    int n_sample = config.stats_.size();
                    assert(n_sample);
                    info.access_stats.Divide(n_sample);
                    info.size /= n_sample;
                    info.link_transfers /= n_sample;
                    info.shape = info.size;
                    info.dataspace_id = (unsigned) pv;
                    info.partition_size = 0;
                    info.distributed_multicast = false;
                    info.content_accesses = info.access_stats.TotalAccesses();
                    info.peer_accesses = 0;
                    info.peer_fills = 0;
                    info.replication_factor = config.logical_x * config.logical_y;
                    info.fanout = config.fanout_x * config.fanout_y; // is this correct?
                    info.SetTensorRepresentation();
                    
                    info.parent_level = std::numeric_limits<unsigned>::max();
                    info.child_level = std::numeric_limits<unsigned>::max();
                    
                    // info.partition_size = info.size;
                }
                ComputeReadUpdate(node);
                ComputeDensityModels(node);
            }
            else {
                ComputeParentShareAccess(tile);
                ComputePeerAccesses(tile);
                ComputeFill(tile);
            }
        }

        void TileStatGenerator::visitOp(const OpNode* node, tiling::CompoundTile& tile) {
            auto& stats_ = analysis_.configs[node].stats_;
            double avg_accesses = 0.0;
            for (auto& kv:stats_) {
                // <space_stamp, >
                avg_accesses += kv.second.compute_info_.accesses;
            }
            avg_accesses /= stats_.size();

            auto & compute_info = tile.compute_info;
            compute_info.accesses += avg_accesses;
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
        }

        void TileStatGenerator::ComputePeerAccesses(tiling::CompoundTile& tile) {
            // TODO 
        }

        void TileStatGenerator::ComputeParentShareAccess(tiling::CompoundTile& tile) {
            auto& info = analysis_.tiles_[node].data_movement_info;
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; 
                ++pv) {
                    info[pv].parent_access_share = 0;
                }
            
            int fanout = 1;

            // find parent;
            const Node* parent = node->get_parent();
            while(parent != nullptr) {
                auto parent_ = static_cast<const TileNode*>(parent);
                if (parent_->get_tile_type() == TileNode::Temporal) break;
                else fanout = analysis_.configs.at(parent_).logical_fanout;
                parent = parent->get_parent();
            }

            if (parent == nullptr) { // root
                return;
            }

            auto& active_tensors = analysis_.configs.at(node).active_fill_tensors;
            auto& parent_info = analysis_.tiles_[parent].data_movement_info;

            for (auto pv: active_tensors) {
                auto & pv_info = info[pv];
                for (auto& x: parent_info[pv].access_stats.stats){
                    auto multicast_factor = x.first.first;
                    auto accesses = x.second.accesses;
                    pv_info.parent_access_share += (accesses * multicast_factor) / (fanout);
                }
            }
        }

        void TileStatGenerator::ComputeReadUpdate(tiling::CompoundTile& tile) {
            auto& info_ = tile.data_movement_info;
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; 
                ++pv) {
                auto& info = info_[pv];
                info.reads = std::round(info.content_accesses + info.peer_accesses);
                
                info.updates = info.temporal_reductions = 0;
                if (problem::GetShape()->IsReadWriteDataSpace.at(pv)) {
                    info.updates = std::round(info.content_accesses);
                    info.temporal_reductions = std::round(info.content_accesses + info.peer_accesses);
                }

                auto& fine_grained_data_accesses = info.fine_grained_data_accesses;
                auto& fine_grained_format_accesses = info.fine_grained_format_accesses;

                fine_grained_format_accesses = {};

                for (unsigned op_id = 0; op_id < tiling::storageOperationTypes.size(); op_id++)
                {
                    auto op_name = tiling::storageOperationTypes[op_id];
                    fine_grained_data_accesses[op_name] = 0;
                    if (op_name.find("metadata") != std::string::npos)
                    {
                        fine_grained_format_accesses[op_name] = {};
                    }
                    else
                    {
                        fine_grained_data_accesses[op_name] = 0;
                    }
                }

                // default to uncompressed without metadata
                fine_grained_data_accesses["random_read"] = info.reads;
                
                fine_grained_data_accesses["random_update"] = info.updates;
            }
        }

        void TileStatGenerator::ComputeFill(tiling::CompoundTile& tile) {
            auto & info_ = tile.data_movement_info;
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces;  
                ++pv){
                auto& info = info_[pv];
                info.fills = info.parent_access_share + info.peer_fills;
                info.fine_grained_data_accesses["random_fill"] = info.fills;
            }
        }

        void TileStatGenerator::ComputeDensityModels(tiling::CompoundTile& tile) {
            auto & info_ = tile.data_movement_info;
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
            {
                auto & info = info_[pv];
                // TODO: might want to have a new data structure for post processed sparse traffic,
                //       now we are carrying both dense and sparse in tile info
                info.expected_density = 1.0;
                info.SetDensityModel(analysis_.common_workload_.GetDensity(pv));
                info.expected_data_occupancy = info.shape;
                info.avg_replication_factor = info.replication_factor;
            }
        } 
        

    } // namespace TileFlow
} // namespaec analysis