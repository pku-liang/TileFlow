#include <stack>
#include <fstream>

#include "tileflow/loop-analysis/nest-analysis.hpp"

using TileFlow::verbose_level;

namespace analysis
{

    namespace TileFlow
    {
        bool gExtrapolateUniformTemporal =
            (getenv("TIMELOOP_DISABLE_TEMPORAL_EXTRAPOLATION") == NULL) ||
            (strcmp(getenv("TIMELOOP_DISABLE_TEMPORAL_EXTRAPOLATION"), "0") == 0);

        NestAnalysis::NestAnalysis(
            const problem::TileFlow::Workloads& workloads_, 
            const mapping::TileFlow::Mapping& mapping_, 
            const model::Engine::Specs& arch_specs_, 
            const model::Topology& topology_)
            : workloads_(workloads_), mapping_(mapping_), arch_specs_(arch_specs_), 
            topology_(topology_), common_workload_(workloads_.get_workload()){
            if (verbose_level > 1) {
                std::cout << "begin analysis..." << std::endl;
                common_workload_.Show();
            }
            problem::Workload::SetCurrShape(common_workload_.GetShape());
        }

        void NestAnalysis::get_loopnest()
        {
            LoopNestConstructor(*this, symb_table_).construct(mapping_.root);
        }

        void NestAnalysis::reset() {
            cycle_ = 0;
            energy_ = 0;
            configs.clear();
            tiles_.clear();
            data_movements_.clear();
        }

        void NestAnalysis::analyze()
        {
            reset();
            get_loopnest();
            get_dimscale();
            get_storage_level();
            get_spatial_offsets(); // depends on get_loopnest
            get_expansion();
            get_datamovement();
        }

        void NestAnalysis::get_storage_level() {
            StorageLevelCalculator pass(*this);
            pass.run(mapping_.root);
        }

        void NestAnalysis::get_spatial_offsets() {
            SpatialOffsetsCalculator pass(*this);
            pass.run(mapping_.root);
        }

        void NestAnalysis::get_expansion(){
            ComputeExpansion pass(*this);
            pass.run(mapping_.root);
        }

        void NestAnalysis::get_datamovement()
        {
            DatamovementCalculator dm(*this);
            RetVal ret = dm.eval(mapping_.root);
            cycle_ = ret.cycle_ + topology_.total_network_latency_;
            energy_ = dm.Energy();
        }

        void NestAnalysis::get_dimscale()
        {
            DimScaleCalculator pass(*this);
            pass.run(mapping_.root);
        }

        void NestAnalysis::Print(std::ostream& o)
        {
            o << "-----------------Nest Analysis----------------" << std::endl;
            Displayer(*this, symb_table_, o).display();
            o << "Cycle: " << cycle_
                    << ", Energy: " << energy_
                     << std::endl;
            o << "--------------END Nest Analysis---------------" << std::endl;
        }

        void NestAnalysis::Report() {
            std::cout << "***TileFlow Result" << std::endl;
            std::cout << ",value" << std::endl;
            std::cout << "Cycle," << cycle_ << std::endl;
            std::cout << "Energy," << energy_ << std::endl;
            std::cout << "***TileFlow Result Ends" << std::endl;
        }

        void NestAnalysis::Export(const std::string& filename) {
            std::ofstream file;
            size_t tmp = filename.find_last_of('.');
            if (tmp == std::string::npos || filename.substr(tmp) != ".csv")
                file.open(filename + ".csv");
            else file.open(filename);
            file << "metric,value" << std::endl;
            file << "Cycle," << cycle_ << std::endl;
            file << "Energy," << energy_ << std::endl;
            if (verbose_level > 1)
                std::cout << "[TileFlow]: result written into " << filename << std::endl;
            file.close();
        }

        void Displayer::visitTile(const TileNode *node)
        {
            auto &configs_ = analysis_.configs;
            if (verbose_level) {
                o_ << prefix_ << node->get_name() << "," << std::endl;
                if (configs_.count(node))
                {
                    auto &config = configs_[node];
                    o_ << prefix_ << "strides,low,high:";
                    for (unsigned dim = 0; dim < problem::GetShape()->NumFactorizedDimensions; dim++){
                        o_ << problem::GetShape()->FactorizedDimensionIDToName.at(dim);
                    }
                    o_ << config.vector_strides_.front() << "," << config.mold_low_.front() << "," << config.mold_high_.front() << std::endl;
                    if (verbose_level > 1) {
                        o_ << prefix_ << "storage:" << node->get_storage_name();
                        o_ << ", fanout:" << config.fanout_x << "," << config.fanout_y << std::endl;
                        o_ << prefix_ << "offset:" << config.spatial_offset_x << "," << config.spatial_offset_y << ",";
                        o_ << prefix_ << "l-fanout" << config.logical_x << "," << config.logical_y << std::endl;
                        o_ << prefix_ << "repFactor:" << config.replication_factor << std::endl;
                        o_ << prefix_ << "strides:" << std::endl;
                        for (unsigned i = 0; i < config.loop_nest.loops.size(); i++) {
                            o_ << prefix_ << config.loop_nest.loops[i].PrintCompact() << ":" 
                                << config.vector_strides_[i] << std::endl;
                        }
                        for (auto& kv: config.stats_) {
                            o_ << prefix_ << "<";
                            for (auto&x: kv.first) o_ << x << ",";
                            o_ << "> max_size, link_transfer, access_stats: " << std::endl;
                            for (int pv = 0; pv < (int) problem::GetShape()->NumDataSpaces; 
                                pv++) {
                                o_ << prefix_ 
                                    << problem::GetShape()->DataSpaceIDToName.at(pv) << ":" 
                                    << kv.second.max_size_[pv] << "," 
                                    << kv.second.link_transfer_[pv] << ","
                                    << kv.second.access_stat_[pv];
                            }
                            o_ << std::endl;
                        }
                    }
                }
                if (analysis_.tiles_.count(node)) {
                    auto& info_ = analysis_.tiles_.at(node).data_movement_info;
                    o_ << info_;
                    // o_ << prefix_ << "size, read, fill, updates:" << std::endl; 
                    // for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces;   
                    //     pv ++) {
                    //     auto& info = info_[pv];
                    //     o_ << prefix_ << "  " << problem::GetShape()->DataSpaceIDToName.at(pv)
                    //         << ": " << info.size << "," << info.reads << "," << info.fills
                    //         << "," << info.updates << std::endl;
                    // }
                }
            }

            node->display(prefix_, false, symbol_table_, o_);
            auto old_prefix = prefix_;
            for (int i = 0; i < (int)node->n_level(); ++i)
                prefix_ += "   ";
            for (auto child : node->get_children())
                const_cast<Node *>(child)->accept(this);
            prefix_ = old_prefix;
        }

        void Displayer::visitScope(const ScopeNode *node)
        {
            node->display(prefix_, false, symbol_table_, o_);
            std::string old_prefix = prefix_;
            o_ << prefix_ << "{" << std::endl;
            prefix_ += "   ";
            for (auto child : node->get_children())
                const_cast<Node *>(child)->accept(this);
            prefix_ = old_prefix;
            o_ << prefix_ << "}" << std::endl;
        }
        void Displayer::visitOp(const OpNode *node)
        {
            node->display(prefix_, false, symbol_table_, o_);
            if (verbose_level && analysis_.tiles_.count(node)) {
                auto& compute_info = analysis_.tiles_.at(node).compute_info;
                o_ << prefix_<< "repFactor:" << compute_info.replication_factor << std::endl;
                o_ << prefix_<< "accesses:" << compute_info.accesses << std::endl;
                o_ << prefix_<< "expanison:" << compute_info.max_x_expansion 
                    << "," << compute_info.max_y_expansion << std::endl;
                if (verbose_level > 1) {
                    o_ << prefix_ << "ComputeInfo:";
                    auto& config = analysis_.configs[node];
                    auto& info = config.stats_;
                    for (auto& kv: info) {
                        o_ << "<";
                        for (auto& idx: kv.first) o_ << idx << ",";
                        o_ << ">:";
                        o_ << kv.second.compute_info_.accesses << "X" 
                            << kv.second.compute_info_.replication_factor << ",";
                    }
                    o_ << std::endl;
                }
            }
        }

        
        void DimScaleCalculator::visitOp(const OpNode *)
        {
            cur_scales.push(std::vector<uint64_t>(n_dim, 1));
        }

        void DimScaleCalculator::visitScope(const ScopeNode *node)
        {
            std::vector<std::uint64_t> cur_scale(n_dim, 1);
            for (auto child : node->get_children())
            {
                child->accept(this);
                auto cur_scale_ = cur_scales.top();
                cur_scales.pop();
                for (int i = 0; i < (int)cur_scale.size(); ++i) {
                    auto x = cur_scale[i];
                    auto y = cur_scale_[i];
                    auto dim = problem::GetShape()->FactorizedDimensionIDToName.at(i);
                    TILEFLOW_ASSERT((1 == x || 1 == y || x == y), dim << " mismatch " << x << " v.s. " << y << std::endl; 
                        node->display("", true, analysis_.symb_table_);
                        child->display("", true, analysis_.symb_table_);
                        std::cerr);
                }
                std::transform(cur_scale.begin(), cur_scale.end(), cur_scale_.begin(),
                               cur_scale.begin(), [this, node, child](uint64_t x, uint64_t y)
                               {
                        assert(1 == x || 1 == y || x == y);
                        return std::max(x,y); });
            }
            cur_scales.push(cur_scale);
        }

        void DimScaleCalculator::visitTile(const TileNode *node)
        {
            node->get_children().front()->accept(this);
            assert(!cur_scales.empty());
            auto cur_scale = cur_scales.top();
            cur_scales.pop();

            auto &config = analysis_.configs[node];
            auto &loops = config.loop_nest.loops;
            config.vector_strides_.resize(loops.size());
            config.mold_low_.resize(loops.size());
            config.mold_high_.resize(loops.size());
            config.mold_high_residual_.resize(loops.size());
            for (int level = loops.size() - 1; level >= 0; --level)
            {
                auto &desc = loops[level];
                for (int dim = 0; dim < n_dim; dim++)
                {
                    config.vector_strides_[level][dim] = cur_scale[dim];
                }
                cur_scale[int(desc.dimension)] *= desc.end - desc.start;
                for (std::uint64_t dim = 0;
                     dim < problem::GetShape()->NumFlattenedDimensions; dim++)
                {
                    config.mold_high_residual_[level][dim] = config.mold_high_[level][dim] = cur_scale[dim] - 1;
                }
            }
            cur_scales.push(cur_scale);
        }

        void StorageLevelCalculator::visitScope(const ScopeNode* node) {
            auto& config = analysis_.configs[node];
            auto& storage_level = config.storage_level = 0;
            for (auto child: node->get_children()) {
                child->accept(this);
                TILEFLOW_ASSERT(!storage_levels_.empty(), "Is there an op node under the Scope Node?");
                storage_level = std::max(storage_level, storage_levels_.top());
                storage_levels_.pop();
            }
            config.fanout_x = analysis_.mapping_.fanoutX_map.at(node->get_storage_level());
            config.fanout_y = analysis_.mapping_.fanoutY_map.at(node->get_storage_level());
            storage_levels_.push(storage_level);
        }

        void StorageLevelCalculator::visitTile(const TileNode* node) {
            for (auto child: node->get_children()) {
                child->accept(this);
                if (!storage_levels_.empty())
                    storage_levels_.pop();
            }
            auto& config = analysis_.configs[node];
            config.storage_level = node->get_storage_level();
            if (node->is_spatial()) {
                config.fanout_x = analysis_.mapping_.fanoutX_map.at(node->get_storage_level());
                config.fanout_y = analysis_.mapping_.fanoutY_map.at(node->get_storage_level());
            }
            else config.fanout_x = config.fanout_y = 1;
            storage_levels_.push(config.storage_level);
        }



        SpatialOffsetsCalculator::offset_t SpatialOffsetsCalculator::merge(
            const SpatialOffsetsCalculator::offset_t& o1,
            const SpatialOffsetsCalculator::offset_t& o2
        ) {
            return {std::max(o1.x, o2.x), std::max(o1.y, o2.y), std::max(o1.max_x, o2.max_x)};
        }
        
        void SpatialOffsetsCalculator::visitScope(const ScopeNode* node) {
            auto type = node->get_scope_type();
            offset_t init_offset = {0,0,0};
            if (!input_offsets.empty()) {
                init_offset = input_offsets.top();
                input_offsets.pop();
            }
            offset_t output_offset;
            if (type == ScopeNode::Sequential || type == ScopeNode::Sharing) {
                output_offset = {0,0,0};
                for (auto child: node->get_children()) {
                    input_offsets.push(init_offset);
                    child->accept(this);
                    assert(!output_offsets.empty());
                    output_offset = merge(output_offset, output_offsets.top());
                    output_offsets.pop();
                }
            }
            else {
                output_offset = init_offset;
                for (auto child: node->get_children()) {
                    input_offsets.push(output_offset);
                    child->accept(this);
                    assert(!output_offsets.empty());
                    output_offset = output_offsets.top();
                    output_offsets.pop();
                }
            }
            output_offsets.push(output_offset);
        } 
        
        void SpatialOffsetsCalculator::visitTile(const TileNode* node) {
            offset_t init_offset = {0,0,0};
            if (!input_offsets.empty()) {
                init_offset = input_offsets.top();
                input_offsets.pop();
            }
            
            auto& config = analysis_.configs[node];
            config.logical_x = config.logical_y = 1;
            if (node->is_spatial()) {
                for (auto loop: config.loop_nest.loops) {
                    int loop_count = 
                        (loop.end - loop.start) / loop.stride;
                    if (loop::IsSpatialX(loop.spacetime_dimension)) {
                        config.logical_x *= loop_count;
                    }
                    else config.logical_y *= loop_count;
                }
            }
            config.logical_fanout = config.logical_x * config.logical_y;
            offset_t output_offset = init_offset;
            if (init_offset.y + config.logical_y <= config.fanout_y) {
                config.spatial_offset_x = init_offset.x;
                config.spatial_offset_y = init_offset.y;
                output_offset.y += config.logical_y;
                output_offset.max_x = std::max(output_offset.max_x, 
                    output_offset.x + (unsigned)config.logical_x);
            }
            else {
                config.spatial_offset_x = output_offset.max_x;
                config.spatial_offset_y = 0;
                output_offset.x = output_offset.max_x;
                output_offset.y = config.logical_y;
                output_offset.max_x = output_offset.x + config.logical_x; 
            }
            if (verbose_level > 1) {
                std::cout << "===========Offset Calculation=========" << std::endl;
                std::cout << "Node:" << node->get_name() << std::endl;
                node->display("\t", false, analysis_.symb_table_);
                std::cout << "input:" << init_offset << std::endl;
                std::cout << "output: " << output_offset << std::endl;
                std::cout << "fanout: " << config.fanout_x << "," << config.fanout_y << std::endl;
                std::cout << "logical:" << config.logical_x << "," << config.logical_y << std::endl;
                std::cout << "storage_level:" << node->get_storage_name() << ", " << node->get_name() << ", " << node->get_storage_level() << std::endl;
                std::cout << "======================================" << std::endl;
            }
            assert(output_offset.y <= config.fanout_y);
            assert(output_offset.max_x <= config.fanout_x);
            output_offsets.push(output_offset);
            // how much resources are used. 

            config.replication_factor = replication_factor;
            if (node->is_spatial()) 
                replication_factor *= config.logical_x * config.logical_y;
            for (auto child: node->get_children()) {
                child->accept(this);
            }
            replication_factor = config.replication_factor;
        }

        void SpatialOffsetsCalculator::visitOp(const OpNode* node) {
            auto& config = analysis_.configs[node];
            config.replication_factor = replication_factor;
        }

        std::ostream& operator << (std::ostream& o, 
            const SpatialOffsetsCalculator::offset_t& offset){
            o << "x,y,maxX:" << offset.x << "," << offset.y << "," << offset.max_x << ","; 
            return o;
        }

        std::ostream& operator<< (
            std::ostream& o, 
            const InputParam& input){
            o << "Input:" << std::endl;
            o << "n_epoch:" << input.num_epochs_ << std::endl;
            o << "t-stamp: ";
            for (auto x:input.time_stamp_) o << x << ",";
            o << std::endl;
            o << "s-stamp: ";
            for (auto x:input.space_stamp_) o << x << ",";
            o << std::endl;
            for (unsigned i = 0; i < input.cur_transform_.Order(); i++)
                std::cout << problem::GetShape()->FactorizedDimensionIDToName.at(i);
            o << ":" << input.cur_transform_ << std::endl;
            o << "init state:" << std::endl;
            input.init_working_set_.show();
            return o;
        }

        std::ostream& operator<< (
            std::ostream& o,
            const RetVal& ret
        ) {
            o << "Ret:" << std::endl;
            o << "\tDeltas:"<<std::endl;
            ret.deltas_.show();
            o << "\tLastWorkingSet:" << std::endl;
            ret.last_working_set_.show();
            o << "\tCycles:" << ret.cycle_ << std::endl;
            o << "\tAccess:" << ret.access_stat_ << std::endl;
            o << "\ttile:" << ret.p_tile_.get() << std::endl;
            return o; 
        }

        void ComputeExpansion::visitTile(const TileNode* node) {
            auto old_expansion_ = expansion_;
            auto& config = analysis_.configs[node];
            expansion_.first *= config.logical_x;
            expansion_.second *= config.logical_y;
            config.max_x_expansion = expansion_.first;
            config.max_y_expansion = expansion_.second;
            for (auto child: node->get_children()) 
                child->accept(this);
            expansion_ = old_expansion_;
        }

        void ComputeExpansion::visitOp(const OpNode* node) {
            auto& info = analysis_.tiles_[node].compute_info;
            info.max_x_expansion = expansion_.first;
            info.max_y_expansion = expansion_.second;
        }
    
        void PerfectLoopnestAnalyzer::init(const problem::Workload *wl, const loop::Nest *nest)
        {
            workload_ = const_cast<problem::Workload*>(wl);

            Reset();
            cached_nest = *nest;

            // Copy over everything we need from the nest.
            storage_tiling_boundaries_ = nest->storage_tiling_boundaries;
            packed_skew_descriptors_ = nest->skew_descriptors;
            no_link_transfer_ = nest->no_link_transfer;
            no_multicast_ = nest->no_multicast;
            no_temporal_reuse_ = nest->no_temporal_reuse;

            // Construct nest_state_.
            int level = nest->loops.size();
            for (auto iter = nest->loops.rbegin(); iter != nest->loops.rend(); ++iter)
            {
                auto& descriptor = *iter;
                analysis::LoopState cur;
                cur.level = --level;
                cur.descriptor = descriptor;
                nest_state_.push_back(cur);    
            }

            num_epochs_ = input_.num_epochs_;
        }

        RetVal PerfectLoopnestAnalyzer::calculateDataMovement()
        {
            auto node = input_.curr_node_;
            assert(node->get_type() == Node::type_t::Tile);
            auto tile_node = static_cast<const TileNode *>(node);
            RetVal ret_val;
            if (nest_state_.size() != 0)
            {
                InitializeNestProperties();
                InitializeLiveState();
                DetectImperfectFactorization();
                // std::cout << "PerfectLoopnestAnalyzer::calculateDataMovement:" << std::endl;
                // for (auto iter = nest_state_.rbegin(); iter!= nest_state_.rend(); iter++) {
                //     std::cout << problem::GetShape()->FlattenedDimensionIDToName.at(iter->descriptor.dimension) << ":";
                //     std::cout << iter->descriptor.start << "," << iter->descriptor.end << "," << iter->descriptor.stride << ",";
                //     std::cout << vector_strides_[iter->level];
                //     std::cout << std::endl;
                // }
                if (!tile_node->is_spatial())
                {
                    ret_val = ComputeTemporalWorkingSet();
                }
                else
                {
                    ret_val = ComputeSpatialWorkingSet();
                }
            }

            if (input_.num_epochs_ && (node->get_parent() == nullptr || 
                node->get_parent()->get_type() != Node::Tile || 
                static_cast<const TileNode*>(node->get_parent())->get_tile_type() != TileNode::Spatial)) {
                assert(ret_val.p_tile_ != nullptr);
                ComputeParentShareAccess(ret_val.access_stat_, *ret_val.p_tile_);
                ComputeFill(*ret_val.p_tile_);
                FinalizeTile(*ret_val.p_tile_);
            }
            // Done.
            working_sets_computed_ = true;

            return ret_val;
        }

        void PerfectLoopnestAnalyzer::InitPerLevelDimScales()
        {
            auto &config = dm_.analysis_.configs[input_.curr_node_];
            vector_strides_ = config.vector_strides_;
            // std::cout << "vector_strides_:" << std::endl;
            // for (auto& vector_stride: vector_strides_) {
            //     vector_stride.Print(std::cout) << std::endl;
            // }
            mold_low_ = config.mold_low_;
            mold_high_ = config.mold_high_;
            mold_high_residual_ = config.mold_high_residual_;
            cur_transform_ = input_.cur_transform_;
            // std::cout << "curr_transform_:";
            // cur_transform_.Print(std::cout) << std::endl;
        }

        void PerfectLoopnestAnalyzer::InitStorageBoundaries()
        {
            storage_boundary_level_.resize(nest_state_.size(), false);
            arch_storage_level_.resize(nest_state_.size());
            disable_temporal_extrapolation_.resize(nest_state_.size(), false);

            unsigned storage_level = static_cast<const TileNode *>(input_.curr_node_)->get_storage_level();
            unsigned loop_level = 0;
            // std::cout << "nest_state_.size(): " << nest_state_.size();
            // std::cout << "storage_tiling_boundaries_: ";
            // for (auto &i: storage_tiling_boundaries_) std::cout << i << ",";
            // std::cout << std::endl;
            for (auto &i : storage_tiling_boundaries_)
            {
                ASSERT(i < storage_boundary_level_.size());
                storage_boundary_level_[i] = true;

                auto skew_it = packed_skew_descriptors_.find(storage_level);
                if (skew_it != packed_skew_descriptors_.end())
                {
                    skew_descriptors_[i] = skew_it->second;

                    // Walk through the skew descriptor and poison all temporal loop
                    // variables it touches.
                    for (auto &term : skew_it->second.terms)
                    {
                        if (term.variable.dimension != problem::GetShape()->NumFlattenedDimensions && !term.variable.is_spatial)
                        {
                            auto dim = term.variable.dimension;
                            // Walk through the loops in this loop block and poison the loop
                            // corresponding to this problem dimension.
                            for (unsigned level = loop_level; level <= i; level++)
                            {
                                if (nest_state_.at(level).descriptor.dimension == dim)
                                {
                                    disable_temporal_extrapolation_.at(level) = true;
                                    // There can only be 1 such loop in each block.
                                    break;
                                }
                            }
                        }
                    }
                }

                // Establish loop level -> storage level map.
                for (; loop_level <= i; loop_level++)
                {
                    arch_storage_level_[loop_level] = storage_level;
                }
                storage_level++;
            }
        }

        RetVal PerfectLoopnestAnalyzer::ComputeTemporalWorkingSet()
        {
            // 1. Get current working set
            problem::OperationSpace point_set = GetCurrentWorkingSet(nest_state_.rbegin());
            
            // We need to mask out the unused tensor from the point_set
            auto &read_tensors = input_.curr_node_->get_active_tensors().read_tensors;
            auto &update_tensors = input_.curr_node_->get_active_tensors().update_tensors;

            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++) {
                if (!read_tensors.count(pv) && !update_tensors.count(pv)) 
                    point_set.GetDataSpace(pv).Reset();
            }
            
            auto sizes = point_set.GetSizes();
            auto & max_size = config_.stats_[input_.space_stamp_].max_size_;
            std::transform(sizes.begin(), sizes.end(),
                           max_size.begin(), max_size.begin(),
                           [](std::size_t x, std::size_t y)
                           { return std::max(x, y); });
            
            // 2. recursive call to simulate the running;
            RetVal ret;
            ret.last_working_set_[0] = point_set;
            if (input_.num_epochs_){
                ret.p_tile_ = std::make_shared<tiling::CompoundTile>();
                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces;
                ++pv) {
                    auto& info = ret.p_tile_->data_movement_info[pv];
                    info.size = info.shape = sizes[pv];
                }   
                problem::PerDataSpace<AccessStatMatrix> child_access_stat;
                ret.cycle_ = SimulateTemporalExecution(child_access_stat);
                InitTile(child_access_stat, *ret.p_tile_);
                ComputeReadUpdate(*ret.p_tile_);
                ComputeDensityModels(*ret.p_tile_);
                problem::OperationSpace delta(workload_);
                if (input_.init_working_set_.getDataSpaces().count(0))
                    ret.deltas_[0] = point_set - input_.init_working_set_.at(0);
                else ret.deltas_[0] = point_set;
                problem::PerDataSpace<std::uint64_t> link_transfers;
                ComputeStats(input_.init_working_set_, ret.deltas_, 
                    ret.access_stat_, link_transfers);
            }
            
            return ret;
        }

        int PerfectLoopnestAnalyzer::SimulateTemporalExecution(
            problem::PerDataSpace<AccessStatMatrix>& access_stat)
        {
            std::uint64_t cycle = 0;
            std::vector<int> dims;
            std::vector<int> scales;
            std::vector<int> loop_counts = {};
            std::vector<int> trip_counts = {1};
            for (auto iter = nest_state_.rbegin(); iter != nest_state_.rend(); ++iter) {
                dims.push_back(iter->descriptor.dimension);
                scales.push_back(vector_strides_[iter->level][dims.back()]);
                loop_counts.push_back((iter->descriptor.end - iter->descriptor.start)
                    / iter->descriptor.stride);
                trip_counts.push_back(trip_counts.back() * loop_counts.back());
            }
            
            
            auto active_tensors = input_.curr_node_->get_active_tensors().read_tensors;
            auto update_tensors = input_.curr_node_->get_active_tensors().update_tensors;
            active_tensors.insert(update_tensors.begin(), update_tensors.end());
            
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces;
                ++pv) {
                access_stat[pv].clear();        
            }
            // the init {} --> (0,0,...,0)
            {
                InputParam input = input_;
                input.init_working_set_ = {};
                input.curr_node_ = input_.curr_node_->get_children().front();
                auto ret_val = dm_.computeDelta(input);
                // std::cout << "<" << input_.curr_node_->get_storage_name() << "::init>" << std::endl;
                // std::cout << "last_index:" << input.cur_transform_ << std::endl;
                // std::cout << input;
                // std::cout << ret_val; 
                // std::cout << "</" << input_.curr_node_->get_storage_name() << "::init>" << std::endl;
                for (auto pv: active_tensors){
                    access_stat[pv].Accumulate(ret_val.access_stat_.at(pv));
                }
                cycle += ret_val.cycle_;
            }

            int all_iter = 1;
            for (auto lc: loop_counts) all_iter *= lc;

            int real_iter = 1;
            /**
             * delta = (00001)
             *          (0009) -> (0010)
            */
            for (int i = 0; i < (int)nest_state_.size(); ++i) {
                if (loop_counts[i] <= 1) continue;
                InputParam input = input_;
                input.curr_node_ = input_.curr_node_->get_children().front();
                input.num_epochs_ = 0;
                for (int j = i+1; j < (int)nest_state_.size(); ++j) {
                    input.cur_transform_[dims[j]] += scales[j] * (loop_counts[j] - 1);
                }
                RetVal ret_val = dm_.computeDelta(input); 
                auto last_index = input.cur_transform_;
                
                input.init_working_set_ = ret_val.last_working_set_;
                input.num_epochs_ = num_epochs_ * (trip_counts[i+1] - trip_counts[i]);
                input.cur_transform_ = input_.cur_transform_;
                input.cur_transform_[dims[i]] += scales[i];
                real_iter += input.num_epochs_;
                ret_val = dm_.computeDelta(input);  
                // std::cout << "<" << input_.curr_node_->get_storage_name() << "::"
                //      << problem::GetShape()->FlattenedDimensionIDToName.at(nest_state_[i].descriptor.dimension)
                //     << ">" << std::endl;
                // std::cout << "last_index:" << last_index << std::endl;
                // std::cout << input;
                // std::cout << ret_val; 
                // std::cout << "</" << input_.curr_node_->get_storage_name() << "::"
                //      << problem::GetShape()->FlattenedDimensionIDToName.at(nest_state_[i].descriptor.dimension)
                //     << ">" << std::endl;
                for (auto pv: active_tensors){
                    access_stat[pv].Accumulate(ret_val.access_stat_.at(pv));
                }
                cycle += ret_val.cycle_;
            }

            auto & acc_access_stat = config_.stats_[input_.space_stamp_].access_stat_;
            for (auto pv: active_tensors){
                acc_access_stat[pv].Accumulate(access_stat.at(pv));
            }

            // utils to create a partial tile
            return cycle;
        }

        RetVal PerfectLoopnestAnalyzer::ComputeSpatialWorkingSet()
        {
            // Only second step: compute the spatial deltas;
            RetVal ret;
            FillSpatialDeltas(nest_state_.rbegin(), 
                              0,                     // base_index,
                              0,                     // extrapolation_stride
                              nest_state_.rbegin(), // extrapolation_level 
                              ret);  // the ret val
            if (input_.num_epochs_) {
                // overwrite the fanout and replication factor
                TILEFLOW_ASSERT(ret.p_tile_, "spatial error" << std::endl; input_.curr_node_->display("", true, dm_.analysis_.symb_table_); std::cerr);
                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; ++pv){
                    auto& info = ret.p_tile_->data_movement_info[pv];
                    info.fanout = config_.fanout_x * config_.fanout_y;
                }
                problem::PerDataSpace<std::uint64_t> link_transfers;
                ComputeStats(input_.init_working_set_, 
                    ret.deltas_, ret.access_stat_, link_transfers);
                assert(ret.p_tile_!=nullptr);
                ComputePeerAccesses(link_transfers, *ret.p_tile_);
            }
            return ret;
        }

        problem::PerDataSpace<Point> PerfectLoopnestAnalyzer::GetCurrentTranslationVectors(
            std::vector<analysis::LoopState>::reverse_iterator cur) {
            problem::PerDataSpace<Point> translation_vectors;
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++) {
                unsigned order = problem::GetShape()->DataSpaceOrder.at(pv);
                auto& point = translation_vectors[pv] = Point(order);
                for (unsigned dim = 0; dim < order; dim++) {
                    point[dim] = 0;
                    for (auto& term: problem::GetShape()->Projections.at(pv).at(dim)){
                        if (term.second == cur->descriptor.dimension) {
                            point[dim] += cur->descriptor.stride 
                                * vector_strides_[cur->level][cur->descriptor.dimension];
                        }
                    }
                }

            }
            return translation_vectors;
        }

        void PerfectLoopnestAnalyzer::FillSpatialDeltas(std::vector<analysis::LoopState>::reverse_iterator cur,
                                                        std::uint64_t base_index,
                                                        int extrapolation_stride,
                                                        std::vector<analysis::LoopState>::reverse_iterator extrapolation_level,
                                                        RetVal& ret)
        {
            if (cur == nest_state_.rend()) {
                if (extrapolation_stride == 0) { 
                    // must be singleton.
                    assert(singleton);
                    singleton = false;
                    std::uint64_t spatial_id = SpatialIDL2P(base_index);
                    InputParam input = input_;
                    input.curr_node_ = input_.curr_node_->get_children().front();
                    input.space_stamp_.push_back(spatial_id); 
                    input.cur_transform_ = cur_transform_;
                    input.init_working_set_ = 
                        MemoryState(0, input.init_working_set_[spatial_id]);
                    RetVal ret_ = dm_.computeDelta(input);
                    ret.deltas_.insert(spatial_id, ret_.deltas_[0]);
                    ret.last_working_set_.insert(spatial_id, ret_.last_working_set_[0]);
                    ret.cycle_ = ret_.cycle_;
                    ret.p_tile_ = std::move(ret_.p_tile_);
                }
                else {
                    auto translation_vectors = GetCurrentTranslationVectors(extrapolation_level);
                    std::uint64_t dst_delta_index = SpatialIDL2P(base_index);
                    std::uint64_t src_delta_index = SpatialIDL2P(base_index - extrapolation_stride);
                    
                    auto& dst_temporal_delta = ret.deltas_[dst_delta_index];
                    auto& src_temporal_delta = ret.deltas_[src_delta_index];
                    auto& dst_last_working_set = ret.last_working_set_[dst_delta_index];
                    auto& src_last_working_set = ret.last_working_set_[src_delta_index];
                    
                    // std::cout << "index:";
                    // for (auto idx: indices_) std::cout << idx << ",";
                    // std::cout << "baseidx: " << base_index  << ",";
                    // std::cout << "translation_vector:" << translation_vectors << std::endl;
                    // std::cout << "dst:" << dst_delta_index << ", src:" << src_delta_index << std::endl;
                    
                    for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
                    {
                        dst_temporal_delta.GetDataSpace(pv) = src_temporal_delta.GetDataSpace(pv);
                        dst_temporal_delta.GetDataSpace(pv).Translate(translation_vectors.at(pv));
                        dst_last_working_set.GetDataSpace(pv) = src_last_working_set.GetDataSpace(pv);
                        dst_last_working_set.GetDataSpace(pv).Translate(translation_vectors.at(pv));
                    }
                }
                return; 
            } 

            int level = cur->level;
            auto dim = cur->descriptor.dimension;

            int end = IsLastGlobalIteration_(level + 1, cur->descriptor.dimension) ? cur->descriptor.residual_end : cur->descriptor.end;

            unsigned num_iterations = 1 + ((end - 1 - cur->descriptor.start) /
                                           cur->descriptor.stride);

            // First, update loop gist. FIXME: handle base!=0, stride!=1.
            ASSERT(cur->descriptor.start == 0);
            ASSERT(cur->descriptor.stride == 1);

            base_index *= end;

            int iterations_run = 0;
            const int iterations_to_run = 1;
            int scale = vector_strides_[level][dim];

            for (indices_[level] = cur->descriptor.start;
                indices_[level] < end;
                indices_[level] += cur->descriptor.stride, iterations_run++)
            {                
                auto next_extrapolation_stride = extrapolation_stride * num_iterations; // * num_iterations?
                auto next_extrapolation_level = extrapolation_level;
                if (iterations_run >= iterations_to_run) // Extrapolate using this level
                {
                    next_extrapolation_stride = cur->descriptor.stride;
                    next_extrapolation_level = cur;
                }

                ++cur;

                FillSpatialDeltas(cur, 
                                base_index + indices_[level],
                                next_extrapolation_stride,
                                next_extrapolation_level, 
                                ret);

                --cur;
                cur_transform_[dim] += scale;
            }
        }

        std::uint64_t PerfectLoopnestAnalyzer::SpatialIDL2P(std::uint64_t logical_id) {
            return ((logical_id % config_.logical_y) + config_.spatial_offset_y) +
                ((logical_id / config_.logical_y) + config_.spatial_offset_x) * config_.fanout_y;
        }

        void PerfectLoopnestAnalyzer::ComputeStats(
            const MemoryState& last_working_set,
            const MemoryState& deltas,
            problem::PerDataSpace<AccessStatMatrix>& access_stats, 
            problem::PerDataSpace<std::uint64_t>& link_transfers
        ) {
            problem::PerDataSpace<std::unordered_set<std::uint64_t>> unaccounted_delta;
            for (auto& delta: deltas.getDataSpaces())
            {
                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
                unaccounted_delta[pv].insert(delta.first);
            }
            const int enable_link_transfer = true;
            if (enable_link_transfer) {
                ComputeLinkTransfer(last_working_set, deltas, link_transfers, unaccounted_delta);
            }
            bool enable_multicast = 
                static_cast<const TileNode*>(input_.curr_node_)->is_multicast_enabled();
            ComputeAccessStat(deltas, unaccounted_delta, access_stats, enable_multicast);
            return;   
        }

       
        void PerfectLoopnestAnalyzer::ComputePeerAccesses(
            const problem::PerDataSpace<std::uint64_t>& link_transfer,
            tiling::CompoundTile& tile){
            auto& infos = tile.data_movement_info;
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces;
            ++pv) {
                auto& info = infos[pv];
                info.link_transfers = link_transfer[pv];
                info.peer_fills = info.peer_accesses 
                    = info.link_transfers / config_.logical_fanout; 
                info.reads += info.peer_accesses;
                if (problem::GetShape()->IsReadWriteDataSpace.at(pv))
                    info.temporal_reductions += info.peer_accesses;
            }
        }

        void PerfectLoopnestAnalyzer::ComputeLinkTransfer(
            const MemoryState& last_working_set, 
            const MemoryState& delta,
            problem::PerDataSpace<std::uint64_t>& link_transfers,
            problem::PerDataSpace<std::unordered_set<std::uint64_t>>& unaccounted_delta
            
        ) {
            const std::vector<std::pair<int, int> > routings = {
                {-1,0}, {0,-1}, {0,1}, {1,0}
            };
            int h_size = config_.fanout_x;
            int v_size = config_.fanout_y;

            auto GetLinearIndex = [&h_size, &v_size](int h_id, int v_id)
                {
                assert(0 <= h_id && h_id < h_size && 0 <= v_id && v_id < v_size);
                std::uint64_t linearIndex = v_id * h_size + h_id;  // row major layout
                return linearIndex;
                };

            auto& cur_spatial_deltas = delta.getDataSpaces();
            auto& prev_spatial_deltas = last_working_set.getDataSpaces();

            int num_spatial_elems = h_size * v_size;

            std::vector<problem::PerDataSpace<bool>> inter_elem_reuse;
            inter_elem_reuse.resize(num_spatial_elems);
            for (int i = 0; i < num_spatial_elems; i++)
            {
                inter_elem_reuse.at(i).fill(false);
            }

            problem::PerDataSpace<bool> no_link_transfer;
            for(unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
            {
                no_link_transfer[pv] = false;
            }

            for (auto& route: routings) {
                for (int h_id = 0; h_id < h_size; h_id++)
                {
                    for (int v_id = 0; v_id < v_size; v_id++)
                    {
                        int h_id_ = h_id + route.first;
                        int v_id_ = v_id + route.second;
                        if (h_id_ < 0 || h_id_ >= h_size || v_id_ < 0 || v_id_ >= v_size)   
                            continue;
                        auto cur_skewed_spatial_index = GetLinearIndex(h_id, v_id);
                        auto prev_skewed_spatial_index = GetLinearIndex(h_id_, v_id_);
                        CompareSpatioTemporalDeltas(cur_spatial_deltas, prev_spatial_deltas,
                                                    cur_skewed_spatial_index, prev_skewed_spatial_index,
                                                    inter_elem_reuse,
                                                    no_link_transfer);
                    }
                }
            }

            // Compute the total number of accesses that can be bypassed
            // by using link transfers
            link_transfers.fill(0);

            problem::PerDataSpace<std::uint64_t>& acc_link_transfers = 
                config_.stats_[input_.space_stamp_].link_transfer_;

            for (auto& delta: cur_spatial_deltas)
            {
                auto& cur_skewed_spatial_index = delta.first;
                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces;
                    pv++)
                {
                    if (inter_elem_reuse.at(cur_skewed_spatial_index)[pv])
                    {
                        acc_link_transfers[pv] += (delta.second.GetSize(pv) * num_epochs_);
                        link_transfers[pv] += (delta.second.GetSize(pv) * num_epochs_);
                        auto unaccounted_it = unaccounted_delta[pv].find(cur_skewed_spatial_index);
                        ASSERT(unaccounted_it != unaccounted_delta[pv].end());
                        unaccounted_delta[pv].erase(unaccounted_it);
                    }
                }
            }
            if (verbose_level > 1) {
                std::cout << "===========BEG LINK===========" << std::endl;
                std::cout << "last working set:" << std::endl;
                last_working_set.show();
                std::cout << "delta:" << std::endl;
                delta.show();
                std::cout << "num_epochs_: " << num_epochs_ << std::endl;
                std::cout << "link transfer:";
                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; ++pv) {
                    std::cout << problem::GetShape()->DataSpaceIDToName.at(pv) 
                        << ":" << link_transfers[pv] << ",";
                }
                std::cout << std::endl;
                std::cout << "===========END LINK===========" << std::endl;
            }
        }

        void PerfectLoopnestAnalyzer::ComputeAccessStat(
            const MemoryState& delta, 
            problem::PerDataSpace<std::unordered_set<std::uint64_t>>& unaccounted_delta, 
            problem::PerDataSpace<AccessStatMatrix>& access_stats,
            bool enable_multicast
        ){
            auto& spatial_deltas = delta.getDataSpaces();

            // For each data type, records the number of unaccounted deltas
            // that the current delta matches with. This will be used
            // to infer the multicast factor for a specific delta.
            // reused across loop iterations to avoid initialization overheads.
            problem::PerDataSpace<uint64_t> num_matches;

            // Prepare a legacy-style multicast/scatter signature to collect the
            // delta data, then populate the new access stats.
            struct TempAccessStats
            {
                double accesses = 0;
                std::uint64_t scatter_factor = 0;
                double hops = 0.0;
            };
            problem::PerDataSpace<std::unordered_map<std::uint64_t, TempAccessStats>> temp_stats;

            // FIXME: we should only be looking at physical dimensions here. The problem
            // is that sparse mappings may appear to exceed the physical dimensions before
            // space-skipping is applied. The very notion of spatial skew and physical
            // location for space-skipping sparse mappings is something we need to figure
            // out.
            auto h_size = config_.fanout_x;
            auto v_size = config_.fanout_y;
            

            for (auto delta_it = spatial_deltas.begin(); delta_it != spatial_deltas.end(); delta_it++)
                //for (std::uint64_t i = 0; i < num_deltas; i++)
            {
                auto& skewed_spatial_index = delta_it->first;
                auto& delta = delta_it->second;

                num_matches.fill(0);
                
                problem::PerDataSpace<std::vector<std::uint64_t>> match_set;

                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
                {
                    auto unaccounted_it = unaccounted_delta[pv].find(skewed_spatial_index);
                    if (unaccounted_it == unaccounted_delta[pv].end())
                        //if (!unaccounted_delta[i][pv])
                    {
                        // this delta was already accounted for,
                        // skip the comparisons.
                        continue;
                    }

                    unaccounted_delta[pv].erase(unaccounted_it);
                    //unaccounted_delta[i][pv] = false;
                    num_matches[pv] = 1;  // we match with ourselves.
                    match_set[pv].push_back(skewed_spatial_index);

                    for (auto delta_other_it = std::next(delta_it); delta_other_it != spatial_deltas.end(); delta_other_it++)
                    //for (std::uint64_t j = i + 1; j < num_deltas; j++)
                    {
                        auto& skewed_other_spatial_index = delta_other_it->first;
                        auto& delta_other = delta_other_it->second;

                        auto unaccounted_other_it = unaccounted_delta[pv].find(skewed_other_spatial_index);
                        if (unaccounted_other_it != unaccounted_delta[pv].end())
                            //if (unaccounted_delta[j][pv])
                        {
                            if (delta.CheckEquality(delta_other, pv))
                            //if (spatial_deltas[i].CheckEquality(spatial_deltas[j], pv))
                            {
                            // We have a match, record it
                            unaccounted_delta[pv].erase(unaccounted_other_it);
                            //unaccounted_delta[j][pv] = false;
                            num_matches[pv]++;
                            match_set[pv].push_back(skewed_other_spatial_index);
                            }
                        }
                    }
                }

                // update the number of accesses at different multicast factors.
                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
                {
                    if (num_matches[pv] > 0 && delta.GetSize(pv) > 0)
                    {
                        auto& temp_struct = temp_stats[pv][num_matches[pv]];
                        temp_struct.accesses += (delta.GetSize(pv) * num_epochs_);
                        temp_struct.scatter_factor++;

                        // Compute the average number of hops from the edge of the array
                        // (at this level) to the nodes in the match set.
                        // Assume injection point is at center of V-axis. Routing algorithm is
                        // to go along H maximally, then drop vertical paths.

                        ASSERT(num_matches[pv] == match_set[pv].size());
                        
                        double hops = 0;
                        
                        std::uint64_t h_max = 0;
                        for (auto& linear_id : match_set[pv])
                        {
                            std::uint64_t h_id = linear_id % h_size;
                            h_max = std::max(h_max, h_id);
                        }
                        hops += double(h_max);
                        
                        double v_center = double(v_size-1) / 2;
                        for (auto& linear_id : match_set[pv])
                        {
                            std::uint64_t v_id = linear_id / h_size;
                            hops += std::abs(double(v_id) - v_center);
                        }

                        // Accumulate this into the running hop count. We'll finally divide this
                        // by the scatter factor to get average hop count.
                        temp_struct.hops += hops;
                    }
                }
            }

            // Populate the actual stats.
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
            {
                for (auto& x: temp_stats[pv])
                {
                    auto multicast = x.first;
                    auto scatter = x.second.scatter_factor;
                    if (enable_multicast)
                        access_stats[pv](multicast, scatter) = { x.second.accesses, x.second.hops };
                    else 
                        access_stats[pv](1,scatter) = {x.second.accesses * multicast, x.second.hops};
                }
            }
        }

        void PerfectLoopnestAnalyzer::InitTile(
            const problem::PerDataSpace<AccessStatMatrix>& access_stat,
            tiling::CompoundTile& tile
        ){
            auto & infos = tile.data_movement_info;
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; pv++){
                auto& info = infos[pv];
                info.subnest = config_.loop_nest.loops;
                info.access_stats = access_stat.at(pv);
                info.link_transfers = 0;
                info.dataspace_id = (unsigned) pv;
                info.partition_size = 0;
                info.distributed_multicast = false;
                info.content_accesses = info.access_stats.TotalAccesses();
                info.peer_accesses = 0;
                info.peer_fills = 0;
                info.replication_factor = config_.replication_factor;
                info.fanout = config_.fanout_x * config_.fanout_y; // is this correct?
                info.SetTensorRepresentation();
                info.max_x_expansion = config_.max_x_expansion;
                info.max_y_expansion = config_.max_y_expansion;

                info.parent_level = std::numeric_limits<unsigned>::max();
                info.child_level = std::numeric_limits<unsigned>::max();
            }
        }

        void PerfectLoopnestAnalyzer::ComputeReadUpdate(
            tiling::CompoundTile& tile) {
            
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; 
                ++pv){
                auto & info = tile.data_movement_info[pv];
                info.reads = info.updates = info.temporal_reductions = 0;
            }

            for (auto pv: input_.curr_node_->get_active_tensors().read_tensors) {
                auto & info = tile.data_movement_info[pv];
                info.reads = std::round(info.content_accesses + info.peer_accesses);
            }

            for (auto pv: input_.curr_node_->get_active_tensors().update_tensors) {
                auto & info = tile.data_movement_info[pv];
                info.updates = std::round(info.content_accesses);
                info.temporal_reductions = std::round(info.content_accesses + info.peer_accesses);
            }
        }

        void PerfectLoopnestAnalyzer::ComputeDensityModels(
            tiling::CompoundTile& tile
            ) {
            auto & infos = tile.data_movement_info;
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
            {
                auto & info = infos[pv];
                // TODO: might want to have a new data structure for post processed sparse traffic,
                //       now we are carrying both dense and sparse in tile info
                info.expected_density = 1.0;
                info.SetDensityModel(dm_.analysis_.common_workload_.GetDensity(pv));
                info.expected_data_occupancy = info.shape;
                info.avg_replication_factor = info.replication_factor;
            }
        } 

        void PerfectLoopnestAnalyzer::ComputeParentShareAccess(
            const problem::PerDataSpace<AccessStatMatrix>& access_stats,
            tiling::CompoundTile& tile) {
            auto& info = tile.data_movement_info;
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; 
                ++pv) {
                    info[pv].parent_access_share = 0;
                }
            
            for (auto pv: input_.curr_node_->get_active_tensors().fill_tensors) {
                auto & pv_info = info[pv];
                for (auto& x: access_stats.at(pv).stats){
                    auto multicast_factor = x.first.first;
                    auto accesses = x.second.accesses;
                    pv_info.parent_access_share += 
                     (accesses * multicast_factor) / (pv_info.fanout);
                }
            }
        }

        void PerfectLoopnestAnalyzer::ComputeFill(tiling::CompoundTile& tile) {
            auto & info_ = tile.data_movement_info;
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces;  
                ++pv){
                auto& info = info_[pv];
                info.fills = info.parent_access_share + info.peer_fills;
            }
        }

        void PerfectLoopnestAnalyzer::FinalizeTile(tiling::CompoundTile& tile) {
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; 
            ++pv) {
                auto& info = tile.data_movement_info[pv];
                
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

                info.fine_grained_data_accesses["random_fill"] = info.fills;
            }
        }

    } // namespace TileFlow

} // namespace analysis

/**
 * ReadWrite 
 * 
*/

