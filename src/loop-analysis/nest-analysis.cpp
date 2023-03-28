#include <stack>

#include "tileflow/loop-analysis/nest-analysis.hpp"

namespace analysis
{

    namespace TileFlow
    {
        bool gExtrapolateUniformTemporal =
            (getenv("TIMELOOP_DISABLE_TEMPORAL_EXTRAPOLATION") == NULL) ||
            (strcmp(getenv("TIMELOOP_DISABLE_TEMPORAL_EXTRAPOLATION"), "0") == 0);

        std::vector<const OpNode *> CollectOpNode::collectOpNodes(Node *root)
        {
            opnodes_.clear();
            root->accept(this);
            return std::move(opnodes_);
        }

        void CollectOpNode::visitOp(const OpNode *node)
        {
            opnodes_.push_back(node);
        }

        NestAnalysis::NestAnalysis(problem::TileFlow::Workloads& workloads_, 
            mapping::TileFlow::Mapping& mapping_, 
            model::Engine::Specs& arch_specs_)
            : workloads_(workloads_), mapping_(mapping_), arch_specs_(arch_specs_),
            common_workload_(workloads_.get_workload()){
            common_workload_.Show();
            problem::Workload::SetCurrShape(common_workload_.GetShape());
        }

        void NestAnalysis::get_loopnest()
        {
            LoopNestConstructor(*this).construct(mapping_.root);
        }

        void NestAnalysis::analyze()
        {
            get_loopnest();
            get_dimscale();
            get_active_tensors();
            get_storage_level();
            get_spatial_offsets();
            get_datamovement();
            collect_info();
        }

        void NestAnalysis::get_storage_level() {
            StorageLevelCalculator pass(*this);
            pass.run(mapping_.root);
        }

        void NestAnalysis::get_spatial_offsets() {
            SpatialOffsetsCalculator pass(*this);
            pass.run(mapping_.root);
        }

        void NestAnalysis::get_datamovement()
        {
            DatamovementCalculator dm(*this);
            dm.run(mapping_.root);
        }

        

        void NestAnalysis::get_active_tensors()
        {
            std::vector<const OpNode *> opnodes = CollectOpNode().collectOpNodes(mapping_.root);
            std::unordered_map<std::string, const Node *> tensor2producer;
            for (auto &t : workloads_.get_ins())
            {
                tensor2producer[t] = mapping_.root;
            }

            std::unordered_map<const Node *, std::vector<problem::Shape::DataSpaceID>> access_pattern;

            for (auto node : opnodes)
            {
                auto &ptr = node->get_workload();
                for (auto &t : ptr->get_ins())
                {
                    TILEFLOW_ASSERT(tensor2producer.count(t), "Op " << node->get_name() << "'s input " << t << " is unclear");
                    auto producer = tensor2producer[t];
                    problem::Shape::DataSpaceID producer_id = producer->get_type() == Node::Op ? 
                        workloads_.get_shape().DataSpaceNameToID.at(static_cast<const OpNode *>(producer)->get_name() + "::" + t) :
                        problem::Shape::DataSpaceID(-1);
                    problem::Shape::DataSpaceID consumer_id = workloads_.get_shape().DataSpaceNameToID.at(node->get_name() + "::" + t);
                    add_access_pattern(producer_id, producer, consumer_id, node);
                }
                tensor2producer[ptr->get_out()] = node;
            }

            for (auto &t : workloads_.get_outs())
            {
                TILEFLOW_ASSERT(tensor2producer.count(t), "Output " << t << " is unclear");
                auto producer_op_name = static_cast<const OpNode *>(tensor2producer[t])->get_name();
                int id = workloads_.get_shape().DataSpaceNameToID.at(producer_op_name + "::" + t);
                add_access_pattern(id, tensor2producer[t], problem::Shape::DataSpaceID(-1), mapping_.root);
            }
        }

        void NestAnalysis::add_access_pattern(
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

            // adds common to read active tensors
            while (!common.empty())
            {
                auto node = common.top();
                common.pop();
                if (consumer_id != problem::Shape::DataSpaceID(-1))
                    configs[node].active_read_tensors.push_back(consumer_id);
                if (producer_id != problem::Shape::DataSpaceID(-1))
                    configs[node].active_read_tensors.push_back(producer_id);
                if (node->get_type() == Node::type_t::Tile && 
                    static_cast<const TileNode *>(node)->get_tile_type() == TileNode::Temporal)
                    break;
            }

            if (consumer_id != problem::Shape::DataSpaceID(-1))
                while (!sc.empty())
                {
                    auto node = sc.top();
                    configs[node].active_read_tensors.push_back(consumer_id);
                    configs[node].active_fill_tensors.push_back(consumer_id);
                    sc.pop();
                }

            if (producer_id != problem::Shape::DataSpaceID(-1))
                while (!sp.empty())
                {
                    auto node = sp.top();
                    configs[node].active_read_tensors.push_back(producer_id);
                    configs[node].active_fill_tensors.push_back(producer_id);
                    sp.pop();
                }
        }

        void NestAnalysis::get_dimscale()
        {
            DimScaleCalculator pass(*this);
            pass.run(mapping_.root);
        }

        void NestAnalysis::collect_info() {
            ComputeAccess pass1(*this);
            pass1.run(mapping_.root);
            ComputeExpansion pass2(*this);
            pass2.run(mapping_.root);
        }

        void NestAnalysis::Print()
        {
            std::cout << "-----------------Nest Analysis----------------" << std::endl;
            std::cout << "ComputeInfo:" << std::endl;
            // for (auto& kv: compute_info_) {
            //     for (auto& spatial_id: kv.first) std::cout << spatial_id << ",";
            //     std::cout << ":";
            //     std::cout << kv.second.accesses << std::endl;
            // }
            Displayer(*this).display();
            std::cout << "--------------END Nest Analysis---------------" << std::endl;
        }

        void Displayer::visitTile(const TileNode *node)
        {
            std::cout << prefix_ << "Tile:" << std::endl;
            auto &configs_ = analysis_.configs;
            if (configs_.count(node))
            {
                auto &config = configs_[node];
                auto& shape = analysis_.workloads_.get_shape();
                std::cout << prefix_ << "read tensors:";
                for (auto id : config.active_read_tensors)
                    std::cout << shape.DataSpaceIDToName.at(id) << ",";
                std::cout << std::endl;
                std::cout << prefix_ << "fill tensors:";
                for (auto id : config.active_fill_tensors)
                    std::cout << shape.DataSpaceIDToName.at(id) << ",";
                std::cout << std::endl;
                std::cout << prefix_ << "storage:" << analysis_.arch_specs_.topology.LevelNames()[config.storage_level];
                std::cout << ", fanout:" << config.fanout_x << "," << config.fanout_y << std::endl;
                std::cout << prefix_ << "offset:" << config.spatial_offset_x << "," << config.spatial_offset_y << ",";
                std::cout << prefix_ << "l-fanout" << config.logical_x << "," << config.logical_y << std::endl;
                std::cout << prefix_ << "repFactor:" << config.replication_factor << std::endl;
                
                for (auto& kv: config.stats_) {
                    std::cout << prefix_ << "<";
                    for (auto&x: kv.first) std::cout << x << ",";
                    std::cout << "> max_size, link_transfer, access_stats: " << std::endl;
                    for (int pv = 0; pv < (int) problem::GetShape()->NumDataSpaces; 
                        pv++) {
                        std::cout << prefix_ 
                            << problem::GetShape()->DataSpaceIDToName.at(pv) << ":" 
                            << kv.second.max_size_[pv] << "," 
                            << kv.second.link_transfer_[pv] << ","
                            << kv.second.access_stat_[pv];
                    }
                    std::cout << std::endl;
                }
            }
            if (analysis_.tiles_.count(node)) {
                auto& info_ = analysis_.tiles_.at(node).data_movement_info;
                std::cout << prefix_ << "size, read, fill, updates:" << std::endl; 
                for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces;   
                    pv ++) {
                    auto& info = info_[pv];
                    std::cout << prefix_ << "  " << problem::GetShape()->DataSpaceIDToName.at(pv)
                        << ": " << info.size << "," << info.reads << "," << info.fills
                        << "," << info.updates << std::endl;
                }
            }

            node->display(prefix_, false);
            auto old_prefix = prefix_;
            for (int i = 0; i < (int)node->n_level(); ++i)
                prefix_ += "   ";
            for (auto child : node->get_children())
                const_cast<Node *>(child)->accept(this);
            prefix_ = old_prefix;
        }
        void Displayer::visitScope(const ScopeNode *node)
        {
            node->display(prefix_, false);
            std::cout << "{" << std::endl;
            std::string old_prefix = prefix_;
            prefix_ += "   ";
            for (auto child : node->get_children())
                const_cast<Node *>(child)->accept(this);
            prefix_ = old_prefix;
            std::cout << prefix_ << "}" << std::endl;
        }
        void Displayer::visitOp(const OpNode *node)
        {
            node->display(prefix_, false);
            std::cout << prefix_ << "ComputeInfo:";
            auto& info = analysis_.configs[node].stats_;
            for (auto& kv: info) {
                std::cout << "<";
                for (auto& idx: kv.first) std::cout << idx << ",";
                std::cout << ">:";
                std::cout << kv.second.compute_info_.accesses << "X" 
                    << kv.second.compute_info_.replication_factor << ",";
            }
            auto& compute_info = analysis_.tiles_.at(node).compute_info;
            std::cout << std::endl;
            std::cout << prefix_<< "repFactor:" << compute_info.replication_factor << std::endl;
            std::cout << prefix_<< "accesses:" << compute_info.accesses << std::endl;
            std::cout << prefix_<< "expanison:" << compute_info.max_x_expansion 
                << "," << compute_info.max_y_expansion << std::endl;
        }

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
            auto& compute_info = analysis_.configs[node].stats_[input.space_stamp_].compute_info_;
            compute_info.accesses += input.num_epochs_;
            compute_info.replication_factor = analysis_.configs[node].replication_factor;
            // Add a mask here.
            RetVal ret;
            auto & active_tensors = analysis_.configs.at(node).active_read_tensors;
            problem::OperationSpace point_set(MemoryState::workload_,
                 input.cur_transform_, input.cur_transform_); // A single point
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++) {
                if (find(active_tensors.begin(), active_tensors.end(), pv) 
                    == active_tensors.end()) 
                    point_set.GetDataSpace(pv).Reset();
            }
            ret.last_working_set_.insert(0, point_set);
            ret.deltas_ = ret.last_working_set_ - input.init_working_set_;

            // std::cout << "--------Op---------" << std::endl;
            // std::cout << input << std::endl;
            // std::cout << ret << std::endl;
            // std::cout << "-------------------" << std::endl;
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
            // Done.
            working_sets_computed_ = true;

            return ret_val;
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
                std::transform(cur_scale.begin(), cur_scale.end(), cur_scale_.begin(),
                               cur_scale.begin(), [](uint64_t x, uint64_t y)
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
            // 1. recursive call to simulate the runnning;
            if (input_.num_epochs_)
                SimulateTemporalExecution();
            // 2. Get current working set
            problem::OperationSpace point_set = GetCurrentWorkingSet(nest_state_.rbegin());
            
            // TODO: add optional temporal reuse here

            // We need to mask out the unused tensor from the point_set
            auto & active_tensors = config_.active_read_tensors;
            for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++) {
                if (find(active_tensors.begin(), active_tensors.end(), pv) == active_tensors.end()) 
                    point_set.GetDataSpace(pv).Reset();
            }
            
            auto sizes = point_set.GetSizes();
            auto & max_size = config_.stats_[input_.space_stamp_].max_size_;
            std::transform(sizes.begin(), sizes.end(),
                           max_size.begin(), max_size.begin(),
                           [](std::size_t x, std::size_t y)
                           { return std::max(x, y); });

            RetVal ret;
            problem::OperationSpace delta(workload_);
            if (input_.init_working_set_.getDataSpaces().count(0))
                ret.deltas_[0] = point_set - input_.init_working_set_.at(0);
            else ret.deltas_[0] = point_set;
            ret.last_working_set_[0] = point_set;
            return ret;
        }

        void PerfectLoopnestAnalyzer::SimulateTemporalExecution()
        {
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
            
            auto& access_stat = config_.stats_[input_.space_stamp_].access_stat_;

            // the init
            {
                InputParam input = input_;
                input.num_epochs_ = 1;
                input.curr_node_ = input_.curr_node_->get_children().front();
                auto ret_val = dm_.computeDelta(input);
                for (auto &kv : ret_val.deltas_.getDataSpaces()){
                    for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++){
                        access_stat[pv](1, 1).accesses += input.num_epochs_ * kv.second.GetSize(pv);
                    }
                }
            }

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
                ret_val = dm_.computeDelta(input);  
                // std::cout << "--------------" << std::endl;
                // std::cout << "last_index:" << last_index << std::endl;
                // std::cout << input;
                // std::cout << ret_val; 
                // std::cout << "--------------" << std::endl;
                for (auto &kv : ret_val.deltas_.getDataSpaces()){
                    for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++){
                        access_stat[pv](1, 1).accesses += input.num_epochs_ * kv.second.GetSize(pv);
                    }
                }
            }
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
                }
                else {
                    auto translation_vectors = GetCurrentTranslationVectors(extrapolation_level);
                    std::uint64_t dst_delta_index = SpatialIDL2P(base_index);
                    std::uint64_t src_delta_index = SpatialIDL2P(base_index - extrapolation_stride);
                    
                    auto& dst_temporal_delta = ret.deltas_[dst_delta_index];
                    auto& src_temporal_delta = ret.deltas_[src_delta_index];
                    auto& dst_last_working_set = ret.last_working_set_[dst_delta_index];
                    auto& src_last_working_set = ret.last_working_set_[src_delta_index];
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
            int iterations_to_run = 1;
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

        void StorageLevelCalculator::visitScope(const ScopeNode* node) {
            auto& config = analysis_.configs[node];
            auto& storage_level = config.storage_level = -1;
            for (auto child: node->get_children()) {
                child->accept(this);
                assert(!storage_levels_.empty());
                if (storage_level == (unsigned)-1) {
                    storage_level = storage_levels_.top();
                }
                else assert(storage_level == storage_levels_.top());
                storage_levels_.pop();
            }
            config.fanout_x = analysis_.mapping_.fanoutX_map.at(config.storage_level);
            config.fanout_y = analysis_.mapping_.fanoutY_map.at(config.storage_level);
            storage_levels_.push(storage_level);
        }

        void StorageLevelCalculator::visitTile(const TileNode* node) {
            for (auto child: node->get_children()) child->accept(this);
            auto& config = analysis_.configs[node];
            config.storage_level = node->get_storage_level();
            if (node->is_spatial()) {
                config.fanout_x = analysis_.mapping_.fanoutX_map.at(config.storage_level);
                config.fanout_y = analysis_.mapping_.fanoutY_map.at(config.storage_level);
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
            offset_t output_offset = init_offset;
            if (type == ScopeNode::Sequential || type == ScopeNode::Sharing) {
                input_offsets.push(init_offset);
                for (auto child: node->get_children()) {
                    child->accept(this);
                    assert(!output_offsets.empty());
                    init_offset = output_offsets.top();
                    output_offsets.pop();
                }
                output_offset = init_offset;
            }
            else {
                for (auto child: node->get_children()) {
                    input_offsets.push(init_offset);
                    child->accept(this);
                    assert(!output_offsets.empty());
                    output_offset = merge(output_offset, output_offsets.top());
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
                for (auto loop: node->get_loops()) {
                    int loop_count = 
                        (loop.end - loop.start) / loop.stride;
                    if (loop::IsSpatialX(loop.spacetime_dimension)) {
                        config.logical_x *= loop_count;
                    }
                    else config.logical_y *= loop_count;
                }
            }
            config.logical_fanout = config.fanout_x * config.fanout_y;
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
            o << "point:" << input.cur_transform_ << std::endl;
            o << "init state:" << std::endl;
            input.init_working_set_.show();
            return o;
        }

        std::ostream& operator<< (
            std::ostream& o,
            const RetVal& ret
        ) {
            o << "Ret:" << std::endl;
            o << "Delats:";
            ret.deltas_.show();
            o << "LastWorkingSet:" << std::endl;
            ret.last_working_set_.show();
            return o; 
        }

        void ComputeAccess::visitTile(const TileNode* node) {
            if (!node->is_spatial()) {
                auto & infos = analysis_.tiles_[node].data_movement_info;
                auto & config = analysis_.configs[node];
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
                }
                ComputeParentShareAccess(node);
                ComputePeerAccesses(node);
                ComputeReadUpdateFill(node);
                ComputeDensityModels(node);
            }

            for (auto child: node->get_children())
                child->accept(this);
        }

        void ComputeAccess::visitOp(const OpNode* node) {
            auto& stats_ = analysis_.configs[node].stats_;
            double avg_accesses = 0.0;
            for (auto& kv:stats_) {
                // <space_stamp, >
                avg_accesses += kv.second.compute_info_.accesses;
            }
            avg_accesses /= stats_.size();
            auto & compute_info = analysis_.tiles_[node].compute_info;
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

        void ComputeAccess::ComputePeerAccesses(const TileNode*) {
            // TODO 
        }

        void ComputeAccess::ComputeParentShareAccess(const TileNode* node) {
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

        void ComputeAccess::ComputeReadUpdateFill(const TileNode* node) {
            auto& info_ = analysis_.tiles_[node].data_movement_info;
            for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; 
                ++pv) {
                auto& info = info_[pv];
                info.reads = std::round(info.content_accesses + info.peer_accesses);
                info.fills = info.parent_access_share + info.peer_fills;
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
                fine_grained_data_accesses["random_fill"] = info.fills;
                fine_grained_data_accesses["random_update"] = info.updates;
            }
        }

        void ComputeAccess::ComputeDensityModels(const TileNode* node) {
            auto & info_ = analysis_.tiles_[node].data_movement_info;
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
        

        void ComputeExpansion::visitTile(const TileNode* node) {
            if (node->get_tile_type() == TileNode::Spatial) {
                auto old_expansion_ = expansion_;
                auto& config = analysis_.configs[node];
                expansion_.first *= config.logical_x;
                expansion_.second *= config.logical_y;
                for (auto child: node->get_children()) 
                    child->accept(this);
                expansion_ = old_expansion_;
            }
            else {
                auto & info = analysis_.tiles_[node].data_movement_info;
                for (int pv = 0; pv < (int)problem::GetShape()->NumDataSpaces; 
                    ++pv) {
                    info[pv].max_x_expansion = expansion_.first;
                    info[pv].max_y_expansion = expansion_.second;
                }
                for (auto child: node->get_children())  
                    child->accept(this);
            }
        }

        void ComputeExpansion::visitOp(const OpNode* node) {
            auto& info = analysis_.tiles_[node].compute_info;
            info.max_x_expansion = expansion_.first;
            info.max_y_expansion = expansion_.second;
        }
    } // namespace TileFlow

} // namespace analysis

/**
 * ReadWrite 
 * 
*/

