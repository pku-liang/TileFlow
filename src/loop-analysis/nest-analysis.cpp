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

        void NestAnalysis::get_loopnest()
        {
            LoopNestConstructor(*this).construct(mapping_.root);
        }

        void NestAnalysis::analyze()
        {
            get_loopnest();
            get_tilewise_workloads();
            get_dimscale();
            get_datamovement();
        }

        void NestAnalysis::get_datamovement()
        {
            problem::Workload::SetCurrShape(&workloads_.get_shape());
            DatamovementCalculator dm(*this);
            dm.run(mapping_.root);
        }

        void NestAnalysis::get_tilewise_workloads()
        {
            std::vector<const OpNode *> opnodes = CollectOpNode().collectOpNodes(mapping_.root);
            std::cout << "Ops: " << std::endl;
            for (auto node : opnodes)
            {
                node->display("", false);
            }
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
                    problem::Shape::DataSpaceID producer_id = producer->get_type() == Node::Op ? workloads_.get_shape().DataSpaceNameToID.at(static_cast<const OpNode *>(producer)->get_name() + "::" + t) : problem::Shape::DataSpaceID(-1);
                    problem::Shape::DataSpaceID consumer_id = workloads_.get_shape().DataSpaceNameToID.at(node->get_name() + "::" + t);
                    add_access_pattern(producer_id, producer, consumer_id, node, access_pattern);
                }
                tensor2producer[ptr->get_out()] = node;
            }

            for (auto &t : workloads_.get_outs())
            {
                TILEFLOW_ASSERT(tensor2producer.count(t), "Output " << t << " is unclear");
                auto producer_op_name = static_cast<const OpNode *>(tensor2producer[t])->get_name();
                int id = workloads_.get_shape().DataSpaceNameToID.at(producer_op_name + "::" + t);
                add_access_pattern(id, tensor2producer[t], problem::Shape::DataSpaceID(-1), mapping_.root, access_pattern);
            }

            problem::Workload::Coefficients coeff_;
            auto & shape_ = workloads_.get_shape();
            coeff_[(unsigned)-1] = 1;
            for (auto& kv: shape_.CoefficientIDToName) {
                coeff_[kv.first] = workloads_.lookup_coeff(kv.second);
            }

            for (auto &kv : access_pattern)
            {
                problem::Shape shape(shape_);
                std::vector<problem::Shape::DataSpaceID> ids;

                for (auto &id_name : shape.DataSpaceIDToName)
                {
                    if (find(kv.second.begin(), kv.second.end(), id_name.first) != kv.second.end())
                        ids.push_back(id_name.first);
                }
                std::map<problem::Shape::FactorizedDimensionID, std::string> DataSpaceIDToName;
                std::map<std::string, problem::Shape::FactorizedDimensionID> DataSpaceNameToID;
                std::vector<problem::Shape::Projection> projections;

                for (int i = 0; i < (int)ids.size(); ++i)
                {
                    int id = ids[i];
                    DataSpaceIDToName[i] = shape.DataSpaceIDToName[id];
                    DataSpaceNameToID[DataSpaceIDToName[i]] = i;
                    projections.push_back(shape.Projections[id]);
                }
                shape.DataSpaceIDToName = std::move(DataSpaceIDToName);
                shape.DataSpaceNameToID = std::move(DataSpaceNameToID);
                shape.Projections = std::move(projections);

                auto &workload = configs[kv.first].workload;
                workload.SetShape(shape);
                workload.SetCoefficients(coeff_);
                workload.SetFactorizedBounds(workloads_.get_factorized_bound());
                workload.DeriveFlattenedBounds();
                configs[kv.first].global_dataspace_ids = std::move(ids);
            }
        }


        void NestAnalysis::add_access_pattern(
            problem::Shape::DataSpaceID producer_id,
            const Node *producer,
            problem::Shape::DataSpaceID consumer_id,
            const Node *consumer,
            std::unordered_map<const Node *, std::vector<problem::Shape::DataSpaceID>> &access_pattern)
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
            while (!common.empty())
            {
                auto node = common.top();
                sc.push(node);
                sp.push(node);
                common.pop();
                if (node->get_type() == Node::type_t::Tile && !static_cast<const TileNode *>(node)->is_spatial())
                    break;
            }

            while (!sc.empty())
            {
                auto node = sc.top();
                access_pattern[node].push_back(consumer_id);
                sc.pop();
            }

            while (!sp.empty())
            {
                auto node = sp.top();
                access_pattern[node].push_back(producer_id);
                sp.pop();
            }
        }

        void NestAnalysis::get_dimscale() {
            DimScaleCalculator pass(*this);
            pass.run(mapping_.root);
        }

        void NestAnalysis::Print()
        {
            std::cout << "-----------------Nest Analysis----------------" << std::endl;
            Displayer(*this).display();
            std::cout << "--------------END Nest Analysis---------------" << std::endl;
        }

        void Displayer::visitTile(const TileNode *node)
        {
            std::cout << prefix_ << "Tile:";
            auto &configs_ = analysis_.configs;
            if (configs_.count(node))
            {
                auto shape = configs_.at(node).workload.GetShape();
                for (auto t : shape->DataSpaceNameToID)
                {
                    std::cout << t.first << ",";
                }
                std::cout << std::endl;
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
        }

        void DatamovementCalculator::run(const Node* root)
        {
            screen_shot_.num_epochs_ = 1;
            screen_shot_.node = root;
            for (unsigned dim = 0; 
                dim < problem::GetShape()->NumFlattenedDimensions; dim++)
            {
                screen_shot_.cur_transform_[dim] = 0;
            }
            root->accept(this);
        }

        void DatamovementCalculator::visitOp(const OpNode *node)
        {
            assert(screen_shot_.node == node);
            compute_info_[screen_shot_.space_stamp_].accesses += screen_shot_.num_epochs_;
            deltas_.push(problem::OperationSpace(node->get_workload().get()));
        }

        void DatamovementCalculator::visitTile(const TileNode *node)
        {
            assert(screen_shot_.node == node);
            auto &config = analysis_.configs[node];
            analysis::TileFlow::PerfectLoopnestAnalyzer analyzer(*this);
            std::map<unsigned, uint64_t> fanoutX_map, fanoutY_map;
            unsigned storage_level = node->get_storage_level();
            fanoutX_map[storage_level] = analysis_.mapping_.fanoutX_map[storage_level];
            fanoutY_map[storage_level] = analysis_.mapping_.fanoutY_map[storage_level];
            analyzer.init(&config.workload, &config.loop_nest, fanoutX_map, fanoutY_map);
            analyzer.calculateDataMovement();
        }

        void DatamovementCalculator::visitScope(const ScopeNode *node)
        {
            assert(screen_shot_.node == node);
            std::vector<problem::OperationSpace> deltas;
            for (auto &child : node->get_children())
            {
                screen_shot_.node = child;
                child->accept(this);
                deltas.emplace_back(deltas_.top());
                deltas_.pop();
            }
            screen_shot_.node = node;
            deltas_.push(combineDeltas(deltas, node->type));
        }

        problem::OperationSpace DatamovementCalculator::combineDeltas(
            const std::vector<problem::OperationSpace> &deltas,
            ScopeNode::type_t type)
        {
            if (type == ScopeNode::Parallel || type == ScopeNode::Pipeline)
            {
                return deltas.front();
            }
            else if (type == ScopeNode::Sequential)
            {
                return deltas.front();
            }
        }

        problem::OperationSpace DatamovementCalculator::computeDelta(const ScreenShot &screen_shot)
        {
            auto old_screen_shot = std::move(screen_shot_);
            screen_shot_ = screen_shot;
            screen_shot_.node->accept(this);
            screen_shot_ = std::move(old_screen_shot);
            assert(!deltas_.empty());
            auto ret = deltas_.top();
            deltas_.pop();
            return ret;
        }

        void PerfectLoopnestAnalyzer::init(problem::Workload *wl, loop::Nest *nest,
                                           std::map<unsigned, std::uint64_t> fanoutX_map,
                                           std::map<unsigned, std::uint64_t> fanoutY_map)
        {
            problem::Workload::SetCurrShape(wl->GetShape());
            Init(wl, nest, fanoutX_map, fanoutY_map);
            num_epochs_ = dm.screen_shot_.num_epochs_;
        }

        void PerfectLoopnestAnalyzer::calculateDataMovement()
        {
            if (nest_state_.size() != 0)
            {
                InitializeNestProperties();
                InitializeLiveState();
                DetectImperfectFactorization();

                // Recursive call starting from the last element of the list.
                dm.deltas_.push(NestAnalysis::ComputeDeltas(nest_state_.rbegin()));
                CollectWorkingSets();
            }

            // Done.
            working_sets_computed_ = true;
        }

        void PerfectLoopnestAnalyzer::ComputeTemporalWorkingSet(
            std::vector<analysis::LoopState>::reverse_iterator cur,
            analysis::ElementState &cur_state)
        {
            // We do two things in this function: (a) calculate the size of the temporal
            // working set for this level, and (b) calculate the number of accesses to
            // this level from the inner level.
            //
            // We used to do both these tasks by accumulating the deltas returned by
            // recursive calls to inner nesting levels. That was problematic for task (a)
            // because inner levels can sometimes buffer data across iterations of *this*
            // level, which sometimes causes the union of the deltas propagated to this
            // level to form a fractured polyhedral space. Note that this fracturing is
            // fine in terms of calculating *accesses* to this level (b), since it
            // reflects filtering.
            //
            // To address this, we first attempted to restrict gradient direction changes
            // during delta computation. However, this only captures a subset of scenarios.
            // It also affects (b), but that is fine because most hardware pattern
            // generators are probably going to be unable to generate patterns that can
            // keep state alive through gradient direction changes.
            //
            // The solution we are now adopting is to use delta accumulation only for (b)
            // and to use an alternative tactic for (a). For certain problem shapes (such
            // as CNN's axis-aligned hyper-rectangles), we can trivially calculate working
            // sets by considering only the corner points in the problem sub-space walked
            // by the subnest starting from this level down. We assume that loops are
            // always ascending (FIXME: check for this during loop construction).

            int level = cur->level;

            bool at_boundary = level == 0;

            bool dump = false; // (level >= 4);

            int end = IsLastGlobalIteration_(level + 1, cur->descriptor.dimension) ? cur->descriptor.residual_end : cur->descriptor.end;

            // First, update loop gist. FIXME: handle base!=0, stride!=1.
            ASSERT(cur->descriptor.start == 0);
            ASSERT(cur->descriptor.stride == 1);
            loop_gists_temporal_[cur->descriptor.dimension] = {0, end};

            //
            // Step II: Compute Accesses by accumulating deltas returned by inner levels.
            //

            std::uint64_t num_iterations = 1 + ((end - 1 - cur->descriptor.start) /
                                                cur->descriptor.stride);

            std::vector<problem::PerDataSpace<std::size_t>> temporal_delta_sizes;
            std::vector<std::uint64_t> temporal_delta_scale;

            bool run_last_iteration = imperfectly_factorized_ || problem::GetShape()->UsesFlattening;
            bool run_second_last_iteration = imperfectly_factorized_ && run_last_iteration;

            if (analysis::TileFlow::gExtrapolateUniformTemporal && !disable_temporal_extrapolation_.at(level))
            {
                // What we would like to do is to *NOT* iterate through the entire loop
                // for this level, but instead fire iterations #0, #1 and #last, and
                // extrapolate the remainder based on the result of iteration #1.

                // Iteration #last is only required for accurate partition size tracking.
                // Otherwise, we reset the point set on any gradient change, and so
                // tracking the point set for the #last iteration is not needed.

                // Note that this entire approach will break if there is any irregularity
                // in working-set movement along the loop (e.g., a modulus in the index
                // expression).

                int dim = int(cur->descriptor.dimension);
                int scale = vector_strides_[level][dim];
                auto saved_transform = cur_transform_[dim];

                // Iteration #0.
                indices_[level] = cur->descriptor.start;
                loop_gists_temporal_.at(dim).index = indices_[level];

                if (num_iterations >= 1)
                {
                    // Invoke next (inner) loop level.
                    ++cur;
                    auto temporal_delta = ComputeDeltas(cur, at_boundary);
                    --cur;

                    temporal_delta_sizes.push_back(temporal_delta.GetSizes());
                    temporal_delta_scale.push_back(1);
                    cur_transform_[dim] += scale;

                    indices_[level] += cur->descriptor.stride;
                    loop_gists_temporal_.at(dim).index = indices_[level];

                    if (at_boundary)
                        time_stamp_.back()++;
                }

                // Iterations #1 through #last-1/last.
                if ((run_second_last_iteration && num_iterations >= 4) ||
                    (run_last_iteration && !run_second_last_iteration && num_iterations >= 3) ||
                    (!run_last_iteration && num_iterations >= 2))
                {
                    // Invoke next (inner) loop level, scaling up the number of epochs
                    // by the number of virtual iterations we want to simulate.
                    std::uint64_t virtual_iterations =
                        run_last_iteration ? num_iterations - 2 : num_iterations - 1;

                    // Run one fewer iteration for imperfect factor support
                    if (run_second_last_iteration)
                        virtual_iterations = virtual_iterations - 1;

                    auto saved_epochs = num_epochs_;
                    num_epochs_ *= virtual_iterations;

                    ++cur;
                    auto temporal_delta = ComputeDeltas(cur, at_boundary);
                    --cur;

                    num_epochs_ = saved_epochs;

                    temporal_delta_sizes.push_back(temporal_delta.GetSizes());
                    temporal_delta_scale.push_back(virtual_iterations);

                    cur_transform_[dim] += (scale * virtual_iterations);

                    indices_[level] += (cur->descriptor.stride * virtual_iterations);
                    loop_gists_temporal_.at(dim).index = indices_[level];

                    if (at_boundary)
                        time_stamp_.back() += virtual_iterations;
                }

                // Iteration # second last to find delta for imperfect factors
                if (run_second_last_iteration && num_iterations >= 3)
                {
                    // Invoke next (inner) loop level.
                    ++cur;
                    auto temporal_delta = ComputeDeltas(cur, at_boundary);
                    --cur;

                    if (num_iterations >= 4)
                    {
                        temporal_delta_scale.back()++;
                    }
                    else
                    {
                        temporal_delta_sizes.push_back(temporal_delta.GetSizes());
                        temporal_delta_scale.push_back(1);
                    }

                    cur_transform_[dim] += scale;

                    indices_[level] += cur->descriptor.stride;
                    loop_gists_temporal_.at(dim).index = indices_[level];

                    if (at_boundary)
                        time_stamp_.back()++;
                }

                // Iteration #last.
                if (run_last_iteration && num_iterations >= 2)
                {
                    // Invoke next (inner) loop level.
                    ++cur;
                    auto temporal_delta = ComputeDeltas(cur, at_boundary);
                    --cur;

                    // If we ran the virtual-iteration logic above, we shouldn't actually
                    // use this returned delta, because we will receive the delta between
                    // iteration #2 and #last. Instead, we just re-use the last delta by
                    // increasing the #virtual iterations (scale) by 1.
                    if (!run_second_last_iteration && num_iterations >= 3)
                    {
                        temporal_delta_scale.back()++;
                    }
                    else
                    {
                        temporal_delta_sizes.push_back(temporal_delta.GetSizes());
                        temporal_delta_scale.push_back(1);
                        cur_transform_[dim] += scale;
                    }

                    indices_[level] += cur->descriptor.stride;
                    loop_gists_temporal_.at(dim).index = indices_[level];

                    if (at_boundary)
                        (time_stamp_.back())++;
                }

                cur_transform_[dim] = saved_transform;
            }
            else // not analysis::TileFlow::gExtrapolateUniformTemporal
            {
                int dim = int(cur->descriptor.dimension);
                int scale = vector_strides_[level][dim];

                auto saved_transform = cur_transform_[dim];

                for (indices_[level] = cur->descriptor.start;
                     indices_[level] < end;
                     indices_[level] += cur->descriptor.stride)
                {
                    loop_gists_temporal_.at(dim).index = indices_[level];

                    // Invoke next (inner) loop level.
                    ++cur;
                    auto temporal_delta = ComputeDeltas(cur, at_boundary);
                    --cur;

                    temporal_delta_sizes.push_back(temporal_delta.GetSizes());
                    temporal_delta_scale.push_back(1);

                    cur_transform_[dim] += scale;

                    if (at_boundary)
                        (time_stamp_.back())++;
                }

                cur_transform_[dim] = saved_transform;
            } // analysis::TileFlow::gExtrapolateUniformTemporal

            if (dump)
            {
                std::cout << "-------\n";
                std::cout << "LEVEL " << level << std::endl;
                std::cout << "-------\n";
            }

            if (at_boundary)
            {
                // Track accesses for only those levels that are relevant
                // in the final analysis after CollapseTiles.
                problem::PerDataSpace<std::size_t> final_delta_sizes;
                final_delta_sizes.fill(0);

                auto num_deltas = temporal_delta_sizes.size();
                for (unsigned i = 0; i < num_deltas; i++)
                {
                    for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
                    {
                        final_delta_sizes[pv] += (temporal_delta_sizes[i][pv] * temporal_delta_scale[i]);
                    }
                }

                for (unsigned pv = 0; pv < problem::GetShape()->NumDataSpaces; pv++)
                {
                    // Set scatter factor (otherwise it will stay at 0 for temporal levels).
                    std::uint64_t scatter_factor = 1;
                    std::uint64_t multicast_factor = 1;

                    auto &access_stats = cur_state.access_stats[pv](multicast_factor, scatter_factor);
                    access_stats.accesses += final_delta_sizes[pv] * num_epochs_;

                    // Set cumulative hops for temporal levels.
                    access_stats.hops = 0.0;

                    // Update delta histogram. Hypothesis is we only need to do this for temporal levels.
                    cur_state.delta_histograms[pv][final_delta_sizes[pv]] += num_epochs_;

                } // for (datatype)
            }     // storage boundary
        }

        problem::OperationSpace PerfectLoopnestAnalyzer::ComputeDeltas(
            std::vector<analysis::LoopState>::reverse_iterator cur,
            bool at_boundary)
        {
            if (!at_boundary)
            {
                return analysis::NestAnalysis::ComputeDeltas(cur);
            }
            ScreenShot screen_shot;
            screen_shot.num_epochs_ = num_epochs_;
            screen_shot.time_stamp_ = time_stamp_;
            screen_shot.space_stamp_ = space_stamp_;
            screen_shot.cur_transform_ = cur_transform_;
            screen_shot.node = dm.screen_shot_.node->get_children().front();
            return dm.computeDelta(screen_shot);
        }

        void DimScaleCalculator::visitOp(const OpNode *)
        {
            cur_scales.push(std::vector<uint64_t>(n_dim, 1));
        }

        void DimScaleCalculator::visitScope(const ScopeNode *node)
        {
            std::vector<std::uint64_t> cur_scale;
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
            auto &config = dm.analysis_.configs[dm.screen_shot_.node];
            vector_strides_ = config.vector_strides_;
            // std::cout << "vector_strides_:" << std::endl;
            // for (auto& vector_stride: vector_strides_) {
            //     vector_stride.Print(std::cout) << std::endl;
            // }
            mold_low_ = config.mold_low_;
            mold_high_ = config.mold_high_;
            mold_high_residual_ = config.mold_high_residual_;
            cur_transform_ = dm.screen_shot_.cur_transform_;
            // std::cout << "curr_transform_:";
            // cur_transform_.Print(std::cout) << std::endl;
        }

        void PerfectLoopnestAnalyzer::InitStorageBoundaries()
        {
            storage_boundary_level_.resize(nest_state_.size(), false);
            arch_storage_level_.resize(nest_state_.size());
            disable_temporal_extrapolation_.resize(nest_state_.size(), false);
            
            unsigned storage_level = static_cast<const TileNode*>(dm.screen_shot_.node)->get_storage_level();
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

    } // namespace TileFlow

} // namespace analysis

// 1. timestamp/spacestamp
// 2. get coefficient
// 3. Flattened v.s.
// 4. no link transfers