#pragma once 

#include <bitset>
#include <stack>

#include "loop-analysis/nest-analysis.hpp"
#include "loop-analysis/loop-state.hpp"

#include "tileflow/mapping/mapping.hpp"
#include "tileflow/common.hpp"
#include "tileflow/problem/problem.hpp"
#include "tileflow/loop-analysis/memory-state.hpp"

using mapping::TileFlow::Node;
using mapping::TileFlow::OpNode;
using mapping::TileFlow::TileNode;
using mapping::TileFlow::ScopeNode;
using mapping::TileFlow::Visitor;


namespace analysis {

namespace TileFlow {

    struct NodeConfig {
        // problem::Workload workload;
        std::vector<problem::Shape::DataSpaceID> active_tensors;
        loop::Nest loop_nest;
        std::vector<problem::OperationPoint> vector_strides_;
        std::vector<problem::OperationPoint> mold_low_;
        std::vector<problem::OperationPoint> mold_high_;
        std::vector<problem::OperationPoint> mold_high_residual_;
        // the state to record the last_point_set and access counting 
        problem::PerDataSpace<AccessStatMatrix> access_stats_;

        // StorageLevelCalculator
        unsigned storage_level;
        std::uint64_t fanout_x, fanout_y;

        // SpatialOffsetsCalculator
        std::uint64_t spatial_offset_x, spatial_offset_y, logical_x, logical_y;

        std::unordered_map<problem::Shape::DataSpaceID, 
            std::vector<analysis::DataMovementInfo> > data_movements;
    };

    class NestAnalysis {
        problem::TileFlow::Workloads& workloads_;
        mapping::TileFlow::Mapping& mapping_;
        model::Engine::Specs& arch_specs_;
        
        void add_access_pattern(
            problem::Shape::DataSpaceID producer_id, 
            const Node* producer, 
            problem::Shape::DataSpaceID consumer_id,
            const Node* consumer);
        /**
         * \brief set the loopnest for tile nodes;
        */
        inline void get_loopnest();
        /**
         * \brief compute datamovement along the mapping;
         * \depends get_dimscale, get_tilewise_workloads
        */
        void get_datamovement();
        /**
         * \brief init the strides for tile nodes;
         * \depends get_loopnest()
        */
        void get_dimscale();
        /**
         * \brief get alive tensors for Tile Nodes
        */
        void get_active_tensors();
        /**
         * \brief set the storage level for tile nodes;
        */
        void get_storage_level();
        /**
         * \brief set the spatial offset for tile nodes;
         * \depends get_storage_level()
         */ 
        void get_spatial_offsets();

        // this is set by get_datamovement
        std::map<std::vector<unsigned>, ComputeInfo> compute_info_;
    public: 
        NestAnalysis(problem::TileFlow::Workloads& workloads_, 
            mapping::TileFlow::Mapping& mapping_, 
            model::Engine::Specs& arch_specs_)
            : workloads_(workloads_), mapping_(mapping_), arch_specs_(arch_specs_){}
        
        std::unordered_map<const Node*, NodeConfig> configs;

        void analyze();
        void Print();
        friend class Displayer;
        friend class DatamovementCalculator;
        friend class DimScaleCalculator;
        friend class LoopNestConstructor;
        friend class StorageLevelCalculator;
    };

    /**
     * \brief the API between pass on nodes. 
    */
    struct InputParam{
        int num_epochs_;
        const Node* curr_node_ = nullptr;
        std::vector<unsigned> time_stamp_;
        std::vector<unsigned> space_stamp_;
        problem::OperationPoint cur_transform_;
        MemoryState init_working_set_;
    };

    std::ostream& operator<< (std::ostream& o, const InputParam& params);

    struct RetVal{
        MemoryState deltas_;
        MemoryState last_working_set_;

        RetVal() {}
        RetVal(const problem::OperationPoint& low_pt, 
            const problem::OperationPoint& high_pt) {
            deltas_.insert(0, low_pt, high_pt);
            last_working_set_ = deltas_;
        }
        RetVal(problem::OperationSpace&& data_space_) {
            deltas_.insert(0, data_space_);
            last_working_set_ = deltas_;
        }
    };

    std::ostream& operator<< (std::ostream& o, const RetVal& params);

    class DatamovementCalculator: public mapping::TileFlow::Visitor {
        NestAnalysis& analysis_;
        problem::Workload& workload_;
        /**
         * \brief the stack to pass parameter between nodes.
        */
        std::stack<InputParam> input_stack_;
        std::stack<RetVal> ret_stack_;
        // const Node* curr_node_;
        RetVal computeDelta(const InputParam& input);
        
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;
    public:
        DatamovementCalculator(NestAnalysis& analysis, problem::Workload& workload): 
            analysis_(analysis), workload_(workload){}
        void run(const Node*) override;
        friend class PerfectLoopnestAnalyzer;
    };

    class PerfectLoopnestAnalyzer: public analysis::NestAnalysis {
        DatamovementCalculator& dm_;
        const InputParam& input_;
        NodeConfig& config_;
        RetVal ComputeTemporalWorkingSet();
        RetVal ComputeSpatialWorkingSet();
        void SimulateTemporalExecution();
        void InitPerLevelDimScales() override;
        void InitStorageBoundaries() override;
        virtual problem::PerDataSpace<Point> GetCurrentTranslationVectors(
            std::vector<analysis::LoopState>::reverse_iterator cur) override;
        void FillSpatialDeltas(std::vector<analysis::LoopState>::reverse_iterator cur,
                         std::uint64_t base_index,
                         int extrapolation_stride,
                         std::vector<analysis::LoopState>::reverse_iterator extrapolation_level,
                         RetVal& ret_val);
        problem::OperationSpace ComputeDeltas(
            std::vector<analysis::LoopState>::reverse_iterator cur,
            bool at_boundary);

        std::uint64_t SpatialIDL2P(std::uint64_t logical_id);
    public: 
        PerfectLoopnestAnalyzer(
            DatamovementCalculator& dm,
            InputParam& input, 
            NodeConfig& config): dm_(dm), input_(input), config_(config){}
        void init(problem::Workload*, loop::Nest*);
        RetVal calculateDataMovement();
        
        
    };

    class CollectOpNode: public mapping::TileFlow::Visitor {
        void visitOp(const OpNode*) override;
        std::vector<const OpNode*> opnodes_;
    public: 
        std::vector<const OpNode*> collectOpNodes(Node*);
    };

    class Displayer: public mapping::TileFlow::Visitor {
        NestAnalysis & analysis_;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;
        std::string prefix_;
    public: 
        Displayer(NestAnalysis& analysis): analysis_(analysis){}
        void display() {analysis_.mapping_.root->accept(this);}
    };

    class DimScaleCalculator: public mapping::TileFlow::Visitor {
        std::stack<std::vector<std::uint64_t> > cur_scales;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;
        NestAnalysis& analysis_;
        int n_dim;
    public: 
        DimScaleCalculator(NestAnalysis& analysis): analysis_(analysis){
            n_dim = analysis_.workloads_.get_shape().NumFactorizedDimensions;
        }
    };

    class LoopNestConstructor: public mapping::TileFlow::Visitor {
        void visitTile(const TileNode* node) override { 
            analysis_.configs[node].loop_nest = node->constructLoopNest(
                analysis_.workloads_.get_shape().FactorizedDimensionNameToID);
            for (auto child: node->get_children()) {child->accept(this);}
        }
        NestAnalysis& analysis_;
    public:
        LoopNestConstructor(NestAnalysis& analysis): analysis_(analysis){}
        void construct(const Node* root) {root->accept(this);}
    };

    class StorageLevelCalculator: public mapping::TileFlow::Visitor {
        std::stack<unsigned> storage_levels_;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        NestAnalysis& analysis_;
    public: 
        StorageLevelCalculator(NestAnalysis& analysis): analysis_(analysis){
        }
    };

    class SpatialOffsetsCalculator: public mapping::TileFlow::Visitor {
        struct offset_t {
            unsigned x, y;
            unsigned max_x;
            
        };
        offset_t merge(const offset_t& o1, const offset_t& o2);
        std::stack<offset_t> input_offsets, output_offsets;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        NestAnalysis& analysis_;
    public: 
        SpatialOffsetsCalculator(NestAnalysis& analysis): analysis_(analysis){
        }
    };
} // namespace TileFlow 

} // namespace analysis 