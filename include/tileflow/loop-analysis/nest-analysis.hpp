#pragma once 

#include <bitset>
#include <stack>

#include "loop-analysis/nest-analysis.hpp"
#include "loop-analysis/loop-state.hpp"

#include "tileflow/mapping/mapping.hpp"
#include "tileflow/common.hpp"
#include "tileflow/problem/problem.hpp"

using mapping::TileFlow::Node;
using mapping::TileFlow::OpNode;
using mapping::TileFlow::TileNode;
using mapping::TileFlow::ScopeNode;
using mapping::TileFlow::Visitor;


namespace analysis {

namespace TileFlow {

    struct NodeConfig {
        problem::Workload workload;
        std::vector<problem::Shape::DataSpaceID> global_dataspace_ids;
        loop::Nest loop_nest;
        std::vector<problem::OperationPoint> vector_strides_;
        std::vector<problem::OperationPoint> mold_low_;
        std::vector<problem::OperationPoint> mold_high_;
        std::vector<problem::OperationPoint> mold_high_residual_;

        std::unordered_map<problem::Shape::DataSpaceID, 
            std::vector<analysis::DataMovementInfo> > data_movements;
        
    };

    class NestAnalysis {
        // opname X tensorname --> AccessPattern
        std::unordered_map<std::string,
            std::unordered_map<std::string, 
            problem::Shape::Projection> >  access_patterns;

        problem::TileFlow::Workloads& workloads_;
        mapping::TileFlow::Mapping& mapping_;
        
        void add_access_pattern(
            problem::Shape::DataSpaceID producer_id, 
            const Node* producer, 
            problem::Shape::DataSpaceID consumer_id,
            const Node* consumer,
            std::unordered_map<const Node*, std::vector<problem::Shape::DataSpaceID> >& access_pattern);
        /**
         * \brief set the workload for tile nodes;
         * \depends get_loopnest
         */
        void get_tilewise_workloads();
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
    public: 
        NestAnalysis(problem::TileFlow::Workloads& workloads_, mapping::TileFlow::Mapping& mapping): workloads_(workloads_), mapping_(mapping){}
        std::unordered_map<const Node*, NodeConfig> configs;

        void analyze();
        void Print();
        friend class Displayer;
        friend class DatamovementCalculator;
        friend class DimScaleCalculator;
        friend class LoopNestConstructor;
    };

    struct ScreenShot{
        int num_epochs_;
        std::vector<unsigned> time_stamp_;
        std::vector<unsigned> space_stamp_;
        problem::OperationPoint cur_transform_;
        const Node* node;
    };

    class DatamovementCalculator: public mapping::TileFlow::Visitor {
        std::map<std::vector<unsigned>, ComputeInfo> compute_info_;

        NestAnalysis& analysis_;
        std::stack<problem::OperationSpace> deltas_;
        ScreenShot screen_shot_; 
        problem::OperationSpace computeDelta(const ScreenShot& screen_shot);
        problem::OperationSpace combineDeltas(
            const std::vector<problem::OperationSpace>& deltas, 
            ScopeNode::type_t type);
        
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;
    public:
        DatamovementCalculator(NestAnalysis& analysis): analysis_(analysis){}
        void run(const Node*) override;
        friend class PerfectLoopnestAnalyzer;
    };

    class PerfectLoopnestAnalyzer: public analysis::NestAnalysis {
        DatamovementCalculator& dm;
        void ComputeTemporalWorkingSet(std::vector<analysis::LoopState>::reverse_iterator cur,
                                 analysis::ElementState& cur_state) override;
        void InitPerLevelDimScales() override;
        void InitStorageBoundaries() override;
        problem::OperationSpace ComputeDeltas(
            std::vector<analysis::LoopState>::reverse_iterator cur,
            bool at_boundary);
    public: 
        PerfectLoopnestAnalyzer(DatamovementCalculator& dm):analysis::NestAnalysis(), dm(dm){}
        void init(problem::Workload*, loop::Nest*, 
            std::map<unsigned, std::uint64_t> fanoutX_map,
            std::map<unsigned, std::uint64_t> fanoutY_map);
        void calculateDataMovement();
        
        
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
} // namespace TileFlow 

} // namespace analysis 