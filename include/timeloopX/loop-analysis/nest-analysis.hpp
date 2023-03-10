#pragma once 

#include <bitset>
#include <stack>

#include "loop-analysis/nest-analysis.hpp"
#include "loop-analysis/loop-state.hpp"

#include "timeloopX/mapping/mapping.hpp"
#include "timeloopX/common.hpp"
#include "timeloopX/problem/problem.hpp"

using mapping::TimeloopX::Node;
using mapping::TimeloopX::OpNode;
using mapping::TimeloopX::TileNode;
using mapping::TimeloopX::ScopeNode;
using mapping::TimeloopX::Visitor;


namespace analysis {

namespace TimeloopX {

    struct NodeConfig {
        problem::Workload workload;
        std::unordered_map<problem::Shape::DataSpaceID, 
            std::vector<analysis::DataMovementInfo> > data_movements;
    };

    class NestAnalysis {
        // opname X tensorname --> AccessPattern
        std::unordered_map<std::string,
            std::unordered_map<std::string, 
            problem::Shape::Projection> >  access_patterns;

        problem::TimeloopX::Workloads& workloads_;
        mapping::TimeloopX::Mapping& mapping_;
        
        void add_access_pattern(
            problem::Shape::DataSpaceID producer_id, 
            const Node* producer, 
            problem::Shape::DataSpaceID consumer_id,
            const Node* consumer,
            std::unordered_map<const Node*, std::vector<problem::Shape::DataSpaceID> >& access_pattern);
    public: 
        NestAnalysis(problem::TimeloopX::Workloads& workloads_, mapping::TimeloopX::Mapping& mapping): workloads_(workloads_), mapping_(mapping){}
        std::unordered_map<const Node*, NodeConfig> configs;
        /**
         * \brief set the alive_consors for Tile Nodes
         */
        void get_tilewise_workloads();
        void Print();
        friend class Displayer;
        friend class DatamovementCalculator;
    };

    struct ScreenShot{
            int num_epochs_;
    };

    class DatamovementCalculator: public mapping::TimeloopX::Visitor {
        NestAnalysis& analysis_;
        std::stack<problem::OperationSpace> deltas_;
        ScreenShot screen_shot_; 
        problem::OperationSpace computeDelta(Node* node);
        problem::OperationSpace combineDeltas(
            const std::vector<problem::OperationSpace>& deltas, 
            ScopeNode::type_t type);

        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;
        friend class PerfectLoopnestAnalyzer;
    };

    class PerfectLoopnestAnalyzer: public analysis::NestAnalysis {
        DatamovementCalculator& dm;
    public: 
        PerfectLoopnestAnalyzer(DatamovementCalculator& dm): dm(dm){}
        void init(problem::Workload*, loop::Nest*, 
            std::map<unsigned, std::uint64_t> fanoutX_map,
            std::map<unsigned, std::uint64_t> fanoutY_map);
        void calculateDataMovement();
        
        
    };

    

    class CollectOpNode: public mapping::TimeloopX::Visitor {
        void visitOp(const OpNode*) override;
        std::vector<const OpNode*> opnodes_;
    public: 
        std::vector<const OpNode*> collectOpNodes(Node*);
    };

    class Displayer: public mapping::TimeloopX::Visitor {
        NestAnalysis & analysis_;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;
        std::string prefix_;
    public: 
        Displayer(NestAnalysis& analysis): analysis_(analysis){}
        void display() {analysis_.mapping_.root->accept(this);}
    };


} // namespace TimeloopX 

} // namespace analysis 