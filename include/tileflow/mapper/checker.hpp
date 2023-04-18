#pragma once 

#include <vector> 

#include "tileflow/mapping/mapping.hpp"

using mapping::TileFlow::Node;
using mapping::TileFlow::OpNode;
using mapping::TileFlow::TileNode;
using mapping::TileFlow::ScopeNode;
using mapping::TileFlow::Visitor;

namespace TileFlow {

    struct Constraint {
        std::shared_ptr<Expr> expr;
        std::string msg;
    };

    class ShapeConstraintParser: public Visitor {
        void visitTile(const TileNode*) override;
        void visitOp(const OpNode*) override;
        std::vector<std::pair<int, std::vector<std::shared_ptr<Expr> > > > 
            exprs_;
        const problem::Workload& workload_;
        std::vector<Constraint > constraints;
    public:
        ShapeConstraintParser(const problem::Workload& workload):
            workload_(workload){}
        std::vector<Constraint > parse(const Node*root);

    };

    class MemoryConstraintParser: public Visitor {
        void visitTile(const TileNode*) override;
        void visitOp(const OpNode*) override;
        void visitScope(const ScopeNode*) override;
        std::shared_ptr<Expr> cal_footprint(unsigned pv);
        std::vector<std::pair<int, std::vector<int> > > factors_;
        const problem::Workload& workload_;
        const model::Topology& topology_;
        std::vector<Constraint> constraints_;
    public:
        MemoryConstraintParser(
            const problem::Workload& workload, 
            const model::Topology& topology): 
            workload_(workload), topology_(topology){}
        std::vector<Constraint > parse(const Node*root);
    };

    class SpatialScopeSwapper: public mapping::TileFlow::Visitor {
        void visitScope(const ScopeNode*) override;
    };

    /**
     * rule1: Tiling fators's multiplication should be equal to the shape;
     * rule2: Spatial TileNode's child must be a Temporal Tile Node
     * rule3: each level should have at most one temporal tile node;
    */
    class SanityChecker: public mapping::TileFlow::Visitor {
        unsigned storage_level_;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override; 
        const model::Topology& topology_;
    public: 
        SanityChecker(const model::Topology& topology): topology_(topology){}
        void run(const Node*) override; 
    };

    class Checker { 
    private: 
        const problem::TileFlow::Workloads& workloads_;
        const mapping::TileFlow::Mapping& mapping_;
        const model::Topology& topology_;
        std::vector<Constraint > constraints;
        void swap_spatial_scope();
        void get_active_tensors();
        void sanity_check();
        void get_shape_constraints();
        void get_memory_constraints();
        void add_access_pattern(
            problem::Shape::DataSpaceID producer_id,
            const Node *producer,
            problem::Shape::DataSpaceID consumer_id,
            const Node *consumer);
    public:
        Checker(const problem::TileFlow::Workloads& workloads,
            const mapping::TileFlow::Mapping& mapping,
            const model::Topology& topology)
            : workloads_(workloads), mapping_(mapping), topology_(topology){}
        const std::vector<Constraint>& get_constraints() const {return constraints;}
        void check();
        void display();
    }; // Mappert 

} // namespace TileFlow 