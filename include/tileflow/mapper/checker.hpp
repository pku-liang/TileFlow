#pragma once 

#include <vector> 

#include "tileflow/mapping/mapping.hpp"

using mapping::TileFlow::Node;
using mapping::TileFlow::OpNode;
using mapping::TileFlow::TileNode;
using mapping::TileFlow::ScopeNode;
using mapping::TileFlow::Visitor;

namespace TileFlow {
    
    class ShapeConstraintParser: public Visitor {
        void visitTile(const TileNode*) override;
        void visitOp(const OpNode*) override;
        std::vector<std::pair<int, std::vector<std::shared_ptr<Expr> > > > 
            exprs_;
        const problem::Workload& workload_;
        bool allow_mismatched_;
        
        std::vector<Constraint > constraints;
    public:
        ShapeConstraintParser(const problem::Workload& workload, bool allow_mismatched = false):
            workload_(workload), allow_mismatched_(allow_mismatched){}
        std::vector<Constraint > parse(const Node*root);

    };

    class MemoryConstraintParser: public Visitor {
        void visitTile(const TileNode*) override;
        void visitOp(const OpNode*) override;
        void visitScope(const ScopeNode*) override;
        std::shared_ptr<Expr> cal_footprint(const Node * node, unsigned pv);
        // std::vector<std::pair<int, std::vector<int> > > factors_;
        std::vector<std::shared_ptr<Expr> > factors_; 
        std::unordered_map<const Node*, std::vector<std::shared_ptr<Expr> > >  node2factors_;
        const problem::Workload& workload_;
        const model::Topology& topology_;
        std::vector<Constraint> constraints_;
        std::vector<const Node*> constraint_nodes_;

        std::unordered_map<
        const Node*, std::unordered_map<unsigned, const Node*> >& init_scope_;
        void add_constraint(const Node*node);
        std::vector<std::shared_ptr<Expr> > factor2expr(
            const std::vector<std::pair<int, std::vector<int> > >& factors);
    public:
        MemoryConstraintParser(
            const problem::Workload& workload, 
            const model::Topology& topology,
            std::unordered_map<
            const Node*, std::unordered_map<unsigned, const Node*> >& init_scope): 
            workload_(workload), topology_(topology), init_scope_(init_scope){}
        std::vector<Constraint> parse(const Node*root);
    };

    class ResourceConstraintParser: public Visitor {
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        std::vector<Constraint> constraints_;
        std::shared_ptr<ResourceExpr> core_usage_;
        const mapping::TileFlow::Mapping& mapping_;
        void add_constraint(const Node* node);
    public:
        ResourceConstraintParser(const mapping::TileFlow::Mapping& mapping):
        mapping_(mapping) {}
        std::vector<Constraint> parse(const Node*root);

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
        bool constraints_parsed_ = false;
        const problem::TileFlow::Workloads& workloads_;
        const mapping::TileFlow::Mapping& mapping_;
        const model::Topology& topology_;
        bool enable_mem_check_;
        bool enable_spatial_check_;
        bool enable_loopcount_check_;

        std::unordered_map<
            const Node*, std::unordered_map<unsigned, const Node*> > init_scope_;

        std::vector<Constraint > constraints;
        void swap_spatial_scope();
        void get_active_tensors();
        void sanity_check();
        void parse_constraints();
        void get_shape_constraints();
        void get_memory_constraints();
        void get_resource_constraints();
        void add_access_pattern(
            problem::Shape::DataSpaceID producer_id,
            const Node *producer,
            problem::Shape::DataSpaceID consumer_id,
            const Node *consumer);
    public:
        Checker(const problem::TileFlow::Workloads& workloads,
            const mapping::TileFlow::Mapping& mapping,
            const model::Topology& topology,
            bool enable_mem_check_ = true,
            bool enable_spatial_check_ = true,
            bool enable_loopcount_check_ = true)
            : workloads_(workloads), mapping_(mapping), topology_(topology), 
            enable_mem_check_(enable_mem_check_), 
            enable_spatial_check_(enable_spatial_check_),
            enable_loopcount_check_(enable_loopcount_check_){}
        const std::vector<Constraint>& get_constraints() const {return constraints;}
        void check(const SymbolTable* symbol_table = nullptr);
        void display(const SymbolTable* symbol_table = nullptr);
    }; // Mappert 

} // namespace TileFlow 