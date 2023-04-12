#pragma once 

#include <bitset>
#include <stack>

#include "loop-analysis/nest-analysis.hpp"
#include "loop-analysis/loop-state.hpp"
#include "model/topology.hpp"


#include "tileflow/mapping/mapping.hpp"
#include "tileflow/common.hpp"
#include "tileflow/problem/problem.hpp"
#include "tileflow/loop-analysis/memory-state.hpp"

using mapping::TileFlow::Node;
using mapping::TileFlow::OpNode;
using mapping::TileFlow::TileNode;
using mapping::TileFlow::ScopeNode;
using mapping::TileFlow::Visitor;

using TileFlow::verbose_level;

namespace analysis {

namespace TileFlow {

    struct stat_t {
        problem::PerDataSpace<AccessStatMatrix> access_stat_;
        problem::PerDataSpace<size_t> max_size_;
        problem::PerDataSpace<std::uint64_t> link_transfer_;
        ComputeInfo compute_info_;
    };

    struct NodeConfig {
        // problem::Workload workload;
        // tensors cause transfer by read (to lower level)
        std::vector<problem::Shape::DataSpaceID> active_read_tensors;
        // tensor that cause transfer by fill (to upper level)
        std::vector<problem::Shape::DataSpaceID> active_fill_tensors;
        loop::Nest loop_nest;
        std::vector<problem::OperationPoint> vector_strides_;
        std::vector<problem::OperationPoint> mold_low_;
        std::vector<problem::OperationPoint> mold_high_;
        std::vector<problem::OperationPoint> mold_high_residual_;
        // the state to record the last_point_set and access counting 

        std::map<std::vector<unsigned>, stat_t> stats_;

        // StorageLevelCalculator
        unsigned storage_level;
        std::uint64_t fanout_x, fanout_y;

        // SpatialOffsetsCalculator
        std::uint64_t spatial_offset_x, spatial_offset_y, logical_x, logical_y, logical_fanout;
        std::uint64_t replication_factor;
        std::uint64_t max_x_expansion, max_y_expansion;
    };

    class NestAnalysis {
        const problem::TileFlow::Workloads& workloads_;
        const mapping::TileFlow::Mapping& mapping_;
        const model::Engine::Specs& arch_specs_;
        const model::Topology& topology_;
        const problem::Workload& common_workload_;

        std::uint64_t cycle_;
        double energy_;
        
        void add_access_pattern(
            problem::Shape::DataSpaceID producer_id, 
            const Node* producer, 
            problem::Shape::DataSpaceID consumer_id,
            const Node* consumer);
        /**
         * \brief Swap the spatial node and its scope child
        */
        inline void swap_spatial_scope();
        /**
         * \brief sanity check
        */
       inline void sanity_check();
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

        /**
         * \brief compute maximum expansion 
        */
        void get_expansion();

        /**
         * \brief Aggregate Datamovement/compute info
        */
        void collect_info();

        std::unordered_map<const Node*, NodeConfig> configs;

        // the interface
        std::unordered_map<const Node*, tiling::CompoundTile> tiles_; 

    public: 
        NestAnalysis(const problem::TileFlow::Workloads& workloads_, 
            const mapping::TileFlow::Mapping& mapping_, 
            const model::Engine::Specs& arch_specs_, 
            const model::Topology& topology_);
        
        const tiling::CompoundTile& get_tile(const Node* node) const 
            {
                if (tiles_.count(node) == 0)    {
                    std::cerr << "ERROR node not found:" << std::endl;
                    node->display("", false);
                }
                return tiles_.at(node);
            }

        void analyze();
        void Print();
        void Report();
        void Export(const std::string& filename);
        friend class Displayer;
        friend class DatamovementCalculator;
        friend class DimScaleCalculator;
        friend class LoopNestConstructor;
        friend class PerfectLoopnestAnalyzer;
        friend class StorageLevelCalculator;
        friend class SpatialOffsetsCalculator;
        friend class ComputeExpansion;
        friend class SanityChecker;
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


    /**
     *                      Temporal        [Spatial]       Collect
     * Tile              new; Read/Write   Fill/Link  Eval -> time = max(read, fill, update)
     * Delta/LastState        new          expand         forward
     * Time                   forward       forward       eval;new 
    */
   /**
    * Sequential: Child -> Sigma_child max(read/write/fill) 
    * Pipeline: time = maxs(max(read/write/update))
    * Parallel: time = max(max(read/write/update))
   */

    struct RetVal{
        // Ret as the delta between temporal tiles
        MemoryState deltas_;
        // Always a ret value
        MemoryState last_working_set_;
        // Ret between Tile node and CollectNode;
        std::shared_ptr<tiling::CompoundTile> p_tile_;

        // Ret between CollectNode and parent; 
        problem::PerDataSpace<AccessStatMatrix> access_stat_;
        std::uint64_t cycle_ = 0;

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
        const problem::Workload& workload_;
        const model::Topology & topology_;
        /**
         * \brief the stack to pass parameter between nodes.
        */
        std::stack<InputParam> input_stack_;
        std::stack<RetVal> ret_stack_;

        double energy_;
        // const Node* curr_node_;
        RetVal computeDelta(const InputParam& input);
        void finalizeStat(unsigned storage_id, RetVal& ret);
        
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;

        bool break_on_failure;
    
    public:
        DatamovementCalculator(NestAnalysis& analysis): 
            analysis_(analysis), workload_(analysis_.common_workload_),
            topology_(analysis.topology_){}
        
        RetVal eval(const Node*);
        double Energy() {return energy_;}
        friend class PerfectLoopnestAnalyzer;
    };

    class PerfectLoopnestAnalyzer: public analysis::NestAnalysis {
        DatamovementCalculator& dm_;
        const InputParam input_;
        NodeConfig& config_;
        RetVal ComputeTemporalWorkingSet();
        RetVal ComputeSpatialWorkingSet();
        int SimulateTemporalExecution(
            problem::PerDataSpace<AccessStatMatrix>& access_stat);
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

        bool singleton = true;
        
        // utils for compute link transfers;

        void ComputeStats(
            const MemoryState& last_working_set,
            const MemoryState& deltas,
            problem::PerDataSpace<AccessStatMatrix>& access_stats, 
            problem::PerDataSpace<std::uint64_t>& link_transfers
        );

        void ComputeLinkTransfer(
            const MemoryState& last_working_set, 
            const MemoryState& delta,
            problem::PerDataSpace<std::uint64_t>& link_transfers,
            problem::PerDataSpace<std::unordered_set<std::uint64_t>>& unaccounted_delta  
        );

        void ComputeAccessStat(
            const MemoryState& delta, 
            problem::PerDataSpace<std::unordered_set<std::uint64_t>>& unaccounted_delta, 
            problem::PerDataSpace<AccessStatMatrix>& access_stat,
            bool enable_multicast = false
        );

        void ComputePeerAccesses(
            const problem::PerDataSpace<std::uint64_t>& link_transfer,
            tiling::CompoundTile& tile);

        void InitTile(
            const problem::PerDataSpace<AccessStatMatrix>& access_stat,
            tiling::CompoundTile& tile
        );

        void ComputeParentShareAccess(
            const problem::PerDataSpace<AccessStatMatrix>& access_stats,
            tiling::CompoundTile& tile
        );
        void ComputeFill(tiling::CompoundTile& tile);

        void FinalizeTile(tiling::CompoundTile& tile);

        void ComputeReadUpdate(tiling::CompoundTile& tile);

        void ComputeDensityModels(tiling::CompoundTile& tile);

    public: 
        PerfectLoopnestAnalyzer(
            DatamovementCalculator& dm,
            InputParam& input, 
            NodeConfig& config): dm_(dm), input_(input), config_(config){}
        void init(const problem::Workload*, const loop::Nest*);
        RetVal calculateDataMovement();
        
        
    };

    class CollectOpNode: public mapping::TileFlow::Visitor {
        void visitOp(const OpNode*) override;
        std::vector<const OpNode*> opnodes_;
    public: 
        std::vector<const OpNode*> collectOpNodes(Node*);
    };
    
    class CollectTileNode: public mapping::TileFlow::Visitor {
        void visitTile(const TileNode* node) override {
            if (node->get_tile_type() == type_)
                nodes_.push_back(node);
            for (auto child: node->get_children())
                child->accept(this);
        }
        std::vector<const TileNode*> nodes_;
        TileNode::type_t type_;
    public: 
        CollectTileNode(TileNode::type_t type = TileNode::Temporal): type_(type){}
        std::vector<const TileNode*> operator() (const Node*root){
            root->accept(this);
            return nodes_;
        }
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
            analysis_.configs[node].loop_nest = node->constructLoopNest();
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
            if (verbose_level) {
                std::cout << "Begin storge level calculation..." << std::endl;
                std::cout << "\tfanoutX: ";
                for (auto& kv: analysis_.mapping_.fanoutX_map)
                    std::cout << kv.first << ":" << kv.second << ",";
                std::cout << std::endl;
                std::cout << "\tfanoutY: ";
                for (auto& kv: analysis_.mapping_.fanoutY_map)
                    std::cout << kv.first << ":" << kv.second << ",";
                std::cout << std::endl;
                
            }
        }
    };

    /**
     * rule1: Tiling fators's multiplication should be equal to the shape;
     * rule2: Spatial TileNode's child must be a Temporal Tile Node
     * rule3: each level should have at most one temporal tile node;
    */
    class SanityChecker: public mapping::TileFlow::Visitor {
        std::unordered_map<int, int> scales_;
        unsigned storage_level_;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override; 
        NestAnalysis& analysis_;
    public: 
        SanityChecker(NestAnalysis& analysis_): analysis_(analysis_){}
        void run(const Node*) override; 
    };

    class SpatialOffsetsCalculator: public mapping::TileFlow::Visitor {
    public:
        struct offset_t {
            unsigned x, y;
            unsigned max_x;
        };
    private:
        offset_t merge(const offset_t& o1, const offset_t& o2);
        int replication_factor = 1;
        std::stack<offset_t> input_offsets, output_offsets;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;
        NestAnalysis& analysis_;
    public: 
        SpatialOffsetsCalculator(NestAnalysis& analysis): analysis_(analysis){
        }
    };    

    std::ostream& operator << (std::ostream& o, const SpatialOffsetsCalculator::offset_t& offset);

    class ComputeExpansion: public mapping::TileFlow::Visitor {
        void visitTile(const TileNode*) override;
        void visitOp(const OpNode*) override;
        NestAnalysis& analysis_;
        std::pair<int, int> expansion_;
    public: 
        ComputeExpansion(NestAnalysis& analysis): analysis_(analysis),
        expansion_(1,1) {}
    };

    class SpatialScopeSwapper: public mapping::TileFlow::Visitor {
        void visitScope(const ScopeNode*) override;
    };

} // namespace TileFlow 

} // namespace analysis 