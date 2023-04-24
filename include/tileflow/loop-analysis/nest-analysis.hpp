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
#include "tileflow/mapper/expr.hpp"

using mapping::TileFlow::Node;
using mapping::TileFlow::OpNode;
using mapping::TileFlow::TileNode;
using mapping::TileFlow::ScopeNode;
using mapping::TileFlow::Visitor;

using TileFlow::verbose_level;

namespace analysis {

namespace TileFlow {
    struct DataMovement {
    private:
        // static const std::map<std::string, std::string> metrics_;
        std::unordered_map<std::string, double> data_;
    public: 
        double& operator[](const std::string& metric) {
            if (data_.count(metric) == 0) data_[metric] = 0.0;
            return data_[metric];
        }
        void report(std::ostream&, const std::string&) const;
    };

    struct DataMovements {
    private: 
        std::unordered_map<std::string, DataMovement> data_;
    public: 
        DataMovement& operator[](const std::string & metric) {
            return data_[metric];
        }
        void clear() {data_.clear();}
        void report(std::ostream&) const;
    };

    struct stat_t {
        problem::PerDataSpace<AccessStatMatrix> access_stat_;
        problem::PerDataSpace<size_t> max_size_;
        problem::PerDataSpace<std::uint64_t> link_transfer_;
        ComputeInfo compute_info_;
    };

    struct NodeConfig {
        // problem::Workload workload;
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
        const SymbolTable* symb_table_ = nullptr;

        std::uint64_t cycle_;
        double energy_;
        DataMovements data_movements_;
        
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
        
        void set_symbol_table(const SymbolTable* symbol_table) {symb_table_ = symbol_table;}
        
        const tiling::CompoundTile& get_tile(const Node* node) const 
        {
            if (tiles_.count(node) == 0)    {
                std::cerr << "ERROR node not found:" << std::endl;
                node->display("", false);
            }
            return tiles_.at(node);
        }
        
        void reset();
        void analyze();
        void Print(std::ostream& o = std::cout);
        void Report();
        void Export(const std::string& filename);
        std::uint64_t get_cycle() const {return cycle_;}
        double get_energy() const {return energy_;}
        const DataMovements& get_data_movements() const {return data_movements_;}
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
        void finalizeStat(unsigned storage_id, RetVal& ret, bool profile);
        
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

    class Displayer: public mapping::TileFlow::Visitor {
        NestAnalysis & analysis_;
        const SymbolTable* symbol_table_;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        void visitOp(const OpNode*) override;
        std::string prefix_;
        std::ostream& o_;
    public: 
        Displayer(NestAnalysis& analysis, const SymbolTable* symbol_table_ = nullptr, std::ostream& o = std::cout)
        : analysis_(analysis), symbol_table_(symbol_table_), o_(o) {}
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
        NestAnalysis& analysis_;
        const SymbolTable* symbol_table_;

        void visitTile(const TileNode* node) override { 
            analysis_.configs[node].loop_nest = node->constructLoopNest(symbol_table_);
            for (auto child: node->get_children()) {child->accept(this);}
        }
    public:
        LoopNestConstructor(NestAnalysis& analysis, const SymbolTable*symbol_table_)
        : analysis_(analysis), symbol_table_(symbol_table_){}
        void construct(const Node* root) {root->accept(this);}
    };

    class StorageLevelCalculator: public mapping::TileFlow::Visitor {
        std::stack<unsigned> storage_levels_;
        void visitTile(const TileNode*) override;
        void visitScope(const ScopeNode*) override;
        NestAnalysis& analysis_;
    public: 
        StorageLevelCalculator(NestAnalysis& analysis): analysis_(analysis){
            if (verbose_level > 1) {
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

} // namespace TileFlow 

} // namespace analysis 