#pragma once 

#include <vector> 
#include <chrono>
#include "tileflow/mapper/checker.hpp"
#include "tileflow/loop-analysis/nest-analysis.hpp"

using mapping::TileFlow::Node;
using mapping::TileFlow::OpNode;
using mapping::TileFlow::TileNode;
using mapping::TileFlow::ScopeNode;
using mapping::TileFlow::Visitor;

namespace TileFlow {

namespace mapper {

enum Objective{
    CYCLE, 
    ENERGY
};

typedef num_t Action;

struct State {
    const int C = 1;

    unsigned n_visit = 0;
    double ave_reward = 0;

    State* parent = nullptr;
    SymbolTable symbol_table_;
    int variable_index;

    std::vector<num_t> candidate_factors_;
    std::unordered_map<num_t, State*> children_; 

    void init_factors(const std::vector<Constraint>& constraints_);

    State(const SymbolTable& table, const std::vector<Constraint>& constraints): 
        symbol_table_(table), variable_index(symbol_table_.get_next_var()) {
            init_factors(constraints);
        }
    State(SymbolTable&& table, const std::vector<Constraint>& constraints): 
        symbol_table_(table), variable_index(symbol_table_.get_next_var()) {
            init_factors(constraints);
        }
    ~State(){for (auto child: children_) delete child.second;}


    std::pair<Action, bool> select_action(bool random);

    State* take_action(const Action& action, const std::vector<Constraint>& constraints_);

    inline bool is_terminated() const {return variable_index > 0;}
    inline bool is_error_out() const {
        return variable_index == ERROR_OUT;
    }
};

std::ostream& operator << (std::ostream& o, const State& s);

/**
 * take action:  -> 
 * state::init(constraint) -> next_variable, candidates :next variable, candidates
 * state::select_candidate()
 * candidates -> candidate 
 * state (symbol table) X <next, variable> --> next global_
*/

struct Env {
private:    
    double punish_ = -10;
    const std::vector<Constraint>& constraints_;
    State* root_ = nullptr;
    const SymbolTable* best_symbol_table_;
    analysis::TileFlow::NestAnalysis& analyzer_;
    Objective obj_;

    State* curr_state_ = nullptr;
    bool terminated_ = false;
    bool expanded_ = false;
    double best_reward_ = -1e9;
    double reward_;

public:
    Env(const std::vector<Constraint>& constraints, 
        const SymbolTable& symbol_table, 
        analysis::TileFlow::NestAnalysis& analyzer,
        Objective obj): constraints_(constraints),
    root_(new State(symbol_table, constraints)), best_symbol_table_(nullptr),
        analyzer_(analyzer), obj_(obj) {reset();
        if (verbose_level > 1) {    
        std::cout << "init env..." << std::endl;
        std::cout << *root_ << std::endl;} }

    ~Env() {delete root_;}
    Action step(bool random = false);
    void reset(State* state = nullptr) {
        curr_state_ = state == nullptr? root_:state;
        assert(curr_state_);
        terminated_ = curr_state_->is_terminated(); 
        expanded_ = false;}
    
    bool is_terminated() const {return terminated_;}
    bool is_expanded() const {return expanded_;}
    double get_reward() const {return reward_;}
    State* get_curr_state() const {return curr_state_;}
    const std::vector<Constraint>& get_constraints() const {return constraints_;}
    const SymbolTable* get_best_symbol_table() const {
        return best_symbol_table_;    
    }
};

std::ostream& operator<< (std::ostream& o, const Env& env);

struct MCTS {
private:
    const int n_iteration = 10000;
    const int n_rollout = 100;
    int n_unexplored;
    Env* env_ = nullptr;

    unsigned timeout_;
    std::chrono::steady_clock::time_point begin_;
    void start_timer() {begin_ = std::chrono::steady_clock::now();}
    unsigned get_elapsed_time() {
        auto cur = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::seconds> (cur - begin_).count();
    }

    void back_prop(State* state, double reward);
    double rollout(State* state);
    State* select_state();

public:
    MCTS(Env* env, unsigned timeout = 120): env_(env), timeout_(timeout){
        n_unexplored = env_->get_curr_state()->candidate_factors_.size();
    }
    void search();
};

class Mapper {
    const std::vector<Constraint>& constraints_;
    const problem::TileFlow::Workloads& workloads_;
    const mapping::TileFlow::Mapping& mapping_;
    const model::Engine::Specs& arch_specs_;
    const model::Topology& topology_;
    SymbolTable optimum_;
    Objective obj_;
    unsigned timeout_;

public:
    Mapper(const std::vector<Constraint>& constraints_,
        const problem::TileFlow::Workloads& workloads_, 
        const mapping::TileFlow::Mapping& mapping_, 
        const model::Engine::Specs& arch_specs_, 
        const model::Topology& topology_, 
        Objective obj, unsigned timeout = 600):
        constraints_(constraints_), workloads_(workloads_), 
        mapping_(mapping_), arch_specs_(arch_specs_), 
        topology_(topology_), obj_(obj), timeout_(timeout) {
            global_symbol_table_.init(constraints_);
            std::cout << "begin mapping..." << std::endl; 
        }
    
    const SymbolTable* search();

    void dump(const std::string& filename);

    void report();

};

} // namespace mapper

} // namespace TileFlow 