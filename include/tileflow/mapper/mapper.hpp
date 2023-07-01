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

struct Reward {
    double value;
    double reward;
    Reward(double reward =0.0, double value = 0.0): 
        value(value), reward(reward) {}
    bool operator>(const Reward & other){
        return reward > other.reward;
    }
    bool operator>(double reward_){
        return reward > reward_;
    }
    Reward operator-(const Reward& other) {
        return {reward - other.reward, value - other.value};
    }
    Reward operator+(const Reward& other) {
        return {reward + other.reward, value + other.value};
    }
    Reward operator/(unsigned n) {
        return {reward / n, value / n};
    }
    void operator+=(const Reward& other) {
        reward += other.reward;
        value += other.value;
    }
};

std::ostream& operator<<(std::ostream& o, const Reward& r);

enum Objective{
    CYCLE, 
    ENERGY
};

typedef num_t Action;

struct State {
    const int C = 1;

    unsigned n_visit = 0;
    Reward ave_reward = 0;

    Reward remembered_reward = 1;

    State* parent = nullptr;
    Action last_action = 0;
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

    void Print(std::ostream& o);
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
public: 

private:    
    double punish_ = -20;
    const std::vector<Constraint>& constraints_;
    State* root_ = nullptr;
    const SymbolTable* best_symbol_table_;
    analysis::TileFlow::NestAnalysis& analyzer_;
    Objective obj_;

    State* curr_state_ = nullptr;
    bool terminated_ = false;
    bool expanded_ = false;
    Reward best_reward_ = {-1e9, 0.0};

    unsigned topk_;
    std::vector<std::pair<const SymbolTable*, double> > topKList_;
    void insert(const SymbolTable*, double);

public:
    Env(const std::vector<Constraint>& constraints, 
        const SymbolTable& symbol_table, 
        analysis::TileFlow::NestAnalysis& analyzer,
        Objective obj, 
        unsigned topK = 1): constraints_(constraints),
    root_(new State(symbol_table, constraints)), 
    best_symbol_table_(nullptr), analyzer_(analyzer), obj_(obj),
    topk_(topK) {reset();
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
    Reward calculate_reward();
    Reward get_reward() const {return best_reward_;}
    State* get_curr_state() const {return curr_state_;}
    const std::vector<Constraint>& get_constraints() const {return constraints_;}
    const SymbolTable* get_best_symbol_table() const {
        return best_symbol_table_;    
    }
    std::ostream& report(std::ostream& o) const;
};

std::ostream& operator<< (std::ostream& o, const Env& env);

struct Algorithm {
    struct log_t{
        unsigned time;
        unsigned iter;
        Reward reward;
        Reward best_reward;
        std::string node;
    };
    enum type_t {
        MCTS,
        RANDOM
    };
protected: 
    unsigned timeout_;
    std::chrono::steady_clock::time_point begin_;
    std::vector<log_t> logs_;
    Env* env_ = nullptr;
    void start_timer() {begin_ = std::chrono::steady_clock::now();}
    unsigned get_elapsed_time() {
        auto cur = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::seconds> (cur - begin_).count();
    }
public:
    Algorithm(unsigned timeout = 120): timeout_(timeout) {}
    void report_csv(std::ostream& o) const;
    virtual void set_env(Env* env) {env_ = env;}
    virtual void search() = 0;
};

struct MCTS: public Algorithm {
private:
    const unsigned n_iteration = 10000;
    const unsigned n_rollout = 100;
    int n_unexplored;
    unsigned iter_;
    bool random_;
    void back_prop(State* state, Reward reward);
    Reward rollout(State* state);
    State* select_state();

public:
    MCTS(unsigned timeout = 120, bool random = false): 
        Algorithm(timeout), random_(random){}
    void search() override;
    void set_env(Env* env){env_ = env;}
};

class Mapper {
    const std::vector<Constraint>& constraints_;
    const problem::TileFlow::Workloads& workloads_;
    const mapping::TileFlow::Mapping& mapping_;
    const model::Engine::Specs& arch_specs_;
    const model::Topology& topology_;
    SymbolTable optimum_;
    Objective obj_;
    std::unique_ptr<Algorithm> alg_;
    unsigned timeout_;
    unsigned topk_;
    void report_csv(std::ostream&);
    void report_mapping(std::ostream&);

public:
    Mapper(const std::vector<Constraint>& constraints_,
        const problem::TileFlow::Workloads& workloads_, 
        const mapping::TileFlow::Mapping& mapping_, 
        const model::Engine::Specs& arch_specs_, 
        const model::Topology& topology_, 
        Objective obj, unsigned timeout = 600,
        const std::string alg = "mcts", unsigned topk = 1):
        constraints_(constraints_), workloads_(workloads_), 
        mapping_(mapping_), arch_specs_(arch_specs_), 
        topology_(topology_), obj_(obj),
        timeout_(timeout), topk_(topk) {
            if (alg == "mcts" || alg == "MCTS")
                alg_ = std::make_unique<MCTS>(timeout);
            else if (alg == "random" || alg == "RANDOM")
                alg_ = std::make_unique<MCTS>(timeout, true);
            global_symbol_table_.init(constraints_);
            std::cout << "begin mapping by " << alg <<  "..." << std::endl; 
        }
    const SymbolTable* search();
    void dump(const std::string& filename);
    void report();
};

} // namespace mapper

} // namespace TileFlow 