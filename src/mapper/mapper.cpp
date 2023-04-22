#include <fstream>
#include <random>

#include "tileflow/mapper/mapper.hpp"
#include "tileflow/mapper/checker.hpp"

namespace TileFlow {

namespace mapper {

const SymbolTable* Mapper::search() {
    std::string obj = obj_ == Objective::CYCLE? "cycle" : "energy";
    TILEFLOW_LOG("Optimize " << obj << "...");
    analysis::TileFlow::NestAnalysis analyzer(workloads_, mapping_, arch_specs_, topology_);
    Env env(constraints_, global_symbol_table_, analyzer, obj_);
    MCTS mcts(&env, timeout_);
    mcts.search();
    auto ret = env.get_best_symbol_table();
    TILEFLOW_ASSERT(ret, "no candidate found");
    TILEFLOW_LOG("best factors: "; ret->show_brief(std::cout); std::cerr);
    optimum_ = *ret;
    return &optimum_;
}

void Mapper::dump(const std::string & filename){
    analysis::TileFlow::NestAnalysis analysis(workloads_, mapping_, arch_specs_, topology_);
    analysis.set_symbol_table(&optimum_);
    analysis.analyze();
    std::ofstream file;
    size_t tmp = filename.find_last_of('.');
    if (tmp == std::string::npos || filename.substr(tmp) != ".csv")
        file.open(filename + ".csv");
    else file.open(filename);
    file << "metric,value" << std::endl;
    file << "Cycle," << analysis.get_cycle() << std::endl;
    file << "Energy," << analysis.get_energy() << std::endl;
    TILEFLOW_LOG("result written into " << filename);        
    file.close();
}

void Mapper::report() {
    analysis::TileFlow::NestAnalysis analysis(workloads_, mapping_, arch_specs_, topology_);
    analysis.set_symbol_table(&optimum_);
    analysis.analyze();
    std::cout << "***TileFlow Result" << std::endl;
    std::cout << ",value" << std::endl;
    std::cout << "Cycle," << analysis.get_cycle() << std::endl;
    std::cout << "Energy," << analysis.get_energy() << std::endl;
    for (auto & constraint: constraints_) {
        if(constraint.type_ == Constraint::MEM) {
            std::cout << "MEM::" << constraint.short_msg << ",";
            auto expr = std::static_pointer_cast<CondExpr>(constraint.expr);
            std::cout << expr->left_->eval(optimum_) / (expr->right_->eval(optimum_) + 0.0) << std::endl;
        }
        else if (constraint.type_ == Constraint::SPATIAL) {
            auto expr = std::static_pointer_cast<PairCondExpr>(constraint.expr);
            auto limit = expr->limit_->eval_pair(optimum_, 0);
            auto usage = expr->expr_->eval_pair(optimum_, limit.second);
            std::cout << "SPATIAL::" << constraint.short_msg << ",";
            std::cout << (usage.first * usage.second + 0.0) / (limit.first * limit.second) << std::endl;
        }
    }
    std::cout << "***TileFlow Result Ends" << std::endl;
}

Action Env::step(bool random) {
    assert(!terminated_);
    assert(curr_state_->candidate_factors_.size());
    Action act; 
    std::tie(act, expanded_) = curr_state_->select_action(random);
    curr_state_ = curr_state_->take_action(act, constraints_);
    terminated_ = curr_state_->is_terminated();
    if (terminated_) {
        if (curr_state_->n_visit > 0) {
            reward_ = curr_state_->ave_reward;
        }
        else if (curr_state_->is_error_out()) {
            reward_ = punish_;
        }
        else 
        {
            analyzer_.set_symbol_table(&curr_state_->symbol_table_);
            if (verbose_level > 1) {
                std::cout << "begin analyze..." << std::endl;
                std::cout << "symbol table: ";
                curr_state_->symbol_table_.show_brief(std::cout);
                std::cout << std::endl;
                analyzer_.Print();
            }
            analyzer_.analyze();
            double value; 
            if (obj_ == Objective::CYCLE) value = (double)(analyzer_.get_cycle());
            else if (obj_ == Objective::ENERGY) value = (double)(analyzer_.get_energy());
            reward_ = -std::log10(value);
            
            if (reward_ > best_reward_) {
                TILEFLOW_LOG("Update best "; curr_state_->symbol_table_.show_brief(std::cerr); std::cerr 
                    << " value: " <<  value);
                best_reward_ = reward_;
                best_symbol_table_ = &curr_state_->symbol_table_;
            }
        }
        if (punish_ + 2 > reward_) 
            punish_ = reward_ - 2;
    }
    return act;
}

State* MCTS::select_state() {
    // std::cout << "---beg state selection---" << std::endl;
    env_->reset();
    // std::cout << *env_->get_curr_state(); 
    while (!env_->is_terminated() && !env_->is_expanded()) {
        env_->step();
        // Action act = env_->step();
        // std::cout << "\tfix to " << act << std::endl;
        // std::cout << *env_->get_curr_state(); 
    }
    auto curr_state = env_->get_curr_state();
    if (env_->is_expanded()) {
        n_unexplored +=
            curr_state->candidate_factors_.size()-1;
        // TILEFLOW_LOG("Expand state has " << curr_state->candidate_factors_.size() << " child.");
    }
    // std::cout << "---end state selection---" << std::endl << std::endl;
    return curr_state;
}

double MCTS::rollout(State* state) {
    double ave_reward = 0.0; 
    for (int i = 0; i < n_rollout; ++i) {
        env_->reset(state);
        while(!env_->is_terminated()) {
            env_->step(true);
        }
        ave_reward += (env_->get_reward() - ave_reward) / (i+1);
    }
    return ave_reward;
}

void MCTS::back_prop(State* state, double reward){
    while (state) {
        state->n_visit ++;
        state->ave_reward += (0.0 + reward - state->ave_reward) / state->n_visit;
        state = state->parent;
    }
}

void MCTS::search() {
    start_timer();
    for (int i = 0; i < n_iteration; i++ ) {
        if (!n_unexplored) {
            if (verbose_level) 
                std::cout << "MCTS: finished searching after " << i << std::endl;
            return;
        }
        auto state = select_state();
        double reward = rollout(state);
        back_prop(state, reward); 
        if (get_elapsed_time() > timeout_) {
            TILEFLOW_LOG("MCTS exit becuase of timeout.");
            return;
        }
    }
}

std::pair<Action, bool> State::select_action(bool random) {
    if (random) {
        auto act = candidate_factors_[std::rand() % candidate_factors_.size()];
        return {act, children_[act] == nullptr};
    }

    double score = -1e9;
    Action act;

    for (auto factor: candidate_factors_) {
        auto& child = children_[factor];
        if (!child || child->n_visit == 0) {
            return {factor, true};            
        }
        double cur_score = child->ave_reward + C * std::sqrt(std::log(n_visit) / child->n_visit);
        if (cur_score > score) {
            act = factor;
            score = cur_score;
        } 
    }
    return {act, false};
}

State* State::take_action(const Action& action,
    const std::vector<Constraint>& constraints_){
    assert(children_.count(action));
    auto& child = children_[action];
    if (!child) {
        SymbolTable table = symbol_table_;
        table.fix_and_update(variable_index, action, constraints_);
        child = new State(std::move(table), constraints_);
        child->parent = this;
    }
    return child;
}

void State::init_factors(const std::vector<Constraint>& constraints){
    TILEFLOW_ASSERT(variable_index > 0 || symbol_table_.count(variable_index),
    "bad variable index " << variable_index);
    if (variable_index > 0) return;

    auto& entry = symbol_table_[variable_index];
    entry.fixed_ = true;
    for (auto it = entry.candidates_.rbegin(); 
        it != entry.candidates_.rend();) {
        entry.value_ = *it;
        bool failed = false;
        for (auto& cons: constraints) {
            if (cons.type_ == Constraint::MEM || cons.type_ == Constraint::SPATIAL){
                if(!cons.expr->eval(symbol_table_)) {
                    failed = true;
                    break;
                }
            }
        }
        if (failed) {
            it = decltype(it)(entry.candidates_.erase( std::next(it).base() ));
        } else {
            break;
        } 
    }
    entry.fixed_ = false;
    for (auto& factor: entry.candidates_) {
        candidate_factors_.push_back(factor);
        children_[factor] = nullptr;
    }
}

std::ostream& operator << (std::ostream& o, const State& s) {
    o << "\t=========State==========" << std::endl;
    o << "\tterminated: " << s.is_terminated() << std::endl;
    o << "\terror out: " << s.is_error_out() << std::endl;
    o << "\tn_visit:" << s.n_visit << std::endl;
    o << "\tave_reward: " << s.ave_reward << std::endl;
    if (!s.is_terminated())  
        o << "\tEntry: " << s.symbol_table_.lookup(s.variable_index) << std::endl;
    o << "\tCandidates: ";
    for (auto factor: s.candidate_factors_) o << factor << ",";
    o << std::endl;
    o << "\t=======End State========" << std::endl;
    return o;
}

std::ostream& operator<< (std::ostream& o, const Env& env) {
    o << "\t============Env===========" << std::endl;
    o << "\tConstraint: ";
    for (auto& cons: env.get_constraints()) {
        o << "\t\t"; 
        cons.expr->display(env.get_curr_state()->symbol_table_);
        o << std::endl;
    }
    o << *env.get_curr_state();
    o << env.get_curr_state()->symbol_table_;
    o << "\t==========End Env=========" << std::endl;
    return o;
}

}

} // namespace TileFlow 