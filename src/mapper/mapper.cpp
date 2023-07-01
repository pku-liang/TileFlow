#include <fstream>
#include <random>

#include "tileflow/mapper/mapper.hpp"
#include "tileflow/mapper/checker.hpp"

namespace TileFlow {

namespace mapper {

void Algorithm::report_csv(std::ostream& o) const {
    o << "index,time,iter,reward,value,best_reward,best_value,node" << std::endl;
    int idx = 0;
    for (auto& log: logs_) {
        o << idx++ << "," << log.time << "," << log.iter 
        << "," << log.reward.reward << "," << log.reward.value
        << "," << log.best_reward.reward << "," << log.best_reward.value
        << "," << log.node << std::endl;
    }
}

const SymbolTable* Mapper::search() {
    std::string obj = obj_ == Objective::CYCLE? "cycle" : "energy";
    TILEFLOW_LOG("Optimize " << obj << "...");
    analysis::TileFlow::NestAnalysis analyzer(workloads_, mapping_, arch_specs_, topology_);
    if (!global_symbol_table_.count_unfixed()){
        optimum_ = global_symbol_table_;
        return &optimum_;
    }
    Env env(constraints_, global_symbol_table_, analyzer, obj_, topk_);
    alg_->set_env(&env);
    alg_->search();
    if(verbose_level) env.report(std::cerr);
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
    std::string prefix = filename.substr(0, filename.find_last_of("."));

    std::ofstream file;
    file.open(prefix + ".csv");
    report_csv(file);
    file.close();
    TILEFLOW_LOG("profiling written into " << prefix << ".csv");        

    file.open(prefix + ".mapping.txt");
    report_mapping(file);
    file.close();

    file.open(prefix + ".tuning.csv");
    alg_->report_csv(file);
    file.close(); 

    TILEFLOW_LOG("mapping written into " << prefix << ".mapping.txt");
}

void Mapper::report() {
    report_mapping(std::cout);
    std::cout << "***TileFlow Result" << std::endl;
    report_csv(std::cout);
    std::cout << "***TileFlow Result Ends" << std::endl;
}

void Mapper::report_mapping(std::ostream& o){
    o << "***Optimal Mapping:" << std::endl;
    analysis::TileFlow::NestAnalysis analysis(workloads_, mapping_, arch_specs_, topology_);
    analysis.set_symbol_table(&optimum_);
    analysis.analyze();
    analysis.Print(o);
    // mapping_.root->display("", true, &optimum_, o);
}

void Mapper::report_csv(std::ostream& o) {
    analysis::TileFlow::NestAnalysis analysis(workloads_, mapping_, arch_specs_, topology_);
    analysis.set_symbol_table(&optimum_);
    analysis.analyze();
    o << "metric,value" << std::endl;
    o << "Cycle," << analysis.get_cycle() << std::endl;
    o << "Energy," << analysis.get_energy() << std::endl;
    analysis.get_data_movements().report(o);
    for (auto & constraint: constraints_) {
        if(constraint.type_ == Constraint::MEM) {
            o << "MEM::" << constraint.short_msg << ",";
            auto expr = std::static_pointer_cast<CondExpr>(constraint.expr);
            o << expr->left_->eval(optimum_) / (expr->right_->eval(optimum_) + 0.0) << std::endl;
        }
        else if (constraint.type_ == Constraint::SPATIAL) {
            auto expr = std::static_pointer_cast<PairCondExpr>(constraint.expr);
            auto limit = expr->limit_->eval_pair(optimum_, 0);
            auto usage = expr->expr_->eval_pair(optimum_, limit.second);
            o << "SPATIAL::" << constraint.short_msg << ",";
            o << (usage.first * usage.second + 0.0) / (limit.first * limit.second) << std::endl;
        }
    }
}

Action Env::step(bool random) {
    assert(!terminated_);
    assert(curr_state_->candidate_factors_.size());
    Action act; 
    std::tie(act, expanded_) = curr_state_->select_action(random);
    curr_state_ = curr_state_->take_action(act, constraints_);
    terminated_ = curr_state_->is_terminated();
    return act;
}

Reward Env::calculate_reward(){
    assert(terminated_ && curr_state_->is_terminated());
    auto& reward = curr_state_->remembered_reward;
    if (!(reward > 0)) {
    } 
    else if (curr_state_->is_error_out()) {
        reward = {punish_, 0.0};
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
        if (obj_ == Objective::CYCLE) reward.value = (double)(analyzer_.get_cycle());
        else if (obj_ == Objective::ENERGY) reward.value = (double)(analyzer_.get_energy());
        reward.reward = -std::log10(reward.value);
        // std::cout << "analyzed reward:" << reward.reward << std::endl;
        
        insert(&curr_state_->symbol_table_, reward.value);

        if (reward > best_reward_) {
            TILEFLOW_LOG("Update best "; if (verbose_level) curr_state_->symbol_table_.show_brief(std::cerr); std::cerr 
                    << " value: " <<  reward.value);
            best_reward_ = reward;
            best_symbol_table_ = &curr_state_->symbol_table_;
        }
        if (punish_ + 2 > reward.reward) 
            punish_ = reward.reward - 2;
    }
    // std::cout << "calculated reward:" << reward.reward << std::endl;
    return reward;
}

void Env::insert(const SymbolTable* table, double value){
    auto iter = topKList_.begin();
    while (iter != topKList_.end() && value > iter -> second) iter ++;
    topKList_.insert(iter, {table, value});
    if (topKList_.size() > topk_) topKList_.pop_back();
}

std::ostream& Env::report(std::ostream& o) const {
    o << "=========Top " << topKList_.size() << " ===========" << std::endl;
    for (unsigned i = 0; i < topKList_.size(); i++) {
        o << i << ":\t";
        auto & item = topKList_.at(i);
        item.first->show_brief(o);
        o << ", value: " << item.second << std::endl;
    }
    o << "=======Top " << topKList_.size() << " End==========" << std::endl;
    return o;
}

State* MCTS::select_state() {
    // std::cout << "---beg state selection---" << std::endl;
    env_->reset();
    // std::cout << *env_->get_curr_state(); 
    while (!env_->is_terminated() && !env_->is_expanded()) {
        env_->step(random_);
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

Reward MCTS::rollout(State* state) {
    // Reward ave_reward = 0.0; 
    Reward max_reward = -100; 
    for (unsigned i = 0; i < n_rollout; ++i) {
        env_->reset(state);
        while(!env_->is_terminated()) {
            env_->step(true);
        }
        // ave_reward += (env_->calculate_reward() - ave_reward) / (i+1);
        auto reward = env_->calculate_reward();
        if (reward > max_reward) 
            max_reward = reward;
    }
    // std::cout << "roll out reward: " << max_reward << std::endl;
    return max_reward;
}

void MCTS::back_prop(State* state, Reward reward){
    while (state) {
        state->n_visit ++;
        state->ave_reward += (reward - state->ave_reward) / state->n_visit;
        state = state->parent;
    }
}

void MCTS::search() {
    assert(env_!=nullptr);
    n_unexplored = env_->get_curr_state()->candidate_factors_.size();
    start_timer();
    for (iter_ = 0; iter_ < n_iteration; iter_++ ) {
        auto state = select_state();
        Reward reward = rollout(state);
        back_prop(state, reward); 
        if (!n_unexplored) {
            TILEFLOW_LOG("MCTS: early exit for exhausting the search space after " << iter_ + 1 << " rounds.");
            return;
        }
        if (get_elapsed_time() > timeout_) {
            TILEFLOW_LOG("MCTS exit becuase of timeout.");
            return;
        }
        std::stringstream s;
        state->Print(s);
        logs_.push_back({get_elapsed_time(), iter_, state->ave_reward, env_->get_reward(), s.str()});
    }
    
}

std::pair<Action, bool> State::select_action(bool random) {
    if (random) {
        auto act = candidate_factors_[std::rand() % candidate_factors_.size()];
        return {act, children_[act] == nullptr || children_[act]->n_visit == 0};
    }

    double score = -1e9;
    Action act;

    for (auto factor: candidate_factors_) {
        auto& child = children_[factor];
        if (!child || child->n_visit == 0) {
            return {factor, true};            
        }
        auto cur_score = child->ave_reward + C * std::sqrt(std::log(n_visit) / child->n_visit);
        if (cur_score > score) {
            act = factor;
            score = cur_score.reward;
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
        child->last_action = action;
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

std::ostream& operator<<(std::ostream& o, const Reward& r) {
    o << r.reward << "," << r.value;
    return o;
}


void State::Print(std::ostream& o) {
    auto node = this;
    while (node && node->parent) {
        auto& entry = node->symbol_table_.lookup(node->parent->variable_index);
        assert(entry.fixed_);
        o << entry.name_ << ":" << entry.value_ << "-";
        node = node->parent;
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