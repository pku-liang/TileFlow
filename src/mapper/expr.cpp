#include <iostream>
#include <algorithm>
#include <math.h>

#include "tileflow/mapper/expr.hpp"
#include "tileflow/common.hpp"
#include "tileflow/mapper/op.hpp"


using namespace TileFlow::Op;

namespace TileFlow {

SymbolTable global_symbol_table_;

num_t SumExpr::eval(const SymbolTable& symb_table) {
    num_t ret = 0;
    for (auto& op: operands_)
        ret += op->eval(symb_table);
    return ret; 
}

num_t MaxExpr::eval(const SymbolTable& symb_table) {
    num_t ret = 0;
    for (auto& op: operands_)
        ret = std::max(ret, op->eval(symb_table));
    return ret; 
}

void SumExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    bool is_first = true;
    if (operands_.size() > 1) o << "(";
    for (auto& op: operands_) {
        if (!is_first)
            o << "+";
        else is_first = false;
        op->display(symb_table, o);
    }
    if (is_first) 
        o << 1;
    if (operands_.size() > 1) o << ")";
}

void MaxExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    bool is_first = true;
    if (operands_.size() > 1) o << "Max(";
    for (auto& op: operands_) {
        if (!is_first)
            o << ",";
        else is_first = false;
        op->display(symb_table, o);
    }
    if (is_first) 
        o << 1;
    if (operands_.size() > 1) o << ")";
}

ProductExpr::ProductExpr(const std::vector<int>& operands) {
    for (auto id: operands) {
        operands_.push_back(variable(id));
    }
}

ProductExpr::ProductExpr(const std::pair<num_t, std::vector<int> >& operands){
    if (operands.first != 1) {
        operands_.push_back(parameter(operands.first));
    }
    for (auto id: operands.second)
        operands_.push_back(std::make_shared<VariableExpr>(id));
}

num_t ProductExpr::eval(const SymbolTable& symb_table) {
    num_t ret = 1;
    for (auto& op: operands_)
        ret *= op->eval(symb_table);
    return ret; 
}

void ProductExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    bool is_first = true;
    for (auto& op: operands_) {
        if (!is_first)
            o << "*";
        else is_first = false;
        op->display(symb_table, o);
    }
    if (is_first) 
        o << 1;
}

num_t VariableExpr::eval(const SymbolTable& symb_table) {
    auto & entry = symb_table.lookup(idx_);
    if (entry.fixed_) return entry.value_;
    return 1;
}

void VariableExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    auto & entry = symb_table.lookup(idx_);
    o << entry.name_;
    if (entry.fixed_) o << "(" << entry.value_ << ")";
}

num_t ParameterExpr::eval(const SymbolTable& ) {
    return value_;
}

void ParameterExpr::display(const SymbolTable& , std::ostream& o) {
    o << value_;
}

num_t CondExpr::eval(const SymbolTable& symb_table) {
    if (op_ == CondExpr::EQU)
        return left_->eval(symb_table) == right_->eval(symb_table);
    else if (op_ == CondExpr::LESS) 
        return left_->eval(symb_table) < right_->eval(symb_table);
    else if (op_ == CondExpr::LEQ) 
        return left_->eval(symb_table) <= right_->eval(symb_table);
    return 0;
}

void CondExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    left_->display(symb_table, o);
    if (op_ == CondExpr::EQU) o << "==";
    else if (op_ == CondExpr::LESS) o << "<";
    else if (op_ == CondExpr::LEQ) o << "<=";
    right_->display(symb_table, o);
}

void PairExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    o << "<";
    x_->display(symb_table, o);
    o << ", ";
    y_->display(symb_table, o);
    o << ">";
}

std::pair<int, int> PairExpr::eval_pair(const SymbolTable& symb_table, int) {
    return {x_->eval(symb_table), y_->eval(symb_table)};
}

void PairSumExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    bool isFirst = true;
    if (operands_.size() > 1) o << "(";
    for (auto& operand: operands_) {
        if (!isFirst) o << "+";
        else isFirst = false;
        operand->display(symb_table, o);
    }
    if (operands_.size() > 1) o << ")";
}

std::pair<int, int> PairSumExpr::eval_pair(const SymbolTable& symb_table, int limit_y) {
    std::vector<std::pair<int, int> > pairs;
    for (auto & operand: operands_) {
        pairs.push_back(operand->eval_pair(symb_table, limit_y));
    }
    std::sort(pairs.begin(), pairs.end(), [](const std::pair<int, int>& p1, const std::pair<int, int>&p2){
        return p1.first > p2.first;
    });
    int max_x = 0, start_y = limit_y;
    int max_y = 0;
    for (auto& pair: pairs) {
        if (start_y + pair.second <= limit_y) {
            start_y += pair.second;
        }
        else {
            max_x += pair.first;
            start_y = pair.second;
        }
        max_y = std::max(max_y, start_y);
    }
    return {max_x, max_y};
}


void PairMaxExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    bool isFirst = true;
    if (operands_.size() > 1) o << "Max(";
    for (auto& operand: operands_) {
        if (!isFirst) o << ",";
        else isFirst = false;
        operand->display(symb_table, o);
    }
    if (operands_.size() > 1) o << ")";
}

std::pair<int, int> PairMaxExpr::eval_pair(const SymbolTable& symb_table, int limit_y) {
    int max_x = 0, max_y = 0;
    for (auto& operand: operands_){
        auto pair = operand->eval_pair(symb_table, limit_y);
        max_x = std::max(max_x, pair.first);
        max_y = std::max(max_y, pair.second);
    }
    return {max_x, max_y};
}

void PairCondExpr::display(const SymbolTable& symb_table, std::ostream& o) {
    expr_->display(symb_table, o);
    if (op_ == PairCondExpr::LEQ) o << "<=";
    limit_->display(symb_table, o);
}

num_t PairCondExpr::eval(const SymbolTable& symb_table) {
    auto limit = limit_->eval_pair(symb_table, 0);
    auto usage = expr_->eval_pair(symb_table, limit.second);
    return limit.first >= usage.first && limit.second >= usage.second;
}

std::pair<int, int> PairCondExpr::eval_pair(const SymbolTable& , int ){
    TILEFLOW_ERROR("NOT IMPLEMENTED ERROR");
    return {0,0};
}

std::set<num_t> get_candidates(num_t num) {
    std::set<num_t> ret;
    int square_root = std::sqrt(num);
    for (int i = 1; i <= square_root; i++) {
        if (num % i == 0) {
            ret.insert(i);
            ret.insert(num / i);
        }
    }
    return ret;
}

void intersect(std::set<num_t>& s, const std::set<num_t>& t) {
    auto is = s.begin(), it = t.begin();
    while(is != s.end() && it != t.end()) {
        while (it != t.end() && *it < *is) {it++;}
        while (is != s.end() && *is < *it) {is = s.erase(is);}
        if (it!= t.end() && is != s.end() && *it == *is) {
            it++; is++;
        }
    }
    while(is != s.end()) is = s.erase(is);
}

num_t gcd(num_t a, num_t b)
{
    if (a == 0)
        return b;
    return gcd(b % a, a);
}

bool SymbolTable::fail_check(const std::vector<Constraint>& constraints_) {
    for (auto& cons: constraints_) {
        if (cons.type_ == Constraint::MEM || cons.type_ == Constraint::SPATIAL) {
            if (!cons.expr->eval(*this)){
                if (verbose_level > 1)
                    TILEFLOW_LOG("Search end because of cond check failure: " << cons.msg);
                return true;
            } 
        }
        else if (cons.type_ == Constraint::LOOPCOUNT) {
            auto cond = std::static_pointer_cast<CondExpr>(cons.expr);
            auto vars = VariableCollector()(cons.expr.get(), [this](int idx){return !lookup(idx).fixed_;});
            auto r = cond->right_->eval(*this);
            auto l = cond->left_->eval(*this);
            if ((vars.empty() && l != r) || (!vars.empty() && r % l != 0)) {
                if (verbose_level > 1)
                    TILEFLOW_LOG("Search end because of cond check failure: " << cons.msg);
                return true;
            }
        }
    }
    return false;
}

void SymbolTable::fix_and_update(int index, num_t value,
        const std::vector<Constraint>& constraints_) {
    assert(value!=0);
    assert(idx2values_.count(index));
    auto& entry = idx2values_[index];
    assert(entry.candidates_.count(value));
    entry.value_ = value;
    entry.fixed_ = true;

    std::vector<std::shared_ptr<CondExpr> > lc_cons;
    for (auto& cons: constraints_) {
        if (cons.type_ == Constraint::LOOPCOUNT) 
            lc_cons.push_back(std::static_pointer_cast<CondExpr>(cons.expr));
    }
    std::unordered_map<int, num_t> fixed;
    for (auto cons: lc_cons) {
        std::set<int> vars = VariableCollector()(cons.get(), [this, index](int idx)->bool
            {return !lookup(idx).fixed_ || (idx == index);});
        if(vars.count(index)) {
            num_t lc = cons->right_->eval(*this) / cons->left_->eval(*this);
            std::set<num_t> candidates = get_candidates(lc);
            int a_index;
            for (auto& idx: vars) {
                if (idx == index) continue;
                a_index = idx;
                auto& entry = idx2values_[idx];
                intersect(entry.candidates_, candidates);
                if (entry.candidates_.size() == 1) {
                    assert(*entry.candidates_.begin() == 1);
                    entry.fixed_ = true;
                    entry.value_ = 1;
                }
            }
            if (vars.size() == 2) {
                fixed[a_index] = lc;
            }
        }
    }
    for(auto& kv: fixed){
        auto& entry = idx2values_[kv.first];
        if ((entry.fixed_ == false && entry.candidates_.count(kv.second)) 
        || (entry.fixed_ && entry.value_ == kv.second)) {
            entry.fixed_ = true; 
            entry.value_ = kv.second;
            entry.candidates_ = {kv.second};
        }
        else {
            entry.fixed_ = false;
            entry.value_ = -1;
            entry.candidates_ = {};
        }
    }
    failed_ = fail_check(constraints_);
}

void SymbolTable::init(const std::vector<Constraint>& constraints_) {
    std::vector<std::shared_ptr<CondExpr> > lc_cons; 
    for (auto& cons: constraints_) {
        if (cons.type_ == Constraint::LOOPCOUNT) 
            lc_cons.push_back(
                std::static_pointer_cast<CondExpr>(cons.expr));
    }
    std::unordered_map<int, std::vector<int> > candidates;
    std::unordered_map<int, num_t> fixed;
    for (auto cons: lc_cons) {
        std::set<int> vars = VariableCollector()(cons.get());
        num_t lc = cons->right_->eval(*this) / cons->left_->eval(*this);
        for (auto& idx: vars) {
            candidates[idx].push_back(lc);
        }
        if (vars.size() == 1) {
            auto v = *vars.begin();
            if (fixed.count(v) && lc != fixed[v])
                fixed[v] = 0; 
            else fixed[v] = lc;
        }
    }
    
    for (auto& kv: candidates) {
        auto& entry = idx2values_[kv.first];
        num_t g = kv.second.front();
        for (int i = 1; i < (int)kv.second.size(); i++) 
            g = gcd(g, kv.second[i]);
        auto factors = get_candidates(g);
        entry.candidates_.insert(factors.begin(), factors.end());
        if (g == 1) {
            assert(*entry.candidates_.begin() == 1);
            entry.fixed_ = true;
            entry.value_ = 1;
        }
    }

    for(auto& kv: fixed){
        auto& entry = idx2values_[kv.first];
        if ((!entry.fixed_ && entry.candidates_.count(kv.second)) 
        || (entry.fixed_ && (entry.value_ == kv.second))) {
            entry.fixed_ = true; 
            entry.value_ = kv.second;
            entry.candidates_ = {kv.second};
        }
        else {
            entry.fixed_ = false;
            entry.value_ = -1;
            entry.candidates_ = {};
        }
    }
    failed_ = fail_check(constraints_);
}

int SymbolTable::get_next_var() const {
    if (failed_) return ERROR_OUT;
    size_t min_candidate = 1e3;
    int var = NORMAL_OUT;
    for (auto& kv: idx2values_) {
        if (!kv.second.fixed_ && (min_candidate > kv.second.candidates_.size())) {
            min_candidate = kv.second.candidates_.size();
            var = kv.first;
            if (min_candidate == 0) {
                TILEFLOW_LOG("ERROR OUT at " << kv.second.name_; std::cerr);
                return ERROR_OUT;
            }
        }
    }
    return var;
}


void ExprVisitor::visitPairExpr(const PairExpr*expr){
    expr->x_->accept(this);
    expr->y_->accept(this);
}
void ExprVisitor::visitPairSumExpr(const PairSumExpr*expr){
    for (auto& operand: expr->operands_) operand->accept(this);
}
void ExprVisitor::visitPairMaxExpr(const PairMaxExpr*expr){
    for (auto& operand: expr->operands_) operand->accept(this);
}
void ExprVisitor::visitPairCondExpr(const PairCondExpr*expr){
    expr->expr_->accept(this); expr->limit_->accept(this);
}
void ExprVisitor::visitProductExpr(const ProductExpr*expr){
    for (auto& operand: expr->operands_) operand->accept(this);
}
void ExprVisitor::visitVariableExpr(const VariableExpr*){}

void ExprVisitor::visitParameterExpr(const ParameterExpr*){}

void ExprVisitor::visitCondExpr(const CondExpr*expr){
    expr->left_->accept(this); expr->right_->accept(this);
}
void ExprVisitor::visitSumExpr(const SumExpr* expr) {
    for (auto& operand: expr->operands_) operand->accept(this);
}
void ExprVisitor::visitMaxExpr(const MaxExpr* expr) {
    for (auto& operand: expr->operands_) operand->accept(this);
}

std::ostream& operator<< (std::ostream& o, const Entry& e) {
    o << "(" << e.name_ << "," << e.fixed_ << "," << e.value_;
    if (e.candidates_.size()) {
        o << "{";
        for (auto x: e.candidates_) 
            o << x << ",";
        o << "}";
    }
    o << ")";
    return o;
}

std::ostream& operator<< (std::ostream& o, const SymbolTable& table) {
    o << "=========Symbol Table=========" << std::endl;
    for (int i = -1; i >= table.idx; i--) {
        o << table.lookup(i) << std::endl;
    }
    o << "======End Symbol Table========" << std::endl;
    return o;
}

void SymbolTable::show_brief(std::ostream& o) const {
    for (int i = -1; i >= idx; i--) {
        auto& entry = idx2values_.at(i);
        o << "<" << entry.name_ << "," << entry.fixed_ << "," << entry.value_ << ">,";
    }
}

} // namespaec TileFlow 