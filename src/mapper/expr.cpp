#include <iostream>
#include <algorithm>
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

void SumExpr::display(const SymbolTable& symb_table) {
    bool is_first = true;
    if (operands_.size() > 1) std::cout << "(";
    for (auto& op: operands_) {
        if (!is_first)
            std::cout << "+";
        else is_first = false;
        op->display(symb_table);
    }
    if (is_first) 
        std::cout << 1;
    if (operands_.size() > 1) std::cout << ")";
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

void ProductExpr::display(const SymbolTable& symb_table) {
    bool is_first = true;
    for (auto& op: operands_) {
        if (!is_first)
            std::cout << "*";
        else is_first = false;
        op->display(symb_table);
    }
    if (is_first) 
        std::cout << 1;
}

num_t VariableExpr::eval(const SymbolTable& symb_table) {
    return symb_table.lookup(idx_).value_;
}

void VariableExpr::display(const SymbolTable& symb_table) {
    std::cout << symb_table.lookup(idx_).name_;
}

num_t ParameterExpr::eval(const SymbolTable& ) {
    return value_;
}

void ParameterExpr::display(const SymbolTable& ) {
    std::cout << value_;
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

void CondExpr::display(const SymbolTable& symb_table) {
    left_->display(symb_table);
    if (op_ == CondExpr::EQU) std::cout << "==";
    else if (op_ == CondExpr::LESS) std::cout << "<";
    else if (op_ == CondExpr::LEQ) std::cout << "<=";
    right_->display(symb_table);
}

void PairExpr::display(const SymbolTable& symb_table) {
    std::cout << "<";
    x_->display(symb_table);
    std::cout << ", ";
    y_->display(symb_table);
    std::cout << ">";
}

std::pair<int, int> PairExpr::eval_pair(const SymbolTable& symb_table, int) {
    return {x_->eval(symb_table), y_->eval(symb_table)};
}

void PairSumExpr::display(const SymbolTable& symb_table) {
    bool isFirst = true;
    if (operands_.size() > 1) std::cout << "(";
    for (auto& operand: operands_) {
        if (!isFirst) std::cout << "+";
        else isFirst = false;
        operand->display(symb_table);
    }
    if (operands_.size() > 1) std::cout << ")";
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


void PairMaxExpr::display(const SymbolTable& symb_table) {
    bool isFirst = true;
    if (operands_.size() > 1) std::cout << "Max(";
    for (auto& operand: operands_) {
        if (!isFirst) std::cout << ",";
        else isFirst = false;
        operand->display(symb_table);
    }
    if (operands_.size() > 1) std::cout << ")";
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

void PairCondExpr::display(const SymbolTable& symb_table) {
    expr_->display(symb_table);
    if (op_ == PairCondExpr::LEQ) std::cout << "<=";
    limit_->display(symb_table);
}

num_t PairCondExpr::eval(const SymbolTable& symb_table) {
    auto limit = limit_->eval_pair(symb_table, 0);
    auto usage = expr_->eval_pair(symb_table, limit.second);
    std::cout << "limit: " << limit.first << ", " << limit.second << ";";
    std::cout << "usage: " << usage.first << ", " << usage.second << std::endl;
    return limit.first >= usage.first && limit.second >= usage.second;
}

std::pair<int, int> PairCondExpr::eval_pair(const SymbolTable& , int ){
    TILEFLOW_ERROR("NOT IMPLEMENTED ERROR");
    return {0,0};
}

} // namespaec TileFlow 