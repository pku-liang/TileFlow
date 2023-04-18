#include <iostream>
#include "tileflow/mapper/expr.hpp"

namespace TileFlow {

SymbolTable global_symbol_table_;

int SumExpr::eval(const SymbolTable& symb_table) {
    int ret = 0;
    for (auto& op: operands_)
        ret += op->eval(symb_table);
    return ret; 
}

void SumExpr::display(const SymbolTable& symb_table) {
    bool is_first = true;
    for (auto& op: operands_) {
        if (!is_first)
            std::cout << "+";
        else is_first = false;
        op->display(symb_table);
    }
    if (is_first) 
        std::cout << 1;
}

ProductExpr::ProductExpr(const std::vector<int>& operands) {
    for (auto id: operands) {
        operands_.push_back(std::make_shared<VariableExpr>(id));
    }
}

int ProductExpr::eval(const SymbolTable& symb_table) {
    int ret = 1;
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

int VariableExpr::eval(const SymbolTable& symb_table) {
    return symb_table.lookup(idx_).value_;
}

void VariableExpr::display(const SymbolTable& symb_table) {
    std::cout << symb_table.lookup(idx_).name_;
}

int ParameterExpr::eval(const SymbolTable& ) {
    return value_;
}

void ParameterExpr::display(const SymbolTable& ) {
    std::cout << value_;
}

int CondExpr::eval(const SymbolTable& symb_table) {
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

} // namespaec TileFlow 