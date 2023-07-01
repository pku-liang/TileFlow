#include <assert.h>

#include "tileflow/mapper/op.hpp"


namespace TileFlow {

namespace Op {


std::shared_ptr<ResourceExpr> max(std::vector<std::shared_ptr<ResourceExpr> > exprs){
    assert(exprs.size());
    if (exprs.size() == 1) return exprs.front();
    return std::make_shared<PairMaxExpr>(exprs);
}

std::shared_ptr<ResourceExpr> sum(std::vector<std::shared_ptr<ResourceExpr> > exprs){
    assert(exprs.size());
    if (exprs.size() == 1) return exprs.front();
    return std::make_shared<PairSumExpr>(exprs);
}

std::shared_ptr<ResourceExpr> operator <= (
    const std::shared_ptr<ResourceExpr>& expr,
    const std::shared_ptr<ResourceExpr>& limit 
){
    return std::make_shared<PairCondExpr> (
        expr, 
        limit, 
        PairCondExpr::type_t::LEQ
    );
}

std::shared_ptr<ResourceExpr> pair(int x, int y) {
    return pair(
        parameter(x),parameter(y)
    );
}

std::shared_ptr<ResourceExpr> pair(
    const std::shared_ptr<Expr>& x, 
    const std::shared_ptr<Expr>& y) {
    return std::make_shared<PairExpr>(x,y);
}

std::shared_ptr<Expr> product(std::vector<std::shared_ptr<Expr>>& exprs) {
    if (!exprs.size()) return parameter(1);
    if (exprs.size() == 1) return exprs.front();
    return std::make_shared<ProductExpr>(exprs);
}

std::shared_ptr<Expr> product(std::initializer_list<std::shared_ptr<Expr>> exprs) {
    return std::make_shared<ProductExpr>(exprs);
}

std::shared_ptr<Expr> product(const std::pair<num_t, std::vector<int> >& exprs){
    return std::make_shared<ProductExpr>(exprs);
}

std::shared_ptr<Expr> sum(std::vector<std::shared_ptr<Expr> >& exprs) {
    if (!exprs.size()) return parameter(0);
    if (exprs.size() == 1) return exprs.front();
    return std::make_shared<SumExpr>(exprs);
}

std::shared_ptr<Expr> max(std::vector<std::shared_ptr<Expr> >& exprs) {
    if (!exprs.size()) return parameter(0);
    if (exprs.size() == 1) return exprs.front();
    return std::make_shared<MaxExpr>(exprs);
}

std::shared_ptr<Expr> variable(int x) {
    return std::make_shared<VariableExpr>(x);
}

std::shared_ptr<Expr> operator == (
    const std::shared_ptr<Expr>& left,
    const std::shared_ptr<Expr>& right
){
    return std::make_shared<CondExpr>(
        left, right, CondExpr::EQU
    );
}

std::shared_ptr<Expr> operator <= (
    const std::shared_ptr<Expr>& left,
    const std::shared_ptr<Expr>& right
){
    return std::make_shared<CondExpr>(
        left, right, CondExpr::LEQ
    );
}

std::shared_ptr<Expr> parameter(num_t x) {
    return std::make_shared<ParameterExpr>(x);
}

} // namespace Op 

} // namespace TileFlow 