#include "tileflow/mapper/expr.hpp"



namespace TileFlow {

namespace Op {

std::shared_ptr<ResourceExpr> max(std::vector<std::shared_ptr<ResourceExpr> > exprs);

std::shared_ptr<ResourceExpr> sum(std::vector<std::shared_ptr<ResourceExpr> > exprs);
 
std::shared_ptr<ResourceExpr> operator <= (const std::shared_ptr<ResourceExpr>&, const std::shared_ptr<ResourceExpr>&);

std::shared_ptr<ResourceExpr> pair(int x, int y);

std::shared_ptr<ResourceExpr> pair(
    const std::shared_ptr<Expr>& x, 
    const std::shared_ptr<Expr>& y);

std::shared_ptr<Expr> product(std::vector<std::shared_ptr<Expr> >& exprs);

std::shared_ptr<Expr> product(std::initializer_list<std::shared_ptr<Expr>> exprs);

std::shared_ptr<Expr> product(const std::pair<num_t, std::vector<int> >& exprs);

std::shared_ptr<Expr> sum(std::vector<std::shared_ptr<Expr> >& exprs);

std::shared_ptr<Expr> max(std::vector<std::shared_ptr<Expr> >& exprs);

std::shared_ptr<Expr> variable(int x);

std::shared_ptr<Expr> operator == (
    const std::shared_ptr<Expr>& left,
    const std::shared_ptr<Expr>& right
);

std::shared_ptr<Expr> operator <= (
    const std::shared_ptr<Expr>& left,
    const std::shared_ptr<Expr>& right
);

std::shared_ptr<Expr> parameter(num_t x);

} // namespace Op 

} // namespace TileFlow 