#pragma once 

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <functional>
#include <iostream>

namespace TileFlow {
    const int ERROR_OUT=1;
    const int NORMAL_OUT=2;

    typedef size_t num_t;

    struct Constraint;

    struct Entry {
        std::string name_;
        num_t value_;
        int idx_;
        std::set<num_t> candidates_;
        bool fixed_ = false; 
    };

    std::ostream& operator<< (std::ostream& o, const Entry&);

    class  SymbolTable {
        std::unordered_map<std::string, int> name2idx_;  
        std::unordered_map<int, Entry> idx2values_;
        bool failed_ = false; 
        bool fail_check(const std::vector<Constraint>& constraints_);

    public:
        int idx = 0;
        const Entry& lookup(const std::string& key) const {return idx2values_.at(name2idx_.at(key));}
        const Entry& lookup(int key) const {return idx2values_.at(key);}
        int count(int key) const {return idx2values_.count(key);}
        bool is_terminated() const {
            for(auto& kv: idx2values_) 
                if(!kv.second.fixed_) return false;
            return true;
        }
        int get_num_variables() const {return -idx;}
        int count_unfixed() const {int ret = idx2values_.size(); for (auto& kv: idx2values_) ret -= kv.second.fixed_; return ret;}
        int get_next_var() const;
        void show_brief(std::ostream& o) const;
        
        
        Entry& operator[](int key) {return idx2values_[key];}
        int insert(const std::string name = "") {
            std::string name_ = name;
            if (name_ == "?" || name_ == "X"){
                for (int i = 0; name2idx_.count(name_); i++){
                    name_ = name + std::to_string(i);
                }
            }
            if (name2idx_.count(name_) == 0) {
                idx--;
                name2idx_[name_] = idx;
                idx2values_[idx] = {name_, 0, idx, {}, false};
            }
            return name2idx_.at(name_);
        }

        void init(const std::vector<Constraint>& constraints_);
        void fix_and_update(int index, num_t value, const std::vector<Constraint>& constraints_);
    };

    std::ostream& operator<< (std::ostream& o, const SymbolTable&);

    extern SymbolTable global_symbol_table_;
    
    struct PairExpr;
    struct PairSumExpr;
    struct PairMaxExpr;
    struct PairCondExpr;
    struct ProductExpr;
    struct VariableExpr;
    struct ParameterExpr;
    struct CondExpr;
    struct SumExpr;
    struct MaxExpr;

    struct ExprVisitor {
        virtual void visitPairExpr(const PairExpr*);
        virtual void visitPairSumExpr(const PairSumExpr*);
        virtual void visitPairMaxExpr(const PairMaxExpr*);
        virtual void visitPairCondExpr(const PairCondExpr*);
        virtual void visitProductExpr(const ProductExpr*);
        virtual void visitSumExpr(const SumExpr*);
        virtual void visitMaxExpr(const MaxExpr*);
        virtual void visitVariableExpr(const VariableExpr*);
        virtual void visitParameterExpr(const ParameterExpr*);
        virtual void visitCondExpr(const CondExpr*);
    };

    struct Expr {
        virtual num_t eval(const SymbolTable&) = 0;
        virtual void display(const SymbolTable&, std::ostream& = std::cout) = 0;
        virtual void accept(ExprVisitor*) const = 0;
    };

    struct ResourceExpr: public Expr {
        virtual num_t eval(const SymbolTable& ) override {return 0;}
        virtual void display(const SymbolTable&, std::ostream&) override {}
        // given y's limit, compute minimum required x
        virtual std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) = 0;
        virtual void accept(ExprVisitor* visitor) const = 0;
    };

    struct PairExpr: public ResourceExpr {
        std::shared_ptr<Expr> x_;
        std::shared_ptr<Expr> y_;
        PairExpr(const std::shared_ptr<Expr>& x, 
            const std::shared_ptr<Expr>& y): x_(x), y_(y) {}
        void display(const SymbolTable& symb_table, std::ostream&) override;
        std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitPairExpr(this);}
    };

    struct PairSumExpr: public ResourceExpr {
        std::vector<std::shared_ptr<ResourceExpr> > operands_;
        PairSumExpr(std::vector<std::shared_ptr<ResourceExpr> >& operands):
            operands_(operands){}
        void display(const SymbolTable& symb_table, std::ostream&) override;
        std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitPairSumExpr(this);}
    };

    struct PairMaxExpr: public ResourceExpr {
        std::vector<std::shared_ptr<ResourceExpr> > operands_;
        PairMaxExpr(std::vector<std::shared_ptr<ResourceExpr> >& operands):
            operands_(operands){}
        void display(const SymbolTable& symb_table, std::ostream&) override;
        std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitPairMaxExpr(this);}
    };

    struct PairCondExpr: public ResourceExpr {
        std::shared_ptr<ResourceExpr> expr_;
        std::shared_ptr<ResourceExpr> limit_;
        enum type_t {
            LEQ
        }op_;
        PairCondExpr(
            std::shared_ptr<ResourceExpr> expr,
            std::shared_ptr<ResourceExpr> limit, 
            type_t op): expr_(expr), limit_(limit), op_(op) {}
        void display(const SymbolTable& symb_table, std::ostream&) override;
        num_t eval(const SymbolTable& symb_table) override;
        std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitPairCondExpr(this);}
    };

    struct SumExpr: public Expr {
        std::vector<std::shared_ptr<Expr> > operands_;
        SumExpr(std::vector<std::shared_ptr<Expr> >& operands):
            operands_(operands){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table, std::ostream&) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitSumExpr(this);}
    };
    
    struct MaxExpr: public Expr {
        std::vector<std::shared_ptr<Expr> > operands_;
        MaxExpr(std::vector<std::shared_ptr<Expr> >& operands):
            operands_(operands){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table, std::ostream&) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitMaxExpr(this);}
    };

    struct ProductExpr: public Expr {
        std::vector<std::shared_ptr<Expr> > operands_;
        ProductExpr(const std::vector<std::shared_ptr<Expr> >& operands):
            operands_(operands){}
        ProductExpr(const std::initializer_list<std::shared_ptr<Expr> >& operands):
            operands_(operands){}
        ProductExpr(const std::vector<int>& operands);
        ProductExpr(const std::pair<num_t, std::vector<int> >& operands);
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table, std::ostream&) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitProductExpr(this);}
    };
    struct VariableExpr: public Expr {
        int idx_;
        VariableExpr(int idx): idx_(idx){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table, std::ostream&) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitVariableExpr(this);}
    };

    struct ParameterExpr: public Expr {
        num_t value_;
        ParameterExpr(num_t value): value_(value){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table, std::ostream&) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitParameterExpr(this);}
    };
    
    struct CondExpr: public Expr {
        enum type_t {
            EQU, 
            LEQ,
            LESS
        }op_;
        std::shared_ptr<Expr> left_;
        std::shared_ptr<Expr> right_;
        CondExpr(std::shared_ptr<Expr> left, 
                std::shared_ptr<Expr> right, 
                CondExpr::type_t op):
        op_(op), left_(left), right_(right){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table, std::ostream&) override;
        void accept(ExprVisitor* visitor) const override {visitor->visitCondExpr(this);}
    };

    struct VariableCollector: public ExprVisitor {
        std::set<int> variables_;
        std::function<bool(int)> _func;

        void visitVariableExpr(const VariableExpr* expr) override {
            if (_func(expr->idx_))
                variables_.insert(expr->idx_);
        }

        std::set<int> operator()(const Expr* root, std::function<bool(int)>func) {
            _func = func;
            root->accept(this);
            return std::move(variables_);
        }
        std::set<int> operator()(const Expr* root) {
            _func = [](int){return true;};
            root->accept(this);
            return std::move(variables_);
        }
    };
/*
    struct Simplifier: public ExprVisitor {
        bool isSum;
        std::vector<std::shared_ptr<Expr> > operands_;
        void visitSumExpr(const SumExpr* expr) override {
            std::vector<std::shared_ptr<Expr> > operands; 
            for (auto operand: expr->operands_) {
                isSum = false;
                operands_.clear();
                operand->accept(this);
                if (isSum) {
                    operands.insert(operands.end(), operands_.begin(), operands_.end());
                }
                else operands.push_back(operand);
            }
            const_cast<SumExpr*>(expr)->operands_ = operands;
            operands_ = operands;
            isSum = true;
        }
        void visitPairExpr(const PairExpr*) {isSum=false; }
        void visitPairSumExpr(const PairSumExpr*) {isSum=false; }
        void visitPairMaxExpr(const PairMaxExpr*) {isSum=false; }
        void visitPairCondExpr(const PairCondExpr*) {isSum=false; }
        void visitProductExpr(const ProductExpr*) {isSum=false; }
        void visitSumExpr(const SumExpr*) {isSum=false; }
        void visitMaxExpr(const MaxExpr*) {isSum=false; }
        void visitVariableExpr(const VariableExpr*) {isSum=false; }
        void visitParameterExpr(const ParameterExpr*) {isSum=false; }
        void visitCondExpr(const CondExpr*) {isSum=false; }
    };
*/
    struct Constraint {
        enum {
            LOOPCOUNT,
            MEM, 
            SPATIAL
        }type_; 
        std::shared_ptr<Expr> expr;
        std::string msg;
        std::string short_msg = "";
    };

} // namespace TileFlow