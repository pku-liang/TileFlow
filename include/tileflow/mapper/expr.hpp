#pragma once 

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>

namespace TileFlow {
    typedef size_t num_t;

    struct Entry {
        std::string name_;
        num_t value_;
        int idx_;
    };

    class  SymbolTable {
        std::unordered_map<std::string, int> name2idx_;  
        std::unordered_map<int, Entry> idx2values_;
        int idx = 0;
    public:
        Entry lookup(const std::string& key) const {return idx2values_.at(name2idx_.at(key));}
        Entry lookup(int key) const {return idx2values_.at(key);}
        int insert(const std::string name = "") {
            std::string name_ = name;
            for (int i = 0; name2idx_.count(name_); i++){
                name_ = name + std::to_string(i);
            }
            idx--;
            name2idx_[name_] = idx;
            idx2values_[idx] = {name_, 0, idx};
            return idx;
        }
        int get_num_variables() const {return -idx;}
        
    };

    extern SymbolTable global_symbol_table_;

    struct Expr {
        virtual num_t eval(const SymbolTable&) = 0;
        virtual void display(const SymbolTable&) = 0;
    };

    struct ResourceExpr: public Expr {
        virtual num_t eval(const SymbolTable& ) override {return 0;}
        virtual void display(const SymbolTable& ) override {}
        // given y's limit, compute minimum required x
        virtual std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) = 0;
    };

    struct PairExpr: public ResourceExpr {
        std::shared_ptr<Expr> x_;
        std::shared_ptr<Expr> y_;
        PairExpr(const std::shared_ptr<Expr>& x, 
            const std::shared_ptr<Expr>& y): x_(x), y_(y) {}
        void display(const SymbolTable& symb_table) override;
        std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) override;
    };

    struct PairSumExpr: public ResourceExpr {
        std::vector<std::shared_ptr<ResourceExpr> > operands_;
        PairSumExpr(std::vector<std::shared_ptr<ResourceExpr> >& operands):
            operands_(operands){}
        void display(const SymbolTable& symb_table) override;
        std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) override;
    };

    struct PairMaxExpr: public ResourceExpr {
        std::vector<std::shared_ptr<ResourceExpr> > operands_;
        PairMaxExpr(std::vector<std::shared_ptr<ResourceExpr> >& operands):
            operands_(operands){}
        void display(const SymbolTable& symb_table) override;
        std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) override;
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
        void display(const SymbolTable& symb_table) override;
        num_t eval(const SymbolTable& symb_table) override;
        std::pair<int, int> eval_pair(const SymbolTable& symb_table, int limit_y) override;
    };

    struct SumExpr: public Expr {
        std::vector<std::shared_ptr<Expr> > operands_;
        SumExpr(std::vector<std::shared_ptr<Expr> >& operands):
            operands_(operands){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table) override;
    };

    template <typename T> 
    struct MaxExpr: public Expr {
        std::vector<std::shared_ptr<T> > operands_;
        MaxExpr(const std::vector<std::shared_ptr<T> >& operands):
            operands_(operands){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table) override;
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
        void display(const SymbolTable& symb_table) override;
    };
    struct VariableExpr: public Expr {
        int idx_;
        VariableExpr(int idx): idx_(idx){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table) override;
    };

    struct ParameterExpr: public Expr {
        num_t value_;
        ParameterExpr(num_t value): value_(value){}
        num_t eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table) override;
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
        void display(const SymbolTable& symb_table) override;
    };

} // namespace TileFlow