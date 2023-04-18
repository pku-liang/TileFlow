#pragma once 

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>

namespace TileFlow {
    struct Entry {
        std::string name_;
        int value_;
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
        virtual int eval(const SymbolTable&) = 0;
        virtual void display(const SymbolTable&) = 0;
    };

    struct SumExpr: public Expr {
        std::vector<std::shared_ptr<Expr> > operands_;
        SumExpr(std::vector<std::shared_ptr<Expr> >& operands):
            operands_(operands){}
        int eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table) override;
    };
    struct ProductExpr: public Expr {
        std::vector<std::shared_ptr<Expr> > operands_;
        ProductExpr(const std::vector<std::shared_ptr<Expr> >& operands):
            operands_(operands){}
        ProductExpr(const std::initializer_list<std::shared_ptr<Expr> >& operands):
            operands_(operands){}
        ProductExpr(const std::vector<int>& operands);
        int eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table) override;
    };
    struct VariableExpr: public Expr {
        int idx_;
        VariableExpr(int idx): idx_(idx){}
        int eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table) override;
    };

    struct ParameterExpr: public Expr {
        int value_;
        ParameterExpr(int value): value_(value){}
        int eval(const SymbolTable& symb_table) override;
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
        int eval(const SymbolTable& symb_table) override;
        void display(const SymbolTable& symb_table) override;
    };
} // namespace TileFlow