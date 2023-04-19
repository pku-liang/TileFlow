#include <fstream>

#include "tileflow/mapper/mapper.hpp"
#include "tileflow/mapper/checker.hpp"

namespace TileFlow {

SymbolTable Mapper::search() {
    if (global_symbol_table_.get_num_variables() == 0) {
        for (auto& cons: constraints_) {
            if (cons.type_ == Constraint::LOOPCOUNT || 
            (cons.type_ == Constraint::MEM && enable_mem_check_) || 
            (cons.type_ == Constraint::SPATIAL && enable_spatial_check_)) {
                TILEFLOW_ASSERT(cons.expr->eval(global_symbol_table_), cons.msg << "("; 
                cons.expr->display(global_symbol_table_); std::cout << ") violated");
            }
        }
        analysis::TileFlow::NestAnalysis analysis(workloads_, mapping_, arch_specs_, topology_);
        analysis.analyze();
        return SymbolTable();
    }
    TILEFLOW_ERROR("NOT Implemented Error");
}

void Mapper::dump(const std::string & filename){
    analysis::TileFlow::NestAnalysis analysis(workloads_, mapping_, arch_specs_, topology_, &optimum_);
    analysis.analyze();
    std::ofstream file;
    size_t tmp = filename.find_last_of('.');
    if (tmp == std::string::npos || filename.substr(tmp) != ".csv")
        file.open(filename + ".csv");
    else file.open(filename);
    file << "metric,value" << std::endl;
    file << "Cycle," << analysis.get_cycle() << std::endl;
    file << "Energy," << analysis.get_energy() << std::endl;
    if (verbose_level)
        std::cout << "[TileFlow]: result written into " << filename << std::endl;
    file.close();
}

void Mapper::report() {
    analysis::TileFlow::NestAnalysis analysis(workloads_, mapping_, arch_specs_, topology_, &optimum_);
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

} // namespace TileFlow 