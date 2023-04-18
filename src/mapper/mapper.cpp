#include <fstream>

#include "tileflow/mapper/mapper.hpp"
#include "tileflow/mapper/checker.hpp"

namespace TileFlow {

SymbolTable Mapper::search() {
    if (global_symbol_table_.get_num_variables() == 0) {
        for (auto& cons: constraints_) {
            TILEFLOW_ASSERT(cons.expr->eval(global_symbol_table_), cons.msg << "("; 
                cons.expr->display(global_symbol_table_); std::cout << ") violated");
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
    analysis.Report();
}

} // namespace TileFlow 