#include "tileflow/loop-analysis/memory-state.hpp"


namespace analysis {

namespace TileFlow {

    MemoryState& MemoryState::Union(const MemoryState& other){
        for (auto& kv: other.data_spaces_) {
            if (data_spaces_.count(kv.first) == 0) {
                data_spaces_.emplace(kv.first, kv.second);
            }
            else {
                data_spaces_.at(kv.first) += kv.second;
            }
        }
        return *this;
    }

    MemoryState& MemoryState::Substract(const MemoryState& other){
        for (auto& kv: data_spaces_) { 
            if (other.data_spaces_.count(kv.first))
                kv.second = kv.second - other.data_spaces_.at(kv.first);
        }
        return *this;
    }

    MemoryState& MemoryState::Add(const MemoryState& other){
        for (auto& kv: other.data_spaces_) {
            if (data_spaces_.count(kv.first) == 0) {
                data_spaces_.emplace(kv.first, kv.second);
            }
            else {
                data_spaces_.at(kv.first) += kv.second;
            }
        }
        return *this;
    }

    MemoryState& MemoryState::Intersect(const MemoryState& ){
        // TODO: realize real intersection logic here.  
        return *this;
    }

    void MemoryState::insert(std::uint64_t spatial_id, 
            const problem::OperationPoint& low_point,
            const problem::OperationPoint& high_point) {
        data_spaces_.emplace(spatial_id, problem::OperationSpace(workload_, low_point, high_point));
    }
    
    void MemoryState::insert(std::uint64_t spatial_id, 
        const problem::OperationSpace& data_space) {
        data_spaces_.emplace(spatial_id, data_space);
    }

    const problem::Workload* MemoryState::workload_;

    void MemoryState::show() const{
        for (auto kv: data_spaces_) {
            std::cout << kv.first << ":";
            kv.second.Print(std::cout);
            std::cout << std::endl;
        }
        std::cout << std::endl;
    } 

    MemoryState MemoryState::operator - (const MemoryState& other) {
        MemoryState ret;
        for (auto& kv: data_spaces_) {
            if (other.getDataSpaces().count(kv.first)) {
                ret.data_spaces_.emplace(kv.first, kv.second - other.getDataSpaces().at(kv.first));
            }
            else {
                ret.data_spaces_.emplace(kv.first, kv.second);
            }
        }
        return ret;
    }

}   // namespace TileFlow

} // namespace analysis