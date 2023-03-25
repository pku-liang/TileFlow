#pragma once 

#include "loop-analysis/loop-state.hpp"

namespace analysis {

namespace TileFlow {


    /**
     * \brief the dataspaces
    */
    struct MemoryState {
        static const problem::Workload* workload_;
        std::unordered_map<std::uint64_t, problem::OperationSpace> data_spaces_;
        
        MemoryState& Union(const MemoryState& other);
        MemoryState& Substract(const MemoryState& other);
        MemoryState& Add(const MemoryState& other);
        MemoryState& Intersect(const MemoryState& other);
        
        // add
        MemoryState& operator += (const MemoryState& other);
        // union
        MemoryState& operator |= (const MemoryState& other);
        // intersect 
        MemoryState& operator &= (const MemoryState& other);
        // substract
        MemoryState& operator -= (const MemoryState& other);

        MemoryState operator - (const MemoryState& other);

        problem::OperationSpace& operator[] (int id) {
            if (data_spaces_.count(id) == 0) {
                data_spaces_.emplace(id, workload_);
            }
            return data_spaces_.at(id);
        }

        inline const problem::OperationSpace& at(std::uint64_t idx) const {
            return data_spaces_.at(idx);
        }

        void insert(std::uint64_t spatial_id, 
            const problem::OperationPoint& low_point,
            const problem::OperationPoint& high_point);

        void insert(std::uint64_t spatial_id, 
        const problem::OperationSpace& data_space);

        MemoryState() = default;
        
        MemoryState(std::uint64_t id, const problem::OperationSpace& data_space){
            data_spaces_.emplace(id, data_space);
        }

        const std::unordered_map<std::uint64_t, problem::OperationSpace>& 
            getDataSpaces() const {return data_spaces_;}


        static void set_workload(const problem::Workload* workload) {
            MemoryState::workload_ = workload;}

        void show() const;
    };

} // namespace TileFlow

} // namespace analysis