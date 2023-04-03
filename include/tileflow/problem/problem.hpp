#pragma once 

#include "mapping/parser.hpp"

#include "tileflow/common.hpp"


namespace problem {

namespace TileFlow {

    typedef unsigned TensorIndex;

    class Workloads;

    class Workload: public problem::Workload {
        std::vector<std::string> ins_;
        std::string out_;
        std::string name_;  
        Workloads& workloads_;
        bool binding_applied = false;
    public: 
        Workload(Workloads& workloads): workloads_(workloads){}
        inline void set_name(const std::string & name){name_ = name;}
        void set_io(const std::vector<std::string>& ins, const std::vector<std::string>& outs);
        inline const std::vector<std::string>& get_ins() const { return ins_; }
        inline const std::string & get_out() const {return out_;}
        inline const std::string & get_name() const {return name_;} 
        void Print();
        friend class Workloads;
        void apply_binding(const std::unordered_map<std::string, std::string>& binding);
    };

    class Workloads{
        std::unordered_map<std::string, std::shared_ptr<problem::TileFlow::Workload> > workloads_;
        std::vector<std::string> ins_;
        std::vector<std::string> outs_;
        config::CompoundConfigNode coeffs_;

        problem::Workload::FactorizedBounds factorized_bounds_;
        problem::Workload::FlattenedBounds flattened_bounds_;
        problem::Workload::Coefficients coefficients_;
        problem::Workload::Densities densities_;
        problem::Shape common_shape_;

        mutable problem::Workload workload_;
        mutable bool workload_constructed_ = false;

    public:
        Workloads() {common_shape_.UsesFlattening = false; coefficients_[-1] = 1;}
        bool add_workload(const std::string& name, std::shared_ptr<problem::TileFlow::Workload>& workload);
        std::shared_ptr<problem::TileFlow::Workload> get_workload(const std::string & op_name) const {
            return workloads_.at(op_name);
        }
        void set_io(const std::vector<std::string>& ins, const std::vector<std::string>& outs);
        void set_coeffs(const config::CompoundConfigNode& coeffs);
        void set_factorized_bound(const std::string& dim, int bound);

        void Print();
        
        const std::vector<std::string>& get_ins() const {return ins_;}
        const std::vector<std::string>& get_outs() const {return outs_;}
        const problem::Workload& get_workload() const;
        const problem::Shape& get_shape() const {return common_shape_;}

        friend class Workload;
    };

    void ParseWorkloads(config::CompoundConfigNode config, Workloads& workloads_);

} // namespace TileFlow 

} // namespace problem