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
        problem::Shape common_shape_;
        Workload::FactorizedBounds factorized_bounds;
        config::CompoundConfigNode coeffs_;

    public:
        Workloads() {common_shape_.UsesFlattening = false;}
        bool add_workload(const std::string& name, std::shared_ptr<problem::TileFlow::Workload>& workload);
        std::shared_ptr<problem::TileFlow::Workload> get_workload(const std::string & op_name) const {
            return workloads_.at(op_name);
        }
        void set_io(const std::vector<std::string>& ins, const std::vector<std::string>& outs);
        const std::vector<std::string>& get_ins() const {return ins_;}
        const std::vector<std::string>& get_outs() const {return outs_;}
        void set_coeff(const config::CompoundConfigNode& coeffs) {coeffs_ = coeffs;}
        int lookup_coeff(const std::string& name) const {int value = 1; coeffs_.lookupValue(name, value); return value;}
        void Print();
        const problem::Shape& get_shape() {return common_shape_;}
        void set_factorized_bound(const std::string& dim, int bound);
        const Workload::FactorizedBounds& get_factorized_bound() {return factorized_bounds;}
        friend class Workload;
        
    };

    void ParseWorkloads(config::CompoundConfigNode config, Workloads& workloads_);

} // namespace TileFlow 

} // namespace problem