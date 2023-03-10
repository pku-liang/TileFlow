#pragma once 

#include "model/engine.hpp"
#include "compound-config/compound-config.hpp"
#include "mapping/mapping.hpp"


#include "timeloopX/problem/problem.hpp"
#include "timeloopX/mapping/loop.hpp"
#include "timeloopX/common.hpp"

namespace mapping {

namespace TimeloopX {

class ScopeNode;
class TileNode;
class OpNode;

class Visitor {
protected:
    virtual void visitScope(const ScopeNode*);
    virtual void visitTile(const TileNode*);
    virtual void visitOp(const OpNode*);
    friend class TileNode;
    friend class ScopeNode;
    friend class OpNode;
};

class Node {
public: 
    enum type_t{
        Tile,
        Op,
        Scope
    };
protected:
    Node::type_t type_;
    Node* parent_ = nullptr;
    std::vector<Node*> children_;
public: 
    
    std::unordered_map<std::string, std::pair<int, int> > ParseFactors(const std::string & factors);
    std::vector<std::string> ParsePermutations(const std::string& buffer);
    unsigned ParseStorageLevel(config::CompoundConfigNode config);

    Node(type_t t_): type_(t_) {}
    type_t get_type() const {return type_;}
    void add_child(Node* child) {assert(child != nullptr); children_.push_back(child); child->set_parent(this);}
    const std::vector<const Node*> get_children() const{
        std::vector<const Node*> retval;
        for (auto child: children_) retval.push_back(child);
        return retval;
    }
    void set_parent(Node* parent) {parent_ = parent;}
    inline const Node* get_parent() const {return parent_;}

    virtual void display(std::string prefix, bool recursive = true) const {std::cout << prefix << std::endl;}
    virtual void accept(Visitor* visitor) const = 0;
    virtual ~Node() {for (auto node: children_) delete node;}
    friend class Visitor;
};

class ScopeNode: public Node {
    public:
    enum type_t {
        Sequential,
        Parallel,
        Pipeline
    }type;

    ScopeNode(config::CompoundConfigNode config);
    void display(std::string prefix, bool recursive) const override;
    void accept(Visitor* visitor) const {visitor->visitScope(this);}
};

class TileNode: public Node {
public:
    enum type_t {
        Temporal,
        Spatial
    };
private:

    // std::pair<int, int> represent the <end, residual end>
    std::vector<loop::TimeloopX::Descriptor> loopnests_;
    unsigned storage_level_;
    TileNode::type_t type;

public:
    TileNode(config::CompoundConfigNode config);
    void display(std::string prefix, bool recursive) const override;
    void accept(Visitor* visitor) const {visitor->visitTile(this);}
    bool is_spatial() const {return type == Spatial;}
    loop::Nest constructLoopNest(
        const std::map<std::string, problem::Shape::FactorizedDimensionID>&) const;
    size_t n_level() const {return loopnests_.size();}
};

class OpNode: public Node {
    std::string name_;
    int op_index_;
    std::shared_ptr<problem::TimeloopX::Workload> p_workload;
    // op dimension --> runtime iteration
    std::unordered_map<std::string, std::string> binding_; 
public:
    OpNode(config::CompoundConfigNode config);
    void display(std::string prefix, bool recursive) const override;
    const std::string & get_name() const {return name_;}
    void accept(Visitor* visitor) const {visitor->visitOp(this);}
    const std::shared_ptr<problem::TimeloopX::Workload>& get_workload() const {return p_workload;}
};

struct Mapping {
    std::map<unsigned, std::uint64_t> fanoutX_map;
    std::map<unsigned, std::uint64_t> fanoutY_map;  
    Node * root = nullptr;
    void Print();
};

Mapping ParseAndConstruct(config::CompoundConfigNode config,
                          model::Engine::Specs& arch_specs,
                          const problem::TimeloopX::Workloads& workload);

}  // namespace TimeloopX  
} // namespace mapping
