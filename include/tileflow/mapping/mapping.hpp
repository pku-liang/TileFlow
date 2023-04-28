#pragma once 

#include "model/engine.hpp"
#include "compound-config/compound-config.hpp"
#include "mapping/mapping.hpp"


#include "tileflow/problem/problem.hpp"
#include "tileflow/mapping/loop.hpp"
#include "tileflow/common.hpp"
#include "tileflow/mapper/expr.hpp"


using TileFlow::SymbolTable;

namespace mapping {

namespace TileFlow {

class Node;
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
public:
    virtual void run (const Node*);
};

struct ActiveTensor {
    std::set<problem::Shape::DataSpaceID> 
    read_tensors, fill_tensors, update_tensors, wb_tensors;
};

class Node {
public: 
    enum type_t{
        Tile,
        Op,
        Scope
    };
protected:
    static const std::unordered_map<type_t, std::string> type2name_; 
    Node::type_t type_;
    std::string name_;
    mutable const Node* parent_ = nullptr;
    std::vector<const Node*> children_;
    mutable std::string storage_level_name_;
    mutable unsigned storage_level_ = unsigned(-1);
    std::vector<std::string> bypassed_;
    bool profile_ = true;

    mutable ActiveTensor active_tensors_;

    void ParseStorageLevel(config::CompoundConfigNode config);
    std::unordered_map<std::string, std::pair<int, int> > ParseFactors(const std::string & factors);
    std::vector<std::string> ParsePermutations(const std::string& buffer);
    void display_active_tensors(std::string prefix, std::ostream& o = std::cout) const;

public: 
    bool is_bypassed(const std::string & tensor) const {
        return std::find(bypassed_.begin(), bypassed_.end(), tensor) != bypassed_.end();}
    bool is_profile() const {return profile_;}
    Node(type_t, config::CompoundConfigNode);
    
    unsigned get_storage_level() const {return storage_level_;}
    void set_storage_level(unsigned storage_level, const std::string& storage_name) const {
        storage_level_ = storage_level; storage_level_name_ = storage_name;}
    std::string get_storage_name() const {return storage_level_name_;}
    std::string get_name() const {return name_;}
    type_t get_type() const {return type_;}

    ActiveTensor& get_active_tensors() const {return active_tensors_;}
    
    void add_child(const Node* child);
    
    void replace_child(const Node* old_child, const Node* new_child){
        auto iter = find(children_.begin(), children_.end(), old_child);
        if (iter != children_.end()) {
            iter = children_.erase(iter);
            children_.insert(iter, new_child);
            new_child->set_parent(this);
        }
    }
    const std::vector<const Node*> get_children() const{
        return children_;
    }
    void set_children(const std::vector<const Node*>& children) {
        children_ = children;
    }

    void reset_children() {children_.clear();}

    const Node* get_first_child() const {
        assert(children_.size());
        return children_.front();
    }
    void set_parent(const Node* parent) const {parent_ = parent;}
    inline const Node* get_parent() const {return parent_;}

    virtual void display(std::string prefix, bool = true, const SymbolTable* = nullptr, std::ostream& o = std::cout) const {o << prefix << std::endl;}
    virtual void accept(Visitor* visitor) const = 0;
    virtual ~Node() {for (auto node: children_) delete node;}
    friend class Visitor;
};

class ScopeNode: public Node {
public:
    enum type_t {
        Sequential,
        Sharing,
        Parallel,
        Pipeline
    };
    ScopeNode(config::CompoundConfigNode config);
    void display(std::string prefix, bool recursive, const SymbolTable* = nullptr, std::ostream& = std::cout) const override;
    void accept(Visitor* visitor) const {visitor->visitScope(this);}
    ScopeNode::type_t get_scope_type() const {return type;}

private: 
    ScopeNode::type_t type;
    
};

class TileNode: public Node {
public:
    enum type_t {
        Temporal,
        Spatial
    };
private:

    // std::pair<int, int> represent the <end, residual end>
    std::vector<loop::TileFlow::Descriptor> loopnests_;
    TileNode::type_t type_;
    bool multicast_enabled_ = true;

public:
    TileNode(config::CompoundConfigNode config);
    void display(std::string prefix, bool recursive, const SymbolTable* = nullptr, std::ostream& = std::cout) const override;
    void accept(Visitor* visitor) const {visitor->visitTile(this);}
    bool is_spatial() const {return type_ == Spatial;}
    bool is_multicast_enabled() const {return multicast_enabled_;}
    TileNode::type_t get_tile_type() const {return type_;}
    
    loop::Nest constructLoopNest(const SymbolTable* symbol_table = nullptr) const;
    size_t n_level() const {return loopnests_.size();}
    const std::vector<loop::TileFlow::Descriptor>& get_loops() const {return loopnests_;}
};

class OpNode: public Node {
    std::string op_name_;
    int op_index_;
    std::shared_ptr<problem::TileFlow::Workload> p_workload;
public:
    OpNode(config::CompoundConfigNode config);
    void display(std::string prefix, bool recursive, const SymbolTable* = nullptr, std::ostream& = std::cout) const override;
    const std::string & get_name() const {return op_name_;}
    void accept(Visitor* visitor) const {visitor->visitOp(this);}
    const std::shared_ptr<problem::TileFlow::Workload>& get_workload() const {return p_workload;}
};

struct Mapping {
    std::map<unsigned, std::uint64_t> fanoutX_map;
    std::map<unsigned, std::uint64_t> fanoutY_map;  
    Node * root = nullptr;
    void Print() const;
};

Mapping ParseAndConstruct(config::CompoundConfigNode config,
                          model::Engine::Specs& arch_specs,
                          const problem::TileFlow::Workloads& workload);

class CollectOpNode: public mapping::TileFlow::Visitor {
    void visitOp(const OpNode* node) override {
        opnodes_.push_back(node);
    }
    std::vector<const OpNode*> opnodes_;
public: 
    std::vector<const OpNode*> collectOpNodes(Node* root){
        opnodes_.clear();
        root->accept(this);
        return std::move(opnodes_);
    }
};
    
class CollectTileNode: public mapping::TileFlow::Visitor {
    void visitTile(const TileNode* node) override {
        if (node->get_tile_type() == type_)
            nodes_.push_back(node);
        for (auto child: node->get_children())
            child->accept(this);
    }
    std::vector<const TileNode*> nodes_;
    TileNode::type_t type_;
public: 
    CollectTileNode(TileNode::type_t type = TileNode::Temporal): type_(type){}
    std::vector<const TileNode*> operator() (const Node*root){
        root->accept(this);
        return nodes_;
    }
};

}  // namespace TileFlow  
} // namespace mapping
