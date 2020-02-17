#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  
  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Type of Node value parameter*/
  using node_value_type = V;
  
  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Type of Edge value parameter. */
  using edge_value_type = E;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
  
  /** Invalid node index value */
  static size_type constexpr invalid_index = size_type(-1);
  
  /** Type of Edge length. */
  using length_type = double;
    
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(){
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    /** Invalid constructor for Node class. */
    Node(): graph_(nullptr), nid_(invalid_index) {
    }
      
    /** Return this node's position. */
    Point& position() {
      return graph_->nodes_[nid_].p_;
    }

    /** Return this node's position as a constant reference to a point. */
    const Point& position() const {
      return graph_->nodes_[nid_].p_;
    }

    /** Return this node's index for users, a number in the range [0, graph_size). */
    size_type index() const {
        return graph_->nodes_[nid_].idx_;
    }

    /** Return the value parameterstored in the node.
     *
     * @tparam node_value_type Type of the value parameter
     * @return The reference to the value stored in this node.
     *
     */
    node_value_type& value() {
        return graph_->nodes_[nid_].v_;
    }
      
    /** Return the constant value parameterstored in the node.
     *
     * @tparam node_value_type Type of the value parameter
     * @return The reference to the constant value stored in this node.
     *
     */
    const node_value_type& value() const {
        return graph_->nodes_[nid_].v_;
    }
    
    /** Return the number of incident edges. */
    size_type degree() const {
        return graph_->adjacent_nodes_[nid_].size();
    }
    
    /** Return the begin of the incident iterator of the node. */
    incident_iterator edge_begin() const {
        return IncidentIterator((*this).graph_, index(), 0);
    }
    
    /** Return the end of the incident iterator of the node. */
    incident_iterator edge_end() const {
        return IncidentIterator((*this).graph_, index(), degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph_== graph_) && (n.nid_ == nid_);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
        
        /** If two nodes are in the same graph, we compare their unique ids.*/
        if(graph_ == n.graph_){
            return nid_ < n.nid_;
        } else {
        /** If two nodes are not in the same graph, we compare their graph pointers.*/
            return std::less<Graph*>()(graph_, n.graph_);
        }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    // Pointer back to the Graph container
    Graph* graph_;
    // This node's unique indentification number
    size_type nid_;
    
    // Private Constructor to construct valid Node objects
    Node(const Graph* graph, size_type nid)
        : graph_(const_cast<Graph*>(graph)), nid_(nid) {
    }
    
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of active nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return active_nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.nid_ == active_nodes_[result_node.index()]
   * @post result_node.index() == nodes_[result_node.nid_].idx_
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
      
      /** Construct new node's information. */
      nodeinfo n;
      n.p_ = position;
      n.v_ = value;
      n.idx_ = active_nodes_.size();
      
      /** Construct new node's information. */
      size_type nid = nodes_.size();
      
      /** Add ths information of the node into the vector that stores all nodes' information in the graph. */
      nodes_.push_back(n);
      /** Add ths unique id of the node into the vector that stores all active nodes' uid in the graph. */
      active_nodes_.push_back(nid);
      
      /** Add the pair {node_nid_, an empty vector} into the map that stores all neighbors of the node.*/
      adjacent_nodes_.insert({nid, std::vector<size_type>()});
    
    return Node(this, nid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently an active Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      
      return n.graph_ == this && n.index() != invalid_index && n.index() < num_nodes();
    
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    
    assert(i < num_nodes());
    assert(i >= 0);
    return Node(this, active_nodes_[i]);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge(): graph_(nullptr), node1_idx_(invalid_index), node2_idx_(invalid_index) {
    }

    /** Return a node of this Edge */
    Node node1() const {
        size_type idx = (*graph_).nodes_[node1_idx_].idx_;
      return (*graph_).node(idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
        size_type idx = (*graph_).nodes_[node2_idx_].idx_;
        return (*graph_).node(idx);
    }
      
    /** Return the value parameterstored in the edge.
     *
     * @tparam edge_value_type Type of the value parameter
     * @return The reference to the value stored in this edge.
     *
     */
    edge_value_type& value() {
        if(node1_idx_ < node2_idx_){
            return graph_->edges_[node1_idx_][node2_idx_].v_;
        } else {
            return graph_->edges_[node2_idx_][node1_idx_].v_;
        }
    }
      
    /** Return the constant value parameterstored in the edge.
     *
     * @tparam node_value_type Type of the value parameter
     * @return The reference to the constant value stored in this node.
     *
     */
    const edge_value_type& value() const {
        if(node1_idx_ < node2_idx_){
            return graph_->edges_[node1_idx_][node2_idx_].v_;
        } else {
            return graph_->edges_[node2_idx_][node1_idx_].v_;
        }
    }
      
    /** Return the initial length of the edge. */
    length_type length() {
        return norm(node1().position() - node2().position());
    }
      
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.graph_ == graph_) && ((e.node1_idx_ == node1_idx_ && e.node2_idx_ == node2_idx_) ||
                                      (e.node2_idx_ == node1_idx_ && e.node1_idx_ == node2_idx_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        
        /** If two edges are in the same graph, compare their unique edge ids. */
        if(graph_ == e.graph_){
            
            size_type first_idx = std::min(node1_idx_, node2_idx_);
            size_type second_idx = std::max(node1_idx_, node2_idx_);
            
            size_type e_first_idx = std::min(e.node1_idx_, e.node2_idx_);
            size_type e_second_idx = std::max(e.node1_idx_, e.node2_idx_);
            
            return (*graph_).edges_[first_idx][second_idx].eid_ < (*graph_).edges_[e_first_idx][e_second_idx].eid_;
            
        } else {
        /** If two edges are not in the same graph, compare their pointers to graph. */
            return std::less<Graph*>()(graph_, e.graph_);
        }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    // Graph reference of the edge
    Graph* graph_;
    // Node1's index of the edge
    size_type node1_idx_;
    // Node2's index of the edge
    size_type node2_idx_;
    // Private constructor of the edge to construct valid Edge object later
    Edge(const Graph* graph, size_type node1_idx, size_type node2_idx)
      : graph_(const_cast<Graph*>(graph)), node1_idx_(node1_idx), node2_idx_(node2_idx) {
     }
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return active_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    
    assert(i >= 0);
    assert(i < num_edges());
    auto search = edges_look_up_.find(active_edges_[i]);
    size_type nid1 = std::get<0>(search->second);
    size_type nid2 = std::get<1>(search->second);
    
    return Edge(this, nid1, nid2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    
    /** Check whether the two nodes are different and they are in this graph. */
    if (!has_node(a) || !has_node(b) || a == b){
        return false;
    }
    
    /** Construct the valid pair of nodes indices of the edge, (smaller index, larger index). */
    std::tuple<size_type, size_type> pair = indices_pair(a, b);
    size_type nid1 = std::get<0>(pair);
    size_type nid2 = std::get<1>(pair);
    
    /** Check whether we have this edge using edges_ quickly. */
    if (edges_.count(nid1)){
        if (edges_.at(nid1).count(nid2)){
            edgeinfo e = edges_.at(nid1).at(nid2);
            if(e.idx_ != invalid_index){
                return true;
            }
        }
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   * @post The newly added egde is active.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    
    /** Check whether two nodes are different and  in this graoh. */
    assert(has_node(a));
    assert(has_node(b));
    assert(a != b);

    /** Construct the valid pair of nodes indices of the edge, (smaller index, larger index). */
    std::tuple<size_type, size_type> pair = indices_pair(a, b);
    size_type nid1 = std::get<0>(pair);
    size_type nid2 = std::get<1>(pair);
    
    /** If doesn't have this edge, add it  edges_map_ and edges_look_up_, as well as strong nodes information in nodes_map.
     All operations below are O(1). */
    if(!has_edge(a,b)){
        
        /** Construct the edge's useful information. */
        edgeinfo e;
        size_type eid = edges_look_up_.size();
        e.eid_ = eid;
        e.v_ = val;
        e.idx_ = active_edges_.size();
        
        /** Add information for this edge. */
        edges_[nid1][nid2] = e;
        /** Update active edges in ths graph. */
        active_edges_.push_back(eid);
        /** Add new neighbors for two nodes. */
        adjacent_nodes_[nid1].push_back(nid2);
        adjacent_nodes_[nid2].push_back(nid1);
        /** Add new edge pairs with their unique id as the key. */
        edges_look_up_.insert({eid, pair});
    } 

    return Edge(this, a.nid_, b.nid_);
        
  }
    
    /** A function to set the value of the node inside a graph with specific index.
     *
     * @param node_idx The index of the traget node.
     * @param value The value we'd like to update for the node.
     * @tparam node_value_type The type of the value parameter of the node.
     *
     * @post node(nodes_[node_idx].idx_).value() = value
     */
    void setNodeValue(size_type node_idx, node_value_type value){
        nodes_[node_idx].v_ = value;
    }
    
    /** A function to get the value of the node inside a graph with specific index.
     *
     * @param node_idx The index of the traget node.
     * @tparam node_value_type The type of the value parameter of the node.
     *
     * @return the value of the node with unique id _node_idx_
     */
    size_type getNodeValue(size_type node_idx){
        return nodes_[node_idx].v_;
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      
      /** Clear all structures for edges. */
      edges_.clear();
      active_edges_.clear();
      edges_look_up_.clear();
      adjacent_nodes_.clear();
      
      /** Clear all structures for nodes. */
      nodes_.clear();
      active_nodes_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator(): graph_(nullptr), cur_(invalid_index) {
    }

    /** The dereference operator for NodeIterator
     *
     * @return The node object that this NodeIterator refers to
     * @post result == (*graph_).node(cur_)
     */
    Node operator*() const {
        return Node(graph_, (*graph_).active_nodes_[cur_]);
    }

    /** The ++ operator for NodeIterator
     *
     * @return The reference to the next NodeIterator object
     * @post new cur_ - 1 = old cur_
     */
    NodeIterator& operator++() {
        cur_++;
        return *this;
    }
    
     /** The == operator for NodeIterator
      *
      * @param i The @a NodeIterator we compare
      * @return whether the two NodeIterators has same pointer to graph and same node iteration index
      */
    bool operator==(const NodeIterator& i) const {
        return (graph_ == i.graph_ && cur_ == i.cur_);
    }

   private:
    friend class Graph;
    
    /** The graph that this NodeIterator points to. */
    Graph* graph_;
    /** The index that this NodeIterator refers to.*/
    size_type cur_;
    
    /** Private Constructor for NodeIterator. */
    NodeIterator(const Graph* ptr, const size_type idx): 
      graph_(const_cast<Graph*>(ptr)), cur_(idx){
    }
    
  };

  /** The begin position of the NodeIterator of a graph. */
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }

  /** The end position of the NodeIterator of a graph. */
  node_iterator node_end() const {
      return NodeIterator(this, num_nodes());
  } 

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    /** The dereference operator for IncidentIterator
     *
     * @return The edge object that this IncidentIterator refers to
     * @post result.node1() == (*graph_).node(idx_)
     * @post (*graph_).has_edge(result.node1(), result.node2()) == true
     */
    Edge operator*() const {
        size_type uid = (*graph_).active_nodes_[idx_];
        size_type neighbor = (*graph_).adjacent_nodes_[uid][cur_];
        return Edge(graph_, uid, neighbor);
    }
    
    /** The ++ operator for IncidentIterator
     *
     * @return The reference to the next IncidentIterator object
     * @post new cur_ - 1 = old cur_
     */
    IncidentIterator& operator++() {
        cur_++;
        return *this;
    }
    
    /** The == operator for IncidentIterator
     *
     * @param ii The @a IncidentIterator we compare
     * @return whether the two IncidentIterator has same pointer to graph and same node iteration index and same neighbor node index
     */
    bool operator==(const IncidentIterator& ii) const {
        return (graph_ == ii.graph_ && idx_ == ii.idx_ && cur_ == ii.cur_);
    }

   private:
    friend class Graph;
    /** The graph pointer of this IncidentIterator. */
    Graph* graph_;
    /** The node index of this IncidentIterator. */
    size_type idx_;
    /** The nighbor node index of this IncidentIterator. */
    size_type cur_;
      
    /** Private constructor for this IncidentIterator. */
    IncidentIterator(const Graph* ptr, const size_type idx, const size_type cur):
        graph_(const_cast<Graph*>(ptr)), idx_(idx), cur_(cur){
    }
    
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** The dereference operator for EdgeIterator
     *
     * @return The edge object that this EdgeIterator refers to
     * @post (*graph_).has_edge(result) == true
     * @post (*graph_).edges_look_up_.at(cur_) == std::tuple<size_type, size_type> value(result.node1_idx, result.node2_idx)
     */
    Edge operator*() const {
        size_type idx1 = std::get<0>((*graph_).edges_look_up_[(*graph_).active_edges_[cur_]]);
        size_type idx2 = std::get<1>((*graph_).edges_look_up_[(*graph_).active_edges_[cur_]]);
        return Edge(graph_, idx1, idx2);
    }
    
    /** The ++ operator for EdgeIterator
     *
     * @return The reference to the next EdgeIterator object
     * @post new cur_ - 1 = old cur_
     */
    EdgeIterator& operator++() {
        cur_++;
        return *this;
    }
    
    /** The == operator for EdgeIterator
     *
     * @param ei The @a EdgeIterator we compare
     * @return whether the two EdgeIterators have same pointer to graph and same edge iteration index
     */
    bool operator==(const EdgeIterator& ei) const {
        return (graph_ == ei.graph_ && cur_ == ei.cur_);
    }

   private:
    friend class Graph;
    
    /** A pointer to the graph that this EdgeIterator belongs to. */
    Graph* graph_;
    
    /** The index of the Edge that this EdgeIterator refers to. */
    size_type cur_;
    
    /** A private constructor EdgeIterator object. */
    EdgeIterator(const Graph* ptr, const size_type idx):
        graph_(const_cast<Graph*>(ptr)), cur_(idx) {}
  };


  /** The begin position of this EdgeIterator. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** The end position of this EdgeIterator. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }
    
    /** Remove an edge of the graph, given two nodes connecting the edge.
     * @return a size_type false value if the graph doean't have this edge;
     *         the edge's user facing id if the graph has this edge
     * @post has_edge(@a a, @a b) == false
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
     *       Else,                        new num_edges() == old num_edges().
     * @post The egde's user facing index is invalid.
     * @post The last active edge will use this edge's user facing id now.
     *
     * Complexity: No more than O(num_nodes() + num_edges()). All operations
     *          are constant here. Except for the loop for neighbors. Since we
     *          assume the graph is sparse, it's definitely less than num_nodes.
     */
    size_type remove_edge(const Node& a, const Node& b){
        
        /** If the graph doesn't have this edge, return a false size_type value.*/
        if(!has_edge(a, b)){
            return (size_type) false;
        }

        /** Construct the valid pair of nodes indices of the edge, (smaller index, larger index). */
        std::tuple<size_type, size_type> pair = indices_pair(a, b);
        size_type nid1 = std::get<0>(pair);
        size_type nid2 = std::get<1>(pair);

        /** Find this edge's unique id and user facing id. */
        size_type eid = edges_.at(nid1).at(nid2).eid_;
        size_type idx = edges_.at(nid1).at(nid2).idx_;
        
        if(num_edges() > 1){
            /** Find last active edge's connecting nodes' nid in order.  */
            size_type last = num_edges() - 1;
            std::tuple<size_type, size_type> target_pair = indices_pair(edge(last).node1(),                 edge(last).node2());
            size_type target_nid1 = std::get<0>(target_pair);
            size_type target_nid2 = std::get<1>(target_pair);
            
            /** Swap the edge to remove and the last active edge,  and update the last active edge's user facing id. */
            std::swap(active_edges_[idx], active_edges_[last]);
            edges_.at(target_nid1).at(target_nid2).idx_ = idx;
        }
        
        /** Invalidate the orginal edge's user facing id. Finally delete the edge from active edges.  */
        edges_.at(nid1).at(nid2).idx_ = invalid_index;
        active_edges_.pop_back();

        /** Remove the edge entry in edge_look_up_. */
        edges_look_up_.erase(eid);
        
        /** Remove  b from the neighbor list of a. */
        for(size_type j = 0; j < a.degree(); ++j){
            if(adjacent_nodes_[a.nid_][j] == b.nid_){
                std::swap(adjacent_nodes_[a.nid_][j], adjacent_nodes_[a.nid_][a.degree()-1]);
                adjacent_nodes_[a.nid_].pop_back();
                break;
            }
        }

        /** Remove  a from the neighbor list of b. */
        for(size_type j = 0; j < b.degree(); ++j){
            if(adjacent_nodes_[b.nid_][j] == a.nid_){
                std::swap(adjacent_nodes_[b.nid_][j], adjacent_nodes_[b.nid_][b.degree()-1]);
                adjacent_nodes_[b.nid_].pop_back();
                break;
            }
        }
        
        return idx;
        
    }
    
    /** Remove an edge of the graph, given an edge.
     * @return a size_type false value if the graph doean't have this edge;
     *         the edge's user facing id if the graph has this edge
     * @post has_edge(@a a, @a b) == false
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
     *       Else,                        new num_edges() == old num_edges().
     * @post The egde's user facing index is invalid.
     * @post The last active edge will use this edge's user facing id now.
     *
     * Complexity: No more than O(num_nodes() + num_edges()). Same
     *          complexity as remove_edge(Node1, Node2).
     */
    size_type remove_edge(const Edge& e){
        
        /** Find two nodes of the edge and call remove_edge(node1, node1). */
        const Node a = e.node1();
        const Node b = e.node2();
        
        return remove_edge(a, b);
        
    }
    
    /** Remove an edge of the graph, given an edge iterator.
     * @return the ending position for edge iterator if passing in the end
     *         the begin position for edge iterator otherwise
     * @post has_edge(@a a, @a b) == false
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
     *       Else,                        new num_edges() == old num_edges().
     * @post The egde's user facing index is invalid.
     * @post The last active edge will use this edge's user facing id now.
     *
     * Complexity: No more than O(num_nodes() + num_edges()). Same
     *          complexity as remove_edge(Node1, Node2).
     */
    edge_iterator remove_edge(edge_iterator e_it){

        if(e_it == edge_end()){
            return edge_end();
        } else {
            Edge e = *e_it;
            Node a = e.node1();
            Node b = e.node2();
            if(has_edge(a, b)){
                remove_edge(a, b);
            }
            return edge_begin();
        }
    }
    
    /** Remove a node of the graph, returning the user facing index of the node
     * if the graph has the node, otherwise return a false size_type value.
     * @post new num_nodes() == old num_nodes() - 1 if the graph has the node
     * @post new num_nodes() == old num_nodes()  if the graph doesn't have the node
     * @post The node n is inactive after the removal.
     *
     * Complexity: No more than O(num_nodes())
     *  Since we assume the graph is sparse, the loop over all incident edges will be less
     *  than the number of nodes by assumption.
     */
    size_type remove_node(const Node& n){
        
        /** If the graph doesn't have this node, return a false size_type value. */
        if(!has_node(n)){
            return (size_type) false;
        }
       
        /** Remove all incident edges to this node. */
        while(n.degree() > 0){
            remove_edge(Edge(this, n.nid_, adjacent_nodes_[n.nid_][0]));
        }

        /** Swap this node with the last active node in this graph, exchange their user facing ids and remove the target node. */
        size_type idx = n.index();
        std::swap(active_nodes_[idx], active_nodes_[num_nodes()-1]);
        nodes_[active_nodes_[idx]].idx_ = idx;
        nodes_[active_nodes_[num_nodes()-1]].idx_ = invalid_index;
        active_nodes_.pop_back();
        
        return idx;
        
    }
    
    /** Remove a node of the graph, returning the end position of the node iterator
     * if given the ending iterator, otherwise return the begin position of the iterator
     * @post new num_nodes() == old num_nodes() - 1 if the graph has the node
     * @post new num_nodes() == old num_nodes()  if the graph doesn't have the node
     * @post The node n is inactive after the removal.
     *
     * Complexity: No more than O(num_nodes())
     *  Since we assume the graph is sparse, the loop over all incident edges will be less
     *  than the number of nodes by assumption.
     */
    node_iterator remove_node(node_iterator n_it){
        
        if(n_it == node_end()){
            return node_end();
        } else {
            Node n = *n_it;
            if(has_node(n)){
                remove_node(n);
            }
            return node_begin();
        }
    }

 private:
    
  // Struct for node
  struct nodeinfo
  {
    // The position of the node
    Point p_;
    // The value of the node
    node_value_type v_;
    // The user facing id of the node
    size_type idx_;
  };
    
  // Struct for edge
  struct edgeinfo
  {
    // The unique id of the edge
    size_type eid_;
    // The value of the edge
    edge_value_type v_;
    // The user facing id of the edge
    size_type idx_;
  };
    
  // A vector attribute to store each node's information in this graph
  std::vector<nodeinfo> nodes_;
    
  // A vector to keep track of each active node's unique id in this graph
  std::vector<size_type> active_nodes_;
    
  // A nested map attribute to store each edge's information in this graph
  std::unordered_map<size_type, std::unordered_map<size_type, edgeinfo>> edges_;
    
  // A vector to keep track of each active edge's unique id in this graph
  std::vector<size_type> active_edges_;

  // A vector attribute to look up edges by their unique ids quickly
  std::unordered_map<size_type, std::tuple<size_type, size_type>> edges_look_up_;
    
  // A map attribute to store all neighbors of each node in this graph
  std::unordered_map<size_type, std::vector<size_type>> adjacent_nodes_;

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** A helper function that constructs valid node indices pair for each edge.
   *
   * @param a The reference to the first @a Node
   * @tparam b The reference to the second @a Node
   * @return A tuple of size two: (smaller node index, larger node index)
      
   * @post std::get<0>(result) < std::get<1>(result)
   */
  std::tuple<size_type, size_type> indices_pair(const Node& a, const Node& b) const {
    
    size_type nid1, nid2;
    if(a.nid_ < b.nid_){
        nid1 = a.nid_;
        nid2 = b.nid_;
    } else {
        nid1 = b.nid_;
        nid2 = a.nid_;
    }
    std::tuple<size_type, size_type> pair{nid1, nid2};
      
    return pair;
  }
  
};

#endif // CME212_GRAPH_HPP
