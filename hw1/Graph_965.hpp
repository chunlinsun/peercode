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
template <typename V>
class Graph {
 private:
  
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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
    const Point& position() const {
      return graph_->nodes_vec_[nid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return nid_;
    }

    /** Return the value parameterstored in the node.
     *
     * @tparam node_value_type Type of the value parameter
     * @return The reference to the value stored in this node.
     *
     */
    node_value_type& value() {
        return graph_->values_vec_[nid_];
    }
      
    /** Return the constant value parameterstored in the node.
     *
     * @tparam node_value_type Type of the value parameter
     * @return The reference to the constant value stored in this node.
     *
     */
    const node_value_type& value() const {
        return graph_->values_vec_[nid_];
    }
    
    /** Return the number of incident edges. */
    size_type degree() const {
        return graph_->nodes_map_[nid_].size();
    }
    
    /** Return the begin of the incident iterator of the node. */
    incident_iterator edge_begin() const {
        return IncidentIterator((*this).graph_, nid_, 0);
    }
    
    /** Return the end of the incident iterator of the node. */
    incident_iterator edge_end() const {
        return IncidentIterator((*this).graph_, nid_, degree());
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
      return nid_ < n.nid_;
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

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_vec_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    
    /** Add ths position of the node into the vector that stores nodes' position of the graph. */
    nodes_vec_.push_back(position);
    /** Add the value of the node into the vector that stores nodes' values of the graph. */
    values_vec_.push_back(value);
    /** Add the pair {node_index, an empty vector} into the map that stores all neighbors of the node.*/
    nodes_map_.insert({num_nodes()-1, std::vector<size_type>()});
    
    return Node(this, num_nodes()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this && n.nid_ < num_nodes();
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
    return Node(this, i);   
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
      return graph_->node(node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_idx_);     
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.graph_ == graph_) && (e.node1_idx_ == node1_idx_ && e.node2_idx_ == node2_idx_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        /** Compare the smaller index of two edges first. */
        if(node1_idx_ < e.node1_idx_){
            return false;
        /** If the smaller index is the same, compare the larger index of two edges. */
        } else if (node1_idx_ > e.node1_idx_) {
            return true;
        } else {
            return node2_idx_ < e.node2_idx_;
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
    return edges_map_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    
    assert(i >= 0);
    assert(i < num_edges());
    auto search = edges_map_.find(i);
    size_type idx1 = std::get<0>(search->second);
    size_type idx2 = std::get<1>(search->second);
    
    return Edge(this, idx1, idx2);       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    
    assert(has_node(a));
    assert(has_node(b));
    if (a.nid_ == b.nid_){
        return false;
    }
    
    /** Construct the valid pair of nodes indices of the edge, (smaller index, larger index). */
    std::tuple<size_type, size_type> value = indices_pair(a, b);
    size_type idx1 = std::get<0>(value);
    size_type idx2 = std::get<1>(value);
    
    /** Check whether we have this edge using edges_loop_up_ quickly. */
    if (edges_look_up_.count(idx1)){
        if (edges_look_up_.at(idx1).count(idx2)){
            return true;
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
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    
    assert(has_node(a));
    assert(has_node(b));
    assert(a.nid_ != b.nid_);

    /** Construct the valid pair of nodes indices of the edge, (smaller index, larger index). */
    std::tuple<size_type, size_type> value = indices_pair(a, b);
    size_type idx1 = std::get<0>(value);
    size_type idx2 = std::get<1>(value);
    
    /** If doesn't have this edge, add it  edges_map_ and edges_look_up_, as well as strong nodes information in nodes_map. */
    if(!has_edge(a,b)){
        size_type eid = num_edges();
        
        edges_map_.insert({eid, value});
        nodes_map_[idx1].push_back(idx2);
        nodes_map_[idx2].push_back(idx1);
        edges_look_up_[idx1][idx2] = eid;
    } 

    return Edge(this, idx1, idx2);
        
  }
    
    /** A function to set the value of the node inside a graph with specific index.
     *
     * @param node_idx The index of the traget node.
     * @param value The value we'd like to update for the node.
     * @tparam node_value_type The type of the value parameter of the node.
     
     * @post node(node_idx).value() = value
     *
     */
    void setNodeValue(size_type node_idx, node_value_type value){
        values_vec_[node_idx] = value;
    }
    
    /** A function to get the value of the node inside a graph with specific index.
     *
     * @param node_idx The index of the traget node.
     * @tparam node_value_type The type of the value parameter of the node.
     *
     * @return the value of the node at index _node_idx_
     */
    size_type getNodeValue(size_type node_idx){
        return values_vec_[node_idx];
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
//--design_0
//--did not clear everything
//--START
  void clear() {
    nodes_vec_.clear();
    edges_map_.clear();
    values_vec_.clear();
    nodes_map_.clear();
  }
//--END

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
    NodeIterator(): nodes_ptr(nullptr), iter_idx(invalid_index) {
    }

    // HW1 #2: YOUR CODE HERE
    /** The dereference operator for NodeIterator
     *
     * @return The node object that this NodeIterator refers to
     * @post result == (*nodes_ptr).node(iter_idx)
     */
    Node operator*() const {
        return Node(nodes_ptr, iter_idx);
    }

    /** The ++ operator for NodeIterator
     *
     * @return The reference to the next NodeIterator object
     * @post new iter_idx - 1 = old iter_idx
     */
    NodeIterator& operator++() {
        iter_idx++;
        return *this;
    }
    
     /** The == operator for NodeIterator
      *
      * @param i The @a NodeIterator we compare
      * @return whether the two NodeIterators has same pointer to graph and same node iteration index
      */
    bool operator==(const NodeIterator& i) const {
        return (nodes_ptr == i.nodes_ptr && iter_idx == i.iter_idx);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    
    /** The graph that this NodeIterator points to. */
    Graph* nodes_ptr;
    /** The index that this NodeIterator refers to.*/
    size_type iter_idx;
    
    /** Private Constructor for NodeIterator. */
    NodeIterator(const Graph* ptr, const size_type idx): 
      nodes_ptr(const_cast<Graph*>(ptr)), iter_idx(idx){
    }
    
  };

  // HW1 #2: YOUR CODE HERE
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
     * @post result.node1() == (*graph_ptr).node(node_idx)
     * @post (*graph_ptr).has_edge(result.node1(), result.node2()) == true
     */
    Edge operator*() const {
        size_type neighbor = (*graph_ptr).nodes_map_[node_idx][neighbor_idx];
        return Edge(graph_ptr, node_idx, neighbor);
    }
    
    /** The ++ operator for IncidentIterator
     *
     * @return The reference to the next IncidentIterator object
     * @post new neighbor_idx - 1 = old neighbor_idx
     */
    IncidentIterator& operator++() {
        neighbor_idx++;
        return *this;
    }
    
    /** The == operator for IncidentIterator
     *
     * @param ii The @a IncidentIterator we compare
     * @return whether the two IncidentIterator has same pointer to graph and same node iteration index and same neighbor node index
     */
    bool operator==(const IncidentIterator& ii) const {
        return (graph_ptr == ii.graph_ptr && node_idx == ii.node_idx && neighbor_idx == ii.neighbor_idx);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    /** The graph pointer of this IncidentIterator. */
    Graph* graph_ptr;
    /** The node index of this IncidentIterator. */
    size_type node_idx;
    /** The nighbor node index of this IncidentIterator. */
    size_type neighbor_idx;
      
    /** Private constructor for this IncidentIterator. */
    IncidentIterator(const Graph* ptr, const size_type no_idx, const size_type ne_idx):
        graph_ptr(const_cast<Graph*>(ptr)), node_idx(no_idx), neighbor_idx(ne_idx){
    }
    
  };

//--documentation_-1
//--good doc strings
//--END

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

    // HW1 #5: YOUR CODE HERE
    /** The dereference operator for EdgeIterator
     *
     * @return The edge object that this EdgeIterator refers to
     * @post (*graph_ptr).has_edge(result) == true
     * @post (*graph_ptr).edges_map.at(edge_idx) == std::tuple<size_type, size_type> value(result.node1_idx, result.node2_idx)
     */
    Edge operator*() const {
        size_type idx1 = std::get<0>((*graph_ptr).edges_map_[edge_idx]);
        size_type idx2 = std::get<1>((*graph_ptr).edges_map_[edge_idx]);
        return Edge(graph_ptr, idx1, idx2);
    }
    
    /** The ++ operator for EdgeIterator
     *
     * @return The reference to the next EdgeIterator object
     * @post new edge_idx - 1 = old edge_idx
     */
    EdgeIterator& operator++() {
        edge_idx++;
        return *this;
    }
    
    /** The == operator for EdgeIterator
     *
     * @param ei The @a EdgeIterator we compare
     * @return whether the two EdgeIterators have same pointer to graph and same edge iteration index
     */
    bool operator==(const EdgeIterator& ei) const {
        return (graph_ptr == ei.graph_ptr && edge_idx == ei.edge_idx);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
    /** A pointer to the graph that this EdgeIterator belongs to. */
    Graph* graph_ptr;
    
    /** The index of the Edge that this EdgeIterator refers to. */
    size_type edge_idx;
    
    /** A private constructor EdgeIterator object. */
    EdgeIterator(const Graph* ptr, const size_type idx):
        graph_ptr(const_cast<Graph*>(ptr)), edge_idx(idx) {}
  };


  // HW1 #5: YOUR CODE HERE
  /** The begin position of this EdgeIterator. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** The end position of this EdgeIterator. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:
  // A vector attribute to store all nodes in this graph
  std::vector<Point> nodes_vec_;

  // A vector attribute to store all edges in this graph
  std::unordered_map<size_type, std::tuple<size_type, size_type>> edges_map_;
 
  // A vector attribute to store values for each node in this graph
  std::vector<node_value_type> values_vec_;
    
  // A map attribute to store all neighbors of each node
  std::unordered_map<size_type, std::vector<size_type>> nodes_map_;
    
  // A nested map attribute to look up edges in a quick manner
  std::unordered_map<size_type, std::unordered_map<size_type,size_type>> edges_look_up_;

  // Private constructor for Graph class
  Graph(std::vector<Point> points, std::unordered_map<size_type, std::tuple<size_type, size_type>> edges,
        std::vector<node_value_type> values, std::unordered_map<size_type, std::vector<size_type>> nodes,
        std::unordered_map<size_type, std::unordered_map<size_type,size_type>> lookup)
      : nodes_vec_(points), edges_map_(edges), values_vec_(values), nodes_map_(nodes), edges_look_up_(lookup) {
  }

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
    
    size_type idx1, idx2;
    if(a.nid_ < b.nid_){
        idx1 = a.nid_;
        idx2 = b.nid_;
    } else {
        idx1 = b.nid_;
        idx2 = a.nid_;
    }
    std::tuple<size_type, size_type> value{idx1, idx2};
      
    return value;
  }
};

#endif // CME212_GRAPH_HPP
