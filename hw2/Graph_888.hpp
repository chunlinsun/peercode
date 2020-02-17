#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

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
 private:

  // HW0: YOUR CODE HERE
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

  typedef V node_value_type;
  typedef E edge_value_type;
  
  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    : nodes_(), i2u_(), adj_list_(), num_edges_(0) {}

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
  class Node : private totally_ordered<Node> {
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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].p_;
    }

    Point& position() {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return a reference to this node's value. */
    node_value_type& value() {
      return graph_->nodes_[uid_].v_;
    }
    
    /** Return a const reference to this node's value. */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].v_;
    }

    /** Return this node's degree, i.e. the number of edges incident on this node. */
    size_type degree() const {
      return graph_->adj_list_[uid_].size();
    }

    /** Return an iterator pointing to the first edge incident on this node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    /** Return an iterator pointing to one past the last edge incident on this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (n.graph_ == graph_) && (n.uid_ == uid_);
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
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_) {
        return uid_ < n.uid_;
      }
      else {
        return graph_ < n.graph_;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer back to the Graph
    graph_type* graph_;
    // This node's unique identification number
    size_type uid_;
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type new_node_id = num_nodes();
    internal_node new_internal_node = internal_node(position, new_node_id, value);
    nodes_.push_back(new_internal_node);
    i2u_.push_back(nodes_.size() - 1);
    adj_list_.push_back(std::vector<std::pair<size_type, edge_value_type>>());
    return Node(this, nodes_.size() - 1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_ == this) && (n.index() < i2u_.size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i2u_[i]);
  }

  /** Remove a node from the graph, if it is valid node of the graph
   * @param[in] n The node to be removed
   * @return the number of nodes removed
   * @post If @a n is a valid node of the graph, new num_nodes() == old num_nodes() - 1
   *       Else,                                 new num_nodes() == old num_nodes()
   * @post If @a n is a valid node of the graph, new num_edges() == old num_edges() - n.degree()
   *       Else,                                 new num_edges() == old num_edges()
   * 
   * Complexity: O((max_degree)^2)
   * 
   * Invalidate: node with uid_ == n.uid_
   * Invalidate: edges with node1() == n or node2() == n
   * Invalidate: node_iterator representing the invalidated node
   * Invalidate: edge_iterators representing the invalidated edges
  */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) {
      return 0;
    }
    size_type idx = n.index();
    auto it = n.edge_begin();
    auto end_it = n.edge_end();
    while (it != end_it) {
      const Edge& e = *it;
      remove_edge(e);
      end_it = n.edge_end();
    }
    i2u_[idx] = i2u_[num_nodes()-1];
    size_type uid = i2u_[idx];
    nodes_[uid].idx_ = idx;
    i2u_.pop_back();
    return 1;
  }

  /** Remove a node from the graph, if it is valid node of the graph
   * @param[in] n_it The iterator to the node to be removed
   * @return The iterator to the next node
   * @post If @a *n_it is a valid node of the graph, new num_nodes() == old num_nodes() - 1
   *       Else,                                     new num_nodes() == old num_nodes()
   * @post If @a *n_it is a valid node of the graph, new num_edges() == old num_edges() - n.degree()
   *       Else,                                     new num_edges() == old num_edges()
   * 
   * Complexity: O((max_degree)^2)
   * 
   * Invalidate: node with uid_ == *n_it.uid_
   * Invalidate: edges with node1() == *n_it or node2() == *n_it
   * Invalidate: node_iterator representing the invalidated node
   * Invalidate: edge_iterators representing the invalidated edges
  */
  node_iterator remove_node(node_iterator n_it) {
    const Node& n = *n_it;
    remove_node(n);
    return n_it;
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type node2_id_ = graph_->adj_list_[node1_id_][node2_adj_id_].first;
      return Node(graph_, node2_id_);
    }

    /** Return node2's index in the adjacency list of node1 */
    size_type node2_adj_id() const {
      return node2_adj_id_;
    }

    /** Return the length of this Edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return a reference to this edge's value. */
    edge_value_type& value() {
      return graph_->adj_list_[node1_id_][node2_adj_id_].second;
    }

    /** Return a const reference to this edge's value. */
    const edge_value_type& value() const {
      return graph_->adj_list_[node1_id_][node2_adj_id_].second;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_) && 
             ((node1_id_ == e.node1_id_ && node2_adj_id_ == e.node2_adj_id_) || 
              (node1_id_ == e.node2_adj_id_ && node2_adj_id_ == e.node1_id_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_) {
        if (node1_id_ == e.node1_id_) {
          return node2_adj_id_ < e.node2_adj_id_;
        }
        else {
          return node1_id_ < e.node1_id_;
        }
      }
      else {
        return graph_ < e.graph_;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph
    graph_type* graph_;
    // Node1's unique identification number
    size_type node1_id_;
    // Node2's index in the adjacency list of node1
    size_type node2_adj_id_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type node1_id, size_type node2_adj_id)
      : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_adj_id_(node2_adj_id) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    for (size_type j = 0; j != adj_list_.size(); ++j) {
      auto neighbors = adj_list_[j];
      if (i < neighbors.size()) {
        return Edge(this, j, i);
      }
      else {
        i -= neighbors.size();
      }
    }
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    auto a_neighbors = adj_list_[a.uid_];
    auto a_deg = a.degree();
    for (size_type j = 0; j < a_deg; ++j) {
      if (a_neighbors[j].first == b.uid_) {
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
    // HW0: YOUR CODE HERE
    auto a_neighbors = adj_list_[a.uid_];
    auto a_deg = a.degree();
    for (size_type j = 0; j < a_deg; ++j) {
      if (a_neighbors[j].first == b.uid_) {
        return Edge(this, a.uid_, j);
      }
    }
    edge_value_type edge_value = edge_value_type();
    std::pair<size_type, edge_value_type> b_pair (b.uid_, edge_value);
    std::pair<size_type, edge_value_type> a_pair (a.uid_, edge_value);
    adj_list_[a.uid_].push_back(b_pair);
    adj_list_[b.uid_].push_back(a_pair);
    ++num_edges_;
    return Edge(this, a.uid_, adj_list_[a.uid_].size() - 1);
  }

  /** Remove an edge from the graph
   * @param[in] a, b Nodes of the edge to be removed
   * @return the number of edges removed
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @post new has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b) == true, new num_edges() == old num_edges() - 1
   *       Else,                                new num_edges() == old num_edges()
   * 
   * Complexity: O(a.degree() + b.degree()) = O(max_degree)
   * 
   * Invalidate: edges with (node1() == a and node2() == b) or (node1() == b and node2() == a)
   * Invalidate: edge iterators representing the invalidated edges
   */
  size_type remove_edge(const Node& a, const Node& b) {
    size_type a_id = a.uid_;
    size_type b_id = b.uid_;
    for (size_type i = 0; i < adj_list_[a_id].size(); ++i) {
      if (adj_list_[a_id][i].first == b_id) {
        adj_list_[a_id][i] = adj_list_[a_id][adj_list_[a_id].size() - 1];
        adj_list_[a_id].pop_back();
      }
    }
    for (size_type i = 0; i < adj_list_[b_id].size(); ++i) {
      if (adj_list_[b_id][i].first == a_id) {
        adj_list_[b_id][i] = adj_list_[b_id][adj_list_[b_id].size() - 1];
        adj_list_[b_id].pop_back();
        num_edges_ -= 1;
        return 1;
      }
    }
    return 0;
  }

  /** Remove an edge from the graph
   * @param[in] e Edge to be removed
   * @return the number of edges removed
   * @pre e.node1() and e.node2() are distinct valid nodes of this graph
   * @post new has_edge(e.node1(), e.node2()) == false
   * @post If old has_edge(e.node1(), e.node2()) == true, new num_edges() == old num_edges() - 1
   *       Else,                                          new num_edges() == old num_edges()
   * 
   * Complexity: O(e.node1().degree() + e.node2().degree()) = O(max_degree)
   * 
   * Invalidate: edges with (node1() == e.node1() and node2() == e.node2()) or (node1() == e.node2() and node2() == e.node1())
   * Invalidate: edge iterators representing the invalidated edges
   */
  size_type remove_edge(const Edge& e) {
    const Node& a = e.node1();
    const Node& b = e.node2();
    return remove_edge(a, b);
  }

  /** Remove an edge from the graph
   * @param[in] e_it Iterator to the edge to be removed
   * @return iterator to the next edge
   * @pre *e_it is a valid edge of this graph
   * @pre *e_it.node1() and *e_it.node2() are distinct valid nodes of this graph
   * @post new has_edge(*e_it.node1(), *e_it.node2()) == false
   * @post If old has_edge(*e_it.node1(), *e_it.node2()) == true, new num_edges() == old num_edges() - 1
   *       Else,                                                  new num_edges() == old num_edges()
   * 
   * Complexity: O(*e_it.node1().degree() + *e_it.node2().degree()) = O(max_degree)
   * 
   * Invalidate: edges with (node1() == *e_it.node1() and node2() == *e_it.node2()) 
   *             or (node1() == *e_it.node2() and node2() == *e_it.node1())
   * Invalidate: edge iterators representing the invalidated edges
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    const Edge& e = *e_it;
    remove_edge(e);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    i2u_.clear();
    adj_list_.clear();
    num_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the Node pointed to by this iterator. 
     * @pre 0 <= node_id_ < graph_->num_nodes()
    */
    Node operator*() const {
      return graph_->node(node_id_);
    }

    /** Return the NodeIterator that points to the next node.
     * @pre node_id_ < graph_->num_nodes()
     * @post (*result).index() == node_id_ + 1
    */
    NodeIterator& operator++() {
      ++node_id_;
      return *this;
    }

    /** Tests whether this NodeIterator and @a node_iter are equal.
     * Equal node iterators have the same graph_ and node_id_.
    */
    bool operator==(const NodeIterator& node_iter) const {
      return node_iter.graph_ == graph_ && node_iter.node_id_ == node_id_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type node_id_;

    /** Private Constructor */
    NodeIterator(const graph_type* graph, size_type node_id)
      : graph_(const_cast<graph_type*>(graph)), node_id_(node_id) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return an iterator pointing to the first node. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return an iterator pointing to one past the last node. */
  node_iterator node_end() const {
    return NodeIterator(this, i2u_.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the Edge pointed to by this operator.
     * @pre 0 <= edge_id_ < graph_->adj_list_[node_id].size()
     * @post result.node1().uid_ == node_id_
     * @post result.node2() is the adjacent node
    */
    Edge operator*() const {
      return Edge(graph_, node_id_, edge_id_);
    }

    /** Returns the IncidentIterator that points to the next edge.
     * @pre edge_id_ < graph_->adj_list_[node_id].size()
     * @post (*result).edge_id_ == edge_id_ + 1
     */
    IncidentIterator& operator++() {
      ++edge_id_;
      return *this;
    }

    /** Tests whether this IncidentIterator and @a inc_iter are equal.
     * Equal incident iterators have the same graph, created by the same node and point to the same edge. 
    */
    bool operator==(const IncidentIterator& inc_iter) const {
      return inc_iter.graph_ == graph_ && inc_iter.node_id_ == node_id_ && inc_iter.edge_id_ == edge_id_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type node_id_;
    size_type edge_id_;

    /** Private constructor */
    IncidentIterator(const graph_type* graph, size_type node_id, size_type edge_id)
      : graph_(const_cast<graph_type*>(graph)), node_id_(node_id), edge_id_(edge_id) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the Edge pointed to by this iterator. */
    Edge operator*() const {
      return *inc_iter_;
    }

    /** Return the EdgeIterator that points to the next edge. 
     * @pre node_iter_ != graph_->node_end()
     * @post If inc_iter_ != *node_iter_.edge_end(), result.node_iter_ == node_iter_ and result.inc_iter_ == inc_iter + 1
     *       Else,                                   result.node_iter_ == node_iter_ + 1 and result.inc_iter_ == *node_iter_.edge_begin()
    */
    EdgeIterator& operator++() {
      while (node_iter_ != graph_->node_end()) {
        node_type curr_node = *node_iter_;
        while (inc_iter_ != curr_node.edge_end()) {
          ++inc_iter_;
          if (inc_iter_ == curr_node.edge_end()) {
            break;
          }
          else if (curr_node.uid_ < ((*inc_iter_).node2()).uid_) {
            return *this;
          }
        }
        ++node_iter_;
        inc_iter_ = (*node_iter_).edge_begin();
      }  
      return *this;    
    }

    /** Tests whether this EdgeIterator and @a edge_iter are equal.
     * Equal edge iterators have the same graph, NodeIterator pointing to the same node,
     * and IncidentIterator pointing to the same edge. 
    */
    bool operator==(const EdgeIterator& edge_iter) const {
      return edge_iter.graph_ == graph_ && edge_iter.node_iter_ == node_iter_ && edge_iter.inc_iter_ == inc_iter_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    node_iterator node_iter_;
    incident_iterator inc_iter_;

    /** Private Constructor */
    EdgeIterator(const graph_type* graph, node_iterator node_iter, incident_iterator inc_iter)
      : graph_(const_cast<graph_type*>(graph)), node_iter_(node_iter), inc_iter_(inc_iter) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return an iterator pointing to the first edge. */
  edge_iterator edge_begin() const {
    node_iterator node_iter = node_begin();
    incident_iterator inc_iter = (*node_iter).edge_begin();
    return EdgeIterator(this, node_iter, inc_iter);
  }

  /** Return an iterator pointing to one past the last edge. */
  edge_iterator edge_end() const {
    node_iterator node_iter = node_end();
    incident_iterator inc_iter = (*node_iter).edge_begin();
    return EdgeIterator(this, node_iter, inc_iter);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Internal type for graph nodes
  struct internal_node {
    Point p_;
    size_type idx_; 
    node_value_type v_;
    internal_node(const Point& p, const size_type idx, const node_value_type& v)
      : p_(p), idx_(idx), v_(v) {}
  };

  // Used to store internal_node for any node added, even if later removed.
  std::vector<internal_node> nodes_;
  // Store the currently "active" set of nodes.
  std::vector<size_type> i2u_;   // Indexed by node idx
  // Adjacency list of graph
  std::vector<std::vector<std::pair<size_type, edge_value_type>>> adj_list_;
  // Number of graph edges
  size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
