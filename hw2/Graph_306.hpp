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
  struct intl_node; // internal struct to store nodes
  struct intl_edge; // internal struct to store edges

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

  typedef V node_value_type;
  typedef E edge_value_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
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
    Node() {
      // HW0: YOUR CODE HERE // FINISH
    }

    /** Return this node's modifiable position. */
    Point& position() {
      return graph_->vec_node_[nidx_].coord_;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE // FINISH
      return graph_->vec_node_[nidx_].coord_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE // FINISH
      return graph_->vec_node_[nidx_].idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return this node's value by reference, which can be changed. */
    node_value_type& value() {
      return graph_->vec_node_[nidx_].val_;
    }

    /** Return this node's value by reference, which cannot be changed. */
    const node_value_type& value() const {
      return graph_->vec_node_[nidx_].val_;
    }

    /** Return this node's degree which is the number of incident edges. */
    size_type degree() const {
      return graph_->adj_node_[nidx_].size();
    }

    /** Return the iterator to the begin of incident edges of current node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, nidx_, 0);
    }
    
    /** Return the iterator to the end of incident edges of current node. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, nidx_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE // FINISH
      return graph_ == n.graph_ && nidx_ == n.nidx_;
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
      // HW0: YOUR CODE HERE // FINISH

      // global order is to compare the graph pointers first
      // and node indices secondly
      return graph_ == n.graph_ ?
             nidx_ < n.nidx_ : graph_ < n.graph_;
    }

    bool valid() const {
      return nidx_ >= 0 && nidx_ < graph_->vec_node_.size()
             && graph_->vec_node_[nidx_].idx_ < graph_->ni2u_.size()
             && graph_->ni2u_[graph_->vec_node_[nidx_].idx_] == nidx_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE // FINISH
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    graph_type* graph_; // pointer to the graph that has this node
    size_type nidx_;     // node index

    Node(const graph_type* graph, size_type node_idx)
    : graph_(const_cast<graph_type*>(graph)), nidx_(node_idx) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE // FINISH
    return ni2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position THe new node's position
   * @param[in] node_value_type THe new node's value type
   * @post new num_node() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations
   */
  Node add_node(const Point& position, const node_value_type& nval = node_value_type()) { 
    vec_node_.push_back(intl_node(ni2u_.size(), nval, position));
    ni2u_.push_back(adj_node_.size());
    adj_node_.push_back(std::vector<std::pair<size_type, size_type>>());
    return Node(this, adj_node_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE // FINISH
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE // FINISH
    return Node(this, ni2u_[i]);
  }

  /** Remove the given node.
   *  @param[in] n The node to be deleted.
   *  return 1 if the node is deleted and 0 otherwise.
   *
   *  @pre has(n) == true.
   *  @post has(n) == false.
   *  @post If deleted, new num_node_ == old num_node_ - 1;
   *        Else,       new num_node_ == old num_node_.
   *  @post The function must not invalidate other indices in the node vector.
   *
   *  COmplexity: at most O(num_nodes()).
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) return 0;
    size_type del_idx = n.index();

    // remove nodes connected with this node from the adjaceny list
    while (n.degree() > 0) {
      remove_edge(*(n.edge_begin()));
    }

    // swap deleted node internal id with the last one in ni2u
    vec_node_[ni2u_.back()].idx_ = del_idx;
    *(ni2u_.begin() + del_idx) = ni2u_.back();
    ni2u_.pop_back();

    return 1;
  }

  /** Remove the node pointed by the given iterator.
   *  @param[in] n_it the iterator pointing to the node to be deleted
   *  return a node_iterator object pointing to the next node or
   *         to the end of nodes if the deleted node is the last one.
   *
   *  @pre has(*n_it) == true.
   *  @post If deleted, new num_node_ == old num_node_ - 1;
   *        Else,       new num_node_ == old num_node_.
   *  @post The function must not invalidate other indices in the node vector.
   *
   *  COmplexity: at most O(num_nodes()).
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE // FINISH
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE // FINISH
      return Node(graph_, node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE // FINISH
      return Node(graph_, node2_idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE // FINISH
      return graph_ == e.graph_ &&
             ((node1_idx_ == e.node1_idx_ && node2_idx_ == e.node2_idx_) ||
              (node1_idx_ == e.node2_idx_ && node2_idx_ == e.node1_idx_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE // FINISH

      // global order is to compare graph points first
      // and edge index secondly
      return graph_ == e.graph_? eidx_ < e.eidx_ : graph_ < e.graph_;
    }

    /** Return the modifiable value of this edge. */
    edge_value_type& value() {
      return graph_->vec_edge_[eidx_].edge_val_;
    }

    /** Return the value of this edge. */
    const edge_value_type& value() const {
      return graph_->vec_edge_[eidx_].edge_val_;
    }

    /** Return the initial length of this edge. */
    double length() const {
      return norm(node1().position() - node2().position());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    graph_type* graph_;   // pointer to the graph that has this node
    size_type eidx_;      // edge index
    size_type node1_idx_;
    size_type node2_idx_;

    // constructor for a valid Edge object
    Edge(const graph_type* graph, size_type eidx, size_type node1, size_type node2)
    : graph_(const_cast<graph_type*>(graph)), eidx_(eidx),
      node1_idx_(node1), node2_idx_(node2) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE // FINISH
    return ei2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE // FINISH
    return Edge(this, ei2u_[i], vec_edge_[ei2u_[i]].node1_, vec_edge_[ei2u_[i]].node2_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE // FINISH
    for (size_type i = 0; i < adj_node_[a.nidx_].size(); i++) {
      if (adj_node_[a.nidx_][i].first == b.nidx_) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& eval = edge_value_type()) {
    // HW0: YOUR CODE HERE // FINISH

    // check if the edge exists
    if (has_edge(a, b)) {
      for (auto e : vec_edge_) {
        if ((e.node1_ == a.nidx_ && e.node2_ == b.nidx_) ||
            (e.node1_ == b.nidx_ && e.node2_ == a.nidx_))
          return Edge(this, ei2u_[e.eidx_], e.node1_, e.node2_);
          //return Edge(this, vec_edge_.size() - 1, a.nidx_, b.nidx_);
      }
    }

    // store new edge, new index, and update the adjaceny list
    ei2u_.push_back(vec_edge_.size());
    adj_node_[a.nidx_].push_back(std::make_pair(b.nidx_, vec_edge_.size()));
    adj_node_[b.nidx_].push_back(std::make_pair(a.nidx_, vec_edge_.size()));
    vec_edge_.push_back(intl_edge(ei2u_.size()-1, a.nidx_, b.nidx_, eval));

    return Edge(this, vec_edge_.size() - 1, a.nidx_, b.nidx_);
  }

  /** Remove the given edge from this graph.
   *  @param[in] n1, n2 two nodes of the edge to be deleted.
   *  @return 1 if the edge is deleted and 0 otherwise.
   *
   *  @pre has_edge(n1, n2) == true.
   *  @post has_edge(n1, n2) == false.
   *  @post If deleted, new @a num_edge_ == old @a num_edge_ - 1
   *        Else        new @a num_edge_ == old @a num_edge_
   *  @post The delete operation must not invalidate edge indices in @a vec_edge_.
   *
   *  Complexity: at most O(num_nodes() + num_edges()).
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!has_edge(n1, n2)) return 0;

    size_type cur_idx = 0;      // Edge id
    size_type cur_intl_idx = 0; // internal edge id

    // remove n2 from n1's adjaceny list
    for (auto it = adj_node_[n1.nidx_].begin(); it != adj_node_[n1.nidx_].end(); ++it) {
      if (it->first == n2.nidx_) {
        cur_idx = it->second;
        cur_intl_idx = vec_edge_[cur_idx].eidx_;
        adj_node_[n1.nidx_].erase(it);
        break;
      }
    }

    // remove n1 from n2's adjaceny list
    for (auto it = adj_node_[n2.nidx_].begin(); it != adj_node_[n2.nidx_].end(); ++it) {
      if (it->first == n1.nidx_) {
        adj_node_[n2.nidx_].erase(it);
        break;
      }
    }

    // swap the deleted edge with the last edge in ei2u
    vec_edge_[ei2u_.back()].eidx_ = cur_intl_idx;
    *(ei2u_.begin() + cur_intl_idx) = ei2u_.back();
    ei2u_.pop_back();

    return 1;
  }

  /** Remove the given edge from this graph.
   *  @param[in] e an edge to be deleted.
   *  @return 1 if the edge is deleted and 0 otherwise.
   *
   *  @pre has_edge(e.node1(), e.node2()) == true.
   *  @post has_edge(e.node1(), e.node2()) == false.
   *  @post If deleted, new @a num_edge_ == old @a num_edge_ - 1
   *        Else        new @a num_edge_ == old @a num_edge_
   *  @post The delete operation must not invalidate edge indices in @a vec_edge_.
   *
   *  Complexity: at most O(num_nodes() + num_edges()).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove the edge indicated by the given iterator from this graph.
   *  @param[in, out] e_it an iterator pointing to the edge for deleting.
   *  @return a edge_iterator type iterator pointing to the next edge or
   *          to the end of edges if the deleted edge is at the end.
   *
   *  @pre e_it points to a valid edge. has_edge(e.node1(), e.node2()) == true.
   *  @post has_edge(e.node1(), e.node2()) == false.
   *  @post If deleted, new @a num_edge_ == old @a num_edge_ - 1
   *        Else        new @a num_edge_ == old @a num_edge_
   *  @post The delete operation must not invalidate edge indices in @a vec_edge_.
   *
   *  Complexity: at most O(num_nodes() + num_edges()).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE //FINISH
    vec_node_.clear();
    vec_edge_.clear();
    adj_node_.clear();
    ni2u_.clear();
    ei2u_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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

    /** Return node by calling node() in class Graph with the node index.
     *  Complexity: O(1).
     */
    Node operator*() const {
        return graph_->node(nitr_idx_);
    }

    /** Return the iterator to the next node or the end of the node vector.
     *  Complexity: O(1).
     */
    NodeIterator& operator++() {
        nitr_idx_++;
        return *this;
    }

    /**
     * @brief Function to check whether two iterators are the same.
     * @param[in] nitr The iterator pointing to the other node in the comparison.
     * @return A boolean value indicating whether two iterators are the same.
     *
     * @post Two iterators are returned as equal if their parent graph pointers
     *       as well as the iterator indices are the same.
     *
     * Complexity: O(1).
     */
    bool operator==(const NodeIterator& nitr) const {
        return graph_ == nitr.graph_ && nitr_idx_ == nitr.nitr_idx_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;  // pointer to the parent graph
    size_type nitr_idx_; // index of the current node being pointed to

    /**
     * @brief Constructor to create an iterator with the graph pointer and
     *        the index of the starting node.
     * @param[in] graph The pointer to the graph that is iterated.
     * @param[in] idx The index of the node which is the start of iteration.
     * @return A NodeIterator type iterator with passed parameters.
     */
    NodeIterator(const graph_type* graph, const size_type idx) :
                 graph_(const_cast<graph_type*>(graph)), nitr_idx_(idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return the iterator to the start of the node vector. */
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }

  /** Return the iterator to the end of the node vector. */
  node_iterator node_end() const {
      return NodeIterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  // class IncidentIterator: private totally_ordered<IncidentIterator> {
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    /** Return node by calling Edge() constructor with two node indices.
     *  Complexity: O(1).
     */
    Edge operator*() const {
      auto cur = graph_->adj_node_[nidx_][cur_edge_];
      return Edge(graph_, cur.second, nidx_, cur.first);
    }

    /** Return the iterator to the next edge or the end of the edge vector.
     *  The iteration range is the incident edges per node.
     *  Complexity: O(1).
     */
    IncidentIterator& operator++() {
      cur_edge_++;
      return *this;
    }
    
    /**
     * @brief Function to check whether two iterators are the same.
     * @param[in] iditr The iterator pointing to the other node in the comparison.
     * @return A boolean value indicating whether two iterators are the same.
     *
     * @post Two iterators are returned as equal if their parent graph pointers
     *       as well as node indices and current edges are the same.
     *
     * Complexity: O(1).
     */
    bool operator==(const IncidentIterator& iditr) const {
      return graph_ == iditr.graph_ &&
             nidx_ == iditr.nidx_ && cur_edge_ == iditr.cur_edge_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_; // pointer to the parent graph
    size_type nidx_;    // the index of the node incident by edges
    size_type cur_edge_;// current edge being iterated

    /**
     * @brief Constructor to create an iterator with the graph pointer and
     *        the index of the starting node.
     * @param[in] graph The pointer to the graph that is iterated.
     * @param[in] idx The index of the node which is the start of iteration.
     * @param[in] edge The starting incident edge being iterated.
     * @return A IncidentIterator type iterator with passed parameters.
     */
    IncidentIterator(const graph_type* graph, size_type idx, size_type edge)
    : graph_(const_cast<graph_type*>(graph)), nidx_(idx), cur_edge_(edge) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /** Return node by calling edge() in class Graph with two node indices.
     *  Complexity: O(1).
     */
    Edge operator*() const {
      return graph_->edge(eidx_);
    }

    /** Return the iterator to the next edge or the end of the edge vector.
     *  The iteration range is the edge vector per graph. So every edge is 
     *  iterated once.
     *  Complexity: O(1). 
     */
    EdgeIterator& operator++() {
      eidx_++;
      return *this;
    }

    /**
     * @brief Function to check whether two iterators are the same.
     * @param[in] eitr The iterator pointing to the other node in the comparison.
     * @return A boolean value indicating whether two iterators are the same.
     *
     * @post Two iterators are returned as equal if their parent graph pointers
     *       as well as the edge indices are the same.
     *
     * Complexity: O(1).
     */
    bool operator==(const EdgeIterator& eitr) const {
      return graph_ == eitr.graph_ && eidx_ == eitr.eidx_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_; // pointer to the parent graph
    size_type eidx_;    // current edge being iterated

    /**
     * @brief Constructor to create an iterator with the graph pointer and
     *        the index of the starting edge.
     * @param[in] graph The pointer to the graph that is iterated.
     * @param[in] idx The index of the edge which is the start of iteration.
     * @return A EdgeIterator type iterator with passed parameters.
     */
    EdgeIterator(const graph_type* graph, size_type idx) 
    : graph_(const_cast<graph_type*>(graph)), eidx_(idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

   /** Return the iterator to the start of the edge vector. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  
   /** Return the iterator to the end of the edge vector. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, ei2u_.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector<intl_node> vec_node_; // vector of internal nodes
  std::vector<intl_edge> vec_edge_; // vector of internal edges

  std::vector<size_type> ni2u_; // internal node id to Node id
  std::vector<size_type> ei2u_; // internal edge id to Edge id

  // adjacency list of nodes
  std::vector<std::vector<std::pair<size_type, size_type>>> adj_node_; 

  // internal struct to store nodes
  struct intl_node {
    size_type idx_;
    node_value_type val_;
    Point coord_;

    intl_node(size_type idx, node_value_type val, const Point& coord)
    : idx_(idx), val_(val), coord_(coord) {}
  };

  // internal struct to store edges
  struct intl_edge {
    size_type eidx_;
    size_type node1_;
    size_type node2_;
    edge_value_type edge_val_;

    intl_edge(size_type idx, size_type node1, size_type node2, edge_value_type val)
    : eidx_(idx), node1_(node1), node2_(node2), edge_val_(val) {}
  };

};

#endif // CME212_GRAPH_HPP