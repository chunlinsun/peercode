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

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  using node_value_type = V;
  using edge_value_type = E;

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
  using pair_type = std::pair<size_type, size_type>;
  using ptr_smaller = std::less<graph_type*>;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    // HW0: YOUR CODE HERE
    : points_(), adjacency_list_() { 
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

    /** Return this node's position. 
     * Made modifiable in HW2.
     */
    Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->points_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return the value of the node */
    node_value_type& value() { return graph_->node_values_[index_]; }

    /** Return the value of the node as constant */
    const node_value_type& value() const { return graph_->node_values_.at(index_); }

    /** Return this node's degree (number of incident edges) */
    size_type degree() const {
      return adjacency_list_[index_].size();
    }

    /** Start of the incident iterator */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, index_, 0);
    }

    /** End of the incident iterator */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, index_, graph_->adjacency_list_[index_].size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (index_ == n.index()) && (graph_ == n.graph_);
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
        return index_ < n.index_;
      }
      else {
        return ptr_smaller{}(graph_, n.graph_);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type index_;

    Node(const graph_type* g, size_type i) 
      : graph_(const_cast<graph_type*>(g)), index_(i) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return points_.size();
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
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    points_.push_back(position);
    node_values_.push_back(node_value);
    edge_values_.push_back(std::vector<edge_value_type>() );
    adjacency_list_.push_back(std::vector<size_type>() );
    return node(points_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_ == this) && (n.index_ < this->points_.size());
  }

  /** Remove a Node from this Graph
   * @brief Removes a Node and all its data from this graph if this graph contains the Node.
   *
   * @param[in] n   Node to be removed.
   * @post          If old has_node(_n_), num_nodes() = old num_nodes() - 1.
   * @post          The Point associated with old _n_ has been removed from this graph.
   * @post          No edges between old _n_ and any other node exist in this graph.
   * @post          No node values or edge values involving old _n_ are stored in this graph.
   * @return        1 if old has_node(_n_), 0 otherwise.
   *
   * The complexity of the removal algorithm is O(max_degree^2). For a graph in
   * which max_degree is O(sqrt(num_nodes)), this is O(num_nodes).
   *
   * If a node is removed, outstanding Nodes, Edges, Node iterators, Edge iterators, and
   * incident iterators will be invalidated.
   */
  size_type remove_node(const Node& n) {

    // if node is in graph, remove it and all data associated with it.
    if (has_node(n)) {

      // (1) remove all edges connected to node in adjacency list. O(max_degree^2)
      auto neighbors = adjacency_list_[n.index()];
      for (size_type i = 0; i < neighbors.size(); ++i) {
        remove_edge(n, node(neighbors[i]));
      }

      // (2) remove node from adjacency list (should be empty vector)
      //     to do this, follow Ex2 and swap node with last node. O(max_degree)
      adjacency_list_[n.index()] = adjacency_list_.back(); // overwrite n with last node
      adjacency_list_.pop_back(); // remove last element of adjacency list

      // (3) set the last node's index within the adjacency list to its new index. O(max_degree^2)
      auto neighbors_last = adjacency_list_[n.index()];
      for (size_type i = 0; i < neighbors_last.size(); ++i) {
        auto neighbors_neighbors = adjacency_list_[neighbors_last[i]];
        for (size_type j = 0; j < neighbors_neighbors.size(); ++j) {
          if (neighbors_neighbors[j] == size() - 1)
            adjacency_list_[neighbors_last[i]][j] = n.index();
        }
      }

      // (4) remove node value from values list. O(1)
      node_values_[n.index()] = node_values_.back();
      node_values_.pop_back();

      // (5) remove node from edge values list. O(max_degree)
      edge_values_[n.index()] = edge_values_.back();
      edge_values_.pop_back();

      // (6) remove node's point from points_. O(1)
      points_[n.index()] = points_.back();
      points_.pop_back();

      return 1;
    }

    // node is not in graph. nothing is removed.
    else {
      return 0;
    }
  }

  /** Remove a Node from this graph
   * @param[in,out] n_it    Node iterator pointing to Node to be removed.
   * @return                Node iterator pointing to the next Node.
   *
   * For full specifications, see remove_node() above.
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return ++n_it;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the sam
   e nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(index1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(index2_);
    }

    /** Return the length of the edge */
    double length() const {
      return norm(graph_->node(index1_).position() - graph_->node(index2_).position());
    }

    /** Return the value of the edge */
    edge_value_type& value() {
      size_type idx_of_node2 = adjacency_list_node2_index();
      return graph_->edge_values_[index1_][idx_of_node2];
    }

    /** Return the value of the edge as constant */
    const edge_value_type& value() const { 
      size_type idx_of_node2 = adjacency_list_node2_index();
      return graph_->edge_values_.at(index1_).at(idx_of_node2); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      bool nodes_are_same = (node1() == e.node1()) && (node2() == e.node2());
      bool reversed_nodes_are_same = (node2() == e.node1()) && (node1() == e.node2());
      return (graph_ == e.graph_) && (nodes_are_same || reversed_nodes_are_same);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (graph_ != e.graph_) {
        return ptr_smaller{}(graph_, e.graph_);
      }
      else {
        // in the same graph, order by nodes
        auto this_minmax = std::minmax(node1(), node2());
        auto e_minmax = std::minmax(e.node1(), e.node2());

        return this_minmax < e_minmax;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type index1_;
    size_type index2_;
    Edge(const graph_type* g, size_type i, size_type j) 
      : graph_(const_cast<graph_type*> (g)), index1_(i), index2_(j) {
    }

    /** Returns node2's index in node1's adjacency vector */
    size_type adjacency_list_node2_index() const {
      size_type idx_of_node2;
      std::vector<size_type> node1_neighbors = graph_->adjacency_list_[index1_];
      for (size_type i = 0; i < node1_neighbors.size(); ++i) {
        if (node1_neighbors[i] == index2_) {
          idx_of_node2 = i;
        }
      }
      return idx_of_node2;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type num_edges{};
    for (size_type i = 0; i < adjacency_list_.size(); ++i) {
      num_edges += adjacency_list_[i].size();
    }
    return num_edges / 2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0 && i < num_edges());

    size_type edge_count{};

    // find the i_th edge in the adjacency list
    for (size_type node = 0; node < adjacency_list_.size(); ++node) {
      for (size_type idx = 0; idx < adjacency_list_[node].size(); ++idx) {

        size_type neighbor = adjacency_list_[node][idx];

        // avoid double counting; only increment count toward the i_th edge for half of edges
        if (neighbor > node) {
          edge_count++;
          if (edge_count > i)
            return Edge(this, node, neighbor);
        }
      }
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // if either Node not in the graph, then their edge can't be in the graph
    if (!has_node(a) || !has_node(b)) {
      return false;
    }

    // look for b in the neighbors of a
    std::vector<size_type> a_neighbors = adjacency_list_[a.index()];

    for(std::size_t i = 0; i < a_neighbors.size(); i++) {
      if (a_neighbors[i] == b.index()) {
        return true;
      }
    }

    // if reached here, did not find edge between a and b
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // a and b have to be in this graph and distinct
    assert(has_node(a) && has_node(b));
    assert(!(a == b));

    if (!has_edge(a, b)) {

      // add edge to adjacency list (of both a and b)
      adjacency_list_[a.index()].push_back(b.index());
      adjacency_list_[b.index()].push_back(a.index());

      // add edge value to edge-value lists
      edge_values_[a.index()].push_back(edge_value);
      edge_values_[b.index()].push_back(edge_value);
    }
    return Edge(this, a.index(), b.index());
  }

  /** Remove an Edge from this Graph
   * @brief Removes an Edge and all its data from this graph if this graph contains the Edge.
   *
   * @param[in] a   node1() of Edge to be removed.
   * @param[in] b   node2() of Edge to be removed.
   * @post          If old has_edge(_a_, _b_), num_edges() = old num_edges() - 1.
   * @post          No edges between _a_ and _b_ exist in this graph.
   * @post          No edge values for old Edge(this, _a_, _b_) are stored in this graph.
   * @return        1 if old has_edge(_a_, _b_), 0 otherwise.
   *
   * The complexity of the removal algorithm is O(max_degree). This satisfies the
   * O(num_nodes + num_edges) requirement.
   *
   * If an edge is removed, outstanding Edges, Edge iterators, and incident iterators
   * will be invalidated.
   */
  size_type remove_edge(const Node& a, const Node& b) {

    // if edge exists, remove all data associated with it
    if (has_edge(a, b)) {

      // remove edge from adjacency list and edge values for a. this step is O(max_degree).
      auto neighbors1 = adjacency_list_[a.index()];
      for (size_type i = 0; i < neighbors1.size(); ++i) {
        if (neighbors1[i] == b.index()) {
          adjacency_list_[a.index()].erase(adjacency_list_[a.index()].begin() + i);
          edge_values_[a.index()].erase(edge_values_[a.index()].begin() + i);
        }
      }

      // remove edge from adjacency list and edge values for b. this step is O(max_degree).
      auto neighbors2 = adjacency_list_[b.index()];
      for (size_type i = 0; i < neighbors2.size(); ++i) {
        if (neighbors2[i] == a.index()) {
          adjacency_list_[b.index()].erase(adjacency_list_[b.index()].begin() + i);
          edge_values_[b.index()].erase(edge_values_[b.index()].begin() + i);
        }
      }

      return 1;
    }

    // edge does not exist in graph, nothing removed
    else {
      return 0;
    }
  }

  /** Remove an Edge from this Graph
   * @param[in] e   Edge to be removed.
   * For additional specifications, see remove_edge() above.
   */
  size_type remove_edge(const Edge& e) {
    auto n1 = e.node1();
    auto n2 = e.node2();

    return remove_edge(n1, n2);
  }

  /** Remove an Edge from this Graph
   * @param[in,out] e_it    Edge iterator pointing to Edge to be removed.
   * @return                Edge iterator pointing to the next Edge.
   * For additional specifications, see remove_edge() above.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return ++e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    points_.clear();
    node_values_.clear();
    edge_values_.clear();
    adjacency_list_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Returns node in graph that iterator is currently at */
    value_type operator*() const { return graph_->node(index_); }

    /** Increments iterator to next node */
    NodeIterator& operator++() { 
      index_++;
      return *this;
    }

    /** Tests whether this iterator is equal to another one */
    bool operator==(const NodeIterator& ni) const {
      return (graph_ == ni.graph_) && (index_ == ni.index_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;

    NodeIterator(const graph_type* g, size_type i) 
      : graph_(const_cast<graph_type*> (g)), index_(i) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns iterator at first node in graph */
  node_iterator node_begin() const { return node_iterator(this, 0); }

  /** Returns iterator 1 past last node in graph */
  node_iterator node_end() const { return node_iterator(this, points_.size()); }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    value_type operator*() const { 
      return Edge(graph_, index_, graph_->adjacency_list_[index_][neighbor_]); 
    }

    IncidentIterator& operator++() { 
      neighbor_++;
      return *this;
    }
    
    bool operator==(const IncidentIterator& ii) const {
      return (graph_ == ii.graph_) && (index_ == ii.index_) & (neighbor_ == ii.neighbor_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;  // node index of node that spawned IncidentIterator
    size_type neighbor_;  // index of neighbor in adjacency list (NOT neighbor's node index)

    IncidentIterator(const graph_type* g, size_type i, size_type j)
      : graph_(const_cast<graph_type*> (g)), index_(i), neighbor_(j) {
    }
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    value_type operator*() const { return graph_->edge(index_); }
    EdgeIterator& operator++() {
      index_++;
      return *this;
    }
    bool operator==(const EdgeIterator& ei) const {
      return (graph_ == ei.graph_) && (index_ == ei.index_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;

    EdgeIterator(const graph_type* g, size_type i)
      : graph_(const_cast<graph_type*> (g)), index_(i) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const { return edge_iterator(this, 0); }
  edge_iterator edge_end() const { return edge_iterator(this, this->num_edges()); }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> points_;
  std::vector<node_value_type> node_values_;
  std::vector<std::vector<edge_value_type>> edge_values_;
  std::vector<std::vector<size_type> > adjacency_list_;

  /** Function for debugging
   * Prints the entire adjacency list in a readable way
   * Keeping this around for future homework debugging
   */
  void print_adjacency() {
    for (size_type i = 0; i < adjacency_list_.size(); ++i) {
      std::cout << "  >>> Node " << i;
      std::cout << ": ";
      for (size_type j = 0; j < adjacency_list_[i].size(); ++j) {
        std::cout << adjacency_list_[i][j] << " ";
      }
      std::cout << " " << std::endl;
    }
  }

};


#endif // CME212_GRAPH_HPP
