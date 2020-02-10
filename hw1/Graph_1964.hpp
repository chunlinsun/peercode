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

//--documentation_-1
//--good doxygen docs, consider putting some information in pre/post conditions
//--END

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : points_(), values_(), node_node_edgemap_(), edges_() {
    // set of edges and nodes is empty
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->points_[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return size_type(uid_);
    }

    /** Returns a reference to the value corresponding to the node
     * @return a reference to the value of the node
     */
    node_value_type& value() {
      return const_cast<graph_type*>(graph_)->values_[uid_];
    }

    /** Returns a const reference corresponding to the node
     * @return a const reference corresponding to the node
     */
    const node_value_type& value() const {
      return graph_->values[uid_];
    }
//--functionality_0
//--Doesn't handle case where degree is 0
//--START
    /** Returns the number of incident edges
     * @return the number of incident edges
     */
    size_type degree() const {
      return graph_->node_node_edgemap_.at(uid_).size();
    }
//--END

    /** Start of the incident iterator
     * @return an incident iterator of the first edge incident to node
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, graph_->node_node_edgemap_.at(uid_).begin());
    }

    /** End of the incident iterator
     * @return an incident iterator of the last edge incident to node
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, graph_->node_node_edgemap_.at(uid_).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && uid_ == n.uid_);
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
      std::less<const graph_type*> fn;
      return fn(graph_, n.graph_) || (graph_ == n.graph_ && uid_ < n.uid_);
    }

   private:
    const graph_type* graph_; // Pointer back to the graph
    size_type uid_; // ID number of the node

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /** Private constructor: */
    Node(const graph_type* g, size_type uid)
      : graph_(g), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return points_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] _position_ The new node's position
   * @param[in] _val_ is the value of the node
   *
   * @post Node position added to points_
   * @post Node value added to values_
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * @return the new node corresponding to the position and value
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    Node n = Node(this, points_.size()); // create a new node
    points_.push_back(position); // push the position of the node onto points_
    values_.push_back(val); // add the value to the value vector
    return n; // return the new node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i); // return the node matching this specification
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return n1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return n2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph_ == e.graph_ && ((n1_ == e.n1_ && n2_ == e.n2_) || (n1_ == e.n2_ && n2_ == e.n1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      std::less<const graph_type*> fn;
      return fn(graph_, e.graph_) || (graph_ == e.graph_ && uid_ < e.uid_);
    }
//--design_1
//--sizeof(Edge) > 32, consider keeping just Node indexes instead of Nodes
//--START
   private:
    const graph_type* graph_; // Pointer to the graph
    size_type uid_; // uid of the edge
    node_type n1_;
    node_type n2_;
    friend class Graph; // Allow Graph to access Edge's private member data and functions.
    /** Private constructor */
    Edge(const graph_type* g, size_type uid, node_type n1, node_type n2)
     : graph_(const_cast<graph_type*>(g)), uid_(uid), n1_(n1), n2_(n2) {
    }
//--END
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      return node_node_edgemap_.count(a.uid_) == 1 && node_node_edgemap_.at(a.uid_).count(b.uid_) == 1;
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
    // check if edge already exists in graph:
    if (node_node_edgemap_.count(a.uid_) == 1 && node_node_edgemap_.at(a.uid_).count(b.uid_) == 1) {
      return Edge(this, node_node_edgemap_.at(a.uid_).at(b.uid_), a, b);
    }
    // create a new node and return otherwise
    Edge e = Edge(this, edges_.size(), a, b);
    node_node_edgemap_[a.uid_][b.uid_] = edges_.size();
    node_node_edgemap_[b.uid_][a.uid_] = edges_.size();
    edges_.push_back(e);
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    values_.clear();
    node_node_edgemap_.clear();
    points_.clear();
    edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
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

    /** Dereferences the Node Iterator and returns the node it corresponds to
     *  @return Node that the iterator corresponds to
     *  Complexity is O(1)
     */
    Node operator*() const {
      return Node(graph_, ind_);
    }

    /** Increments the Node Iterator
     *  @return a reference to the incremented iterator
     *  Complexity is O(1)
     */
    NodeIterator& operator++() {
      ind_++;
      return *this;
    }

    /** Tests whether two iterators are the same
     *  @param[in] _x_ is the NodeIterator to compare against
     *  @returns true/false if the NodeIterators represent the same graph and node
     *  Complexity is O(1)
     */
    bool operator==(const NodeIterator& x) const {
      return graph_ == x.graph_ && ind_ == x.ind_;
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type ind_;

    // Private constructor for the Graph class:
    NodeIterator(const graph_type* graph, size_type ind) {
      graph_ = graph;
      ind_ = ind;
    }
  };

  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  node_iterator node_end() const {
    return NodeIterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    /** Dereferences the Incident Iterator and returns the edge that is currently
     *  being iterated on
     *  @return the edge, corresponding to the incident iterator, with node1() being the source node
     *  Complexity is O(1)
     */
    Edge operator*() const {
      size_type n2_ = m_->first;
      size_type eid_ = m_->second;
      return Edge(graph_, eid_, Node(graph_, n1_), Node(graph_, n2_));
    }

    /** Increments the incident iterator to go to the next edge
     *  @return a reference to the incremented iterator
     *  Complexity is O(1)
     */
    IncidentIterator& operator++() {
      ++m_;
      return *this;
    }

    /** Tests whether two incident iterators are the same
     *  @param[in] _x_ corresponds to the incident iterator being tested for equality against
     *  @return true/false if the _x_ and _this_ are the same
     *  Complexity is O(1)
     */
    bool operator==(const IncidentIterator& x) const {
      return graph_ == x.graph_ && n1_ == x.n1_ && m_ == x.m_;
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type n1_;
    std::unordered_map<size_type, size_type>::const_iterator m_;

    IncidentIterator(const graph_type* graph, size_type n1, std::unordered_map<size_type, size_type>::const_iterator m) {
      graph_ = graph;
      n1_ = n1;
      m_ = m;
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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

    /** Dereferences the EdgeIterator being operated on
     *  and returns the edge corresponding to the iterator's current position
     *  @return the edge corresponding to the iterator's current position
     *  Complexity is O(1)
     */
    Edge operator*() const {
      return graph_->edges_[eid_];
    }

    /** Increments the EdgeIterator to move onto the next edge
     *  @return a reference to the incremented EdgeIterator
     *  Complexity is O(1)
     */
    EdgeIterator& operator++() {
      eid_++;
      return *this;
    }

    /** Tests whether two EdgeIterators are the same
     *  @param[in] _x_ is an EdgeIterator to test equality against
     *  @return true if _x_ is the same as _this_
     *  Complexity is O(1)
     */
    bool operator==(const EdgeIterator& x) const {
      return graph_ == x.graph_ && eid_ == x.eid_;
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type eid_;

    EdgeIterator(const graph_type* graph, size_type eid) {
      graph_ = graph;
      eid_ = eid;
    }
  };

  /** Returns an EdgeIterator corresponding to the first edge for the graph
   *  @return an EdgeIterator corresponding to the first edge for the graph
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Returns an EdgeIterator corresponding to the last edge for the graph
   *  @return an EdgeIterator corresponding to the last edge for the graph
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.size());
  }

 private:
  std::vector<Point> points_;
  std::vector<V> values_;
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> node_node_edgemap_;
  std::vector<Edge> edges_;

};

#endif // CME212_HPP
