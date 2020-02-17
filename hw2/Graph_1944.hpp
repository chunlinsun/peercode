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
  using node_value_type = V;
  using edge_value_type = E;
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
    : nodes_(), edges_(), i2u_node_(), i2u_edge_(), node_node_edgemap_() {
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

    /** Returns a reference to this node's position */
    Point& position() {
      return const_cast<graph_type*>(graph_)->nodes_[uid_].p_;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].idx_;
    }

    /** Returns a reference to the value corresponding to the node
     * @return a reference to the value of the node
     */
    node_value_type& value() {
      return const_cast<graph_type*>(graph_)->nodes_[uid_].v_;
    }

    /** Returns a const reference corresponding to the node
     * @return a const reference corresponding to the node
     */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].v_;
    }

    /** Returns the number of incident edges
     * @return the number of incident edges
     */
    size_type degree() const {
      // if statement is redundant due to add_node, but is here to placate Kyle Shan :)
      if (graph_->node_node_edgemap_.count(uid_) != 0) {
        return graph_->node_node_edgemap_.at(uid_).size();
      }
      return 0;
    }

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
    return i2u_node_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] _position_ The new node's position
   * @param[in] _val_ is the value of the node
   *
   * @post Node data added to nodes_
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * @return the new node corresponding to the position and value
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    size_type newuid = nodes_.size();
    size_type newidx = i2u_node_.size();
    Node n = Node(this, newuid); // create a new node with uid based on total nodes ever created
    nodes_.push_back(nodeinfo(position, val, newidx)); // add new information to nodes_ and give next idx
    i2u_node_.push_back(newuid); // new node's idx maps to new uid
    node_node_edgemap_[newuid]; // initialize the map
    return n; // return the new node
  }

  /** Remove a node form a graph and all its incident edges.
   * Invalidates node_iterator node_end()
   * @param[in] _n_ is the node to remove
   *
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - number of incident edges to _n_
   * @return 1 if node is removed, 0 otherwise
   * Complexity: O(n.degree())
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) {
      return 0;
    }
    // remove edges incident to node:
    incident_iterator it = n.edge_begin();
    while (it != n.edge_end()) {
      remove_edge(*it);
      it = n.edge_begin();
    }
    // get necessary data to do the swap:
    size_type idx_old = nodes_[n.uid_].idx_;
    size_type uid_new = i2u_node_.back();

    // swap!
    i2u_node_[idx_old] = uid_new;
    i2u_node_.pop_back();
    nodes_[uid_new].idx_ = idx_old;

    return 1;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this && n.uid_ < nodes_.size() && n.index() < i2u_node_.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u_node_.at(i)); // return the node matching this specification
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
      return Node(graph_, n1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, n2_);
    }

    /** Returns a reference to the value corresponding to the edge
     * @return a reference to the value of the edge
     */
    edge_value_type& value() {
      return const_cast<graph_type*>(graph_)->edges_[uid_].v_;
    }

    /** Returns a const reference for the value corresponding to the edge
     * @return a const reference for the value corresponding to the edge
     */
    const edge_value_type& value() const {
      return graph_->edges_[uid_].v_;
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
   private:
    const graph_type* graph_; // Pointer to the graph
    size_type uid_; // uid of the edge
    size_type n1_;
    size_type n2_;
    friend class Graph; // Allow Graph to access Edge's private member data and functions.
    /** Private constructor */
    Edge(const graph_type* g, size_type uid, size_type n1, size_type n2)
     : graph_(const_cast<graph_type*>(g)), uid_(uid), n1_(n1), n2_(n2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return i2u_edge_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i2u_edge_[i], edges_[i2u_edge_[i]].nuid1_, edges_[i2u_edge_[i]].nuid2_);
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // check if edge already exists in graph:
    if (node_node_edgemap_.count(a.uid_) == 1 && node_node_edgemap_.at(a.uid_).count(b.uid_) == 1) {
      return Edge(this, node_node_edgemap_.at(a.uid_).at(b.uid_), a.uid_, b.uid_);
    }
    // create a new node and return otherwise
    size_type newuid = edges_.size();
    size_type newidx = i2u_edge_.size();
    Edge e = Edge(this, newuid, a.uid_, b.uid_);
    node_node_edgemap_[a.uid_][b.uid_] = newuid;
    node_node_edgemap_[b.uid_][a.uid_] = newuid;
    edges_.push_back(edgeinfo(val, newidx, a.uid_, b.uid_));
    i2u_edge_.push_back(newuid);
    return e;
  }
  /** Removes edge corresponding to the two input nodes if it's part of graph
   * @param _a_ is one node of the undirected edge
   * @param _b_ is the other node of the undirected edge
   *
   * Invalidates incident_iterator edge_end() and edge_iterator edge_end()
   * @post new num_edges() == old num_edges() - 1
   * @return 1 if edge is removed and 0 otherwise
   * Complexity: O(1)
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // Check if edge exists:
    if (!has_edge(a, b)) {
      return 0;
    }
    // find the edge uid_:
    size_type euid = node_node_edgemap_.at(a.uid_).at(b.uid_);
    return remove_edge(Edge(this, euid, a.uid_, b.uid_));

  }

  /** Removes edge _e_
   * @param _e_ is the undirected edge sought to be removed
   *
   * Invalidates incident_iterator edge_end() and edge_iterator edge_end()
   * @post new num_edges() == old num_edges() - 1
   * @return 1 if edge is removed and 0 otherwise
   * Complexity: O(1)
   */
  size_type remove_edge(const Edge& e) {
    if (!has_edge(e.node1(), e.node2())) {
      return 0;
    }
    // get necessary data to do the swap:
    size_type idx_old = edges_[e.uid_].idx_; // idx corresponding to edge to remove
    size_type a_uid = edges_[e.uid_].nuid1_; // uid corresponding to node A of edge to remove
    size_type b_uid = edges_[e.uid_].nuid2_; // uid corresponding to node B of edge to remove
    size_type uid_new = i2u_edge_.back(); // uid that idx will map to next

    // swap!
    i2u_edge_[idx_old] = uid_new;
    i2u_edge_.pop_back();
    edges_[uid_new].idx_ = idx_old;

    // update adjacency map:
    node_node_edgemap_.at(a_uid).erase(b_uid); // erase edge from adjacency map
    node_node_edgemap_.at(b_uid).erase(a_uid);
    return 1;

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    i2u_node_.clear();
    i2u_edge_.clear();
    node_node_edgemap_.clear();
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
      return graph_->node(ind_);
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

  /** Removes node corresponding to input n_it
   * @param _n_it_ is the node_iterator whose node should be removed
   * Invalidates node_iterator node_end()
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - number of incident edges to node
   * @return node_iterator corresponding to the next node after the node corresponding to _n_it_
   * Complexity: O((*n_it).degree())
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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
      return Edge(graph_, eid_, n1_, n2_);
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
      return graph_->edge(eid_);
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
    return EdgeIterator(this, num_edges());
  }

  /** Removes edge corresponding to input e_it
   * @param _e_it_ is the node_iterator whose node should be removed
   * Invalidates incident_iterator edge_end() and edge_iterator edge_end()
   * @post new num_edges() == old num_edges() - 1
   * @return edge_iterator corresponding to the next edge after the edge corresponding to _e_it_
   * Complexity: O(1)
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    remove_edge(e);
    return e_it;
  }

  struct nodeinfo {
    Point p_;
    node_value_type v_;
    size_type idx_;
    nodeinfo(Point p, node_value_type v, size_type idx) : p_(p), v_(v), idx_(idx) {}
  };

  struct edgeinfo {
    edge_value_type v_;
    size_type idx_;
    size_type nuid1_;
    size_type nuid2_;
    edgeinfo(edge_value_type v, size_type idx, size_type nuid1, size_type nuid2)
        : v_(v), idx_(idx), nuid1_(nuid1), nuid2_(nuid2) {}
  };

 private:
  std::vector<nodeinfo> nodes_;
  std::vector<edgeinfo> edges_;
  std::vector<size_type> i2u_node_;
  std::vector<size_type> i2u_edge_;
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> node_node_edgemap_;

};

#endif // CME212_HPP
