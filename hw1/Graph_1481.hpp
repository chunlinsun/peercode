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
template <typename V>
class Graph {
 private:

  // Internal classes
  struct internal_node;
  struct internal_edge;
  struct incident_edges;

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
  /** Synonym for some user-specified Node value. */
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    size_ = 0;
    num_edges_ = 0;
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
      return graph_->nodes_[idx_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx_;
    }

    /** Return this node's value. */
    node_value_type& value() {
      return graph_->nodes_[idx_].val_;
    }

    /** Return this node's value that can not be changed. */
    const node_value_type& value() const {
      return graph_->nodes_[idx_].val_;
    }

    /** Return this node's degree. */
    size_type degree() const {
      return graph_->nodes_[idx_].inciedges_.size();
    }

    /** Return the start of the incident iterator. */
    incident_iterator edge_begin() const {
      return IncidentIterator {graph_, this, 0};
    }

    /** Return the end of the incident iterator. */
    incident_iterator edge_end() const {
      return IncidentIterator {graph_, this, this->degree()};
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph_ == graph_) && (n.idx_ == idx_);
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
      // If 2 nodes are in the same graph,
      // return whether the index of x < y.
      if (n.graph_ == graph_)
        return idx_ < n.idx_;

      // Otherwise, return whether the graph id of x < y.
      return graph_ < n.graph_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Node has 2 attributes:
    //    -- a unique graph id to which it belongs
    //    -- a unique node id indicating its index in the nodes list of the graph
    graph_type* graph_{};
    size_type idx_{};

    Node(const graph_type* graph, const size_type idx)
    : graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size_;
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // Here the newly added node is always the last one in the nodes list,
    // so its index is the current size_.
    // After adding a node, size_ += 1.
    std::vector<incident_edges> inciedges{};
    nodes_.emplace_back(size_, position, val, inciedges);
    size_ += 1;

    return Node {this, size_-1};        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // A node belongs to this graph when:
    //    -- the graph id matches
    //    -- the index is within the range
    return (n.graph_ == this) && (n.idx_ < size_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node {this, i};
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node {graph_, nidx1_};
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node {graph_, nidx2_};
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // 2 undirected edges are equal when:
      //    -- their graph id matches
      //    -- their edge id matches
      //    -- they have the same 2 nodes
      bool case1 = (node1() == e.node1()) && (node2() == e.node2());
      bool case2 = (node2() == e.node1()) && (node1() == e.node2());

      return (case1 || case2) && (graph_ == e.graph_) && (eidx_ == e.eidx_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    // If 2 edges are in the same graph,
    //    -- compare the index of their node1
    //    -- compare the index of their node2
    // Otherwise, compare their graph id.
    bool operator<(const Edge& e) const {
        //HW0: YOUR CODE HERE
      if (graph_ == e.graph_) {
          if (nidx1_ == e.nidx1_)
              return nidx2_ < e.nidx2_;
          return nidx1_ < e.nidx1_;
      }
      return graph_ < e.graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Edge has 4 attributes:
    //    -- a unique graph id to which it belongs
    //    -- a unique edge id for its index in the edges list
    //    -- a unique node id for node1
    //    -- a unique node id for node2
    graph_type* graph_;
    size_type eidx_;
    size_type nidx1_;
    size_type nidx2_;

    Edge(const graph_type* graph, const size_type eidx, const size_type nidx1, const size_type nidx2)
    : graph_(const_cast<graph_type*>(graph)), eidx_(eidx), nidx1_(nidx1), nidx2_(nidx2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge {this, i, edges_[i].nidx1_, edges_[i].nidx2_};
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // For all the incident edges of node a,
    // check whether the other node of this edge is node b.
    for (auto inciedge : nodes_[a.idx_].inciedges_)
      if (inciedge.nidx_ == b.idx_)
        return true;

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
//--functionality_1
//--off by one error: when you add to inciedges_ the edge id should just be
//--num_edges_ since you haven't incremented it yet.
//--START
  Edge add_edge(const Node& a, const Node& b) {
    // If edge (a, b) is not in the graph,
    //    -- add this edge to the edges list
    //    -- add this edge to both the incident edges of a and b
    //    -- the edge index of newly added edge is num_edges_-1
    //    -- after adding the edge, num_edges += 1
    // If edge (a, b) is already in the graph, return it
    uint eidx = this->num_edges();

    for (auto inciedge : nodes_[a.idx_].inciedges_)
      if (inciedge.nidx_ == b.idx_)
        eidx = inciedge.eidx_;

    if (eidx == num_edges_) {
      edges_.emplace_back(num_edges_, a.idx_, b.idx_);
      nodes_[a.idx_].inciedges_.emplace_back(num_edges_-1, b.idx_);
      nodes_[b.idx_].inciedges_.emplace_back(num_edges_-1, a.idx_);
      num_edges_ += 1;
      return Edge {this, num_edges_-1, a.idx_, b.idx_};
    } else {
      return Edge {this, eidx, a.idx_, b.idx_};
    }
  }
//--END

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Set size and number of edges to 0, clear the vector of nodes and edges.
    size_ = 0;
    num_edges_ = 0;
    nodes_.clear();
    edges_.clear();
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
//--documentation_1
//--docs should be in doxygen style
//--START
    /** Return the node specified by @a nidx_. */
    Node operator*() const {
      return graph_->node(nidx_);
    }

    /** Return the iterator pointing to the next node. */
    NodeIterator& operator++() {
        nidx_++;
      return *this;
    }

    /** Test whether this NodeIterator and @a niter are equal. */
    bool operator==(const NodeIterator& niter) const {
      return (graph_ == niter.graph_) && (nidx_ == niter.nidx_);
    }
//--END

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type * graph_;
    size_type nidx_;

    NodeIterator(const graph_type* graph, const size_type nidx)
    : graph_(const_cast<graph_type*>(graph)), nidx_(nidx) {}
  };

  /** Return the start of the node iterator. */
  node_iterator node_begin() const {
    return NodeIterator {this, 0};
  }

  /** Return the end of the node iterator. */
  node_iterator node_end() const {
    return NodeIterator {this, this->num_nodes()};
  }

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

    /** Return the incident edge specified by @a node_ and @a iidx_. */
    Edge operator*() const {
      size_type nidx = node_->index();
      return Edge {graph_, graph_->nodes_[nidx].inciedges_[iidx_].eidx_,
        nidx, graph_->nodes_[nidx].inciedges_[iidx_].nidx_};
    }

    /** Return the iterator pointing to the next incident edge. */
    IncidentIterator& operator++() {
      iidx_++;
      return *this;
    }

    /** Test whether this IncidentIterator and @a iiter are equal. */
    bool operator==(const IncidentIterator iiter) const {
      return (graph_ == iiter.graph_) && (iidx_ == iiter.iidx_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    node_type* node_;
    size_type iidx_;

    IncidentIterator(const graph_type* graph, const node_type* node, const size_type iidx)
    : graph_(const_cast<graph_type*>(graph)), node_(const_cast<node_type*>(node)), iidx_(iidx) {}
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

    /** Return the edge specified by @a eidx_. */
    Edge operator*() const {
      return graph_->edge(eidx_);
    }

    /** Return the iterator pointing to the next edge. */
    EdgeIterator& operator++() {
      eidx_++;
      return *this;
    }

    /** Test whether this EdgeIterator and @a eiter are equal. */
    bool operator==(const EdgeIterator eiter) const {
      return (graph_ == eiter.graph_) && (eidx_ == eiter.eidx_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type eidx_;

    EdgeIterator(const graph_type* graph, const size_type eidx)
    : graph_(const_cast<graph_type*>(graph)), eidx_(eidx) {}

  };

  /** Return the start of the edge iterator. */
  edge_iterator edge_begin() const {
    return EdgeIterator {this, 0};
  }

  /** Return the end of the edge iterator. */
  edge_iterator edge_end() const {
    return EdgeIterator{this, this->num_edges()};
  }

 private:
  // Store the number of nodes and edges,
  // change them only when adding node or edge.
  size_type size_;
  size_type num_edges_;

  // Struct for storing the incident edges of each node with 2 attributes:
  //    -- edge index
  //    -- the other node of this incident edge
  struct incident_edges {
    size_type eidx_;
    size_type nidx_;

    incident_edges(const size_type eidx, const size_type nidx)
    : eidx_(eidx), nidx_(nidx) {}
  };

  // Struct of internal node class with 3 attributes:
  //    -- node index
  //    -- position of this node
  //    -- incident edges of this node (using incident_edges struct)
  struct internal_node {
    size_type nidx_;
    Point position_;
    node_value_type val_;
    std::vector<incident_edges> inciedges_;

    internal_node(const size_type nidx, const Point& position,
      const node_value_type val, const std::vector<incident_edges> inciedges)
    : nidx_(nidx), position_(position), val_(val), inciedges_(inciedges) {}
  };

  // Struct of internal edge class with 3 attributes:
  //    -- edge index
  //    -- node index of node1
  //    -- node index of node2
  struct internal_edge {
    size_type eidx_;
    size_type nidx1_;
    size_type nidx2_;

    internal_edge(const size_type eidx, const size_type nidx1, const size_type nidx2)
    : eidx_(eidx), nidx1_(nidx1), nidx2_(nidx2) {}
  };

  // Vectors storing the nodes and edges in the graph.
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

};

#endif // CME212_GRAPH_HPP
