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

  // Internal classes
  struct internal_node;
  struct internal_edge;

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
  /** Synonym for some user-specified Edge value. */
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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

    /** Return the reference of this node's position. */
    Point& position() {
      assert(valid());
      return graph_->nodes_[idx_].position_;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[idx_].position_;
    }

    /** Return this node's index, a number in the range [0, i2u_.size()). */
    size_type index() const {
      return graph_->nodes_[idx_].nidx_;
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
      return graph_->adj_mat_[idx_].size();
    }

    /** Return the start of the incident iterator. */
    incident_iterator edge_begin() const {
      return IncidentIterator {graph_, idx_, 0};
    }

    /** Return the end of the incident iterator. */
    incident_iterator edge_end() const {
      return IncidentIterator {graph_, idx_, this->degree()};
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
    graph_type* graph_;
    size_type idx_;

    /** Funtcion to check whether a node is valid */
    bool valid() const {
      return (0 <= idx_) && (idx_ < graph_->nodes_.size()) && (graph_->nodes_[idx_].nidx_ >= 0) &&
             (graph_->nodes_[idx_].nidx_ < graph_->i2u_.size()) && (graph_->i2u_[graph_->nodes_[idx_].nidx_] == idx_);
    }

    Node(const graph_type* graph, const size_type idx)
    : graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // Here the newly added node is always the last one in the nodes list,
    // so its index is the current size_.
    // After adding a node, size_ += 1.
    size_type uidx = nodes_.size();
    nodes_.push_back(internal_node(num_nodes(), position, val));
    i2u_.push_back(uidx);
    std::vector<internal_edge> inciedges{};
    adj_mat_.push_back(inciedges);

    return Node {this, uidx};        // Invalid node
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
    return (n.graph_ == this) && (i2u_[n.index()] < nodes_.size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node {this, i2u_[i]};
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

    /** Return this edge's value. */
    edge_value_type& value() {
      size_type min_idx = std::min(nidx1_, nidx2_);
      size_type max_idx = std::max(nidx1_, nidx2_);
      size_type iidx = 0;

      for (size_type i = 0; i < graph_->adj_mat_[min_idx].size(); i++)
        if (graph_->adj_mat_[min_idx][i].nidx_ == max_idx)
          iidx = i;

      return graph_->adj_mat_[min_idx][iidx].val_;;
    }

    /** Return this edge's value that can not be changed. */
    const edge_value_type& value() const {
      return value();
    }

    /** Return the length of this edge. */
    double length() const {
      return norm(node1().position() - node2().position());
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

      return (case1 || case2) && (graph_ == e.graph_);
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
    size_type nidx1_;
    size_type nidx2_;

    Edge(const graph_type* graph, const size_type nidx1, const size_type nidx2)
    : graph_(const_cast<graph_type*>(graph)), nidx1_(nidx1), nidx2_(nidx2) {
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
    size_type eidx = 0;

    for (auto it = edge_begin(); it != this->edge_end(); ++it)
      if (eidx == i)
        return *it;
      else
        eidx++;

    return *edge_begin();
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
    for (auto inciedge : adj_mat_[a.idx_])
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // If edge (a, b) is not in the graph,
    //    -- add this edge to the edge list
    //    -- add this edge to both the incident edges of a and b
    //    -- the edge index of newly added edge is num_edges_-1
    //    -- after adding the edge, num_edges += 1
    // If edge (a, b) is already in the graph, return it
    if (has_edge(a, b))
      return Edge {this, a.idx_, b.idx_};
    else {
      size_type min_idx = std::min(a.idx_, b.idx_);
      size_type max_idx = std::max(a.idx_, b.idx_);

      adj_mat_[min_idx].emplace_back(max_idx, val);
      adj_mat_[max_idx].emplace_back(min_idx, val);
      num_edges_++;
      return Edge {this, a.idx_, b.idx_};
    }
  }

  /** Given two nodes, remove edge of these two nodes from the graph,
   *  
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return 1 if this edge has been removed, else return 0
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Complexity: O(d) < O(m), where d is the maximum degree of nodes in the graph,
   *             since we only go through the incident edges of @a a and @a b
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // If edge (a, b) is in the graph
    if (has_edge(a, b)) {
      // Remove (a, b) from the incident edge list of a
      for (size_type i = 0; i < a.degree(); i++)
      	if (adj_mat_[a.idx_][i].nidx_ == b.idx_) {
          adj_mat_[a.idx_].erase(adj_mat_[a.idx_].begin() + i);
          break;
      	}

      // Remove (a, b) from the incident edge list of a
      for (size_type i = 0; i < b.degree(); i++)
        if (adj_mat_[b.idx_][i].nidx_ == a.idx_) {
          adj_mat_[b.idx_].erase(adj_mat_[b.idx_].begin() + i);
          break;
        }

  	  // Update @a num_edges_
  	  num_edges_ -= 1;

  	  return 1;
    }

    return 0;
  }

  /** Remove a edge by using its nodes as input */
  size_type remove_edge(const Edge& e) {
  	return remove_edge(e.node1(), e.node2());
  }

  /** Remove a edge in the iterator by using its reference as input */
  edge_iterator remove_edge(edge_iterator e_it) {
  	remove_edge(*e_it);
  	return e_it;
  }

  /** Remove a node in the graph, return 1 if success, else return 0
   * 
   * @post has_node(@a n) == false
   * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
   *       Else,                  new num_nodes() == old num_nodes().
   *
   * Complexity: O(n), since we will go through the index of remaing nodes
   */
  size_type remove_node(const Node& n) {
  	if (has_node(n)) {
  	  // Remove all incidente edges of @a n
      while (n.degree() > 0) {
      	remove_edge(n, Node(this, adj_mat_[n.idx_][0].nidx_));
      }
      // Remove @a n from the node list
      i2u_.erase(i2u_.begin() + n.index());

      // Update the index of nodes after @a n in the node list
      for (size_type i = n.index(); i < num_nodes(); i++)
      	nodes_[i2u_[i]].nidx_ = i;

      return 1;
  	}

  	return 0;
  }

  /** Remove a node in the iterator by using its reference as input */
  node_iterator remove_node(node_iterator n_it) {
  	remove_node(*n_it);
  	return n_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Set number of edges to 0, clear the vector of nodes and edges.
    num_edges_ = 0;
    nodes_.clear();
    i2u_.clear();
    adj_mat_.clear();
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

   private:
    friend class Graph;
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
      return Edge {graph_, nidx_, graph_->adj_mat_[nidx_][iidx_].nidx_};
    }

    /** Return the iterator pointing to the next incident edge. */
    IncidentIterator& operator++() {
      iidx_++;
      return *this;
    }

    /** Test whether this IncidentIterator and @a iiter are equal. */
    bool operator==(const IncidentIterator iiter) const {
      return (graph_ == iiter.graph_) && (nidx_ == iiter.nidx_) && (iidx_ == iiter.iidx_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type nidx_;
    size_type iidx_;

    IncidentIterator(const graph_type* graph, const size_type nidx, const size_type iidx)
    : graph_(const_cast<graph_type*>(graph)), nidx_(nidx), iidx_(iidx) {}
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
      return (*iiter_);
    }

    /** Return the iterator pointing to the next edge. */
    EdgeIterator& operator++() {
      while (niter_ != graph_-> node_end()) {
        while (iiter_ != (*niter_).edge_end()) {
          ++iiter_;

          if (iiter_ == (*niter_).edge_end())
            break;
          else if ((*niter_).idx_ <= (*iiter_).node2().idx_)
            return *this;
        }

        ++niter_;

        if (niter_ != graph_->node_end())
          iiter_ = (*niter_).edge_begin();
      }

      return *this;
    }

    /** Test whether this EdgeIterator and @a eiter are equal. */
    bool operator==(const EdgeIterator eiter) const {
      return (graph_ == eiter.graph_) && (niter_ == eiter.niter_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    NodeIterator niter_;
    IncidentIterator iiter_;

    EdgeIterator(const graph_type* graph, const NodeIterator niter, const IncidentIterator iiter)
    : graph_(const_cast<graph_type*>(graph)), niter_(niter), iiter_(iiter) {}

  };

  /** Return the start of the edge iterator. */
  edge_iterator edge_begin() const {
    return EdgeIterator {this, node_begin(), (*node_begin()).edge_begin()};
  }

  /** Return the end of the edge iterator. */
  edge_iterator edge_end() const {
    return EdgeIterator {this, node_end(), (*node_end()).edge_begin()};
  }

 private:
  // Store the number of edges
  size_type num_edges_;

  // Struct of internal node class with 3 attributes:
  //    -- node index
  //    -- position of this node
  //    -- ndde value
  struct internal_node {
    size_type nidx_;
    Point position_;
    node_value_type val_;

    internal_node(const size_type nidx, const Point& position, const node_value_type val)
    : nidx_(nidx), position_(position), val_(val) {}
  };

  // Struct of internal edge class with 2 attributes:
  //    -- node index of the other node in this edge
  //    -- edge value
  struct internal_edge {
    size_type nidx_;
    edge_value_type val_;

    internal_edge(const size_type nidx, const edge_value_type val)
    : nidx_(nidx), val_(val) {}
  };

  // Vectors storing the nodes and edges in the graph.
  std::vector<internal_node> nodes_;

  // Vectors storing the index of current nodes
  std::vector<size_type> i2u_;

  // Adjacency matrix
  std::vector<std::vector<internal_edge>> adj_mat_;
};

#endif // CME212_GRAPH_HPP
