#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * Note: discussed graph design (node, edges) with Lauren Pendo
 */

#include <algorithm>
#include <functional>
#include <map>
#include <set>
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
class Graph {
 private:

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    : nodes_(), node_size_(0), edges_(), node_map_(), edge_size_(0) {
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
  class Node {
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
      return (Point&) graph_->nodes_[idx_];
      assert(false);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (idx_ == n.idx_ && graph_ == n.graph_) return true;
      return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * compare graphs using std::less template
     * return some comparison if graphs are not equal, otherwise return
     * the less than indices
     */
    bool operator<(const Node& n) const {
      if (graph_ != n.graph_) return std::less<Graph*>()(graph_, n.graph_);
      else return idx_ < n.idx_; 
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    // Node members a pointer to the actual graph and the index of the node
    Graph* graph_;
    size_type idx_;

    // create a valid constructor of the node
    Node(const Graph* graph, size_type idx) 
      : graph_(const_cast<Graph*>(graph)), idx_(idx) {
      }

    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_size_;
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
  Node add_node(const Point& position) {
    // get the size of the node, which will be the current index
    size_type idx = size();
    // add the point position to nodes
    nodes_.push_back(position);
    // increase node size
    node_size_++;
    return Node(this, idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */

  bool has_node(const Node& n) const {
    // check if the current graph is the same graph
    if (this == n.graph_) return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // return node with pointer to a graph and index
    assert(i < size());
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
<<<<<<< HEAD
      return edge_idx_ == e.edge_idx_ and graph_ == e.graph_;
=======
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return false;
>>>>>>> CME212-2020/master
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * compare graphs using std::less template
     */
    bool operator<(const Edge& e) const {
<<<<<<< HEAD
      if (graph_ != e.graph_) return std::less<Graph*>()(graph_, e.graph_);
      else return edge_idx_ < e.edge_idx_;
=======
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return false;
>>>>>>> CME212-2020/master
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // graph member types include graph pointer and indices to nodes and edges
    Graph* graph_;
    size_type edge_idx_;
    size_type node1_idx_;
    size_type node2_idx_;

    // Edge constructor
    Edge(const Graph* graph, 
      size_type edge_idx, 
      size_type node1_idx, 
      size_type node2_idx) 
      : graph_(const_cast<Graph*>(graph)), 
      edge_idx_(edge_idx), 
      node1_idx_(node1_idx), 
      node2_idx_(node2_idx) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // check if index exists 
    assert(i < edge_size_);
    // translate the set type into a vector for easy indexing
    std::vector<size_type> e(edges_[i].begin(), edges_[i].end());
    return Edge(this, i, e[0], e[1]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check if a vector exists at all in this map index
    if (node_map_.find(a.index()) == node_map_.end()) return false;

    // check if node b exists in the vector of node a, otherwise return true
    std::vector<std::vector<size_type>> np = node_map_.at(a.index());
    for (std::size_t i = 0; i < np.size(); i++) {
      if (np[i][0] == b.index()) return true;
    }

    // otherwise return false
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
    // if there exists an edge, return an edge
    // node_map_ is a map from node a to a vector of <node b, edge> vectors
    // iterate through node_map_ at a to get indices of already existing edges
    if (has_edge(a, b)) {
      size_type idx = 0;
      std::vector<std::vector<size_type>> np = node_map_.at(a.index());
      for (std::size_t i = 0; i < np.size(); i++) {
        if (np[i][0] == b.index()) idx = i;
      }
      return Edge(this, idx, a.index(), b.index());
    }
    // create a new edge
    else {
      // get the index of the edge
      size_type idx = num_edges();

      // if a vector at node a already exists, add a new vector
      if(node_map_.find(a.index()) == node_map_.end()) {
        std::vector<size_type> node_edge {b.index(), idx};
        std::vector<std::vector<size_type>> node_edge_set = {node_edge};
        node_map_[a.index()] = node_edge_set;
      }
      // otherwise push a vector of <node b, current edge index> to node a
      else {
        std::vector<size_type> node_edge {b.index(), idx};
        node_map_.at(a.index()).push_back(node_edge);
      }

      // if a vector at node b already exists, add a new vector
      // this repetition is to allow for easy look up of nodes from 
      // a -> b and b -> a
      if(node_map_.find(b.index()) == node_map_.end()) {
        std::vector<size_type> node_edge {a.index(), idx};
        std::vector<std::vector<size_type>> node_edge_set = {node_edge};
        node_map_[b.index()] = node_edge_set;
      }
      // otherwise push a node and current edge index 
      else {
        std::vector<size_type> node_edge {a.index(), idx};
        node_map_.at(b.index()).push_back(node_edge);
      }


      // push back a and b nodes as a set
      std::set<size_type> nodes = {a.index(), b.index()};
      edges_.push_back(nodes);
      // increase the size of the edge
      edge_size_++;
      return Edge(this, idx, a.index(), b.index());
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   * Use stl functions to clear everything in graph.
   */
  void clear() {
    node_size_ = 0;
    edge_size_ = 0;
    nodes_.clear();
    edges_.clear();
    node_map_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
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

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // get a vector of points that are nodes
  std::vector<Point> nodes_;

  // node size
  size_type node_size_;

  // get a mapping of node sets to its indices
  std::vector<std::set<size_type>> edges_;

  // create a map between a node and a vector of <node 2, edge> pairs that
  // connect between a and b
  std::map<size_type, std::vector<std::vector<size_type>>> node_map_;

  // edge size
  size_type edge_size_;

};

#endif // CME212_GRAPH_HPP
