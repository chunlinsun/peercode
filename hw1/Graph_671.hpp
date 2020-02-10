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

  using node_value_type = V;
  
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
    : nodes_(), num_nodes_(0), adj_list_(), num_edges_(0) {}

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
      return graph_->nodes_[node_id_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return a reference to this node's value. */
    node_value_type& value() {
      return graph_->nodes_[node_id_].value;
    }
    
    /** Return a const reference to this node's value. */
    const node_value_type& value() const {
      return graph_->nodes_[node_id_].value;
    }

    /** Return this node's degree, i.e. the number of edges incident on this node. */
    size_type degree() const {
      return graph_->adj_list_[node_id_].size();
    }

    /** Return an iterator pointing to the first edge incident on this node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, node_id_, 0);
    }

    /** Return an iterator pointing to one past the last edge incident on this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, node_id_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (n.graph_ == graph_) && (n.node_id_ == node_id_);
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
        return node_id_ < n.node_id_;
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
    size_type node_id_;
    /** Private Constructor */
    Node(const Graph* graph, size_type node_id)
      : graph_(const_cast<Graph*>(graph)), node_id_(node_id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return num_nodes_;
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
    size_type new_node_id = num_nodes_;
    ++num_nodes_;
    internal_node new_internal_node = internal_node(position, new_node_id, value);
    nodes_.push_back(new_internal_node);
    adj_list_.push_back(std::vector<size_type>());
    return Node(this, new_node_id);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_ == this) && (n.node_id_ < num_nodes_);
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
      return Node(graph_, node2_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_) && 
             ((node1_id_ == e.node1_id_ && node2_id_ == e.node2_id_) || 
              (node1_id_ == e.node2_id_ && node2_id_ == e.node1_id_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_) {
        if (node1_id_ == e.node1_id_) {
          return node2_id_ < e.node2_id_;
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
    // Node2's unique identification number
    size_type node2_id_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type node1_id, size_type node2_id)
      : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_id_(node2_id) {}
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
      std::vector<size_type> neighbors = adj_list_[j];
      if (i < neighbors.size()) {
        return Edge(this, j, neighbors[i]);
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
    std::vector<size_type> a_neighbors = adj_list_[a.node_id_];
    for (auto j = a_neighbors.begin(); j != a_neighbors.end(); ++j) {
      if (*j == b.node_id_) {
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
    if (has_edge(a, b) == false) {
      adj_list_[a.node_id_].push_back(b.node_id_);
      adj_list_[b.node_id_].push_back(a.node_id_);
      ++num_edges_;
    }
    return Edge(this, a.node_id_, b.node_id_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    num_nodes_ = 0;
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
    return NodeIterator(this, num_nodes_);
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
     * @pre 0 <= edge_id_ < graph_->adj_list_[node1_id].size()
     * @post result.node1().index() == node_id_
     * @post result.node2() is the adjacent node
    */
    Edge operator*() const {
      size_type node1_id = node_id_;
      size_type node2_id = graph_->adj_list_[node1_id][edge_id_];
      return Edge(graph_, node1_id, node2_id);
    }

    /** Returns the IncidentIterator that points to the next edge.
     * @pre edge_id_ < graph_->adj_list_[node1_id].size()
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
          else if (curr_node.index() < ((*inc_iter_).node2()).index()) {
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
    Point position;
    size_type node_id; 
    node_value_type value;
    internal_node(const Point& position_, const size_type node_id_, const node_value_type& value_)
      : position(position_), node_id(node_id_), value(value_) {}
  };

  // Vector of graph nodes
  std::vector<internal_node> nodes_;
  // Number of graph nodes
  size_type num_nodes_;
  // Adjacency list of graph
  //--design_0
  //--Consider using different structures like maps to speed up your implementation.
  //--START
  std::vector<std::vector<size_type>> adj_list_;
  //--END
  // Number of graph edges
  size_type num_edges_;

};

//--functionality_0
//--Great job!
//--END
#endif // CME212_GRAPH_HPP
