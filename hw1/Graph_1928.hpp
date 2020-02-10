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
  struct internal_node;
  struct internal_edge;
  struct incident_edge;
  // Proxy nodes and edges in this graph.
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_; 

  unsigned size_ = 0; // Number of nodes.
  unsigned num_edges_ = 0;  // Number of edges.

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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[uid_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // assert(uid_ == graph_->nodes_[uid_].uid);
      return graph_->nodes_[uid_].uid;
    }

    // HW1: YOUR CODE HERE
    /** return the value on the node.
     */
    node_value_type& value() {
      return graph_->nodes_[uid_].node_val;
    }

    /** return the value on the node.
     */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].node_val;
    }

    /** return the number of incident edges.
     */
    size_type degree() const {
      return graph_->nodes_[uid_].inc_edges.size();
    }

    /** Start of the incident iterator.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, this, 0);
    }

    /** End of incident iterator.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, this, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((graph_==n.graph_) && (uid_ == n.uid_)) {
        return true;
      } else {
        return false;
      }
      
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
      if (uid_ < n.uid_) {
        return true;
      } else {
        return false;
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
    size_type uid_;
    // private constructor of Node.
    Node(const graph_type* graph, size_type uid): graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& node_val=node_value_type()) {
    // create a empty vector for holding the incident edges later.
    std::vector<incident_edge> inc_edges{};
    nodes_.push_back(internal_node(position, size_, node_val, inc_edges));
    ++size_;
    return Node(this, size_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    size_type idx = n.index();
    if (idx >= size())
      return false;
    else {
      Node node = Node(this, idx);
      if (node == n)
        return true;
      else
        return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < size()) {
      return Node(this, i);  
    }
    return Node();
    
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1() == e.node1() && node2() == e.node2()) || (node1() == e.node2() && node2() == e.node1());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (edge_uid_ < e.edge_uid_)
        return true;
      else
        return false;

    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // Private Edge constructor
    graph_type* graph_;
    size_type node1_uid_;
    size_type node2_uid_;
    size_type edge_uid_;
    Edge(const graph_type* graph, size_type uid1, size_type uid2, size_type euid)
      :graph_(const_cast<graph_type*>(graph)), node1_uid_(uid1), node2_uid_(uid2), edge_uid_(euid) {
    };

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
    // HW0: YOUR CODE HERE
    if (i < num_edges_) 
      return Edge(this, edges_[i].n1_uid, edges_[i].n2_uid, i); 
    else
      return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (unsigned i = 0; i < nodes_[a.uid_].inc_edges.size(); ++i) {
      if (nodes_[a.uid_].inc_edges[i].id2 == b.uid_) {
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
    if (!(a == b)) {
		  if (!has_edge(a, b)) {
			  edges_.push_back(internal_edge(a.uid_, b.uid_, num_edges_));
        // add an incident edge to both node a and b.
			  nodes_[a.uid_].inc_edges.push_back(incident_edge(num_edges_, b.uid_));
			  nodes_[b.uid_].inc_edges.push_back(incident_edge(num_edges_, a.uid_));
			  num_edges_++;
			  return Edge(this, a.uid_, b.uid_, num_edges_-1);
		  }
		  else {
			  return Edge(this, a.uid_, b.uid_, num_edges_-1); // return the current edge
      }
	  }
	else
	 	return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
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
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    // HW1 #2: YOUR CODE HERE
    // Returns the node that corresponds to the iter_idx_.
    Node operator*() const {
      return graph_->node(iter_idx_);
    }

    // Increment the index of the ndoe that this iterator is at.
    NodeIterator& operator++() {
      ++iter_idx_;
      return *this;
    }

    /** Test whether this node iterator and @a node_iterator are equal.
     *
     * Equal node iterators have the same graph and the same index.
     */
    bool operator==(const NodeIterator& node_iterator) const {
      // True if the two iterators are on the same graph and have the same idx.
      return (node_iterator.graph_ == graph_ && node_iterator.iter_idx_ == iter_idx_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type iter_idx_;
    // HW1 #2: YOUR CODE HERE
    // Private constructor for NodeIterator
    NodeIterator(const graph_type* graph, unsigned idx): graph_(const_cast<graph_type*>(graph)), iter_idx_(idx) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return an iterator pointing at the start of the nodes.
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return an iterator pointing at the end of the nodes.
   */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
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
    /** Return the edge in the graph
     * The first node of the edge is @a node_.
     * The second node of the edge is identitified by @a iter_idx_;
     */
    Edge operator*() const {
      size_type node_id = node_->uid_;
      return Edge(
        graph_,
        node_id, graph_->nodes_[node_id].inc_edges[iter_idx_].id2,
        graph_->nodes_[node_id].inc_edges[iter_idx_].edge_uid
      );
    }

    // Move on to the next inc edge by incrementing @a iter_idx_.
    IncidentIterator& operator++() {
      ++iter_idx_;
      return *this;
    }

    /** Test whether this iterator and @a inc_iterator are equal.
     *
     * Equal incident iterators have the same graph, the same index,
     * and the same node.
     */
    bool operator==(const IncidentIterator& inc_iterator) const {
      return ((inc_iterator.graph_ == graph_) && (inc_iterator.iter_idx_ == iter_idx_) 
               && (inc_iterator.node_ == node_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    node_type* node_;
    size_type iter_idx_;
    // Private constructor of IncidentIterator
    IncidentIterator(const graph_type* graph, const node_type* node, size_type iter_idx):
      graph_(const_cast<graph_type*>(graph)),node_(const_cast<node_type*>(node)),iter_idx_(iter_idx) {}
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
    /** Returns the edge in the graph that is at index @a iter_idx_.
     */
    Edge operator*() const {
      return graph_->edge(iter_idx_);
    }

    /** Move on to the next edge by incrementing @a iter_idx_.
     */
    EdgeIterator& operator++() {
      ++iter_idx_;
      return *this;
    }

    /** Test whether this iterator and @a edge_iter are equal.
     *
     * Equal edge iterators have the same graph and the same index.
     */
    bool operator==(const EdgeIterator& edge_iter) const {
      return ((edge_iter.graph_ == graph_) && (edge_iter.iter_idx_ == iter_idx_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type iter_idx_;
    // Private construtor of EdgeIterator
    EdgeIterator(const graph_type* graph, size_type iter_idx): 
      graph_(const_cast<graph_type*>(graph)), iter_idx_(iter_idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return an iterator pointing at the start of the edges.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  /** Return an iterator pointing at the end of the edges.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct incident_edge {
    size_type edge_uid;    // The unique identifcation for this edge
    size_type id2;        // The unique identifcation for the other node.
    incident_edge(size_type edge_uid_, size_type id2_) :edge_uid(edge_uid_), id2(id2_) {};
  };

  struct internal_node {
    Point p;   // The point held by this node
    size_type uid;      // The unique identifcation for this node
    node_value_type node_val;   // The value that this node stores.
    std::vector<incident_edge> inc_edges;   // The edges incident to this node. 
    internal_node(Point point, size_type id, node_value_type node_value, std::vector<incident_edge> inc_e): p(point), uid(id), node_val(node_value), inc_edges(inc_e) {};
  };

  struct internal_edge {
    size_type n1_uid;     // The unique identifcation for one of the nodes
    size_type n2_uid;     // The unique identifcation for one of the nodes
    size_type uid;    // The unique identifcation for this edge
    internal_edge(size_type a, size_type b, size_type id): n1_uid(a), n2_uid(b), uid(id) {};
  };
};
//--functionality_0
//--Good job!
//--END
#endif // CME212_GRAPH_HPP
