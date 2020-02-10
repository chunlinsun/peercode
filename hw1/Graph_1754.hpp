#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
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

  std::vector<Point> node_vector;
  std::vector<V> node_value;
  std::map<unsigned, std::vector<unsigned>> edge_map;

  std::vector<unsigned> n1, n2;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
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
      return (graph_)->node_vector[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    // HW1: YOUR CODE HERE

    /* @brief return the value of a node
     * @pre method is called by a valid Node
     */
    node_value_type& value() { return (graph_)->node_value[uid_]; };

    /* @brief return the value of a node
     * @pre method is called by a valid Node
     */
    const node_value_type& value() const { return (graph_)->node_value[uid_]; };

    /* @brief count the degree of node
     * @pre method is called by a valid Node
     * @post result >= 0 and result == deg(Node)
     */
    size_type degree() const { 
      if ((graph_)->edge_map.count(uid_) > 0)
        return (graph_)->edge_map[uid_].size();
      return 0;
    };

    /* @brief the start point of an iterator for all edge incident to Node
     * @pre Node has at least one adjacent node (one edge incident to it)
     * @post has_edge((*result), Node)
     */
    incident_iterator edge_begin() const { return IncidentIterator(graph_, uid_, 0); };
    /* @brief the start point of an iterator for all edge incident to Node
     * @pre Node has at least one adjacent node (one edge incident to it)
     * @post has_edge((*result), Node)
     */
    incident_iterator edge_end() const { return incident_iterator(graph_, uid_, degree()); };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) and (uid_ == n.uid_);
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
      return (graph_ == n.graph_) and (uid_ < n.uid_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.

    Graph* graph_;  // pointer back to graph address 
    size_type uid_; // the id of this node in the graph

    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {}

    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return (size_type)(node_vector.size());
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
  //--functionality_0
  //--For this assignment add_node is supposed to take in a value as argument
  //--START
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    node_vector.push_back(position);
    node_value.push_back(value);
    return Node(this, node_vector.size() - 1);
  }	
  //--END

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
  class Edge : private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const { return Node(graph_ptr, node1_uid); }

    /** Return the other node of this Edge */
    Node node2() const { return Node(graph_ptr, node2_uid); }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ptr == e.graph_ptr)
        if (((node1_uid == e.node1_uid) and (node2_uid == e.node2_uid)) or ((node2_uid == e.node1_uid) and (node1_uid == e.node2_uid)))
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ptr == e.graph_ptr)
        if ((node1_uid < e.node1_uid) or ((node1_uid == e.node1_uid) and (node2_uid < e.node2_uid)))
          return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.

    Graph* graph_ptr;  // pointer back to graph address 
    size_type node1_uid, node2_uid;

    Edge(const Graph* graph, const Node& a, const Node& b)
      : graph_ptr(const_cast<Graph*>(graph)), node1_uid(a.uid_), node2_uid(b.uid_) {}

    friend class Graph;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return n1.size(); }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { return Edge(this, Node(this, n1[i]), Node(this, n2[i])); }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (edge_map.count(a.uid_) > 0) {
      std::vector<unsigned> v = edge_map.at(a.uid_);
      if (std::find(v.begin(), v.end(), b.uid_) != v.end())
        return true;
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

    if (has_edge(a, b))
      return Edge(this, a, b);

    edge_map[a.uid_].push_back(b.uid_);
    edge_map[b.uid_].push_back(a.uid_);
    n1.push_back(a.uid_);
    n2.push_back(b.uid_);
    return Edge(this, a, b);        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_vector.clear();
    node_value.clear();
    n1.clear();
    n2.clear();
    edge_map.clear();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /* @brief return the node pointed by the current iterator
     * @pre a valid node iterator
     * @post 0 <= (*result).index < num_nodes()
     */
    Node operator*() const { return Node(graph_ptr, id); }
    /* @brief move the iterator one position forward
     * @pre a valid node iterator
     * @post 0 < (*result).index <= num_nodes()
     */
    NodeIterator& operator++() { id += 1; return *this; }
    /* @brief compare whether two iterators point to the same node
     * @pre a valid node iterator
     * @post 0 < (*result).index <= num_nodes()
     */
    bool operator==(const NodeIterator& n) const { return ((graph_ptr == n.graph_ptr) and (id == n.id)); }

   //private:
    friend class Graph;
    Graph* graph_ptr;
    size_type id;
    NodeIterator(const Graph* graph, size_type id_) : graph_ptr{const_cast<Graph*>(graph)}, id{id_} {}
  };

  // HW1 #2: YOUR CODE HERE
  /* @brief the start point of an iterator for all edge incident to Node
   * @post (*result).index() == 0   
   */
  node_iterator node_begin() const { return NodeIterator(this, (size_type) 0); }
  /* @brief the start point of an iterator for all edge incident to Node
   * @post (*result).index() == num_nodes
   */
  node_iterator node_end() const { return NodeIterator(this, (size_type) node_vector.size()); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator  : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    // HW1 #3: YOUR CODE HERE
    /* @brief return the edge pointed by the current iterator
     * @pre a valid source node
     * @post 0 <= (*result).index() < deg(node)
     */
    Edge operator*() const { return Edge(graph_ptr, Node(graph_ptr, source), Node(graph_ptr, (graph_ptr->edge_map)[source][neighbor_id])); }
    /* @brief move the iterator one position forward
     * @pre a valid source node 
     * @post 0 < (*result).index <= deg(node)
     */
    IncidentIterator& operator++() { neighbor_id += 1; return *this; }
    /* @brief check whether two edges are the same edge
     * @pre a valid source node 
     * @post true if and only if edge.node1 == _i_.node1 and edge.node2 == _i_.node2
     */
    bool operator==(const IncidentIterator& i) const { 
      return ((graph_ptr == i.graph_ptr) and (source == i.source) and (neighbor_id == i.neighbor_id));
    }

   private:
    friend class Graph;
    Graph* graph_ptr;
    size_type source, neighbor_id;
    IncidentIterator(const Graph* graph, size_type s, size_type n) : 
      graph_ptr{const_cast<Graph*>(graph)}, source{s}, neighbor_id{n} {}
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    // HW1 #5: YOUR CODE HERE
    /* @brief return the edge pointed by the current iterator
     * @pre a valid source node
     * @post 0 <= (*result).index() < num_edges()
     */
    Edge operator*() const { return Edge(graph_ptr, Node(graph_ptr, graph_ptr->n1[id]), Node(graph_ptr, graph_ptr->n2[id])); }
    /* @brief move the iterator one position forward
     * @pre a valid source node 
     * @post 0 < (*result).index() <= num_edges()
     */
    EdgeIterator& operator++() { id += 1; return *this; }
    /* @brief check whether two edges are the same edge
     * @pre a valid source node 
     * @post true if and only if two edges have same end nodes
     */
    bool operator==(const EdgeIterator& e) const { return ((graph_ptr == e.graph_ptr) and (id == e.id)); }

   private:
    friend class Graph;
    Graph* graph_ptr;
    size_type id;
    EdgeIterator(const Graph* graph, size_type id_) : graph_ptr{const_cast<Graph*>(graph)}, id{id_} {}
  };

  // HW1 #5: YOUR CODE HERE
  /* @brief the start point of an iterator for all edges
   * @post (*result).index() == 0   
   */
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }
  /* @brief the start point of an iterator for all edges
   * @post (*result).index() == num_edges()
   */
  edge_iterator edge_end() const { return EdgeIterator(this, n1.size()); }

 private:
};

#endif // CME212_GRAPH_HPP
