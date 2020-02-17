#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <iterator>
#include <utility>
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

  struct node_element {
    Point p; V v;
    node_element(Point point, V value) : p{point}, v{value} {}
  };
  struct edge_element {
    unsigned n1, n2; E v;
    edge_element(unsigned node1, unsigned node2, E value) : n1{node1}, n2{node2}, v{value} {}
  };

  std::vector<node_element> nodes;
  std::vector<edge_element> edges;
  std::map<unsigned, std::map<unsigned, unsigned>> adjacency;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
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
  using mapiterator = typename std::map<unsigned, unsigned>::iterator;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {}

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
    Node() {}

    /** Return this node's position. */
    const Point& position() const { return (graph_ptr)->nodes[nid].p; }

    Point& position() { return (graph_ptr)->nodes[nid].p; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return nid; }

    /* @brief return the value of a node
     * @pre method is called by a valid Node
     */
    node_value_type& value() { return (graph_ptr)->nodes[nid].v; }

    /* @brief return the value of a node
     * @pre method is called by a valid Node
     */
    const node_value_type& value() const { return (graph_ptr)->nodes[nid].v; }

    /* @brief count the degree of node
     * @pre method is called by a valid Node
     * @post result >= 0 and result == deg(Node)
     */
    size_type degree() const { 
      if ((graph_ptr)->adjacency.count(nid) > 0) return (graph_ptr)->adjacency.at(nid).size();
      else return 0;
    }

    /* @brief the start point of an iterator for all edge incident to Node
     * @pre Node has at least one adjacent node (one edge incident to it)
     * @post has_edge((*result), Node)
     */
    incident_iterator edge_begin() const { return incident_iterator(graph_ptr, nid, ((graph_ptr->adjacency).at(nid)).begin()); }
    /* @brief the start point of an iterator for all edge incident to Node
     * @pre Node has at least one adjacent node (one edge incident to it)
     * @post has_edge((*result), Node)
     */
    incident_iterator edge_end() const { return incident_iterator(graph_ptr, nid, ((graph_ptr->adjacency).at(nid)).end()); }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const { return (graph_ptr == n.graph_ptr) and (nid == n.nid); }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const { return (nid < n.nid); }

   private:
    // Allow Graph to access Node's private member data and functions.

    Graph* graph_ptr;  // pointer back to graph address 
    size_type nid; // the id of this node in the graph

    Node(const Graph* ptr, size_type uid) : graph_ptr{const_cast<Graph*>(ptr)}, nid{uid} {}

    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { return (size_type)(nodes.size()); }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    nodes.push_back(node_element(position, value));
    return Node(this, nodes.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const { return (this == n.graph_ptr and n.nid < num_nodes()); }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { return Node(this, i); }

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
    Node node1() const { return Node(graph_ptr, node1_id); }

    /** Return the other node of this Edge */
    Node node2() const { return Node(graph_ptr, node2_id); }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ptr == e.graph_ptr)
        if (((node1_id == e.node1_id) and (node2_id == e.node2_id)) or ((node2_id == e.node1_id) and (node1_id == e.node2_id)))
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ptr == e.graph_ptr) return (node1_id * node2_id < e.node1_id * e.node2_id);
      else return (graph_ptr < e.graph_ptr);
    }

    edge_value_type& value() { return (graph_ptr->edges).at((graph_ptr->adjacency).at(node1_id).at(node2_id)).v; }

    const edge_value_type& value() const { return (graph_ptr->edges).at((graph_ptr->adjacency).at(node1_id).at(node2_id)).v; }

   private:
    // Allow Graph to access Edge's private member data and functions.

    Graph* graph_ptr;  // pointer back to graph address 
    size_type node1_id, node2_id;

    Edge(const Graph* graph, const Node& a, const Node& b)
      : graph_ptr(const_cast<Graph*>(graph)), node1_id(a.nid), node2_id(b.nid) {}

    friend class Graph;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return edges.size(); }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { return Edge(this, Node(this, edges[i].n1), Node(this, edges[i].n2)); }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const { return (adjacency.count(a.nid) > 0 and adjacency.at(a.nid).count(b.nid) > 0); }

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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    if (has_edge(a, b)) return Edge(this, a, b);
    adjacency[a.nid][b.nid] = num_edges(); adjacency[b.nid][a.nid] = num_edges();
    edges.push_back(edge_element(a.nid, b.nid, value));
    return Edge(this, a, b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    adjacency.clear();
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

   private:
    Graph* graph_ptr;
    size_type id;

    NodeIterator(const Graph* graph, size_type id_) : graph_ptr{const_cast<Graph*>(graph)}, id{id_} {}

    friend class Graph;
  };

  /* @brief the start point of an iterator for all edge incident to Node
   * @post (*result).index() == 0   
   */
  node_iterator node_begin() const { return node_iterator(this, 0); }
  /* @brief the start point of an iterator for all edge incident to Node
   * @post (*result).index() == num_nodes
   */
  node_iterator node_end() const { return node_iterator(this, nodes.size()); }

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

    /* @brief return the edge pointed by the current iterator
     * @pre a valid source node
     * @post 0 <= (*result).index() < deg(node)
     */
    Edge operator*() const { return Edge(graph_ptr, Node(graph_ptr, source), Node(graph_ptr, (*it).first)); }
    /* @brief move the iterator one position forward
     * @pre a valid source node 
     * @post 0 < (*result).index <= deg(node)
     */
    IncidentIterator& operator++() { ++it; return *this; }
    /* @brief check whether two edges are the same edge
     * @pre a valid source node 
     * @post true if and only if edge.node1 == _i_.node1 and edge.node2 == _i_.node2
     */
    bool operator==(const IncidentIterator& i) const { 
      return ((graph_ptr == i.graph_ptr) and (source == i.source) and (it == i.it));
    }

   private:
    Graph* graph_ptr;
    size_type source;
    mapiterator it;

    IncidentIterator(const Graph* graph, size_type s, mapiterator it_) : 
      graph_ptr{const_cast<Graph*>(graph)}, source{s}, it{it_} {}

    friend class Graph;
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
    Edge operator*() const { 
      return Edge(graph_ptr, Node(graph_ptr, (graph_ptr->edges)[id].n1), Node(graph_ptr, (graph_ptr->edges)[id].n2));
    }
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
    Graph* graph_ptr;
    size_type id;

    EdgeIterator(const Graph* graph, size_type id_) : graph_ptr{const_cast<Graph*>(graph)}, id{id_} {}

    friend class Graph;
  };

  /* @brief the start point of an iterator for all edges
   * @post (*result).index() == 0   
   */
  edge_iterator edge_begin() const { return edge_iterator(this, 0); }
  /* @brief the start point of an iterator for all edges
   * @post (*result).index() == num_edges()
   */
  edge_iterator edge_end() const { return edge_iterator(this, edges.size()); }

  /** Remove a node from the graph, and return 1 if successful, 0 otherwise.
   * @pre @a n is a node of this graph
   * @return 0 or 1 indicating whether the removal is successful
   * @post has_node(@a n) == false
   * @post For any valid node @a m, has_edge(@a n, @a m) == false
   * @post For any valid node @a m, has_edge(@a m, @a n) == false
   *
   * Can invalidate node indices -- in other words, old node(@a i) might not
   * equal new node(@a i).
   * Can invalidate node iterators -- in other words, old node_iterator(@ it)
   * might point to a new node (@a n).
   *
   * Complexity: No more than O(num_nodes())
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) return 0;

    Node last_node = Node(this, num_nodes() - 1);

    if (n.degree() > 0)
      for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
        size_type neighbor = (*e).node2().nid;
        size_type id = adjacency[n.nid][neighbor];

        edges[id] = edges.back(); edges.pop_back();

        adjacency[edges[id].n1][edges[id].n2] = id;
        adjacency[edges[id].n2][edges[id].n1] = id;
        adjacency[neighbor].erase(n.nid);
      }

    if (last_node.degree() > 0)
      for (auto e = last_node.edge_begin(); e != last_node.edge_end(); ++e) {
        size_type neighbor = (*e).node2().nid;
        size_type id = adjacency[last_node.nid][neighbor];

        if (edges[id].n1 == last_node.nid) edges[id].n1 = n.nid;
        if (edges[id].n2 == last_node.nid) edges[id].n2 = n.nid;

        adjacency[neighbor][n.nid] = adjacency[neighbor][last_node.nid];
        adjacency[neighbor].erase(last_node.nid);
      }
    adjacency[n.nid] = adjacency[last_node.nid]; adjacency[last_node.nid].clear();

    nodes[n.nid] = nodes.back(); nodes.pop_back();
    return 1;
  }

  /** Remove a node from the graph, and return the iterator to another node or invalid node.
   * @pre @a n is a node iterator of this graph
   * @return the iterator to another node if there exists one, or an invalid iterator otherwise.
   * @post has_node(*(@a n)) == false
   * @post For any valid node @a m, has_edge((*@a n), @a m) == false 
   * @post For any valid node @a m, has_edge(@a m, (*@a n)) == false
   *
   * Can invalidate node indices -- in other words, old node(@a i) might not
   * equal new node(@a i). Does not invalidate iterator objects.
   * Can invalidate node iterators -- in other words, old node_iterator(@ it)
   * might point to a new node (@a n).
   *
   * Complexity: No more than O(num_nodes())
   */
  node_iterator remove_node(node_iterator n_it) { remove_node(*n_it); return n_it; }

  /** Remove an edge from the graph, and return 1 if successful or 0 otherwise.
   * @pre @a a, @a b are valid nodes of this graph
   * @return 0 or 1 indicating whether the edge removal is successful
   * @post has_edge(@a a, @a b) == false
   * @post has_edge(@a b, @a a) == false
   *
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   * Can invalidate edge iterators -- in other words, old edge_iterator(@ it)
   * might point to a new edge (@a e).
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Node& a, const Node& b) {
    for (auto it = edge_begin(); it != edge_end(); ++it)
      if (((*it).node1() == a and (*it).node2() == b) or ((*it).node1() == b and (*it).node2() == a)) {

        size_type last_n1 = edges.back().n1; size_type last_n2 = edges.back().n2;
        adjacency[last_n1][last_n2] = it.id; adjacency[last_n2][last_n1] = it.id;

        edges[it.id] = edges.back(); edges.pop_back();
        adjacency[a.nid].erase(b.nid); adjacency[b.nid].erase(a.nid);
        return 1;
      }
    return 0;
  }

  /** Remove an edge from the graph, and return 1 if successful or 0 otherwise.
   * @pre @a e is a valid edge of this graph
   * @return 0 or 1 indicating whether the edge removal is successful
   * @post has_edge(@a e.node1(), @a e.node2()) == false
   * @post has_edge(@a e.node1(), @a e.node2()) == false
   *
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   * Can invalidate edge iterators -- in other words, old edge_iterator(@ it)
   * might point to a new edge (@a e).
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Edge& e) { return remove_edge(e.node1(), e.node2()); }

  /** Remove an edge from the graph, and return 1 if successful or 0 otherwise.
   * @pre @a e_it is a valid edge iterator of this graph
   * @return 0 or 1 indicating whether the edge removal is successful
   * @post has_edge(@a (*e_it).node1(), @a (*e_it).node2()) == false
   * @post has_edge(@a (*e_it).node1(), @a (*e_it).node2()) == false
   *
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   * Can invalidate edge iterators -- in other words, old edge_iterator(@ e_it)
   * might point to a new edge (@a e_it).
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it) { remove_edge(*e_it); return e_it; }

 private:
  void graph_test() {
    for (auto e = edge_begin(); e != edge_end(); ++e) {
      Edge e1 = *e;
      Node n1 = e1.node1(); Node n2 = e1.node2();
      assert(has_edge(n1, n2)); assert(has_edge(n2, n1));
      size_type id = adjacency.at(n1.nid).at(n2.nid);
      assert(adjacency.at(n2.nid).at(n1.nid) == id);
      assert(edges[id].n1 == n1.nid or edges[id].n2 == n1.nid);
      assert(edges[id].n1 == n2.nid or edges[id].n2 == n2.nid);
      assert(edges[id].n1 != edges[id].n2);
    }
  }

};

#endif // CME212_GRAPH_HPP
