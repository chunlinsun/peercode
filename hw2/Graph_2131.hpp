#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * 
 * The code from HW0 has been copied and modified from Graph_1152.hpp from
 * the peer code
 */

#include <algorithm>
#include <vector>
#include <map>
#include <functional>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// spring constant to be added to all edges in the graph
static constexpr double spring_const = 100;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {

  /** Private predeclaration of the struct internal_node */
  // NEEDS TO BE PRIVATE SINCE IT IS DECLARED AS PRIVATE IN GRAPH CLASS
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

  /** Templatizing the Graph class for nodes (HW1) */
  using node_value_type = V;

  /** Templatizing the Graph class for edges (HW2) */
  using edge_value_type = E;

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
    Point & position() const {
      return const_cast<internal_node&>(node_ptr->nodes_vec.at(node_uid)).coords;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_uid;
    }

    /** @brief Returns a node_value
    * 
    * @return the value for a given node
    *
    * This function returns the node_value variable for each node
    * This may be an int like the distance, as measured in 
    * number of edges, from a given node 
    */
    node_value_type& value() {
      return const_cast<internal_node&>(node_ptr->nodes_vec.at(node_uid)).property;
    }

    /** @brief Returns a const node_value
    *
    * @return the value for a given node
    *
    * This function returns the node_value variable for each node
    * This may be an int like the distance, as measured in 
    * number of edges, from a given node 
    */
    const node_value_type& value() const {
      return node_ptr->nodes_vec.at(node_uid).property;
    }

    /** @brief Returns the degree of the node
    * 
    * Complexity: O(number of incident iterators)
    *
    * This function returns the degree for each node
    * The degree of a node is the number of edges connected to the node
    */
    size_type degree() const {
        size_type degree = 0;

        for(auto it = edge_begin(); it != edge_end(); ++it) {
            degree++;
        }

        return degree;
    }

    /** @brief Returns a begin iterator to an edge that connects to the given node
    * 
    * This function returns the begin incident iterator for the node
    * This can be used to iterate over all the edges that connect to the node
    * Such as when finding the distance of the node from a root node
    * (see shortest_path.cpp)
    */
    incident_iterator edge_begin() const {
      return incident_iterator(node_uid, node_ptr, node_ptr->nested_map.at(node_uid).begin());
    }


    /** @brief Returns an end iterator to an edge that connects to the given node
    * 
    * This function returns the end incident iterator for the node
    * This can be used to iterate over all the edges that connect to the node
    * Such as when finding the distance of the node from a root node
    * (see shortest_path.cpp)
    */
    incident_iterator edge_end() const {
      return incident_iterator(node_uid, node_ptr, node_ptr->nested_map.at(node_uid).end());
    }


    /** @brief Test whether this node and _n_ are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node & n) const {
      return node_ptr == n.node_ptr && node_uid == n.node_uid;
    }

    /** @brief Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node & n) const {
      return std::less<const Graph *>{}(node_ptr, n.node_ptr) || 
             (node_ptr == n.node_ptr && node_uid < n.node_uid);
    }

    /** @brief constructor used by Graph and NodeIterator to create Node objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] index of the node in owning Graph
     */
    Node(const graph_type* graph, const size_type index)
      : node_ptr(graph),
        node_uid(index) {}

   private:
    friend class Graph;         /** Allow Graph to access Node's private member data and functions. */
    const graph_type* node_ptr; /** Pointer to the containing Graph */
    size_type node_uid;         /** Index of the node in the containing Graph */
  };

  /** @brief Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_i2u.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** @brief Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val The new node's node_value_type
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point & position, const node_value_type& val = node_value_type()) {
    size_type index = nodes_vec.size();
    nodes_vec.emplace_back(position, val, true, index);
    node_i2u.push_back(index);

    return Node(this, index);
  }

  /** @brief Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node & n) const {
    return n.node_ptr == this && n.index() < num_nodes();
  }

  /** @brief Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
  }

  /** @brief Remove node _n_ from the Graph
  * @param[in] n    The Graph::Node to be removed
  *
  * @post The new num_nodes = old num_nodes - 1
  *
  * Complexity : O(num_edges)
  * Remove the node _n_ from the Graph as well as the edges incident to the node _n_
  */
  size_type remove_node(const Node& n) {
    if(has_node(n)) {
        size_type n_uid = n.index();
        auto it = std::find(node_i2u.begin(), node_i2u.end(), n_uid);

        size_type temp = *it;
        *it = node_i2u.back();
        node_i2u.back() = temp;

        node_i2u.pop_back();

        nodes_vec.at(n_uid).valid = false;
        nodes_vec.at(*it).node_idx = n_uid;

        // removing incident edges
        for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
            remove_edge(*it);
        }

        return 1; // remove was successful
    }

    else {
        return 0; // remove was unsuccessful
    }
  }

  /** @brief Remove node that @ n_it points to 
  * @param[in] n_it     Graph::NodeIterator that points to the node
  *                     that needs to be removed
  * 
  * @post The new num_nodes = old num_nodes - 1
  * Complexity : O(num_edges)
  * Remove the node that @ n_it points to, as well as the edges
  * incident to the node
  */
  node_iterator remove_node(node_iterator n_it) {
    node_type n = *n_it;
    remove_node(n);

    return n_it;
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const { 
      return Node(edge_ptr, first_node);
    }

    /** Return the other node of this Edge */
    Node node2() const {  
      return Node(edge_ptr, second_node);
    }

    /** @brief Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge & e) const
    { 
      // Graph class maintains the invariant first_node <= second_node.
      // Therefore we don't have to check the other combination.
      return edge_ptr == e.edge_ptr &&
             first_node == e.first_node &&
             second_node == e.second_node;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge & e) const
    {  
      return std::less<const Graph *>{}(edge_ptr, e.edge_ptr) ||
             (edge_ptr == e.edge_ptr && (first_node < e.first_node || 
                                      (first_node == e.first_node && second_node < e.second_node)));
    }

    /** Return this edge's index, a number in the range [0, graph_size). */
    size_type index() {  
      return edge_uid;
    }

    /** Return this edge's length. */
    double length() const {  
      return norm(node1().position() - node2().position());  
    }

    /** Return this edge's value. */
    edge_value_type& value() {  
        return const_cast<internal_edge&>(edge_ptr->edges_vec.at(edge_uid)).property;
    }

    /** Return this edge's value. */
    const edge_value_type& value() const {  
        return edge_ptr->edges_vec.at(edge_uid).property;
    }

    /** constructor used by Graph, IncidentIterator, EdgeIterator to create Edge objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] node1 index of first node in owning Graph
     * @param[in] node2 index of second node in owning Graph
     */
    Edge(const graph_type* graph, const size_type index,
         const size_type node1,   const size_type node2)
    : edge_ptr(graph),
      edge_uid(index),
      first_node(node1),
      second_node(node2) {}

   private:
    friend class Graph;         // Allow Graph to access Edge's private member data and functions.
    const graph_type* edge_ptr; /** Pointer to the containing Graph */
    size_type edge_uid;         /** Index of the edge */
    size_type first_node;       /** Index of first node */
    size_type second_node;      /** Index of second node */
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edge_i2u.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    internal_edge p = edges_vec.at(i);
    return Edge(this, i, p.node1_index, p.node2_index);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(log num_edges)
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(nested_map.empty()) {
        return false;
    }

    std::map<size_type, size_type> mp = nested_map.at(a.index());
    return mp.find(b.index()) != mp.end();
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
  Edge add_edge(const Node & a, const Node & b, edge_value_type val = edge_value_type()) {
    const size_type index = edges_vec.size();
    std::pair<size_type, size_type> nodes = std::minmax(a.index(), b.index());

    if(!has_edge(a, b)) { // if this edge doesn't exist
      size_type count = edges_vec.size();

      val.set_values(norm(a.position() - b.position()), spring_const);
      edges_vec.emplace_back(nodes.first, nodes.second, val, true, count);
      edge_i2u.push_back(count);

      nested_map.insert(std::make_pair(nodes.first, std::map<size_type, size_type>()));
      nested_map.at(nodes.first).insert(std::make_pair(nodes.second, index));

      nested_map.insert(std::make_pair(nodes.second, std::map<size_type, size_type>()));
      nested_map.at(nodes.second).insert(std::make_pair(nodes.first, index));
    }

    else {
      nested_map.at(nodes.first).insert(std::make_pair(nodes.second, index));
      nested_map.at(nodes.second).insert(std::make_pair(nodes.first, index));
    }
    
    return Edge(this, index, nodes.first, nodes.second);
  }

  /** @brief Remove the edge connecting nodes @ a and @ b
  *
  * @return If has_edge(@ a, @ b), 1.
            Else,                  0.
  * @post The edge connecting @ a and @ b is invalidated
  * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
  *       Else,                        new num_edges() == old num_edges().
  *
  * Complexity : O(log num_edges)
  */
  size_type remove_edge(const Node& a, const Node& b) {
    if(has_edge(a, b)) { // if the edge does exist, then remove it
        size_type e_uid = nested_map.at(a.index()).at(b.index());

        auto it = std::find(edge_i2u.begin(), edge_i2u.end(), e_uid);  

        size_type temp = *it;  
        *it = edge_i2u.back();  // *it is the uid of the edge being swapped with
        edge_i2u.back() = temp; // swapping the uid's of the edges (to be removed and the last one)

        edge_i2u.pop_back(); // remove the uid of the edge to be removed

        edges_vec.at(e_uid).valid = false; // make the edge invalid
        nested_map.at(a.index()).erase(b.index());
        nested_map.at(b.index()).erase(a.index()); 

        edges_vec.at(*it).edge_idx = e_uid;

        return 1;
    }

    else { // the edge doesn't exist (is invalidated), so nothing to do
        return 0;
    }
  }

  /** @brief Remove the edge @ e
  *
  * @return If has_edge(@ e.node1(), @ e.node2()), 1.
            Else,                                  0.
  * @post If old has_edge(@ e.node1(), @ e.node2()), new num_edges() == old num_edges() - 1.
  *       Else,                                      new num_edges() == old num_edges().
  *
  * @post The edge @ e is invalidated
  * Complexity : O(log num_edges)
  */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** @brief Remove the edge @ *e_it
  *
  * @return If has_edge(@ (*e_it).node1(), @ (*e_it).node2()), 1.
            Else,                                  0.
  * @post If old has_edge(@ (*e_it).node1(), @ (*e_it).node2()), new num_edges() == old num_edges() - 1.
  *       Else,                                                  new num_edges() == old num_edges().
  * @post The edge @ *e_it is invalidated
  *
  * Complexity : O(log num_edges)
  */
  edge_iterator remove_edge(edge_iterator& e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_vec.clear();
    nodes_vec.clear();
    edge_i2u.clear();
    nested_map.clear();
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
      // no code added here
    }

    /** @brief Returns the node that the node iterator points to
    * 
    * This function returns the Node that the node iterator points to
    * This can be used to obtain information about the node of interest
    */
    node_type operator*() const {
        return Node(nodeiter_ptr, *node_iter);
    }

    /** @brief Increments the node iterator
    * 
    * This function returns the Node iterator that points to the next Node
    * in the graph. This can be useful when one needs to iterate over all 
    * the nodes in a given graph.
    */
    node_iterator& operator++() {
        node_iter++;
        return *this;
    }

    /** Tests the equality of the current instance of the node iterator to
    * another instance
    * 
    * Two node iterators are the same when the node_iter's are equal
    * the index of the node, being pointed to, are also equal
    */
    bool operator==(const NodeIterator& NodeIt) const {
      return node_iter == NodeIt.node_iter;
    }

    /** Constructor of a node iterator
    * 
    * allows for the construction of an node iterator
    */
    NodeIterator(std::vector<size_type>::const_iterator p, const graph_type* g)
    : node_iter(p),
      nodeiter_ptr(g) {}

   private:
    friend class Graph;
    std::vector<size_type>::const_iterator node_iter;
    const graph_type*                      nodeiter_ptr;
  };

  /** Returns a node iterator that points to the beginning node within a graph
  * 
  * This can be used to iterate over all the nodes in a graph.
  */
  node_iterator node_begin() const {
    return NodeIterator(node_i2u.begin(), this);
  }

  /** Returns a node iterator that points to the end node within a graph
  * 
  * This can be used to iterate over all the nodes in a graph.
  */
  node_iterator node_end() const {
    return NodeIterator(node_i2u.end(), this);
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
    IncidentIterator() {}

    /** Returns an edge that connects to the node of interest
    * 
    * When dereferencing the iterator, the user obtains an Edge object
    * @post The first node of the edge is always the node of interest
    */
    edge_type operator*() {
      return Edge(incid_ptr, map_iter->second, first_idx, map_iter->first);
    }

    /** Increments the incident iterator
    * 
    * This results in the incident iterator pointing to the next edge that
    * connects to the node of interest
    */
    incident_iterator& operator++() {
      ++map_iter;
      return *this;
    }

    /** Tests equality of two incident iterators
    * 
    * Two incident iterators are equal when the two iterators point
    * to the same edge 
    */
    bool operator==(const incident_iterator& it) const {
      return map_iter == it.map_iter;
    }

    /** Constructor for the incident iterator
    * 
    */
    IncidentIterator(size_type index,
                     const graph_type* g,
                     std::map<size_type, size_type>::const_iterator p)
    : first_idx(index),
      incid_ptr(g),
      map_iter(p) {}

   private:
    friend class Graph;
    size_type                                      first_idx;
    const graph_type*                              incid_ptr;
    std::map<size_type, size_type>::const_iterator map_iter;
    // 1st node idx = first_idx
    // 2nd node idx = map_iter->first
    // edge idx     = map_iter->second
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
    EdgeIterator() {}

    /** Returns an edge that the edge iterator points to
    * 
    * @post Dereferencing the end iterator returns an invalid Edge
    */
    Edge operator*() const {  
      size_type index = *edgeiter; // obtain edge's uid
      internal_edge e = edgeiter_ptr->edges_vec.at(index);

      return Edge(edgeiter_ptr, index, e.node1_index, e.node2_index);
    }

    /** Increments the edge iterator
    * 
    * The iterator now points to the next edge in the graph.
    */
    edge_iterator& operator++() {  
      edgeiter++;
      return *this;
    }

    /** Tests the equality of two edge iterators
    * 
    * Two edge iterators are equal when they point to the same graph,
    * the same edge 
    */
    bool operator==(const EdgeIterator& it) const {  
      bool pred2 = edgeiter_ptr == it.edgeiter_ptr;
      bool pred3 = edgeiter     == it.edgeiter;

      return pred2 && pred3;
    }

   private:
    friend class Graph;
    const graph_type*                      edgeiter_ptr;
    std::vector<size_type>::const_iterator edgeiter;

    /** Private constructor 
    * 
    */
    EdgeIterator(std::vector<size_type>::const_iterator it)
    : edgeiter(it) {}
  };

  /** Returns an edge iterator that points to the first edge of the graph
  * 
  * This can be used to iterate over all edges within a graph
  */
  edge_iterator edge_begin() const {  
    EdgeIterator it(edge_i2u.begin());
    it.edgeiter_ptr = this;

    return it;
  }

  /** Returns an edge iterator that points to one past the last edge of the graph
  * 
  * This can be used to iterate over all edges within a graph
  */
  edge_iterator edge_end() const {  
    EdgeIterator it(edge_i2u.end());
    it.edgeiter_ptr = this;
     
    return it;
  }

 private:

  // Enforce the size constraint
  static_assert(sizeof(Node) <= 16, "Node class exceeds size limit");
  static_assert(sizeof(Edge) <= 32, "Edge class exceeds size limit");

  struct internal_node {
    Point           coords;
    node_value_type property;
    bool            valid;
    size_type       node_idx;

    internal_node(const Point p, const node_value_type val, 
                  const bool t,  const size_type i)
    : coords(p), property(val), valid(t), node_idx(i) {}
  };

  struct internal_edge {
    size_type       node1_index;
    size_type       node2_index;
    edge_value_type property;
    bool            valid;
    size_type       edge_idx;

    internal_edge(const size_type n1, const size_type n2, const edge_value_type v, 
                  const bool t, const size_type i)
    : node1_index(n1), node2_index(n2), property(v), valid(t), edge_idx(i) {}
  };
  
  /** Storage for node points */
  std::vector<internal_node> nodes_vec;
  
  /** Storage for edges, smaller node number first */
  std::vector<internal_edge> edges_vec;

  /** Fast iteration over all edges */
  std::map<size_type, std::map<size_type, size_type>> nested_map;
  // <1st node index, <2nd node index, edge index>>

  /** vector to store the edge_uid; indexed by edge_idx */
  std::vector<size_type> edge_i2u;

  /** vector to store the node_uid, indexed by node_idx */
  std::vector<size_type> node_i2u;
};

#endif
