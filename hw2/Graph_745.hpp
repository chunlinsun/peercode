#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <map>
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
    /** Construct an invalid node. */
    Node() {
    }

    /** Change this node's position. */
    Point& position() {
      return g->nodes.at(uid);
    }

    /** Return this node's position. */
    const Point& position() const {
      return g->nodes.at(uid);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return g->node_uid_to_id.at(uid);
    }
    
    /** Set a Node's value.
     * @param[in] val of the type V, the type of values in the graph.
     * @post The Node's value has been updated to val.
     *
     * Complexity: O(1) amortized operations.
     */
    node_value_type& value() {
      return g->nodeValues.at(uid);
    }


    node_value_type& value_set(V val){
      g->nodeValues.at(uid) = val;
      return g->nodeValues.at(uid);
    }

    /** Return the value of a node.
     * Complexity: O(1) amortized operations.
     */
    const node_value_type& value() const {
      return g->nodeValues.at(uid);
    }
    
    /** Return the degree of a node
     */
    size_type degree() const {
      return g->edgeMap.at(uid).size();
    }

    /** Return an IncidentIterator that points to the beginning of the node's
    list of adjacent edges. */
    incident_iterator edge_begin() const {
      size_type my_id = g->node_uid_to_id.at(uid);
      return IncidentIterator(g, my_id, 0);
    }

    /** Return an IncidentIterator that points to the end of the node's
    list of adjacent edges. */
    incident_iterator edge_end() const {
      size_type my_id = g->node_uid_to_id.at(uid);
      return IncidentIterator(g, my_id, g->edgeMap.at(uid).size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (g == n.g && uid == n.uid);
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
      return ((g == n.g && uid < n.uid) || std::less<Graph *>{}(this->g, n.g));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /**Node attributes*/
    Graph* g;
    size_type uid;
    
    /**Valid Node constructor*/
    Node(const Graph* g, size_type id) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->uid = g->node_id_to_uid.at(id);
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_id_to_uid.size();
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
    nodes.push_back(position);
    nodeValues.push_back(val);

    //Update node_id_to_uid and node_uid_to_id
    node_id_to_uid.push_back(nodes.size() - 1);
    node_uid_to_id.push_back(node_id_to_uid.size() - 1);

    std::vector <std::tuple <size_type, size_type>> vec;
    edgeMap.insert(std::pair<size_type, std::vector<std::tuple<size_type, size_type>>>(node_id_to_uid.at(num_nodes() - 1), vec));
    return Node(this, num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    size_type my_id = node_uid_to_id.at(n.uid);
    return (num_nodes() > my_id);
  }

   /** @post All Node, Edge, NodeIterator, IncidentIterator,
     * and EdgeIterator objects are invalidated.
     *
     * @post has_edge on any edge with this node as node1()
     * or node2() returns false.
     *
     * @post has_edge on any edge without this node as node1()
     * or node2()returns the same result as before the function was called.
     * 
     * @post The node index corresponding to any node may change.
     *
     * @post The edge index corresponding to any pair of nodes may change.
     *
     * @runtime O(n^2) for a dense graph.
     * If each node has a small number of neighbors, then O(1).
    **/
  size_type remove_node(const Node& n) {
    if (!has_node(n)) {
      return 0;
    }
    //remove all edges adjacent to the node
    if (edgeMap.count(n.uid) != 0) {
      while (edgeMap.at(n.uid).size() > 0) {
        size_type edge_uid = std::get<1>(edgeMap.at(n.uid).at(0));
        size_type edge_id = edge_uid_to_id.at(edge_uid);
        Edge e = Edge(this, edge_id, false);
        remove_edge(e);
      }
    }
    size_type largest = 0;
    for (unsigned int i = 0; i < node_uid_to_id.size(); i++) {
      if (node_uid_to_id.at(i) > largest)
        largest = node_uid_to_id.at(i);
    }
    size_type largest_uid = 0;
    for (unsigned int i = 0; i < node_id_to_uid.size(); i++) {
      if (node_id_to_uid.at(i) > largest_uid)
        largest_uid = node_id_to_uid.at(i);
    }
    size_type uid_of_highest_id = node_id_to_uid.back();
    node_id_to_uid.at(node_uid_to_id.at(n.uid)) = uid_of_highest_id; //swap
    node_uid_to_id.at(uid_of_highest_id) = node_uid_to_id.at(n.uid); // adjust id's
    node_id_to_uid.pop_back(); //pop

    largest = 0;
    for (unsigned int i = 0; i < node_uid_to_id.size(); i++) {
      if (node_uid_to_id.at(i) > largest)
        largest = node_uid_to_id.at(i);
    }
    
    return 1;
  }

   /** @pre n_it must be a valid NodeIterator that points to
     * a node in the graph.
     *
     * @post All Node, Edge, NodeIterator, IncidentIterator,
     * and EdgeIterator objects are invalidated.
     *
     * @post has_edge on any edge with this node as node1()
     * or node2() returns false.
     *
     * @post has_edge on any edge without this node as node1()
     * or node2()returns the same result as before the function was called.
     * 
     * @post The node index corresponding to any node may change.
     *
     * @post The edge index corresponding to any pair of nodes may change.
     *
     * The node_iterator returned will be valid if num_nodes() is
     * > 1 before the function is called.
     *
     * @runtime O(n^2) for a dense graph.
     * If each node has a small number of neighbors, then O(1).
    **/
  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    remove_node(n);
    return node_begin();
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    edge_value_type& value() {
      return g->edgeValues.at(uid);
    }

    /** Return the value of a node.
     * Complexity: O(1) amortized operations.
     */
    const edge_value_type& value() const {
      return g->edgeValues.at(uid);
    }

    size_type index() const {
      return g->edge_uid_to_id.at(uid);
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (backwards == false) {
        return Node(g, g->node_uid_to_id.at(std::get<0>(g->edgeVector.at(uid))));
      }
      else {
        return Node(g, g->node_uid_to_id.at(std::get<1>(g->edgeVector.at(uid))));
      }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (backwards == false) {
        return Node(g, g->node_uid_to_id.at(std::get<1>(g->edgeVector.at(uid))));
      }
      else {
        return Node(g, g->node_uid_to_id.at(std::get<0>(g->edgeVector.at(uid))));
      }
    }

    double length() const {
      Point diff = node2().position() - node1().position();
      return sqrt((diff.x)*(diff.x) + (diff.y)*(diff.y) + (diff.z)*(diff.z));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (g == e.g && uid == e.uid);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((g == e.g && uid < e.uid) || std::less<Graph *>{}(this->g, e.g));
    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    /** Edge attributes */
    Graph* g;
    size_type uid;
    bool backwards;

    /** Valid Edge constructor */
    Edge(const Graph* g, size_type id, bool backwards) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->uid = g->edge_id_to_uid.at(id);
      this->backwards = backwards;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_id_to_uid.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i, false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (a == b){
      return false; //Cannot have self-loops
    }
    else {
      //Check whether a is in the list of nodes adjacent to b.
      if (edgeMap.count(b.uid) == 0)
        return false;
      if (edgeMap.at(b.uid).size() == 0)
        return false;
      else {
        for (unsigned int i = 0; i < edgeMap.at(b.uid).size(); i++) {
          if (std::get<0>(edgeMap.at(b.uid).at(i)) == a.uid)
            return true;
        }
      }
    }
    return false; //This means we did not find it.
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

  Edge add_edge(const Node& a, const Node& b, const edge_value_type& eval = edge_value_type()) { 
    assert((a == b) == false); //Cannot create self-loops.
    //If we already have the edge, return that edge.
    if (has_edge(a, b)){
      for (unsigned int i = 0; i < edgeMap.at(a.uid).size(); i++) {
        if (std::get<0>(edgeMap.at(a.uid).at(i)) == b.uid) {
          size_type my_id = edge_uid_to_id.at(std::get<1>(edgeMap.at(a.uid).at(i)));
          return Edge(this, my_id, (b < a));
          }
        }
      }
    //Add a new edge to the edgeVector.
    if (a < b) {
      std::tuple <size_type, size_type> tupV = std::make_tuple(a.uid, b.uid);
      edgeVector.push_back(tupV);
    }
    else {
      std::tuple <size_type, size_type> tupV = std::make_tuple(b.uid, a.uid);
      edgeVector.push_back(tupV);
    }

    //Add the edge's value to the edgeValues vector.
    edgeValues.push_back(eval);

    //Update edge_id_to_uid and edge_uid_to_id
    edge_id_to_uid.push_back(edgeVector.size() - 1);
    edge_uid_to_id.push_back(edge_id_to_uid.size() - 1);

    //Add a new edge to the edgeMap at index a.
    std::tuple <size_type, size_type> tupM = std::make_tuple(b.uid, edgeVector.size() - 1);
    edgeMap.at(a.uid).push_back(tupM);
    
    //Add a new edge to the edgeMap at index b.
    tupM = std::make_tuple(a.uid, edgeVector.size() - 1);
    edgeMap.at(b.uid).push_back(tupM);

    return Edge(this, num_edges() - 1, (b < a)); //Return the edge we added.
  }

   /** @pre e must have index less than num_edges.
     *
     * @post All Edge objects, IncidentIterator objects,
     * and EdgeIterator objects are invalidated.
     *
     * @post has_edge on the nodes in this edge will return false.
     *
     * @post has_edge on any other pair of nodes will return the same result
     * as before the function was called.
     *
     * @post The edge index corresponding to any pair of nodes may change.
     *
     * @runtime O(n) for a dense graph.
     * If each node has a small number of neighbors, then O(1).
     */
  size_type remove_edge(const Edge& e) {
    if (!has_edge(e.node1(), e.node2())) {
      return 0;
    }
    Node my_node1 = e.node1();
    Node my_node2 = e.node2();
    //remove node2 from node1's adjacency list.
    for (unsigned int i = 0; i < edgeMap.at(my_node1.uid).size(); i++) {
      if (std::get<0>(edgeMap.at(my_node1.uid).at(i)) == my_node2.uid) {
        edgeMap.at(my_node1.uid).at(i) = edgeMap.at(my_node1.uid).back();
        edgeMap.at(my_node1.uid).pop_back();
      }
    }
    //remove node1 from node2's adjacency list.
    for (unsigned int i = 0; i < edgeMap.at(my_node2.uid).size(); i++) {
      if (std::get<0>(edgeMap.at(my_node2.uid).at(i)) == my_node1.uid) {
        edgeMap.at(my_node2.uid).at(i) = edgeMap.at(my_node2.uid).back();
        edgeMap.at(my_node2.uid).pop_back();
      }
    }
    //update edge_id_to_uid and edge_uid_to_id
    size_type last_uid = edge_id_to_uid.back();
    edge_id_to_uid.at(edge_uid_to_id.at(e.uid)) = last_uid; //swap
    edge_uid_to_id.at(last_uid) = edge_uid_to_id.at(e.uid); // adjust id's
    edge_id_to_uid.pop_back(); //pop
    return 1;
  }

   /** @pre Both nodes must have index less than num_nodes().
     *
     * @post All Edge objects, IncidentIterator objects,
     * and EdgeIterator objects are invalidated.
     *
     * @post has_edge on these nodes returns false.
     *
     * @post has_edge on any other pair of nodes returns the same result
     * as before the function was called.
     *
     * @post If the nodes had an edge between them before the function was called,
     * the function returns 1. Else, the function returns 0.
     *
     * @post The edge index corresponding to any pair of nodes may change.
     *
     * @runtime O(n) for a dense graph.
     * If each node has a small number of neighbors, then O(1).
     */
  size_type remove_edge(const Node& my_node1, const Node& my_node2) {
    //find the edge uid.
    if (!has_edge(my_node1, my_node2))
      return 0;
    Edge e = Edge();
    for (unsigned int i = 0; i < edgeMap.at(my_node1.uid).size(); i++) {
      if (std::get<0>(edgeMap.at(my_node1.uid).at(i)) == my_node2.uid) {
        size_type my_uid = std::get<1>(edgeMap.at(my_node1.uid).at(i));
        e = edge(edge_uid_to_id.at(my_uid));
      }
    }
    return remove_edge(e);
  }

   /** @pre e_it must be a valid edge_iterator that points to an edge.
     *
     * @post All Edge objects, IncidentIterator objects,
     * and EdgeIterator objects are invalidated.
     *
     * @post The edge_iterator returned will be valid if num_edges is
     * > 1 before the function is called.
     *
     * @post has_edge on the nodes in this edge returns false.     
     *
     * @post has_edge on any other pair of nodes will return the same result
     * as before the function was called.
     *
     * @post The edge index corresponding to any pair of nodes may change.
     *
     * @runtime O(n) for a dense graph.
     * If each node has a small number of neighbors, then O(1).
     */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    remove_edge(e);
    return edge_begin();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    std::vector<Point> n;
    nodes = n;
    std::vector<node_value_type> nval;
    nodeValues = nval;
    std::vector<edge_value_type> eval;
    edgeValues = eval;
    std::map<size_type, std::vector<std::tuple <size_type, size_type>>> em;
    edgeMap = em;
    std::vector<std::tuple<size_type, size_type>> ev;
    edgeVector = ev;
    std::vector<size_type> nullvec;
    edge_uid_to_id = nullvec;
    edge_id_to_uid = nullvec;
    node_uid_to_id = nullvec;
    node_id_to_uid = nullvec;
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
    NodeIterator() {
    }

    /** Return the node that a NodeIterator is pointing to.
     * @pre idx is a valid index in the nodes list in Graph g.
     * @post Graph g has not been changed.
     *
     * Complexity: O(1) amortized operations.
     */  
    Node operator*() const {
      return Node(g, g->node_uid_to_id.at(uid));
    }

    /** Increment the NodeIterator and return it.
     * @pre idx is a valid index in the nodes list in Graph g.
     *
     * Complexity: O(1) amortized operations.
     */
    NodeIterator& operator++() {
      uid++;
      return *this;
    }

    /* Return true if this NodeIterator is equal to NodeIterator other. 
     */
    bool operator==(const NodeIterator& other) const {
      return (g == other.g && uid == other.uid);
    }

   private:
    friend class Graph;
    Graph* g;
    size_type uid;

    /* Valid NodeIterator Constructor
     */
    NodeIterator(const Graph* g, size_type uid) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->uid = uid;
    }
  };

  /** Return an NodeIterator that points to the first Node.
   */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }


  /** Return an NodeIterator that points to the last Node.
   */
  NodeIterator node_end() const {
    return NodeIterator(this, num_nodes());
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Return the Edge that the IncidentIterator is pointing to.
     */
    Edge operator*() const {
      size_type node2_uid = std::get<0>(g->edgeMap.at(node_uid).at(edge_num));
      size_type edge_uid = std::get<1>(g->edgeMap.at(node_uid).at(edge_num));
      size_type edge_id = g->edge_uid_to_id.at(edge_uid);
      return Edge(g, edge_id, (node2_uid < node_uid));
    }

    /** Increment the IncidentIterator and return it.
     */
    IncidentIterator& operator++() {
      edge_num++;
      return *this;
    }

    /** Return true if this IncidentIterator is equal to IncidentIterator other.
     */
    bool operator==(const IncidentIterator& other) const {
      return (g == other.g && node_uid == other.node_uid && edge_num == other.edge_num);
    }

   private:
    friend class Graph;
    Graph* g;
    size_type node_uid;
    size_type edge_num;

    /* Valid IncidentIterator Constructor
     */
    IncidentIterator(const Graph* g, size_type node_id, size_type edge_num) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->node_uid = g->node_id_to_uid.at(node_id);
      this->edge_num = edge_num;
    }
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
    EdgeIterator() {
    }
    
    /** Return the Edge that the EdgeIterator is pointing to.
     */
    Edge operator*() const {;
      return Edge(g, id, false);
    }

    /** Increment the EdgeIterator and return it.
     */
    EdgeIterator& operator++() {
      id++;
      return *this;
    }

    /** Return true if the EdgeIterator is equal to EdgeIterator other.
     */
    bool operator==(const EdgeIterator& other) const {
      return (g == other.g && id == other.id);
    }

   private:
    friend class Graph;
    Graph* g;
    size_type id;

    /* Valid EdgeIterator Constructor
     */
    EdgeIterator(const Graph* g, size_type id) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->id = id;
    }
  };
  
  /* Return an EdgeIterator that points to the first Edge in the graph.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  
  /* Return an EdgeIterator that points to the last Edge in the graph.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  std::vector<Point> nodes;
  std::vector<node_value_type> nodeValues;
  std::map<size_type, std::vector<std::tuple <size_type, size_type>>> edgeMap;
  std::vector<std::tuple<size_type, size_type>> edgeVector;
  std::vector<edge_value_type> edgeValues;
  std::vector<size_type> edge_id_to_uid;
  std::vector<size_type> edge_uid_to_id;
  std::vector<size_type> node_id_to_uid;
  std::vector<size_type> node_uid_to_id;
};

#endif // CME212_GRAPH_HPP
