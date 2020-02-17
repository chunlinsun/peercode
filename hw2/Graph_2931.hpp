#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>
#include <map>
#include <functional>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph{
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
  Graph() {}

  /** Default destructor */
  ~Graph() = default;

  //
  // nodes
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node>{
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
      this->graph = nullptr;
      this->idx = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      return (const_cast<Graph*>(this->graph))->\
          nodes[(const_cast<Graph*>(this->graph))->n_i2u[this->idx]].first;
    }

    /** Return a reference to this node's position. */
    Point& position(){
      return (const_cast<Graph*>(this->graph))->\
          nodes[(const_cast<Graph*>(this->graph))->n_i2u[this->idx]].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    /**@return A reference to the value associated with this node
     */
    node_value_type& value() {
      return (const_cast<Graph*>(this->graph))->\
          values[(const_cast<Graph*>(this->graph))->n_i2u[this->idx]];
    }

    /**@return A const reference to the value associated with this node
     */
    const node_value_type& value() const {
      return (const_cast<Graph*>(this->graph))->\
          values[(const_cast<Graph*>(this->graph))->n_i2u[this->idx]];
    }

    /**@return an unsigned int that is equal to the number of edges incident
     * to this node
     */
    size_type degree() const {
      return this->graph->adj.at(this->graph->n_i2u[idx]).size();
    }

    /**@return an incident iterator that points to the first incident edge
     */
    incident_iterator edge_begin() const{
      typename std::map<unsigned int, unsigned int>::iterator ai = \
            (const_cast<Graph*>(this->graph))->adj.at((const_cast<Graph*>(this->graph))->\
                    n_i2u[idx]).begin();
      return IncidentIterator(this->graph, ai, this->index_g());
    }

    /**@return an incident iterator that represents the end of the collection
     * of incident edges.
     */
    incident_iterator edge_end() const{
      typename std::map<unsigned int, unsigned int>::iterator ai =\
            (const_cast<Graph*>(this->graph))->adj.at((const_cast<Graph*>(this->graph))->\
                    n_i2u[idx]).end();
      return IncidentIterator(this->graph, ai, this->index_g());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (this->graph == n.graph && this->idx == n.idx){
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
      std::set<Graph*> graphs = {const_cast<Graph*>(this->graph), const_cast<Graph*>(n.graph)};
      if (const_cast<Graph*>(this->graph) == const_cast<Graph*>(n.graph)){
        if (this->index_g() < n.index_g()){
          return true;
        } else {
          return false;
        }
      } else if (*graphs.begin()== const_cast<Graph*>(this->graph)){
        return true;
      } else {
        return false;
      }
    }

    private:
    unsigned int index_g() const{
      return this->graph->n_i2u[this->idx];
    }

    bool valid() const{
    return graph->n_i2u[idx] >= 0 && graph->n_i2u[idx] < graph->nodes.size()         // uid in range.
           && idx < graph->n_i2u.size()  // idx in range.
           && graph->n_i2u[idx] == graph->nodes[graph->n_i2u[idx]].second; // uid in sync.
    }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    const Graph* graph;
    unsigned int idx;

    Node(const Graph* g, unsigned int idx_){
      this->graph = g;
      this->idx = idx_;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->n_i2u.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    unsigned int index_global = this->nodes.size();
    unsigned int index_user = this->n_i2u.size();
    this->nodes[index_global].first = position;
    this->nodes[index_global].second = index_user;
    this->values.push_back(value);
    this->n_i2u.push_back(index_global);
    this->adj[index_global] = {};

    return Node(this, index_user);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this == n.graph){
      return true;
    } else {
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
    // HW0: YOUR CODE HERE
    if (i < this->num_nodes()){
      return Node(this, i);
    } else {
      std::cout << "index chosen >= number of nodes in graph, invalid node \
        returned" << std::endl;
      return Node();
    }
  }

  //
  // edges
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      this->n1_gid = 0;
      this->n2_gid = 0;
      this->idx = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      Node node1 = Node(this->graph, const_cast<Graph*>(this->graph)->nodes[\
              this->n1_gid].second);
      return node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      Node node2 = Node(this->graph, const_cast<Graph*>(this->graph)->nodes[\
              this->n2_gid].second);
      return node2;
    }

    /**@return A reference to the value associated with this edge
     */
    edge_value_type& value() {
      return (const_cast<Graph*>(this->graph))->\
          e_values[const_cast<Graph*>(this->graph)->e_i2u[this->idx]];
    }

    /**@return A const reference to the value associated with this edge
     */
    const edge_value_type& value() const {
      return (const_cast<Graph*>(this->graph))->\
          e_values[const_cast<Graph*>(this->graph)->e_i2u[this->idx]];
    }

    /**@return the current length of the edge */
    double length() {
      Node node1 = Node(this->graph, const_cast<Graph*>(this->graph)->nodes[\
              this->n1_gid].second);
      Node node2 = Node(this->graph, const_cast<Graph*>(this->graph)->nodes[\
              this->n2_gid].second);
      Point p1 = node1.position();
      Point p2 = node2.position();
      return norm(p1-p2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->idx == e.idx && this->graph == e.graph &&\
              ((this->n1_gid == e.n1_gid && this->n2_gid == e.n2_gid) ||\
              (this->n1_gid == e.n2_gid && this->n2_gid == e.n1_gid))){
        return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      std::set<Graph*> graphs = {const_cast<Graph*>(this->graph), const_cast<Graph*>(e.graph)};
      if (const_cast<Graph*>(this->graph) == const_cast<Graph*>(e.graph)){
        if (this->index_g() < e.index_g()){
          return true;
        } else {
          return false;
        }
      } else if (*graphs.begin()== const_cast<Graph*>(this->graph)){
        return true;
      } else {
        return false;
      }
    }

    private:
    unsigned int index_g() const{
      return const_cast<Graph*>(this->graph)->e_i2u[this->idx];
    }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const Graph* graph;
    unsigned int n1_gid;
    unsigned int n2_gid;
    unsigned int idx;
      
    Edge(const Graph* g, unsigned int n1, unsigned int n2, unsigned int idx_){
      this->graph = g;
      this->n1_gid = n1;
      this->n2_gid = n2;
      this->idx = idx_;
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->e_i2u.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    //if (this->edges.find(const_cast<Graph*>(this)->e_i2u[i]) == this->edges.end()) {
    //  std::cout << "index not found; outputting invalid Edge." << std::endl;
    //  return Edge();
    //}
    unsigned int n1_gi = (this->edges.at(const_cast<Graph*>(this)->e_i2u[i])).first[0];
    unsigned int n2_gi = (this->edges.at(const_cast<Graph*>(this)->e_i2u[i])).first[1];
    return Edge(this, n1_gi, n2_gi, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (this->adj.size() == 0){
      return false;
    } else if (this->adj.find(this->n_i2u[a.index()]) == this->adj.end()) {
      return false;
    } else if (this->adj.at(this->n_i2u[a.index()]).size() == 0){
      return false;
    } else if (this->adj.at(this->n_i2u[a.index()]).find(this->n_i2u[b.index()])\
            == this->adj.at(this->n_i2u[a.index()]).end() ){
      return false;
    } else {
      return true;
    }
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
    if(this->has_edge(a, b)){
      unsigned int idx_global = this->adj[this->n_i2u[a.index()]][this->n_i2u[b.index()]];
      unsigned int idx_user = this->edges[idx_global].second;
      return this->edge(idx_user);
    } else if (a != b) {
      // Add entry to edges
      unsigned int idx_global = this->edges.size();
      unsigned int idx_user = this->e_i2u.size();
      unsigned int n1_ui = a.index();
      unsigned int n2_ui = b.index();
      unsigned int n1_gi = this->n_i2u[n1_ui];
      unsigned int n2_gi = this->n_i2u[n2_ui];
      std::vector<unsigned int> nds = {n1_gi, n2_gi};
      this->edges[idx_global].first = nds;
      this->edges[idx_global].second = idx_user;
      this->e_i2u.push_back(idx_global);
      // Add nodes to adjacency matrix
      this->adj[this->n_i2u[a.index()]][this->n_i2u[b.index()]] = idx_global;
      this->adj[this->n_i2u[b.index()]][this->n_i2u[a.index()]] = idx_global;
      // Add default-constructed value
      this->e_values.push_back(edge_value_type());

      return Edge(this, n1_gi, n2_gi, idx_user);
    } else {
      std::cout << "Tried to add edge where node is attached to itself, outputting invalid edge" << std::endl;
      return Edge();
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    this->nodes = {};
    this->edges = {};
    this->adj = {};
    this->values = {};
    this->e_values = {};
    this->n_i2u = {};
    this->e_i2u = {};
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
    /** @return the node associated with this node iterator
     */
    Node operator*() {
      return Node(graph, idx);
    }

    /** Increment the node iterator by one position
     * @return the node iterator post-increment.
     */
    NodeIterator& operator++(){
      this->idx = this->idx + 1;
      return *this;
    }

    /**@param[in] _ni_ a node iterator to check for equivalence
     * @return true if the two node iterators are equivalent
     */
    bool operator==(const NodeIterator& ni) const {
      if (this->graph == ni.graph && this->idx == ni.idx){
        return true;
      } else {
        return false;
      }
    }

    /**@param[in] _ni_ a node iterator to check for non-equivalence
     * @return true if the two node iterators are not equivalent
     */
    bool operator!=(const NodeIterator& ni) const {
      if (this->graph == ni.graph && this->idx == ni.idx){
        return false;
      } else {
        return true;
      }
    }

  private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph;
    unsigned int idx;
    
    NodeIterator(Graph* g, unsigned int idx_){
      graph = g;
      idx = idx_;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /**@return a node iterator that points to the first node in the graph
   */
  node_iterator node_begin() const{
    return NodeIterator(const_cast<Graph*>(this), 0);
  }

  /**@return a node iterator that represents the end of the graph nodes
   */
  node_iterator node_end() const{
    unsigned int end = this->n_i2u.size();
    return NodeIterator(const_cast<Graph*>(this), end);
  }

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
    /** @return the edge associated with this incident iterator
     */
    Edge operator*() const{
      unsigned int idx = (this->adj_iter)->second;
      unsigned int idx_user = const_cast<Graph*>(this->graph)->edges[idx].second;
      if (const_cast<Graph*>(this->graph)->edges[idx].first[0] == this->base_node_gi){
        unsigned int n2_gi = const_cast<Graph*>(this->graph)->edges[idx].first[1];
        return Edge(this->graph, this->base_node_gi, n2_gi, idx_user);
      } else {
        unsigned int n2_gi = const_cast<Graph*>(this->graph)->edges[idx].first[0];
        return Edge(this->graph, this->base_node_gi, n2_gi, idx_user);
      }
    }

    /** Increment the incident iterator by one position
     * @return the incident iterator post-increment.
     */
    IncidentIterator& operator++(){
      ++(this->adj_iter);
      return *this;
    }

    /**@param[in] _ii_ an incident iterator to check for equivalence
     * @return true if the two incident iterators are equivalent
     */
    bool operator==(const IncidentIterator& ii) const{
      if (this->graph == ii.graph && this->adj.iter == ii.adj_iter &&\
              this->base_node_gi == ii.base_node_gi){
        return true;
      } else {
        return false;
      }
    }

    /**@param[in] _ii_ an incident iterator to check for non-equivalence
     * @return true if the two incident iterators are not equivalent
     */
    bool operator!=(const IncidentIterator& ii) const{
      if (this->graph == ii.graph && this->adj_iter == ii.adj_iter &&\
              this->base_node_gi == ii.base_node_gi){
        return false;
      } else {
        return true;
      }
    }
      

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph;
    typename std::map<unsigned int, unsigned int>::iterator adj_iter;
    unsigned int base_node_gi;

    IncidentIterator(const Graph* g, typename std::map<unsigned int, unsigned int>::iterator ai,\
            unsigned int n1){
      this->graph = g;
      this->adj_iter = ai;
      this->base_node_gi = n1;
    }
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
    /** @return the edge associated with this edge iterator
     */
    Edge operator*() const{
      unsigned int n1_gi = this->graph->edges.at(this->graph->e_i2u[this->idx]).first[0];
      unsigned int n2_gi = this->graph->edges.at(this->graph->e_i2u[this->idx]).first[1];
      return Edge(this->graph, n1_gi, n2_gi, this->idx);
    }

    /** Increment the edge iterator by one position
     * @return the edge iterator post-increment.
     */
    EdgeIterator& operator++(){
      this->idx = this->idx + 1;
      return *this;
    }

    /**@param[in] _ni_ an edge iterator to check for equivalence
     * @return true if the two edge iterators are equivalent
     */
    bool operator==(const EdgeIterator& ei) const{
      if (this->graph == ei.graph && this->idx == ei.idx){
        return true;
      } else {
        return false;
      }
    }

    /**@param[in] _ni_ an edge iterator to check for non-equivalence
     * @return true if the two edge iterators are not equivalent
     */
    bool operator!=(const EdgeIterator& ei) const{
      if (this->graph == ei.graph && this->idx == ei.idx){
        return false;
      } else {
        return true;
      }
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph;
    unsigned int idx;

    EdgeIterator(const Graph* g, unsigned int i){
      this->graph = g;
      this->idx = i;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  /**@return an edge iterator that points to the first edge in the graph.
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /**@return an edge iterator that represents the end of the graph edges.
   */
  edge_iterator edge_end() const{
    unsigned int num_edges = this->e_i2u.size();
    return EdgeIterator(this, num_edges);
  }

/**
 * @brief Edge removal given two nodes
 *
 * @param[in] a   A node object
 * @param[in] b   A different node object
 * @return idx_user   The active user id of the edge that replaced the removed edge
 *
 * @pre e_i2u.size() > 0 (there are edges in the graph)
 * @pre _a_ and _b_ are both valid nodes, belong to the same graph, and have an edge between them
 * @post There is one less edge in the graph: this amounts to e_i2u losing one element
 */
  unsigned int remove_edge(const Node& a, const Node& b){
    if (this->e_i2u.size() > 0){
      unsigned int idx_global = this->adj[this->n_i2u[a.index()]][this->n_i2u[b.index()]];
      unsigned int idx_user = this->edges[idx_global].second;
      // swap & pop
      this->e_i2u[idx_user] = this->e_i2u[this->e_i2u.size()-1];
      this->e_i2u.pop_back();
      // update active user id
      this->edges[this->e_i2u[idx_user]].second = idx_user;
      // remove from adj list
      this->adj[this->n_i2u[a.index()]].erase(this->n_i2u[b.index()]);
      this->adj[this->n_i2u[b.index()]].erase(this->n_i2u[a.index()]);

      return idx_user;
    } else {
      return 0;
    }
  }

/**
 * @brief Edge removal given an edge object
 *
 * @param[in] e   An edge object
 * @return idx_user   The active user id of the edge that replaced the removed edge
 *
 * @pre e_i2u.size() > 0 (there are edges in the graph)
 * @pre _e_ is a valid edge object
 * @post There is one less edge in the graph: this amounts to e_i2u losing one element
 */
  unsigned int remove_edge(const Edge& e){
    Node node1 = Node(this, this->nodes[e.n1_gid].second);
    Node node2 = Node(this, this->nodes[e.n2_gid].second);
    return remove_edge(node1, node2);
  }

/**
 * @brief Edge removal given an edge iterator
 *
 * @param[in] e_it   An edge iterator
 * @return The edge iterator corresponding to the first edge in the user-facing edge order
 *
 * @pre e_i2u.size() > 0 (there are edges in the graph)
 * @pre *_e_it_ is a valid edge object
 * @post There is one less edge in the graph: this amounts to e_i2u losing one element
 */
  edge_iterator remove_edge(edge_iterator e_it){
    unsigned int i = remove_edge(*e_it);
    (void) i;
    return this->edge_begin();
  }

/**
 * @brief Node removal given node reference
 *
 * @param[in] n   A node object
 * @return idx_user   The active user id of the node that replaced the removed node
 *
 * @pre n_i2u.size() > 0 (there are edges in the graph)
 * @pre _n_ is a valid node object
 * @post There is one less node in the graph: this amounts to n_i2u losing one element
 */
  unsigned int remove_node(const Node& n){
    unsigned int idx_user = n.index();
    unsigned int idx_global = this->n_i2u[idx_user];
    if (this->num_nodes() > 0){
      // remove edges attached to this node
      IncidentIterator ii = n.edge_begin();
      bool cond = true;
      if (this->adj[idx_global].size() == 0){
        cond = false;
      }
      while (cond){
        remove_edge(*ii);
        ii = n.edge_begin();
        if (this->adj[idx_global].size() == 0){
          break;
        }
      }
      // swap & pop
      this->n_i2u[idx_user] = this->n_i2u[this->num_nodes()-1];
      this->n_i2u.pop_back();
      // update active user id
      this->nodes[this->n_i2u[idx_user]].second = idx_user;
      //this->adj.erase(idx_global);

      return idx_user;
    } else {
      return 0;
    }
  }
  
/**
 * @brief Node removal given node iterator
 *
 * @param[in] n_it   A node iterator
 * @return The node iterator that corresponds to the first node in the user-facing order
 *
 * @pre n_i2u.size() > 0 (there are nodes in the graph)
 * @pre *_n_it_ is a valid node object
 * @post There is one less node in the graph: this amounts to n_i2u losing one element
 */
  node_iterator  remove_node(node_iterator  n_it){
    remove_node(*n_it);
    return this->node_begin();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::map<unsigned int, std::pair<Point, unsigned int>> nodes;
  std::vector<V> values;
  std::vector<E> e_values;
  std::map<unsigned int, std::pair<std::vector<unsigned int>, unsigned int>> edges;
  std::map<unsigned int, std::map<unsigned int, unsigned int>> adj;
  std::vector<unsigned int> n_i2u;
  std::vector<unsigned int> e_i2u;
};

#endif // CME212_GRAPH_HPP