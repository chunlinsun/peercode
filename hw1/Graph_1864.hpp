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
      return this->graph->nodes[this->idx];
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
      return (const_cast<Graph*>(this->graph))->values[this->idx];
    }

    /**@return A const reference to the value associated with this node
     */
    const node_value_type& value() const {
      return (const_cast<Graph*>(this->graph))->values[this->idx];
    }

    /**@return an unsigned int that is equal to the number of edges incident
     * to this node
     */
    size_type degree() const {
      //--functionality_1
      //--Avoid using objects as keys in your maps. 
      //--This breaks your degree function since different instances of
      //--objects do not correspond to the same keys.
      //--Consider using node indexes as the key instead.
      //--START
      return this->graph->adj.at(*this).size();
      //--END
    }

    /**@return an incident iterator that points to the first incident edge
     */
    incident_iterator edge_begin() const{
      typename std::map<Node, unsigned int>::iterator ai = \
            (const_cast<Graph*>(this->graph))->adj.at(*this).begin();
      return IncidentIterator(this->graph, ai, this);
    }

    /**@return an incident iterator that represents the end of the collection
     * of incident edges.
     */
    incident_iterator edge_end() const{
      typename std::map<Node, unsigned int>::iterator ai =\
            (const_cast<Graph*>(this->graph))->adj.at(*this).end();
      return IncidentIterator(this->graph, ai, this);
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
      if (this->idx < n.idx){
        return true;
      } else {
        return false;
      }
    }

    private:
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
    return this->nodes.size();
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
  //--The node is supposed to take in a value as argument. 
  //--No points deducted since you may have been confused by the missing argument name.
  //--Please incorporte the following fix.
  //--START
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    this->nodes.push_back(position);
    this->values.push_back(value);
    unsigned int index = this->size()-1;
    Graph* graph = this;
    return Node(graph, index);
  }
  //--END

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
      this->node_1 = nullptr;
      this->node_2 = nullptr;
      this->idx = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return *(this->node_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return *(this->node_2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      const Graph* g1 = (this->node_1)->graph;
      const Graph* g2 = (e.node_1)->graph;
      if (this->idx == e.idx && g1 == g2){
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
      if (this->idx < e.idx){
        return true;
      } else {
        return false;
      }
    }

    private:
    // Allow Graph to access Edge's private member data and functions.
      friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
      const Node* node_1;
      const Node* node_2;
      unsigned int idx;
      
      Edge(const Node* n1, const Node* n2, unsigned int idx_){
        this->node_1 = n1;
        this->node_2 = n2;
        this->idx = idx_;
      }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (this->edges.find(i) == this->edges.end()) {
      std::cout << "index not found; outputting invalid Edge." << std::endl;
      return Edge();
    }
    const Node* n1 = &((this->edges.at(i))[0]);
    const Node* n2 = &((this->edges.at(i))[1]);
    return Edge(n1, n2, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (this->adj.find(a) == this->adj.end()) {
      return false;
    }
    if (this->adj.at(a).find(b) == this->adj.at(a).end() ){
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
      unsigned int idx = this->adj[a][b];
      return this->edge(idx);
    } else {
      // Add entry to edges
      unsigned int idx = this->edges.size();
      std::vector<Node> nodes = {a, b};
      this->edges[idx] = nodes;
      // Add nodes to adjacency matrix
      this->adj[a][b] = idx;
      this->adj[b][a] = idx;

      return Edge(&a, &b, idx);
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
    friend class Graph<V>;
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
    unsigned int end = this->nodes.size();
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
      const Node* n_2 = &((this->adj_iter)->first);
      unsigned int idx = (this->adj_iter)->second;
      return Edge(this->base_node, n_2, idx);
    }

    /** Increment the incident iterator by one position
     * @return the incident iterator post-increment.
     */
    IncidentIterator& operator++(){
      ++(this->adj_iter);
      return *this;
    }

    /**@param[in] _ni_ an incident iterator to check for equivalence
     * @return true if the two incident iterators are equivalent
     */
    bool operator==(const IncidentIterator& ii) const{
      if (this->graph == ii.graph && this->adj.iter == ii.adj_iter && this->base_node == ii.base_node){
        return true;
      } else {
        return false;
      }
    }

    /**@param[in] _ni_ an incident iterator to check for non-equivalence
     * @return true if the two incident iterators are not equivalent
     */
    bool operator!=(const IncidentIterator& ii) const{
      if (this->graph == ii.graph && this->adj_iter == ii.adj_iter && this->base_node == ii.base_node){
        return false;
      } else {
        return true;
      }
    }
      

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph;
    typename std::map<Node, unsigned int>::iterator adj_iter;
    const Node* base_node;

    IncidentIterator(const Graph* g, typename std::map<Node, unsigned int>::iterator ai,\
            const Node* n1){
      this->graph = g;
      this->adj_iter = ai;
      this->base_node = n1;
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
      const Node* n_1 = &(this->graph->edges.at(this->idx)[0]);
      const Node* n_2 = &(this->graph->edges.at(this->idx)[1]);
      return Edge(n_1, n_2, this->idx);
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
    unsigned int num_edges = this->edges.size();
    return EdgeIterator(this, num_edges);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> nodes;
  std::vector<V> values;
  std::map<unsigned int, std::vector<Node>> edges;
  std::map<Node, std::map<Node, unsigned int>> adj;
};

#endif // CME212_GRAPH_HPP
