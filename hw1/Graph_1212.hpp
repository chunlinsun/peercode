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
template<typename V>
class Graph{
 private:
    struct internal_point;
  // internal struct to store point information
    
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
  Graph() {
      points = std::vector<internal_point>();
      size_ = 0;
      edges = std::vector<std::vector<size_type>>();
      adjs = std::vector<std::vector<size_type>>();
      num_edges_ = 0;
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
        graph_ = nullptr;
      // HW0: YOUR CODE HERE
    }
     
    /** Return this node's position. */
    const Point& position() const {
      assert (graph_ != nullptr);
      return graph_->points[pid_].position;
      // HW0: YOUR CODE HERE
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert (graph_ != nullptr);
      return pid_;
      
    }

    /** @return the value of the node, also allow modification of the value*/
    node_value_type& value(){
        return graph_->points[pid_].value_;
    }

    /** @return the value of the node, does not allow modification of the value*/
    const node_value_type& value() const{
        return graph_->points[pid_].value_;
    }

    /** @return the number of incident edges*/
    size_type degree() const{
        return graph_->adjs[pid_].size();
    }

    /** @return the start of the incident iterator*/
    incident_iterator edge_begin() const{
        return IncidentIterator(graph_, 0, pid_);
    }

    /** @return the end of the incident iterator*/
    incident_iterator edge_end() const{
        return IncidentIterator(graph_, degree(), pid_);
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
      if (this->graph_ == n.graph_ && this->pid_ == n.pid_){
          return true;
      }
      // HW0: YOUR CODE HERE
      return false;
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
      std::less <void*> compare;
      // a function to compare pointers
      if (compare(this->graph_, n.graph_)){
          return true;
      }
      if (this->graph_ == n.graph_ && this->pid_ < n.pid_){
          return true;
      }
      // HW0: YOUR CODE HERE
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type pid_;
    /** Construct a valid node. */
    Node(const Graph* graph, size_type pid){
        graph_ = const_cast<Graph*>(graph);
        pid_ = pid;
    }
    // private constructor method of node
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
    // HW0: YOUR CODE HERE
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
    internal_point my_point(position, value);
    points.push_back(my_point);
    size_++;
    std::vector<size_type>nd{};
    adjs.push_back(nd);
    return Node(this, size_ - 1);
    // HW0: YOUR CODE HERE
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this) return true;
    // HW0: YOUR CODE HERE
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < size_) return Node(this, i);
    // HW0: YOUR CODE HERE
    return Node();        // Invalid node
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
      graph_ = nullptr;
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (graph_ != nullptr){
          return graph_->node(nd1_);
      }
      // HW0: YOUR CODE HERE
      return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (graph_ != nullptr){
          return graph_->node(nd2_);
      }
      // HW0: YOUR CODE HERE
      return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ == e.graph_ ){
          if (nd1_ == e.nd1_ && nd2_ == e.nd2_)
              return true;
          if (nd1_ == e.nd2_ && nd2_ == e.nd1_)
              return true;
      }
      //HW0: YOUR CODE HERE
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      std::less <void*> compare;
      // a function to compare pointers
      if (compare(this->graph_, e.graph_)){
          return true;
      }
      if (this->graph_ == e.graph_ ){
          if (std::min(nd1_, nd2_) < std::min(e.nd1_, e.nd2_))
              return true;
          if (std::min(nd1_, nd2_) == std::min(e.nd1_, e.nd2_) &&
                   std::max(nd1_, nd2_) < std::max(e.nd1_, e.nd2_))
              return true;
      }
      //HW0: YOUR CODE HERE
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type nd1_;
    size_type nd2_;

    /** Construct a valid Edge. */
    Edge(const Graph* graph, size_type nd1, size_type nd2){
        graph_ = const_cast<Graph*>(graph);
        nd1_ = nd1;
        nd2_ = nd2;
    }
    // private constructor method of edge
    
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
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
    if (i < num_edges_){
        return Edge(this, edges[i][0], edges[i][1]);
    }
    // HW0: YOUR CODE HERE
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check if a and b are distinct valid nodes of this graph
    if (!(a == b) && a.graph_ == this && b.graph_ == this){
          if (std::find(adjs[a.index()].begin(), adjs[a.index()].end(), b.index()) != adjs[a.index()].end())
              return true;
    }
    
    // HW0: YOUR CODE HERE
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
    // check if a and b are distinct valid nodes of this graph
    if (!(a == b) && a.graph_ == this && b.graph_ == this){
         if (has_edge(a, b))
             return Edge(this, a.index(), b.index());
         // has_edge(a,b) == false, need to create new edge
         adjs[a.index()].push_back(b.index());
         adjs[b.index()].push_back(a.index());
         std::vector<size_type> my_edge;
         my_edge.push_back(a.index());
         my_edge.push_back(b.index());
         edges.push_back(my_edge);
         // add new edge to the internal struct of edges
         num_edges_++;
         // new num_edges() == old num_edges() + 1.
         return Edge(this, a.index(), b.index());
    }
    // HW0: YOUR CODE HERE
    return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      size_ = 0;
      num_edges_ = 0;
      points = std::vector<internal_point>();
      edges = std::vector<std::vector<size_type>>();
      adjs = std::vector<std::vector<size_type>>();
    // HW0: YOUR CODE HERE
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** The dereference operator
     *@return the node the iterator refers to*/
    Node operator*() const{
        return grh->node(i);
    }

    /** the increment operator
     *@return the iterator which refers to the subsequent node*/
    NodeIterator& operator++(){
        i++;
        return *this;
    }

    /** the equality comparision operator
     *@param[in] ndit the iterator which is compared with the present iterator
     *@return true if they refer to the same node; false otherwise */
    bool operator==(const NodeIterator& ndit) const{
        return (grh == ndit.grh && i == ndit.i);
    }

    /** the inequality comparision operator
     *@param[in] ndit the iterator which is compared with the present iterator
     *@return false if they refer to the same node; true otherwise */
    bool operator!=(const NodeIterator& ndit) const{
        return (grh != ndit.grh || i != ndit.i);
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    Graph* grh;
    size_type i;

    /** Construct a valid NodeIterator. */
    NodeIterator(const Graph* g, size_type j) {
        grh = const_cast<Graph*>(g);
        i = j;
    }
    // HW1 #2: YOUR CODE HERE
  };

  /** @return the start of the node iterator*/
  node_iterator node_begin() const{
      return NodeIterator(this, 0);
  }

  /** @return the end of the node iterator*/
  node_iterator node_end() const{
      return NodeIterator(this, size_);
  }

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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }
    /** The dereference operator
     *@return the edge the iterator refers to*/
    Edge operator*() const{
        return grh->add_edge(grh->node(nd), grh->node(grh->adjs[nd][i]));
    }

    /** the increment operator
     *@return the iterator which refers to the subsequent edge*/
    IncidentIterator& operator++(){
        i++;
        return *this;
    }

    /** the equality comparision operator
     *@param[in] idit the iterator which is compared with the present iterator
     *@return true if they are based on the same node and refer to the
                 same edge; false otherwise */
    bool operator==(const IncidentIterator& idit) const{
        return (grh == idit.grh && nd == idit.nd && i == idit.i);
    }

    /** the inequality comparision operator
     *@param[in] idit the iterator which is compared with the present iterator
     *@return false if they are based on the same node and refer to the
                 same edge; true otherwise */
    bool operator!=(const IncidentIterator& idit) const{
        return (grh != idit.grh || nd != idit.nd || i != idit.i);
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;

    Graph* grh;
    size_type i;
    size_type nd;
    /** Construct a valid IncidentIterator. */
    IncidentIterator(Graph* g, size_type j, size_type nd_){
        grh = g;
        i = j;
        nd = nd_;
    }
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

    /** The dereference operator
     *@return the edge the iterator refers to*/
    Edge operator*() const{
        return grh->edge(i);
    }

    /** the increment operator
     *@return the iterator which refers to the subsequent edge*/
    EdgeIterator& operator++(){
        i++;
        return *this;
    }

    /** the equality comparision operator
     *@param[in] egit the iterator which is compared with the present iterator
     *@return true if they refer to the same edge; false otherwise */
    bool operator==(const EdgeIterator& egit) const{
        return (grh == egit.grh && i == egit.i);
    }

    /** the inequality comparision operator
     *@param[in] egit the iterator which is compared with the present iterator
     *@return false if they refer to the same edge; true otherwise */
    bool operator!=(const EdgeIterator& egit) const{
        return (grh != egit.grh || i != egit.i);
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    Graph* grh;
    size_type i;
    /** Construct a valid EdgeIterator. */
    EdgeIterator(const Graph* g, size_type j) {
        grh = const_cast<Graph*>(g);
        i = j;
    }
    // HW1 #5: YOUR CODE HERE
  };

  /** @return the start of the edge iterator*/
  edge_iterator edge_begin() const{
      return EdgeIterator(this, 0);
  }

  /** @return the end of the edge iterator*/
  edge_iterator edge_end() const{
      return EdgeIterator(this, num_edges_);
  }

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:
    struct internal_point{
        Point position;
        node_value_type value_;
        /** Construct internal_point*/
        internal_point(const Point& position, node_value_type value){
            this->position = position;
            this->value_ = value;
        }
    };
    // internal struct to store point information
    
    std::vector<internal_point> points;
    std::vector<std::vector<size_type>> adjs;
    std::vector<std::vector<size_type>> edges;
    // internal struct to store edge information
    size_type size_;
    size_type num_edges_;

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP


