#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
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
      return graph_->points_[idx_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /* @return This node's value. */
    node_value_type& value(){
      return graph_->values_[idx_];
    }
    /* @return This node's value. */
    const node_value_type& value() const{
      return graph_->values_[idx_];
    }
    
    /* @return The number of incident edges.*/
    size_type degree() const{
      return graph_->ma_edges_[idx_].size();    
    }

    /* @return An iterator pointing 
     * at the start of the edges incident to the node 
     */
    incident_iterator edge_begin() const{
      return incident_iterator(graph_,idx_,true);
    }

    /* @return An iterator pointing 
     * at the end of the edges incident to the node
     */
    incident_iterator edge_end() const{
      return incident_iterator(graph_,idx_,false);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ((n.graph_ == graph_) && (n.idx_ == idx_))
        return true;
      else
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
      if (idx_ < n.idx_)
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type idx_;
    Node(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx){
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return points_.size();
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
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()){
    points_.push_back(position);
    values_.push_back(value);
    return Node(this, points_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if ((n.idx_ < points_.size()) && (n.graph_ == this)) 
      return true;
    else 
      return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i >= 0);
    assert(i < num_nodes());
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graphE_->node(idx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graphE_->node(idx2_);
    }
//--functionality_0
//--in an undirected graph, the same edge is determined by having both of the same nodes. Your checking direction. You should also check the graph they belong to. 
//--START
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      if ((e.node1() == node1()) && (e.node2() == node2()))
        return true;
      else 
        return false;

    }
//--END

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (idxE_ < e.idxE_)
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graphE_;
    size_type idx1_;
    size_type idx2_;
    size_type idxE_;
    Edge(const Graph* graph, size_type idx1, size_type idx2, size_type idx)
      : graphE_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2),idxE_(idx){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >= 0);
    assert(i < edges_.size());
    return Edge(this, edges_[i][0], edges_[i][1], i);      
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const{
    assert(has_node(a));
    assert(has_node(b));

    auto iter = ma_edges_.find(a.index());

    if (iter == ma_edges_.end()){
      return false;
    }
    else if (iter->second.find(b.index())== iter->second.end()){
      return false;
    }
    return true;
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
    assert(has_node(a));
    assert(has_node(b));
    assert(a.index()!=b.index());

    if (has_edge(a,b)){
      size_type temp_idx = ma_edges_[a.index()][b.index()];
      return Edge(this, a.index(), b.index(), temp_idx);
    }

    std::vector<size_type> temp = {a.index(),b.index()};
    edges_.push_back(temp);

    ma_edges_[a.index()][b.index()]=edges_.size()-1;
    ma_edges_[b.index()][a.index()]=edges_.size()-1;
    
    return Edge(this, a.index(), b.index(), edges_.size()-1);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
//--functionality_1
//--need to update your clear function
//--START
  void clear() {
    points_.clear();
    edges_.clear();
  }
//--END
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

    /* @brief Dereference operator
     * @return An node pointed by the iterator
     */
    Node operator*() const{
      return itr_graph_->node(itr_idx_);  
    }

    /* @brief Increment to the next node
     * @return A node iterator pointing to the next node
     * @post The iterator points to the next node
     */
    NodeIterator& operator++(){  
      itr_idx_++;
      return *this;
    }
    
    /* @brief Defines equality between two iterators
     * @return True if two iterators points to the same node
     */
    bool operator==(const NodeIterator& node_itr) const{
      return ((itr_graph_==node_itr.itr_graph_)
               &&(itr_idx_ == node_itr.itr_idx_));
    }

    /* @brief Defines inequality between two iterators
     * @return True if two iterators don't point to the same node
     */
    bool operator!=(const NodeIterator& node_itr) const{
      return ((itr_graph_!=node_itr.itr_graph_)
              ||(itr_idx_ != node_itr.itr_idx_));
    }

   private:
    friend class Graph;
    Graph* itr_graph_;
    size_type itr_idx_;
    NodeIterator(const Graph* itr_graph, size_type itr_idx)
      : itr_graph_(const_cast<Graph*>(itr_graph)),itr_idx_(itr_idx){
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /* 
   * @return An iterator that points at the start of the nodes
   */
  node_iterator node_begin() const{
    return node_iterator(this, 0); 
  }

  /*
   * @return An iterator that points at end end of the nodes
   */
  node_iterator node_end() const{
    return node_iterator(this, num_nodes());
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

    /* @brief Dereference operator
     * @return An edge pointed by the iterator
     */
    Edge operator*() const{
      size_type node2_idx_ = map_itr_->first;
      size_type edge_idx_ = map_itr_->second;
      return Edge(ii_graph_,node_idx_,node2_idx_,edge_idx_);
    }

    /* @brief Increment to the next edge incident to the node
     * @return An edge iterator pointing to the next edge incident to the node
     * @post The iterator points to the next edge incident to the node
     */
    IncidentIterator& operator++(){
      map_itr_++;
      return *this;
    }

    /* @brief Defines equality between two iterators
     * @return True if two iterators points to the same edge 
     */
    bool operator==(const IncidentIterator& inc_itr) const{
      return ((ii_graph_ == inc_itr.ii_graph_)
	      &&(node_idx_ == inc_itr.node_idx_)
	      &&(map_itr_ == inc_itr.map_itr_));
    }

    /* @brief Defines inequality between two iterators
     * @return True if two iterators don't point to the same edge
     */
    bool operator!=(const IncidentIterator& inc_itr) const{
      return ((ii_graph_ != inc_itr.ii_graph_)
              ||(node_idx_ != inc_itr.node_idx_)
              ||(map_itr_ != inc_itr.map_itr_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* ii_graph_;
    size_type node_idx_;
    std::map<size_type,size_type>::iterator map_itr_; 
    IncidentIterator(const Graph* ii_graph, size_type node_idx, bool begin_)
      : ii_graph_(const_cast<Graph*>(ii_graph)),node_idx_(node_idx){
      if (begin_ == true){
	map_itr_ = ii_graph_->ma_edges_[node_idx].begin();
      }
      else{
      map_itr_ = ii_graph_->ma_edges_[node_idx].end();
      } 
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

    /* @brief Dereference operator
     * @return An edge pointed by the iterator
     */
    Edge operator*() const{
      return ei_graph_->edge(edge_idx_);  
    }

    /* @brief Increment to the next edge
     * @return An edge iterator pointing to the next edge
     * @post The iterator points to the next edge 
     */
    EdgeIterator& operator++(){
      edge_idx_++;
      return *this; 
    }

    /* @brief Defines equality between two iterators 
     * @return True if two iterators points to the same edge
     */
    bool operator==(const EdgeIterator& ei) const{
      return ((ei_graph_== ei.ei_graph_)
              &&(edge_idx_ == ei.edge_idx_));
    }

    /* @brief Defines inequality between two iterators 
     * @return True if two iterators don't point to the same edge 
     */
    bool operator!=(const EdgeIterator& ei) const{
      return ((ei_graph_!= ei.ei_graph_)
              ||(edge_idx_ != ei.edge_idx_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* ei_graph_;
    size_type edge_idx_;
    EdgeIterator(const Graph* ei_graph, size_type edge_idx)
      : ei_graph_(const_cast<Graph*>(ei_graph)),edge_idx_(edge_idx){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /* @return An iterator pointing at the start of the edges */
  edge_iterator edge_begin() const{
    return edge_iterator(this, 0);
  }

  /* @return An iterator pointing at the end of the edges */
  edge_iterator edge_end() const{
    return edge_iterator(this, num_edges());
  }


 private:
  std::vector<Point> points_;
  std::vector<std::vector<size_type>> edges_;
  std::map<size_type,std::map<size_type,size_type>> ma_edges_;
  std::vector<node_value_type> values_;
};

#endif // CME212_GRAPH_HPP
