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
template <typename V, typename E>
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
  using pair_type = std::pair<size_type,size_type>;
  /** Construct an empty graph. */
  // RK: Double check this initalization 
  Graph() : points_(), edges_() {
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
  class Node: private totally_ordered<Node>{
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
      return graph_->points_[node_index_];
    }

    Point& position() {
      // HW2: Modifiable position! for HW2
      return graph_->points_[node_index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_index_;
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
      // HW0: YOUR CODE HERE
      // Same graph and same index 
      return (graph_ == n.graph_ && n.node_index_== this->node_index_);
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
      // Same graph, check the node index
      if (graph_ == n.graph_){
        return this->node_index_ < n.node_index_;
      }
      // check the graph ordering
      else{
        return std::less<Graph*>{}(graph_, n.graph_);
      }
    }
    // Node value HW1:
    node_value_type& value(){
      return this->graph_->node_values_[node_index_];
    }

    // Const node value HW1:
    const node_value_type& value() const {
      return const_cast<node_value_type&>(this->graph_->node_values_[node_index_]);
    }

    // Return the number of edges incident on the node
    size_type degree() const {
      return graph_->adjacency_[node_index_].size();
    }

    /** start  of the  incident  iterator.
    * @return IncidentIterator to the first edge.
    *
    * Complexity: O(1).
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(node_index_, 0, graph_);
    }

    /** End of  incident  iterator.
    * @return IncidentIterator to the one after the last edge.
    *
    * Complexity: O(1).
    */
    incident_iterator edge_end() const {
      return IncidentIterator(node_index_, degree(), graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type node_index_;
    graph_type* graph_;
    // Constructor
    Node(const graph_type* graph, size_type node_index_)
        : node_index_(node_index_), graph_(const_cast<graph_type*>(graph)) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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

  Node add_node(const  Point& position, const node_value_type& value = node_value_type()){
    //HW 1: Add node with values
    points_.push_back(position);
    node_values_.push_back(value);

    std::vector<pair_type> empty_list;
    adjacency_.push_back(empty_list);
    return node(points_.size()-1); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_ == this && n.node_index_ < this->size()){
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //std::cout << "i: " << i << "num: " << num_nodes() << std::endl;
    //assert(i>=0 && i < num_nodes()); //pre condition
    return Node(this, i);        // Invalid node
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
      //std::cout << "node1 idx: " << node_index_1_ << "\n" << std::endl;
      return graph_->node(node_index_1_);   // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node_index_2_);      // Invalid Node
    }

    // HW2: Edge values:
    edge_value_type& value() {
      return this->graph_->edge_values_[edge_idx_];
    } 

    const edge_value_type& value() const{
      return const_cast<edge_value_type&>(this->graph_->edge_values_[edge_idx_]);
    }

    //length
    double length() const {
      return norm_2(node1().position() - node2().position());
    } 
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (this->graph_ == e.graph_) && (this->edge_idx_ == e.edge_idx_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // graph odering
      if (this->graph_ == e.graph_){
        return this->edge_idx_ < e.edge_idx_;
      }
      // check the graph ordering
      else{
        return std::less<Graph*>{}(this->graph_, e.graph_);
      }
    }

   size_type index() const{
      return edge_idx_;
   }
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type edge_idx_;

    size_type node_index_1_;
    size_type node_index_2_;
    // Constructor
    Edge(const graph_type* graph, size_type edge_idx, size_type node_1, size_type node_2)
      : graph_(const_cast<graph_type*>(graph)), edge_idx_(edge_idx), node_index_1_(node_1), node_index_2_(node_2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i>=0 && i < num_edges()); // pre condition
    return Edge(this, i, this->edges_[i].first, this->edges_[i].second);  
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  pair_type sort_nodes(const Node& a, const Node& b) const {
    size_type min_idx = a.node_index_;
    size_type max_idx = b.node_index_;
    if (min_idx > max_idx){ // swap if wrong
      min_idx = b.node_index_;
      max_idx = a.node_index_;
    }
    return pair_type(min_idx, max_idx);
  }

  bool has_edge(const Node& a, const Node& b) const {
    if (!(has_node(a) && has_node(b))){
      return false;
    }
    // HW0: YOUR CODE HERE
    for(auto incidence : adjacency_.at(a.index())){
        if (incidence.first == b.index()){
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // Both nodes be part og the graph
    assert(has_node(a) && has_node(b));
    assert(!(a==b));
    // sort the nodes
    pair_type this_edge = sort_nodes(a, b);
    // check for the existing edge
    for (std::size_t i=0; i < num_edges(); i++){
      if (this_edge == this->edges_[i]) return edge(i);
    }

    edges_.push_back(this_edge);
    edge_values_.push_back(value); // initaliza with the lenght
    // Update the adjacency list
    adjacency_.at(this_edge.first).push_back(pair_type(this_edge.second, edges_.size()-1));
    adjacency_.at(this_edge.second).push_back(pair_type(this_edge.first, edges_.size()-1));
    return edge(edges_.size()-1);
  }

  /** Remove an edge
   * @param[in] edge   The edge
   * @return           1 if the edge was removed, 0 otherwise. 
   *
   * @pre  @a edge valid. 
   * @post has_edge(@a edge.node1, @a edge.node2) == false 
   * @post EdgeIterators are invalidated.
   *
   * Complexity: O(max_degree). 
   */
  size_type remove_edge(const Edge& edge){
    // First check if the edge exist
    if (!has_edge(edge.node1(), edge.node2())){
          return 0;
    }
    size_type index = edge.index();

    // If the edge exist, pop the last one
    // replace removed edge with last one
    size_type n_edge = num_edges();
    this->edges_[index] = this->edges_[n_edge-1];
    this->edges_.pop_back();
    // Edge values:
    this->edge_values_[index] = this->edge_values_[n_edge-1];
    this->edge_values_.pop_back();
   
    // delete the index from adjcnecy list
    std::vector<pair_type>& incidence3 = adjacency_.at(edge.node1().index());
    for(auto it = incidence3.begin(); it != incidence3.end(); ++it){
      if (it->second == index){
        *it = incidence3.back();
        break;
      }
    }
    incidence3.pop_back();

    // delete the index from adjcnecy list
    std::vector<pair_type>& incidence4 = adjacency_.at(edge.node2().index());
    for(auto it = incidence4.begin(); it != incidence4.end(); ++it){
      if (it->second == index){
        *it = incidence4.back();
        break;
      }
    }
    incidence4.pop_back();

    // Update index of n_edge - 1 to index
    std::vector<pair_type>& incidence = adjacency_.at(edge.node1().index());
    for (auto it = incidence.begin(); it != incidence.end(); ++it){
      if (it->second == n_edge-1){
        it->second = index;
      }
    }
    // Same for second node
    std::vector<pair_type>& incidence2 = adjacency_.at(edge.node2().index());
    for (auto it = incidence2.begin(); it != incidence2.end(); ++it){
      if (it->second == n_edge-1){
        it->second = index;
      }
    }
    return 1;
    }

  /** Remove an edge
   * @param[in] node1   The first node
   * @param[in] node2   The second node
   * @return       1 if the edge was removed, 0 otherwise. 
   *
   * @pre  @a node1 and @a node2 are valid. 
   * @post has_edge(@a node1, @a node2) == false 
   * @post EdgeIterators are invalidated.
   *
   * Complexity: O(max_degree). 
   */
  size_type remove_edge(const Node& node1, const Node& node2){
    if (!has_edge(node1, node2)){
      return 0; 
    }

    pair_type this_edge = sort_nodes(node1, node2);
    std::vector<pair_type>& incidence = adjacency_.at(this_edge.first);
    size_type edge_idx = 0;
    for(auto it = incidence.begin(); it != incidence.end(); ++it){
      if (it->first == this_edge.second){
        edge_idx = it->second;
        break;
      }
    }
   
   this->remove_edge(Edge(this, edge_idx, this_edge.first, this_edge.second)); 
   return 0;
  }

  edge_iterator remove_edge(edge_iterator e_it) {
    this->remove_edge(*e_it);
    return e_it;
  }

  /** Remove a node
   * @param[in] node   The node to remove.
   * @return           1 if @a node was removed, 0 otherwise. 
   *
   * @pre  @a node is a valid.
   * @post has_node(@a n) == false 
   * @post All edges connecting @a n to any other node are removed.
   * @post NodeIterators and EdgeIterators are invalidated.
   * @post Incident iterator associated to @a n is invalidated. 
   *
   * Complexity: O(num_nodes()) assuming the graph is sparse. 
   */
  size_type remove_node(const Node& node){

      // First remove edges
      auto end_it = node.edge_end();
      for (auto it=node.edge_begin(); it != end_it; ++it){
         this->remove_edge(*it);
      }

      // size_type node_idx = node.index();
      // while (adjacency_[node_idx].size() > 0) {
      //    remove_edge(node, this->node(adjacency_[node_idx][0].first));
      // }

      //Delete node info
      size_type n_nodes = num_nodes();
      points_[node.index()] = points_[n_nodes-1];
      points_.pop_back();

      node_values_[node.index()] = node_values_[n_nodes-1];
      node_values_.pop_back();

      adjacency_[node.index()] = adjacency_[n_nodes-1];
      adjacency_.pop_back();

      // Updated indexes:
      size_type index = node.index();

      std::vector<pair_type>& incidence2 = adjacency_[index];
      for (auto it = incidence2.begin(); it != incidence2.end(); ++it){
         for(auto in_it = adjacency_[it->first].begin(); in_it != adjacency_[it->first].end(); ++in_it){
            if (in_it->first == n_nodes-1){
               in_it->first = index;
            }
            auto edge = edges_[in_it->second];
            if (edge.first == n_nodes-1){
               edge = sort_nodes(this->node(index), this->node(edge.second));
            }
            else if (edge.second == n_nodes-1){
               edge = sort_nodes(this->node(index), this->node(edge.first));
            }
         }
         
         auto edge = edges_[it->second];
         if (edge.first == n_nodes-1){
            edge = sort_nodes(this->node(index), this->node(edge.second));
         }
         else if (edge.second == n_nodes-1){
            edge = sort_nodes(this->node(index), this->node(edge.first));
         }
      }
      return 0;
  }

  node_iterator remove_node(node_iterator n_it){
     remove_node(*n_it);
     return n_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    edges_.clear();
    points_.clear();
    node_values_.clear();
    edge_values_.clear();
    adjacency_.clear();
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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return the element iterator is pointing to
     * @return node the iterator is pointing to
     * @pre valid iterator (node_idx_ is valid)
     * @post return a valid node
     * O(1) complexity
     */
    value_type operator*() const{
         //std::cout << "* operatopr idx: " << node_index_ << "\n" << std::endl;
        return this->graph_->node(node_index_);
    }

    /** Increement the pointer and return the iterator
    *
    * Complexity: O(1).
    */
    NodeIterator& operator++(){
      ++node_index_;
      return *this;
    }

    /** Check equality of this iterator and @a iter
    * @return boolean true/ false
    * @pre valid iterator @a iter
    * Complexity: O(1).
    */
    bool operator==(const NodeIterator& iter) const{
      return (this->graph_ == iter.graph_) && (this->node_index_ == iter.node_index_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_; 
    size_type node_index_;
    // private constructor
    NodeIterator(const graph_type* graph, size_type node_index) :
      graph_(const_cast<graph_type*>(graph)), node_index_(node_index) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** return the begin iterator 
  * @return  pointer to first node.
  *
  * Complexity: O(1).
  */
  NodeIterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** return the one after end iterator
  * @return NodeIterator to one after end
  *
  * Complexity: O(1).
  */
  NodeIterator node_end() const{
    return NodeIterator(this, this->num_nodes());
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return the iterator element
     *
     * @return node the iterator is pointing to
     * @pre valid iterator
     * @post return a valid value type
     *
     * Complexity: O(1).
     */
    value_type operator*() const {
      pair_type edge = graph_->adjacency_[incident_node_index_][incident_index_];
      return Edge(graph_, edge.second, incident_node_index_, edge.first);
    }

    /** Increment the pointer
    * @return iterator
    * @post valid iterator
    * Complexity: O(1).
    */
    IncidentIterator& operator++(){
      ++incident_index_;
      return *this;
    }

    /** Check equality with @a iter
    * @return bool
    * @pre valid iterator @a iter
    * Complexity: O(1).
    */
    bool operator==(const IncidentIterator& iter) const{
      bool cond_1 = this->graph_ == iter.graph_;
      bool cond_2 = this->incident_node_index_ == iter.incident_node_index_;
      bool cond_3 = this->incident_index_ == iter.incident_index_;
      return cond_1 && cond_2 && cond_3;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type incident_node_index_;
    size_type incident_index_;
    graph_type* graph_;
    IncidentIterator(size_type node_index, size_type index, const graph_type* graph): 
          incident_node_index_(node_index), incident_index_(index), graph_(const_cast<graph_type*>(graph)){

          }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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

    /** Return the edge that iterator points to.
    *
    * Complexity: O(1).
    */ 
    Edge operator*() const {
      pair_type this_edge = graph_->edges_[index_];
      return Edge(graph_, index_, this_edge.first, this_edge.second);
    }

    /** Increment and Return the edge iterator 
    *
    * Complexity: O(1).
    */
    edge_iterator& operator++() {
      ++index_;
      return *this;
    }

    /** Test whether edge iterator and @a itr are equal.
    *
    * Equal edge iterators have the same graph and index.
    *
    * Complexity: O(1).
    */
    bool operator==(const edge_iterator& itr) const {
      return ((this->graph_ == itr.graph_) && (this->index_ == itr.index_));
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;
    // constructor
    EdgeIterator(size_type index, const graph_type* graph): graph_(const_cast<graph_type*>(graph)), index_(index){
          }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 /** Return  an iterator to the first edge
   * @return an iterator to the first edge
   *
   * Complexity: O(1).
   */
  EdgeIterator edge_begin() const{
     return EdgeIterator(0, this);
  }
  
   /** Return the end of an iterator
     * @return an iterator pointing to after last edge
     *
     * Complexity: O(1).
     */
  EdgeIterator edge_end() const{
     return EdgeIterator(edges_.size(), this);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> points_;
  std::vector<pair_type> edges_;
  std::vector<std::vector<pair_type>> adjacency_;
  std::vector<node_value_type> node_values_;
  std::vector<edge_value_type> edge_values_;  
};

#endif // CME212_GRAPH_HPP
