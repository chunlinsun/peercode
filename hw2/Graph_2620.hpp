#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>
#include <set>
#include <iterator>

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

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;
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

  using P = Point;
  using pair_size_type = std::pair<size_type,size_type>;
  using ptr_smaller = std::less<const graph_type*>;


  /** Structure that stores everything that characterizes a node
   * position_ : position of the node in space
   * value_ : value of the node
   * adj_ adjacency list from the node
   * idx user-side index of the node
   */ 
  struct NodeData{
    P position_;
    node_value_type value_;
    std::vector<std::pair<size_type,size_type>> adj_;
    size_type idx;
  };

  /** Structure that stores everything that characterizes an edge
   * nodes pair with head and tail nodes
   * value_ value of the edge
   */
  struct EdgeData{
    pair_size_type nodes;
    edge_value_type value_;
  };
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
    Node() {}

    /** Return this node's position. 
     *  @pre Undefined behaviour if this is call from an invalid node 
     */
    Point& position(){
      return graph_->elements_[uid_].position_;
    }

    /** Return this node's position as a constant. 
     *  @pre Undefined behaviour if this is call from an invalid node 
     */
    const Point& position() const {
      return graph_->elements_[uid_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). 
     *  @pre Undefined behaviour if this is call from an invalid node 
     */
    size_type index() const {
      return graph_->elements_[uid_].idx;
    }

    /**
     * @brief Allows to access and modify the value of a node
     * @return The value of the node by reference
     * @pre This is called from a valid node
     */
    node_value_type& value(){
      return graph_->elements_[uid_].value_;
    }

    /**
     * @brief Allows to access the value of a node without
     * modifying it
     * @return The value of the node as a constant reference
     * @pre This is called from a valid node
     */
    const node_value_type& value() const{
      return graph_->elements_[uid_].value_;
    }


    /**
     * @brief Getter function for the degree of a node
     * @return The degree of the node
     * @pre This is called from a valid node
     */
    size_type degree() const{
      return graph_->elements_[uid_].adj_.size();
    }

    /**
     * @return The begin iterator that allows iteration over the edges incident
     * to the node
     * @pre This is called from a valid node
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,uid_,graph_->elements_[uid_].adj_.cbegin());
    }

    /**
     * @return The end iterator that allows iteration over the edges incident
     * to the node
     * @pre This is called from a valid node
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_,uid_,graph_->elements_[uid_].adj_.cend());
    }

    /** Test whether this node and @a n are equal.
     *  @pre Both nodes are valid
     *  Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((this->graph_ == n.graph_) && (this->uid_ == n.uid_));
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * @pre Both nodes are valid
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      return (this->uid_< n.uid_) || ( (this->uid_ == n.uid_) && (ptr_smaller{}(this->graph_, n.graph_)));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    graph_type* graph_;
    size_type uid_;

    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

    // Checks if this node is valid
    bool valid() const
    {
      return graph_->has_node(*this);
    }


  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
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
    i2u_.push_back(elements_.size());
    elements_.push_back({position,value,{},size()-1});
    return Node(this,elements_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return ( n.graph_ == this && (n.index() < this->i2u_.size()) 
            && (n.uid_ < this->elements_.size())) && (i2u_[elements_[n.uid_].idx] == n.uid_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i2u_[i]);
  }

  /**
   * @brief Removes node _n_, 
   * Can detect nodes from other graphs or already deleted nodes
   * @param _n_ node to be removed
   * @pre _n_ is a valid node
   * @post Invalidates any existing NodeIterator at the end position
   * @post Might invalidate existing IncidentIterators
   * @post Might invalidate existing EdgeIterators
   * @post The node _n_ is invalidated, as is every edge incident to that node
   * @post Either return value is 0 and the graph didn't have that node
   * or return value is 1 and the node has been deleted
   *
   * Complexity: O(max_degree^2)
   */
  size_type remove_node ( const Node & n){
  
    if (!(this->has_node(n))) return 0;


    while (elements_[n.uid_].adj_.size() > 0){
      remove_edge(*(n.edge_begin()));
    }

    size_type index = n.index();

    std::swap(i2u_[index],i2u_[i2u_.size()-1]);
    elements_[i2u_[index]].idx = index;
    
    i2u_.pop_back();

    assert(!has_node(n));
    return 1;
  }

  /**
   * @brief Removes node pointed at by _n_it_
   * @pre _n_it_ doesn't point to the end
   * @post Invalidates any existing NodeIterator at the end position, therefore if iterating
   * you need to reinitialize the end iterator
   * @post Might invalidate existing IncidentIterators
   * @post Might invalidate existing EdgeIterators
   * @post The node pointed at is invalidated, as is every edge incident to that node
   * @post The node is deleted
   * @post _n_it_ points to an element unseen by the iterator
   *
   * Complexity: O(max_degree^2)
   */
  node_iterator remove_node(node_iterator n_it){
    remove_node(*n_it);
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Access value of edge 
     * @pre This edge is valid (hasn't been deleted)
     */
    edge_value_type& value(){
      return graph_->edges_[edge_index()].value_;
    }

    /** Access value of edge without modifying 
     * @pre This edge is valid (hasn't been deleted)
     */
    const edge_value_type& value() const{
      return graph_->edges_[edge_index()].value_;
    }


    /** Return a node of this Edge
     * @pre This edge is valid (hasn't been deleted)
     */
    Node node1() const {
      return Node(graph_,n1id_);   
    }

    /** Return the other node of this Edge
     * @pre This edge is valid (hasn't been deleted)
     */
    Node node2() const {
      return Node(graph_,n2id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->graph_ == e.graph_) 
              && (std::minmax(this->n1id_,this->n2id_) == std::minmax(e.n1id_,e.n2id_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (std::minmax(this->n1id_,this->n2id_) < std::minmax(e.n1id_,e.n2id_))
              || ( (std::minmax(this->n1id_,this->n2id_) == std::minmax(e.n1id_,e.n2id_)) 
                    && (ptr_smaller{}(this->graph_, e.graph_)));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    graph_type* graph_;
    size_type n1id_;
    size_type n2id_;

    Edge(const graph_type* graph, size_type n1id, size_type n2id)
      : graph_(const_cast<graph_type*>(graph)), n1id_(n1id), n2id_(n2id){
    }

    /** Access to index of edge 
     * @pre This edge is valid (hasn't been deleted)
     */
    size_type& edge_index(){
      auto it = graph_->elements_[node1().uid_].adj_.begin();
      while ( (node2().uid_!= (*it).first) && 
            (it != graph_->elements_[node1().uid_].adj_.end()) ){
        ++it;
      }

      assert(it != graph_->elements_[node1().uid_].adj_.end());
      return (*it).second;
    }

    /** Access to index of edge 
     * @pre This edge is valid (hasn't been deleted)
     */
    const size_type& edge_index() const{
      auto it = graph_->elements_[node1().uid_].adj_.begin();
      while ( (node2().uid_!= (*it).first) && 
            (it != graph_->elements_[node1().uid_].adj_.end()) ){
        ++it;
      }

      assert(it != graph_->elements_[node1().uid_].adj_.end());
      return (*it).second;
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
    return Edge(this,edges_[i].nodes.first,edges_[i].nodes.second); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto it = elements_[a.uid_].adj_.begin() ; it != elements_[a.uid_].adj_.end() ; ++it){
      if (b.uid_== (*it).first){
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type() ) {
    if (! has_edge(a,b) ){

      elements_[a.uid_].adj_.push_back({b.uid_,this->num_edges()});
      elements_[b.uid_].adj_.push_back({a.uid_,this->num_edges()});
      edges_.push_back({std::minmax(a.uid_,b.uid_),value});   

    }
    return Edge(this,a.uid_,b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_.clear();
    elements_.clear();
    i2u_.clear();

  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>{
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

    /**
     * @brief Dereference operator for NodeIterator
     * @pre This is called from a valid NodeIterator
     * @pre this iterator is not equal to the end iterator
     * @return The node where the iterator is pointing
     */
    value_type operator*() const{
      return graph_->node(index_);
    } 


    /**
     * @brief ++ operator for NodeIterator
     * @pre This is called from a valid NodeIterator
     * @pre this iterator is not equal to the end iterator
     * @return The updated iterator
     */
    NodeIterator& operator++(){
      index_++;
      return *this;
    }

    /**
     * @brief Test whether @ it and this iterator are equal
     * Two NodeIterators are equal if they have the same underlying graph
     * and are at the same index
     */
    bool operator==(const NodeIterator& it) const{
      return (this->graph_== it.graph_) && (this->index_ == it.index_);
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type index_;

    NodeIterator(const graph_type* graph, size_type index) :
      graph_(graph), index_(index){}


  };

  /**
   * @brief Provide begin NodeIterator to iterate over the nodes
   * of this graph
   */
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }

  /**
   * @brief Provide end NodeIterator to iterate over the nodes
   * of this graph
   */
  node_iterator node_end() const{
    return NodeIterator(this,this->num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>{
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

    /**
     * @brief Dereference operator for IncidentIterator
     * @pre This is called from a valid IncidentIterator
     * @pre this iterator is not equal to the end iterator
     * @post node1() of the returned Edge is always node
     * over whose incident Edges we are iterating
     * @return The node where the iterator is pointing
     */
    Edge operator*() const{
      return Edge(graph_, n1id_, (*it_).first);
    }

    /**
     * @brief ++ operator for IncidentIterator
     * @pre This is called from a valid IncidentIterator
     * @pre this iterator is not equal to the end iterator
     * @return The updated iterator
     */
    IncidentIterator& operator++(){
      it_++;
      return *this;
    }

    /**
     * @brief Test whether @ it2 and this iterator are equal
     * Two IncidentIterator are equal if they have the same underlying graph,
     * the same "origin" node and they are pointing to the same edge
     */
    bool operator==(const IncidentIterator& it2) const{
      return (this->graph_== it2.graph_) && (this->it_ == it2.it_)
             && (this->n1id_ == it2.n1id_);
    }


   private:
    friend class Graph;
    const graph_type* graph_;
    size_type n1id_;
    std::vector<std::pair<size_type,size_type>>::const_iterator it_;
    
    IncidentIterator(const graph_type* graph, size_type n1id,
                    std::vector<std::pair<size_type,size_type>>::const_iterator it) :
                    graph_(graph),n1id_(n1id),it_(it){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
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


    /**
     * @brief Dereference operator for EdgeIterator
     * @pre This is called from a valid EdgeIterator
     * @pre this iterator is not equal to the end iterator
     * @return The edge where the iterator is pointing
     */
    Edge operator*() const{
      return Edge(graph_, (graph_->edges_[index_]).nodes.first, (graph_->edges_[index_]).nodes.second);
    }

    /**
     * @brief ++ operator for EdgeIterator
     * @pre This is called from a valid EdgeIterator
     * @pre this iterator is not equal to the end iterator
     * @return The updated iterator
     */
    EdgeIterator& operator++(){
      index_++;
      return *this;
    }

    /**
     * @brief Test whether @ it2 and this iterator are equal
     * Two EdgeIterator are equal if they have the same underlying graph,
     * and they are pointing to the same edge
     */
    bool operator==(const EdgeIterator& it2) const{
      return (this->graph_ == it2.graph_) && (this->index_ == it2.index_);
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type index_;

    EdgeIterator(const graph_type* graph, size_type index):
                 graph_(graph),index_(index){}
  };

  /**
   * @brief Provide begin EdgeIterator to iterate over the edges
   * of this graph
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /**
   * @brief Provide end EdgeIterator to iterate over the edges
   * of this graph
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, edges_.size());
  }


  /**
   * @brief Removes edge between _a_ and _b_
   * @pre _a_ and _b_ are valid nodes from this graph
   * @post Might invalidate existing IncidentIterators
   * @post Invalidates any "end" EdgeIterator
   * @post Invalidates any copy of this edge
   * @post Either the return value is 0 and the edge didn't exist in the graph
   * or the return value is 1 and the edge has been deleted
   *
   * Complexity : O(max_degree)
   */

  size_type remove_edge(const Node& a, const Node& b){

    bool has_edge = false;
    size_type index{};

    // std::cout << a.uid_ << " " << b.uid_ << std::endl;


    for (auto it = elements_[a.uid_].adj_.begin();
         it != elements_[a.uid_].adj_.end(); ++it){
      if ((*it).first == b.uid_){
        index = (*it).second;
        std::iter_swap(it, elements_[a.uid_].adj_.end() - 1);
        elements_[a.uid_].adj_.pop_back();
        has_edge=true;
        break;
      }
    }

    if (!has_edge) return 0;


    for (auto it = elements_[b.uid_].adj_.begin();
         it != elements_[b.uid_].adj_.end(); ++it){
      if ((*it).first == a.uid_){
        std::iter_swap(it, elements_[b.uid_].adj_.end() - 1);
        elements_[b.uid_].adj_.pop_back();
        break;
      }
    }

    if (index != num_edges() - 1){
      std::swap(edges_[index],edges_[num_edges()-1]);
      pair_size_type ns = edges_[index].nodes;
      Edge(this,ns.first,ns.second).edge_index() = index;
      Edge(this,ns.second,ns.first).edge_index() = index;
    }

    edges_.pop_back();
    


    return 1;
  }


  /**
   * @brief Removes edge _e_
   * @pre _e_ is a valid edge
   * @post Might invalidate existing IncidentIterators
   * @post Invalidates any "end" EdgeIterator
   * @post Invalidates any copy of this edge
   * @post Either the return value is 0 and the edge didn't exist in the graph
   * or the return value is 1 and the edge has been deleted
   *
   * Complexity : O(max_degree)
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(),e.node2());
  }

  /**
   * @brief Removes edge pointed at by _e_it_ and returns an iterator pointing 
   * to an unseen element (or end iterator)
   * @pre _e_ is not the end iterator
   * @post Might invalidate existing IncidentIterators
   * @post Invalidates any "end" EdgeIterator
   * @post Invalidates the edge pointed at
   * @post The edge is deleted from the graph
   * @post The iterator points at an unseen element
   *
   * Complexity : O(max_degree)
   */
  edge_iterator remove_edge( edge_iterator e_it){
    remove_edge(*e_it);
    return e_it;
  }

 private:


  /**
   * elements_ stores each node as a NodeData structure that contains a position,
   * a value and an adjacency list for the node
   * edges_ stores edges without duplicates for easy iteration and access to i-th
   * i2u_ mapping between user-side index and unique id of nodes
   */
  std::vector<NodeData> elements_;
  std::vector<EdgeData> edges_;
  std::vector<size_type> i2u_;
};

#endif // CME212_GRAPH_HPP
