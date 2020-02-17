#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * @collab Yue Li
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
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
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  // HW0: YOUR CODE HERE
  //Private Graph Constructor
    : nodes_(),  i2u_(), adj_list_(), edges_(), i2u_edges_() {}

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
      //Invalid node constructor takes no arguments
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_.at(uid_).position_;
    }

    /** Return this node's position. */
    Point& position() {
      return graph_->nodes_.at(uid_).position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      return graph_->nodes_.at(uid_).idx_;
    }


    /** Return this node's value, an argument of the Graph Template
    * This value can be updated
    * @return node value type of node
    */
    node_value_type& value(){
      return graph_->nodes_.at(uid_).node_value_;
    };

    /** Return this node's value as a const, an argument of the Graph Template
    * This value cannot be updated
    * @return node value type of node
    */
    const node_value_type& value() const{
      return graph_->nodes_.at(uid_).node_value_;
    };

    /** Return this node's degree,
    * The degree of a node is defined as the number of edges incident to a node
    * @return size type indicating degree of node
    */
    size_type degree() const{
      return graph_->adj_list_.at(uid_).size();
    };

    /** Return an incident iterator object indicating the beginning of an
    * iterator
    * @pre  edge_begin != edge_end
    * @return and incident iterator object
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, uid_,
        graph_->adj_list_.at(uid_).begin(),
        graph_->adj_list_.at(uid_).end());
    };

    /** Return an incident iterator object indicating the end of an
    * iterator
    * @return and incident iterator object
    */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_,
        graph_->adj_list_.at(uid_).end(),
        graph_->adj_list_.at(uid_).end());
    };

    /** Test whether this node and @a n are equal.
     * @param _n_ node to compare this node with
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      return graph_ == n.graph_ && uid_ == n.uid_
      && graph_->nodes_.at(uid_).idx_ == graph_->nodes_.at(n.uid_).idx_;
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
      if (graph_ == n.graph_){ //check if graphs are equal
        // return index ordering
        return graph_->nodes_.at(uid_).idx_ < graph_->nodes_.at(n.uid_).idx_;
      }
      else {
        //compare pointers to maintain trichotomy
        return std::less<Graph*>()(graph_, n.graph_);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_; //pointer to graph object
    size_type uid_; //unique identifier for node object


    /**Private Constructor */
    Node(const Graph* grph, size_type uid)
      : graph_(const_cast<Graph*>(grph)), uid_(uid){}


  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return i2u_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
    const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE

    size_type uid = nodes_.size();
    size_type idx = i2u_.size();
    i2u_.emplace_back(uid);
    nodes_.emplace_back(internal_node(idx, position, node_value));

    std::map<size_type, size_type> empty;
    adj_list_[uid] = empty;
    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    assert(n.uid_ < nodes_.size());
    return this == n.graph_ && nodes_.at(n.uid_).idx_ < i2u_.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < i2u_.size());
    return Node(this, i2u_.at(i));
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
      //Empty Constructor
    }

    /** Return a node of this Edge */
    Node node1() const {
      //Get first node from edge object
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //Get second node from edge object
      return Node(graph_, node2_uid_);
    }

    /** Return this node's value, an argument of the Graph Template
    * This value can be updated
    * @return node value type of node
    */
    edge_value_type& value(){
      return graph_->edges_.at(edge_uid_).edge_value_;
    };

    /** Return this node's value as a const, an argument of the Graph Template
    * This value cannot be updated
    * @return node value type of node
    */
    const edge_value_type& value() const{

      return graph_->edges_.at(edge_uid_).edge_value_;
    };

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((graph_ == e.graph_) &&
              (edge_uid_ == e.edge_uid_) &&
              ((node1_uid_ == e.node1_uid_
              && node2_uid_ == e.node2_uid_) ||
              (node1_uid_ == e.node2_uid_ &&
              node2_uid_ == e.node1_uid_)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_ ){ //Check if graphs are equal
        return graph_->edges_.at(edge_uid_).edge_idx_
        < graph_->edges_.at(e.edge_uid_).edge_idx_;
      }
      else {
        return std::less<Graph*>()(graph_, e.graph_); //compare pointers
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;

    size_type edge_uid_;
    size_type node1_uid_;
    size_type node2_uid_;


    /**Private Constructor */
    Edge(const Graph* grph, size_type edge_uid, size_type node1_uid, size_type
        node2_uid)
      : graph_(const_cast<Graph*>(grph)),
        edge_uid_(edge_uid),
        node1_uid_(node1_uid),
        node2_uid_(node2_uid)
        {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return i2u_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < i2u_edges_.size());
    size_type edge_uid = i2u_edges_.at(i);
    return Edge(this, edge_uid, edges_.at(edge_uid).node1_uid_,
                edges_.at(edge_uid).node2_uid_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (adj_list_.find(a.uid_) == adj_list_.end()) {
      //Did not find node A, return false
      return false;
    } else {
      //Found node A, return bool if finds node B
      return adj_list_.at(a.uid_).count(b.uid_);
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
  Edge add_edge(const Node& a, const Node& b,
     const edge_value_type& edge_value = edge_value_type()) {

    if (has_edge(a, b)){
      //If edge exists, call constructor with edge index
      return Edge(this, adj_list_.at(a.uid_).find(b.uid_)->second,
                  a.uid_, b.uid_);

    }

    else {
    //Add Pairing to vector of tuples
    size_type edge_uid = edges_.size();
    size_type edge_idx = i2u_edges_.size();
    edges_.push_back(internal_edge(edge_idx, a.uid_, b.uid_, edge_value));
    i2u_edges_.push_back(edge_uid);


    //Add to adjacency list
    adj_list_[a.uid_][b.uid_] = edge_uid;
    adj_list_[b.uid_][a.uid_] = edge_uid;

    return Edge(this, edge_uid, a.uid_, b.uid_);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    i2u_.clear();
    edges_.clear();
    i2u_edges_.clear();
    adj_list_.clear();

  }

/** Removes a given node from the graph
* @param[in] _n_ node to be removed from graph
* @post if node exists, num_nodes() < num_nodes() Before Removal - 1
* @return number of nodes removed from Graph
*/
size_type remove_node(const Node& n){
  size_type curr_node_idx = nodes_.at(n.uid_).idx_; //Get curr idx
  //Check if graph has node
  if (has_node(n)){
    //Extract adjacency list for node and iterate through
    std::map<size_type, size_type> curr_node_adj = adj_list_.at(n.uid_);                                                                                   adj_list_.find(n.uid_);
    for (auto it = curr_node_adj.begin(); it != curr_node_adj.end(); ++it){
        size_type node2_uid = (*it).first; //Get node 2
        Node n2 = Node(this, node2_uid);
        //remove edge for incident edge
        remove_edge(n, n2);
    }
    adj_list_.erase(adj_list_.find(n.uid_)); //Erase node from adj_list_

    //Swap and Pop to remove node from node data structures
    size_type end_uid = i2u_.back();
    i2u_.at(curr_node_idx) = end_uid;
    nodes_.at(end_uid).idx_ = curr_node_idx;
    i2u_.pop_back();
    return 1;
  }
  //If node not in graph, do not remove; return 0
  else{
  return 0;
}
}

/** Removes a given node from the graph
*
* @pre node iterator is valid
* @param[in] _n_it_ node iterator of a to be removed from graph
* @post If node exists,
                num_nodes() After Remove = num_nodes() Before Removal - 1
* @post iterator is valid                 
* @return node iterator pointing to new object n_its place
*/
node_iterator remove_node(node_iterator n_it){
//Dereference Iterator
Node curr_node = (*n_it);
//Call Remove Node
remove_node(curr_node);
//Return Iterator
return n_it;
}


/** Removes an edge from the graph based on two nodes
*
* @param[in] Node's _a_ and _b_ to remove the edge between
* @post If edge exists,
              num_edges() after remove < num_edges() Before Removal - 1
* @return number of edges removed {0,1}
*/
size_type remove_edge(const Node& a , const Node& b ){
  //Check if graph has edge
  if (has_edge(a,b)){
    //Construct edge object given node a and b
    Edge e =  Edge(this, adj_list_.at(a.uid_).find(b.uid_)->second,
                a.uid_, b.uid_);
    //Remove Edge object
    return remove_edge(e);
    }
    //If edge doesn't exist, return 0
    else {
      return 0;
    }
  }

/** Removes an edge from the graph based on an edge object
*
* @param[in] edge e to remove from graph
* @post If edge exists,
                    num_edges() after removal < num_edges() Before Removal - 1
* @return number of edges removed
*/
size_type remove_edge(const Edge& e){
  //Get curr edge idx
  size_type curr_idx = edges_.at(e.edge_uid_).edge_idx_;
  //Check if curr_idx is a valid idx in current graph
  if (has_edge(Node(this, e.node1_uid_), Node(this, e.node2_uid_))){

  //Get last edge
  size_type end_uid = i2u_edges_.back();

  //Remove edge from adjacency list
  (adj_list_.at(e.node1_uid_)).erase((adj_list_.at(e.node1_uid_)).
                                                find(e.node2_uid_));
  (adj_list_.at(e.node2_uid_)).erase((adj_list_.at(e.node2_uid_)).
                                                find(e.node1_uid_));
  //Swap with last edge
  i2u_edges_.at(curr_idx) = end_uid;
  //Update index of previous last edge
  edges_.at(end_uid).edge_idx_ = curr_idx;
  //Remove edge
  i2u_edges_.pop_back();
  return 1;
  }
  //If not valid edge, do not remove, return 0
  else {
    return 0;
  }
}

/** Removes an edge from the graph based on an edge iterator
*
* @pre iterator is valid
* @param[in] edge _e_it_ iterator pointing to edge to remove from graph
* @post If edge exists,
*            num_edges() After Remove < num_edges() Before Removal - 1
* @post iterator is valid
* @return number of edges removed
*/
edge_iterator remove_edge(edge_iterator e_it){
  //Dereference iterator to get edge object
  Edge e = *e_it;
  //Call Remove Edge on edge object
  remove_edge(e);
  //Return iterator pointing to current object
  return e_it;
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

    /** Dereference Operator
    * @pre iterator != end
    * @return a node object
    */
    Node operator*() const{
      //Get uid of node object
      size_type uid = graph_->i2u_.at(curr_idx_);
      //return constructed node object
      return Node(graph_, uid);
    }

    /** Increment Operator
    * @pre iterator not at end
    * @return a node iterator
    */
    NodeIterator& operator++(){
        ++curr_idx_;
        return (*this);
    }

    /** Equality Operator
    * @param _node_iter_ other node iterator to compare this with
    * @return boolean indicating equality
    */
    bool operator==(const NodeIterator& node_iter) const{
      return (curr_idx_ == node_iter.curr_idx_
      && graph_ == node_iter.graph_);
    }

    /** Inequality Iterator
    * @param _node_iter_ other node iterator
    * @return boolean indicating inequality
    */
    bool operator!=(const NodeIterator& node_iter) const{
      return !(*this == node_iter);
    }

   private:
    friend class Graph;

    Graph* graph_;
    size_type curr_idx_;

    //Node Iterator Constructor
    NodeIterator(const Graph* grph,
      size_type idx)
    				: graph_(const_cast<Graph*>(grph)), curr_idx_(idx){};
  };



  /** Return the beginning of a node iterator
  * @pre nodes_.size() >= 0 and i2u_.size >= 0
  * @return the beginning of a node iterator
  */
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }

  /**Return the end of a node iterator
  * @pre nodes.size() >= 0 and i2u_.size >= 0
  * @return the end of a node iterator based on number of active nodes
  */
  node_iterator node_end() const{
    return NodeIterator(this, i2u_.size());
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


    /** Dereference Operator for Incident Iterators
    * @pre iterator != end
    * @return an edge object
    */
    Edge operator*() const{
      return Edge(graph_, incident_iter_begin_->second,
                  node1_uid_, incident_iter_begin_->first);
    };

    /** Increment Operator for Incident Iterators
    * @pre iterator != end
    * @return an incident iterator object
    */
    IncidentIterator& operator++(){
      ++incident_iter_begin_;
      return(*this);
    };

    /** Equality Operator for Incident Iterators
    * @param _inc_iter_ other incident iterator
    * @return boolean indicating equality
    */
    bool operator==(const IncidentIterator& inc_iter) const{
      return (graph_ == inc_iter.graph_ &&
              node1_uid_ == inc_iter.node1_uid_ &&
              incident_iter_begin_ == inc_iter.incident_iter_begin_ &&
              incident_iter_end_ == inc_iter.incident_iter_end_);
    };

    /** Inequality Operator for Incident Iterators
    * @param _inc_iter_ other incident iterator
    * @return boolean indicating inequality
    */
    bool operator!=(const IncidentIterator& inc_iter) const{
      return !(*this == inc_iter);
    };

   private:
    friend class Graph;
    Graph* graph_;
    size_type node1_uid_;
    std::map<size_type, size_type>::iterator incident_iter_begin_;
    std::map<size_type, size_type>::iterator incident_iter_end_;

    //Incident Iterator Constructor
    IncidentIterator( const Graph* grph,
      size_type node_uid,
      std::map<size_type, size_type>::iterator incident_iter_begin,
      std::map<size_type, size_type>::iterator incident_iter_end)
						: graph_(const_cast<Graph*>(grph)),
              node1_uid_(node_uid),
              incident_iter_begin_(incident_iter_begin),
						  incident_iter_end_(incident_iter_end){};
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


    /** Dereference Operator for Edge Iterators
    * @pre iterator not at end
    * @return an edge object
    */
    Edge operator*() const{
      size_type edge_uid = graph_->i2u_edges_.at(curr_edge_idx_);
      return Edge(graph_, edge_uid,
                 graph_->edges_.at(edge_uid).node1_uid_,
                 graph_->edges_.at(edge_uid).node2_uid_);
    };


    /** Increment Operator for Edge Iterators
    * @pre iterator != end
    * @return a reference to an edge iterator object
    */
    EdgeIterator& operator++(){
      ++curr_edge_idx_;
      return (*this);
    };


    /** Equality Operator for Edge Iterators
    * @pre edge_iter valid iterator
    * @param _edge_iter_ other edge iterator
    * @return boolean indicating equality between iterators
    */
    bool operator==(const EdgeIterator& edge_iter) const{
      return (graph_ == edge_iter.graph_ &&
              curr_edge_idx_ == edge_iter.curr_edge_idx_);
    };


    /** Inequality Operator for Edge Iterators
    * @param _edge_iter_ other edge iterator
    * @pre _edge_iter_ valid iterator
    * @return boolean indicating inquality between iterators
    */
    bool operator!=(const EdgeIterator& edge_iter) const{
      return !(*this == edge_iter);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type curr_edge_idx_;


    //Edge iterator constructor
    EdgeIterator(const Graph* grph, size_type curr_edge_idx)
						: graph_(const_cast<Graph*>(grph)),
            curr_edge_idx_(curr_edge_idx){};
  };



  /** Returns beginning of an edge iterator
  * @pre num _edges_.size() >= 0 && i2u_edges_.size() >= 0
  * @return edge iterator at beginning of iterator
  */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** Returns end of an edge iterator
  * @pre _edges_.size() >= 0 && i2u_edges_.size() >= 0
  * @return edge iterator at end of iterator
  */
  edge_iterator edge_end() const{
    return EdgeIterator(this, i2u_edges_.size());
  }


 private:


  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.

  //Internal type for node elements
  struct internal_node{
    size_type idx_; //cur valid idx of node
    Point position_; // Position
    node_value_type node_value_; // node value


  //Internal Node Struct representation
  internal_node(const int idx, const Point& position,
     const node_value_type& node_value = node_value_type())
      : idx_(idx), position_(position), node_value_(node_value) {};

  };

  //Internal type for node elements
  struct internal_edge{
    size_type edge_idx_; //cur valid idx of edge
    size_type node1_uid_; //uid of node 1
    size_type node2_uid_; // uid of node 2
    edge_value_type edge_value_; // edge value



  //Internal Node Struct representation
  internal_edge(const int edge_idx, const int node1_uid,
    const int node2_uid,
     const edge_value_type& edge_value = edge_value_type())
      : edge_idx_(edge_idx), node1_uid_(node1_uid), node2_uid_(node2_uid),
      edge_value_(edge_value){};
  };

  std::vector<internal_node> nodes_; //indexed by uid
  std::vector<size_type> i2u_; // contains uid, indexed by idx


  //map of maps representing node relations in an adjacency list with
  //an additional map to edge uid. All indexing based on internal node
  // and edge uids
  std::map<size_type, std::map<size_type, size_type>> adj_list_;

  //vector of node pairs in an edge
  // std::vector<std::vector<size_type>> edges_;
  std::vector<internal_edge> edges_; //indexed by edge_uid_
  std::vector<size_type> i2u_edges_; //contained edge_uid_, indexed by edge_idx

};

#endif // CME212_GRAPH_HPP
