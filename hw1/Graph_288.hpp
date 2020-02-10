#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template < typename V, typename E = int >
class Graph {

 public:
  // using node_value_type = V ;
 	typedef V node_value_type;
 	typedef E edge_value_type;

  //
  // PUBLIC TYPE DEFINITIONS
  //

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
private:
//--style_0
//--internal classes should be declared up here, but defined at the bottom
//--of the file.
//--START
  /** Key internal class for storing the identity (a unique 
  random number that is different for any two nodes) to 
  uniquely identity a node, and a Point class object p
  for storing the location of the node. 
  * p: the loation of the node
  * value: the value of a node, for example the 
  * degree: the number of nodes that share an edge with this
  node, but the index is smaller. This avoids a lot of double counting.*/
  class Node_set{
    friend class Graph;
    Point p;
    node_value_type value;
    size_type degree;
    Node_set(Point p, node_value_type value, size_type degree){
      this->p = p;
      this->value = value;
      this->degree = degree;
    }
  };
//--END

  /**
  Key internal class for storing the index of the jth neighbor of the
  ith node. I use a proxy class because proxy class can be easily
  extended for more complex purposes. 
  */
  class Edge_set{  
  public:
    Edge_set(){}
    Edge_set(size_type neighbor_index, edge_value_type value = edge_value_type ()){
      this->neighbor_index = neighbor_index;
      this->value = value;
    }
    bool operator<(const Edge_set& another_edge_set) const{
    	return this->neighbor_index< another_edge_set.neighbor_index;
    }
    bool operator==(const Edge_set& another_edge_set) const{
    	return this->neighbor_index == another_edge_set.neighbor_index;
    }
  private:
    friend class Graph;
    size_type neighbor_index;
    edge_value_type value;
  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
public:
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
    }
//--style_0
//--these constructors should be private, since the only way for a user to
//--create a valid node should be with Graph::node (or dereferincing a NodeIterator).
//--START
    Node(size_type index_val, const graph_type* Graph_node){
      this->index_val = index_val;
      this->Graph_node = const_cast<graph_type*>(Graph_node);
    }
//--END

    /** Return this node's position. */
    const Point& position() const {
      return Graph_node->Node_sets[index_val].p;
    }

    Point& position(){
    	return Graph_node->Node_sets[index_val].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_val;
    }

    node_value_type & value (){
      return Graph_node->Node_sets[index_val].value;
    }
    const node_value_type & value () const {
      return Graph_node->Node_sets[index_val].value;
    }
    size_type degree() const{
      return Graph_node->Edge_sets[index_val].size();
    }
    incident_iterator edge_begin() const{
      return IncidentIterator(this->index(), 0, this->Graph_node);
    }
    incident_iterator edge_end() const{
      return IncidentIterator(this->index(), degree(), this->Graph_node);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return this->Graph_node==n.Graph_node && this->index_val==n.index_val;
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
      return this->index_val<n.index_val;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;
    // index_val is the index of the node on the graph.
    size_type index_val;
    // The graph that the node is on.
    graph_type* Graph_node;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return Node_sets.size();
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
  Node add_node ( const Point & position, const node_value_type & node_value = node_value_type ()){
    size_type index = this->size();
    this->Node_sets.push_back(Node_set(position, node_value, 0));
    this->Edge_sets.push_back(std::vector<Edge_set>(0));
    return Node(index, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   * A node is on the graph if and only the graph of the node
   * is this graph, and the index of the node is smaller than
   * the size of the graph.
   */
  bool has_node(const Node& n) const {
    return n.Graph_node==this && n.index()< size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(i, this);
  }

  size_type remove_node(const Node& node){
  	size_type return_type = 1;
  	size_type index_1 = node.index();
  	if (has_node(node)){
		for (auto edge = node.edge_begin(); edge!= node.edge_end(); ++edge){
			size_type index_2 = (*edge).node2().index();
			auto location_2 = std::lower_bound(Edge_sets[index_2].begin(), Edge_sets[index_2].end(), Edge_set(index_1));
    		Edge_sets[index_2].erase(location_2);
    		number_edges--;
  		}
  		// Right now this looks too complicated. Will revise in the future. 
  		for (size_type i = 0; i<Edge_sets.size(); i++){
  			for (size_type j = 0; j<Edge_sets[i].size(); j++){
  				if (Edge_sets[i][j].neighbor_index>index_1){
  					Edge_sets[i][j].neighbor_index--;
  				}
  			}
  		}
  		Edge_sets.erase(Edge_sets.begin()+node.index());
  		Node_sets.erase(Node_sets.begin()+node.index());
  		return return_type;
  	}
  	return 0;
  }

  node_iterator remove_node(node_iterator n_it){
  	size_type return_type = remove_node(*n_it);
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    Edge(const graph_type* Graph_edge, size_type index_1, size_type index_2){
      this->index_1 = index_1;
      this->index_2 = index_2;
      this->Graph_edge = const_cast<graph_type*>(Graph_edge);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(index_1, this->Graph_edge);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(index_2, this->Graph_edge);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
//--functionality_1
//--needs to consider that index1 and index2 could be in opposite order.
//--for example your add_edge function does not check the order when the
//--edge to add already exists.
//--START
    bool operator==(const Edge& e) const {
      return (this->index_1== e.index_1) && (this->index_2== e.index_2)\
       && (this->Graph_edge== e.Graph_edge);
    }
//--END

    // double length() const{
    // 	return norm(node1().position()- node2().position());
    // }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
    	if (this->Graph_edge== e.Graph_edge){
    		if (this->index_1!=e.index_1){
        		return this->index_1<e.index_1;
      		}
      		else {
        		return this->index_2<e.index_2;
      		}
    	}
    	else{
    		return this->Graph_edge<e.Graph_edge;
    	}
      
    }

    edge_value_type& value(){
    	auto location = std::lower_bound(Graph_edge->Edge_sets[index_1].begin(),Graph_edge->Edge_sets[index_1].end(), 
    		Edge_set(index_2));
    	return (*location).value;
    }

    const edge_value_type& value() const{
    	auto location = std::lower_bound(Graph_edge->Edge_sets[index_1].begin(),Graph_edge->Edge_sets[index_1].end(), 
    		Edge_set(index_2));
    	return (*location).value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type index_1;
    size_type index_2;
    graph_type* Graph_edge;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return number_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
   * First of all, find the first node index j of the ith edge, and then 
   * find which neighbor the ith edge corresponds to the jth node. After
   * find the index, return the edge. 
   */
  Edge edge(size_type i) const {
    size_type index_1 = 0;
    while (i>=Node_sets[index_1].degree){
      i-= Node_sets[index_1].degree;
      index_1++;
    }
    return Edge(this, index_1,Edge_sets[index_1][i].neighbor_index);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
   * 
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type index_1 = a.index();
    size_type index_2 = b.index();
    auto location_1= std::lower_bound(Edge_sets[index_1].begin(), Edge_sets[index_1].end(), Edge_set(index_2));
    auto location_2= std::upper_bound(Edge_sets[index_1].begin(), Edge_sets[index_1].end(), Edge_set(index_2));
    return location_1!= location_2;
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
   * When I add edges, I make sure that I did not double count the edges. 
   * In Edge_sets, I make sure that Edge_sets[i], which is the neighbors of 
   * the ith node, the index of all of its neighbors are larger than i. 
   *
   * Take Kyle's suggestion into consideration: Has to determine whether 
   * has_edge first.
   */
  Edge add_edge(const Node& a, const Node& b) {
    if(has_edge(a, b)){
      return Edge(this, a.index(), b.index());
    }
    size_type index_1 = a.index();
    size_type index_2 = b.index();
    edge_value_type edge_value = norm(a.position()- b.position());
    /** I want to make sure that the wat that I insert the edge are
    in order. **/
    auto location_1 = std::lower_bound(Edge_sets[index_1].begin(), Edge_sets[index_1].end(), Edge_set(index_2));
    auto location_2 = std::lower_bound(Edge_sets[index_2].begin(), Edge_sets[index_2].end(), Edge_set(index_1));
    Edge_sets[index_1].insert(location_1, Edge_set(index_2, edge_value));
    Edge_sets[index_2].insert(location_2, Edge_set(index_1, edge_value));
    if (index_1>index_2){
      Node_sets[index_1].degree++;
    }
    else{
      Node_sets[index_2].degree++;
    }   
    number_edges++;
    return Edge(this, index_1, index_2);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Node_sets.clear();
    Edge_sets.clear();
    number_edges = 0;
  }

  size_type remove_edge(const Node& a, const Node& b){
  	if (has_edge(a, b)){
  		size_type index_1 = a.index();
  		size_type index_2 = b.index();
  		auto location_1 = std::lower_bound(Edge_sets[index_1].begin(), Edge_sets[index_1].end(), Edge_set(index_2));
    	auto location_2 = std::lower_bound(Edge_sets[index_2].begin(), Edge_sets[index_2].end(), Edge_set(index_1));
    	Edge_sets[index_1].erase(location_1);
    	Edge_sets[index_2].erase(location_2);
    	if (index_1>index_2){
    		Node_sets[index_1].degree--;
    	}
    	else{
    		Node_sets[index_2].degree--;
    	}
    	number_edges--;
    	return 1;
  	}
  	std::cout<<"The edge between a and b does not exist"<<std::endl;
  	return 0; 
  }

  size_type remove_edge(const Edge& edge){
  	Node node1 = edge.node1();
  	Node node2 = edge.node2();
  	return remove_edge(node1, node2);
  }

  edge_iterator remove_edge(edge_iterator e_it){
  	size_type edge_bool = remove_edge(*e_it);
  	return e_it;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
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

    NodeIterator(size_type index, const graph_type* Graph_node) {
      this->index = index;
      this->Graph_node = const_cast<graph_type*>(Graph_node);
    }
//--documentation_2
//--almost no documentation for HW1 Graph methods, and not in doxygen style
//--START
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const{
      return Node(index, this->Graph_node);
    }
    NodeIterator& operator++(){
      if (index<Graph_node->size()){
        index++;
      }
      return *this;
    }
    bool operator==(const NodeIterator& another_iterator) const {
      return this->index == another_iterator.index;
    }
//--END

   private:
    friend class Graph;
    size_type index;
    graph_type* Graph_node;
  };

  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const{
    return NodeIterator(0, this);
  }
  node_iterator node_end() const{
    return NodeIterator(size(), this);
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

    IncidentIterator(size_type index_1, size_type adj_2, const graph_type* Graph_edge) {
      this->index_1 = index_1;
      this->adj_2 = adj_2;
      this->Graph_edge = const_cast<graph_type*>(Graph_edge);
    }

    Edge operator*() const{
      return Edge(Graph_edge, index_1, Graph_edge->Edge_sets[index_1][adj_2].neighbor_index);
    }
    IncidentIterator& operator++(){
      if (adj_2<Graph_edge->Edge_sets[index_1].size()){
        adj_2++;
      }
      return *this;
    }
    bool operator==(const IncidentIterator& another_incident_iterator) const {
      return this->index_1== another_incident_iterator.index_1 && this->adj_2== another_incident_iterator.adj_2;
    }
//--style_1
//--It would make more sense if this function and index_1 were a part of EdgeIterator. 
//--But as is, this function should not be public since the user should not be using it.
//--START
    /** This helper function is not needed in NodeIterator, but rather heavily used in EdgeIterator.
    The usage is that when the edges of one nodes have already got traversed, then we are moving to the 
    next node.**/
    void update_next(){
      while (adj_2==Graph_edge->Node_sets[index_1].degree && index_1<Graph_edge->size()){
        index_1++;
        adj_2=0;
      }
    }
//--END

   private:
    friend class Graph;
    size_type index_1;
    size_type adj_2;
    graph_type* Graph_edge;
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

    EdgeIterator(size_type index_1, size_type adj_2, const graph_type* Graph_edge) {
      current_incident = IncidentIterator(index_1, adj_2, Graph_edge);
      // Find the first valid edge.
      current_incident.update_next();
    }

    Edge operator*() const{
      return *current_incident;
    }

    /** How do I guarantee that each edge only get iterated once?
    I guarantee this because for each index_1, Node_sets[index_1].degree stores
    the number of edges incident to index_1 whose index_2 (the index of node2)
    is smaller than index_1. This way, I only iterate up until some edges of index_1
    and then I will still iterating over the incident edges of index_1.**/
    EdgeIterator& operator++(){
      ++current_incident;
      current_incident.update_next();
      return *this;
    }
    /**Two edge iterators are equal if and only if their incident iterators are equal **/
    bool operator==(const EdgeIterator& another_edge_iterator) const{
        return this->current_incident== another_edge_iterator.current_incident;
    }

   private:
    friend class Graph;
    friend class IncidentIterator; 
    /** I implement the edge_iterator based on the incident_iterator. IncidentIterator provides me 
    with an easy way to implement the EdgeIterator. I don't need other data variables at all!**/
    IncidentIterator current_incident;
  };

  edge_iterator edge_begin() const{
  return EdgeIterator(0, 0, this);
  }
  edge_iterator edge_end() const{
    return EdgeIterator(this->size(),0,this);
  }

 private:
  /** Node_sets are a vector of Node_set, where Node_set is 
  an internal class that stores the identity and position of the node.
  Node_sets[i], i is the ith node in terms of index on the graph.**/
  std::vector<Node_set> Node_sets;
  /**
  Edge_sets are a vector of Edge_set, where Edge_set is an internal class
  that stores the neighbor index of a node. Edge_sets[i][j] is the ith 
  node in terms of index on the graph, and j is the jth neighbor of the 
  ith node. 
  **/
  std::vector<std::vector<Edge_set>> Edge_sets;

  /** Take Kyle's suggestions to make num_edges() faster to run.**/
  size_type number_edges = 0;
};

#endif // CME212_GRAPH_HPP



