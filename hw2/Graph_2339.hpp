#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>
#include <unordered_map>

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

  // Internal Node class that stores value and position
  class Node_set;

  //Internal Node class that stores index and the number of nodes
  // with smaller index than it. (degree)
  class Node_helper;


  // Internal edge class that stores value and neighbor uid.
  class Edge_set;

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
  private:
    Node() {
    }

    Node(size_type uid, const graph_type* Graph_node){
      this->uid = uid;
      this->Graph_node = const_cast<graph_type*>(Graph_node);
    }
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
    

    /** Return this node's position. */
    const Point& position() const {
      return Graph_node->Node_sets[uid].p;
    }

    Point& position(){
    	return Graph_node->Node_sets[uid].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return Graph_node->uid_index[uid].index_val;
    }
    /**Return the node's unique uid.**/
    size_type get_uid() const{
    	return uid;
    }
    /**Return the value of the node**/
    node_value_type & value (){
      return Graph_node->Node_sets[uid].value;
    }
    /**Return the value of the node**/
    const node_value_type & value () const {
      return Graph_node->Node_sets[uid].value;
    }

    /**Return the degree of the node**/
    size_type degree() const{
      return Graph_node->adjacency_list[uid].size();
    }

    /**Return an incident iterator that points to the first edge of the graph.**/
    incident_iterator edge_begin() const{
      return IncidentIterator(this->index(), 0, this->Graph_node);
    }

    /**Return an incident iterator that points to the last edge of the graph.**/
    incident_iterator edge_end() const{
      return IncidentIterator(this->index(), degree(), this->Graph_node);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return this->Graph_node==n.Graph_node && this->uid==n.uid;
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
      return this->uid<n.uid;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;
    // index_val is the index of the node on the graph.
    size_type uid;
    // The graph that the node is on.
    graph_type* Graph_node;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return index_uid.size();
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
  Node add_node (const Point & position, const node_value_type & node_value = node_value_type ()){
    size_type uid = Node_sets.size();
    size_type index = size();
    this->Node_sets.push_back(Node_set(position, node_value));
    this->uid_index[uid] = Node_helper(index, 0);
    this->index_uid[index] = uid;
    this->Edge_sets.push_back(std::vector<Edge_set>(0));
    this->adjacency_list.insert({index, std::vector<size_type>(0)});
    return Node(uid, this);
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
    return n.Graph_node==this && uid_index.count(n.get_uid());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(index_uid.at(i), this);
  }


  /**Remove a node, first of all I remove the edges from adjacency_list
  Then renumber the indices of the other nodes because the uid never gets changed,
  but the index of nodes get changed. **/
  size_type remove_node(const Node& node){
  	size_type uid_1 = node.get_uid();
  	size_type index_1 = node.index();
  	if (has_node(node)){
		  for (auto edge = node.edge_begin(); edge!= node.edge_end(); ++edge){
        size_type uid_2 = (*edge).node2().get_uid();
        auto location_1= std::lower_bound(adjacency_list[uid_2].begin(), adjacency_list[uid_2].end(), uid_1);
    		adjacency_list[uid_2].erase(location_1);
    		if (uid_2>uid_1){
    			uid_index[uid_2].degree--;
    		}    		
    		number_edges--;
  		}
  		adjacency_list.erase(uid_1);

  		index_uid.clear();
  		uid_index.erase(uid_1);
  		/**Change the index of some nodes**/
  		for (auto it = uid_index.begin(); it!= uid_index.end(); ++it){
  			if (it->second.index_val>index_1){
  				it->second.index_val--;
  			}
  			index_uid.insert({it->second.index_val, it->first});
  		}
  		return 1;
  	}
  	return 0;
  }

  /**Remove a node, take a node iterator and make use of the function above.**/
  node_iterator remove_node(node_iterator n_it){
  	size_type return_type = remove_node(*n_it);
    // ++n_it;
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

    Edge(const graph_type* Graph_edge, size_type uid_1, size_type uid_2){
      this->uid_1 = uid_1;
      this->uid_2 = uid_2;
      this->Graph_edge = const_cast<graph_type*>(Graph_edge);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(uid_1, this->Graph_edge);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(uid_2, this->Graph_edge);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((this->uid_1== e.uid_1) && (this->uid_2== e.uid_2)\
       && (this->Graph_edge== e.Graph_edge)){
        return true;
      }
      if ((this->uid_1== e.uid_2) && (this->uid_2== e.uid_1)\
       && (this->Graph_edge== e.Graph_edge)){
        return true;
      }
      return false;
    }


    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
    	if (this->Graph_edge== e.Graph_edge){
    		if (this->uid_1!=e.uid_1){
        		return this->uid_1<e.uid_1;
      		}
      		else {
        		return this->uid_2<e.uid_2;
      		}
    	}
    	else{
    		return this->Graph_edge<e.Graph_edge;
    	}
      
    }

    /** Return the value of an edge**/
    edge_value_type& value(){
    	auto location = std::lower_bound(Graph_edge->Edge_sets[uid_1].begin(),\
    		Graph_edge->Edge_sets[uid_1].end(), \
    		Edge_set(uid_2));
    	return (*location).value;
    }
    /** Return the value of an edge**/
    const edge_value_type& value() const{
    	auto location = std::lower_bound(Graph_edge->Edge_sets[uid_1].begin(),
    		Graph_edge->Edge_sets[uid_1].end(), \
    		Edge_set(uid_2));
    	return (*location).value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type uid_1;
    size_type uid_2;
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
    while (i>=uid_index.at(index_uid.at(index_1)).degree){
      i-= uid_index.at(index_uid.at(index_1)).degree;
      index_1++;
    }
    return Edge(this, index_uid.at(index_1), adjacency_list.at(index_uid.at(index_1))[i]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
   * 
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type uid_1 = a.get_uid();
    size_type uid_2 = b.get_uid();
    auto location_1= std::lower_bound(adjacency_list.at(uid_1).begin(), adjacency_list.at(uid_1).end(), uid_2);
    auto location_2= std::upper_bound(adjacency_list.at(uid_1).begin(), adjacency_list.at(uid_1).end(), uid_2);
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
      return Edge(this, a.get_uid(), b.get_uid());
    }

    size_type uid_1 = a.get_uid();
    size_type uid_2 = b.get_uid();

    // Add neighbor uid to the adjacency_list. Adjacency_list get removed from time to time.
    //

    auto location_adj1 = std::lower_bound(adjacency_list[uid_1].begin(), adjacency_list[uid_1].end(), uid_2);
    auto location_adj2 = std::lower_bound(adjacency_list[uid_2].begin(), adjacency_list[uid_2].end(), uid_1);
    adjacency_list[uid_1].insert(location_adj1, uid_2);
    adjacency_list[uid_2].insert(location_adj2, uid_1);

    // Add edge neighbor uid and the edge value to Edge_sets. Edge_sets will not get 
    //removed after they get initialized.

    edge_value_type edge_value = norm(a.position()- b.position());
    auto location_edge1 = std::lower_bound(Edge_sets[uid_1].begin(), Edge_sets[uid_1].end(), Edge_set(uid_2));
    auto location_edge2 = std::lower_bound(Edge_sets[uid_2].begin(), Edge_sets[uid_2].end(), Edge_set(uid_1));
    Edge_sets[uid_1].insert(location_edge1, Edge_set(uid_2, edge_value));
    Edge_sets[uid_2].insert(location_edge2, Edge_set(uid_1, edge_value));

    if (uid_1>uid_2){
      uid_index[uid_1].degree++;
    }
    else{
      uid_index[uid_2].degree++;
    }   
    number_edges++;
    return Edge(this, uid_1, uid_2);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Node_sets.clear();
    Edge_sets.clear();
    index_uid.clear();
    adjacency_list.clear();
    number_edges = 0;
  }


  /** When I try to remove an edge, I only remove the index
  * from the adjacency list**/
  size_type remove_edge(const Node& a, const Node& b){
  	if (has_edge(a, b)){
  		size_type uid_1 = a.get_uid();
  		size_type uid_2 = b.get_uid();
  		auto location_1 = std::lower_bound(adjacency_list[uid_1].begin(), adjacency_list[uid_1].end(), uid_2);
    	auto location_2 = std::lower_bound(adjacency_list[uid_2].begin(), adjacency_list[uid_2].end(), uid_1);
    	adjacency_list[uid_1].erase(location_1);
    	adjacency_list[uid_2].erase(location_2);
    	if (uid_1>uid_2){
    		uid_index[uid_1].degree--;
    	}
    	else{
    		uid_index[uid_2].degree--;
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
    // ++e_it;
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
    // A valid node iterator constructor
    NodeIterator(size_type index, const graph_type* Graph_node) {
      this->index = index;
      this->Graph_node = const_cast<graph_type*>(Graph_node);
    }

    // Dereference operator and return the node that the iterator is pointing to.
    Node operator*() const{
    	// std::cout<<index<<" "<<Graph_node->index_uid.at(index)<<std::endl;
      return Node(Graph_node->index_uid[index], this->Graph_node);
    }

    // ++ operator and return the iterator pointing to the next node.
    NodeIterator& operator++(){
      if (index<Graph_node->size()){
        index++;
      }
      return *this;
    }

    // Check if two node iterators are the same by checking their indexes.
    bool operator==(const NodeIterator& another_iterator) const {
      return this->index == another_iterator.index;
    }

   private:
    friend class Graph;
    size_type index;
    graph_type* Graph_node;
  };

  // An iterator that points to the first node of the graph
  node_iterator node_begin() const{
    return NodeIterator(0, this);
  }

  // An end iterator that points to the place just after the last 
  // node of the graph.
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

    // An valid incident iterator constructor.
    IncidentIterator(size_type index_1, size_type adj_2, const graph_type* Graph_edge) {
      this->index_1 = index_1;
      this->adj_2 = adj_2;
      // graph_size = Graph_edge->size();
      this->Graph_edge = const_cast<graph_type*>(Graph_edge);
    }

    // Dereference operator and return the Edge that the incident iterator is pointing to 
    Edge operator*() const{
      return Edge(Graph_edge, Graph_edge->index_uid[index_1], \
      	Graph_edge->adjacency_list[Graph_edge->index_uid[index_1]][adj_2]);
    }

    // ++ operator and returns the next incident iterator.
    IncidentIterator& operator++(){
      if (adj_2<Graph_edge->adjacency_list[Graph_edge->index_uid[index_1]].size()){
        adj_2++;
      }
      return *this;
    }

    // check if the two incident iterator are the same
    bool operator==(const IncidentIterator& another_incident_iterator) const {
      return this->index_1== another_incident_iterator.index_1 && this->adj_2== another_incident_iterator.adj_2;
    }


   private:
    /** This helper function is not needed in NodeIterator, but rather heavily used in EdgeIterator.
    The usage is that when the edges of one nodes have already got traversed, then we are moving to the 
    next node.**/
    
    void update_next(){
      while (Graph_edge->index_uid.count(index_1) && adj_2==Graph_edge->uid_index[Graph_edge->index_uid[index_1]].degree ){
        index_1++;
        adj_2=0;
      }
    }
    friend class Graph;
    size_type index_1;
    size_type adj_2;
    graph_type* Graph_edge;
    // const size_type graph_size;
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

    // Constructor for edge_iterator and update the internal incident iterator
    // until I arrive at the first valid edge.
    EdgeIterator(size_type index_1, size_type adj_2, const graph_type* Graph_edge) {
      current_incident = IncidentIterator(index_1, adj_2, Graph_edge);
      // Find the first valid edge.
      current_incident.update_next();
    }

    // Use incident iterator's de reference operator to return the edge 
    // the incident iterator is pointing to.
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

  // Return the first edge of the graph.
  edge_iterator edge_begin() const{
  return EdgeIterator(0, 0, this);
  }

  // Return the pointer beyond the last edge of the graph.
  edge_iterator edge_end() const{
    return EdgeIterator(this->size(),0,this);
  }

 private:

  /** Key internal class for storing the identity (a unique 
  random number that is different for any two nodes) to 
  uniquely identity a node, and a Point class object p
  for storing the location of the node. 
  * p: the loation of the node
  * value: the value of a node, for example the 
  * degree: the number of nodes that share an edge with this
  node, but the index is smaller. This avoids a lot of double counting.*/
  class Node_set{
    Node_set(Point p, node_value_type value){
      this->p = p;
      this->value = value;
    }
    friend class Graph;
    Point p;
    node_value_type value;
  };


  /**A key internal class that stores the index and degree of a node. This Node_helper 
  class can get removed during the process.**/
  class Node_helper{
  public:
    Node_helper(){}
    Node_helper(size_type index_val, size_type degree){
      this->index_val = index_val;
      this->degree= degree;
    }
    friend class Graph;
    size_type index_val;
    size_type degree;
  };



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
  /** Node_sets are a vector of Node_set, where Node_set is 
  an internal class that stores the identity and position of the node.
  Node_sets[i], i is the ith node in terms of index on the graph.**/
  std::vector<Node_set> Node_sets;

  /** The map from uid to index**/
  std::unordered_map<size_type, Node_helper> uid_index;

  /** The map from index to uid.*/
  std::unordered_map<size_type, size_type> index_uid;

  /**adjacency_list stores the edges that are still on the graph.
	That is, Edge_sets store all of the available edges on the graph
	 and will never get deleted, but the edges in adjacency_list will
	 not get deleted. However, adjacency_list only store the neighbor
	node uid and not the edge value.
  **/
  std::unordered_map<size_type, std::vector<size_type>> adjacency_list;
  /**
  Edge_sets are a vector of Edge_set, where Edge_set is an internal class
  that stores the neighbor uid of neighbor index and the value of the edge. Edge_sets[i][j] is the ith 
  node in terms of index on the graph, and j is the jth neighbor of the 
  ith node. 
  **/
  std::vector<std::vector<Edge_set>> Edge_sets;

  /** Take Kyle's suggestions to make num_edges() faster to run.**/
  size_type number_edges = 0;
};

#endif // CME212_GRAPH_HPP



