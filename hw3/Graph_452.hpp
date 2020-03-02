#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
#include <functional>
#include<tuple>
#include <unordered_map>
#include <unordered_set>

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
  // PUBLIC TYPE DEFINITIONS -- added in HW0
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  // HW1
  using node_value_type = V;
  class Node;
  // HW1
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  // HW2
  using edge_value_type = E;
  class Edge;
  // HW2
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
	: internal_nodes_(), internal_edges_(), size_(0), next_node_idx_(0), next_edge_idx_(0), neighbors_(), 
	internal_node_val_() , internal_edge_val_(), idx_(), i2u_(), degrees_(){
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
	  uid_ = -1;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
	  return * (*graph_).get_node_position(uid_);
    }
	
	// HW2 - modifiable Node position
	Point& position() {
	  return * (*graph_).get_node_position(uid_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {								//added the 1st const
      // HW0: YOUR CODE HERE
	  return (*graph_).get_node_idx(uid_);				//HW2
    }
	

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
	
	/** Return the node's value. */
    node_value_type& value() {
		return  (*graph_).get_node_val(uid_);			//HW2
	}
	
	/** Return the node's value. */
    const node_value_type& value() const {
		return graph_.get_node_val(uid_);
	}
	
	/** Test whether  this node is active on the graph
	 * 
	 * Complexity O(1)
	 */
	 bool valid() const {
		 return uid_>= 0 && graph_->real_size()				// uid_ in range
					&& index() < graph_->active_size()			// idx in range
					&& graph_->i2u_[index()] == uid_ ;			//uid_ in sync
	 }

	//test
	void uid() const {
		std::cout << "uid: " << uid_  << std::endl;
		std::cout << "ind: " << index() << std::endl;
	}
	size_type uid_val() const {
		return uid_;
	}

	/** Return the number of incident edges to this node. */
	size_type degree () const {
		/*std::unordered_set<size_type>  incident_edges = ((*graph_).get_incident_edges(uid_));
		std::cout << "HERE" << std::endl;
		if (incident_edges.empty()){
			std::cout << "HERE" << std::endl;
			return 0;
		}
		else
			return ((*graph_).get_incident_edges(uid_)).size();
		*/
		return (*graph_).get_degree(uid_);
	}

	/** Return a incident_iterator object pointing at the first  Edge object incident to this node.*/
	incident_iterator edge_begin() const {
		graph_type* gr  =  const_cast<graph_type*>(graph_);
		if (degree() == 0)
			return IncidentIterator(nullptr, -1, std::unordered_set<size_type>{});
		else{
			std::unordered_set<size_type>  incident_edges   =  (*graph_).get_incident_edges(uid_);
			return IncidentIterator(gr, uid_, incident_edges);
		}
	}

	/** Return a incident_iterator object denoting there are no other incident edges to visit.
	 * @post the iterator returned has all its members set to invalid values, namely nullptr, -1 and empty set, according to their type.
	 */
	incident_iterator edge_end() const {
		return IncidentIterator(nullptr, -1, std::unordered_set<size_type>{});		
	}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
	  if (graph_ == n.graph_ and uid_ == n.uid_)
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
      // HW0: YOUR CODE HERE
	  if (uid_ <  n.uid_)
		  return true;
	  if (uid_ > n.uid_)																	//
		  return false;
	  if (uid_ == n.uid_ and std::less<graph_type*>{}(graph_, n.graph_))					
		  return true;
      return false;
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
	
	// Pointer back to the Graph
	graph_type* graph_;
	// Node's index
	size_type uid_;
	
	/** Private Constructor */
	Node(const graph_type* graph, size_type uid)
		: graph_(const_cast<graph_type*>(graph)),  uid_(uid){
	}
  };
  
  // HW1
   /**
	* Return a pointer to the collection of incident edge indices to a node
	*/
	std::unordered_set<size_type> get_incident_edges(size_type node_ind) const{
		return (neighbors_.find(node_ind))->second;
	}
	
	size_type get_degree(size_type node_ind) const{
		return (degrees_.find(node_ind))->second;
	}
	
  /**
	* Return a pointer to the collection of all edges in the graph
	*/
	std::tuple<size_type, size_type> get_edge(size_type edge_ind) const {
		return (internal_edges_.find(edge_ind))->second;
	}

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return active_size();														
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return active_size();														
  }
  
  /** Return the number of all nodes ever added to the graph.
   *
   * Complexity: O(1).
   */
  size_type real_size() const {
	  return size_;						
  }  
  
  /** Return the number of active nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type active_size() const {
	  return i2u_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
	internal_nodes_[next_node_idx_] = position;
	internal_node_val_[next_node_idx_] = val;
	degrees_[next_node_idx_] = 0;
	std::unordered_set<size_type> empty_set {};
	neighbors_[next_edge_idx_] = empty_set;
	idx_[next_node_idx_] = i2u_.size();										// HW2
	i2u_.push_back(next_node_idx_);
	++size_;
	++next_node_idx_;
    return Node(this, next_node_idx_-1);     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n)  const{ 
    // HW0: YOUR CODE HERE
	if (internal_nodes_.find(n.uid_) != internal_nodes_.end() and
		(internal_nodes_.find(n.uid_)->second).x == n.position() .x and			
		(internal_nodes_.find(n.uid_)->second).y == n.position().y and
		(internal_nodes_.find(n.uid_)->second).z == n.position().z)
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
  Node node(size_type i) const{
    // HW0: YOUR CODE HERE
	assert(i < active_size());
	Node result_node = Node(this, i2u_[i]);
	return result_node;
  }
  
  /** Return the user facing idx of node with index @a i
 	* @pre 0 <= @a i < num_nodes()
	*
	* Complexity: O(1).
	*/
  size_type get_node_idx(size_type i) const {						// HW2
	  assert(i < real_size());
	  assert(idx_.count(i) > 0);
	  return (idx_.find(i))->second;
  }
  
  /** Return the value of node with index @a i
	* @pre 0 <= @a i < num_nodes()
	*
	* Complexity: O(1).
	*/
	const node_value_type& get_node_val(size_type i) const {
		assert (i < real_size());
		assert(internal_node_val_.count(i) > 0 );
		return (internal_node_val_.find(i))->second;
	}
	node_value_type& get_node_val(size_type i) {
		assert (i < real_size());
		assert(internal_node_val_.count(i) > 0 );
		return (internal_node_val_.find(i))->second;
	}
	
	/** Set the value of node with index @a i to @a v
	* @pre 0 <= @a i < num_nodes()
	*
	* Complexity: O(1).
	*/
	void set_node_val(size_type i, node_value_type v) {
		assert (i < real_size());
		internal_node_val_[i] = v;
	}
	
  /** Return the position of node with index @a i
	* @pre 0 <= @a i < num_nodes()
	*
	* Complexity: O(1).
	*/
	Point * get_node_position(size_type i) const{
		assert (i < real_size());
		return const_cast<Point*>(&(internal_nodes_.find(i))->second);
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
      // HW0: YOUR CODE HERE
	  index_ = -1;
    }

    /** Return a node of this Edge 
	*    Arbitrarily the node with the smallest index
	*    is returned as the first one
	*/
    Node node1() const {
      // HW0: YOUR CODE HERE
	  size_type n_uid = std::get<0>(nodes_ind_);
	  size_type n_idx = (*graph_).get_node_idx(n_uid);
	  Node result_node = (*graph_).node(n_idx);
	  return result_node;
    }

    /** Return the other node of this Edge
	*    Arbitrarily the node with the largest index
	*    is returned as the second one
	*/
    Node node2() const {
      // HW0: YOUR CODE HERE
	  size_type n_uid = std::get<1>(nodes_ind_);
	  size_type n_idx = (*graph_).get_node_idx(n_uid);
	  Node result_node = (*graph_).node(n_idx);
	  return result_node;
    }
	
	 /** Return this edge's index */
    size_type index() const {
      // HW0: YOUR CODE HERE
	  return index_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes in the same graph.
     */
    bool operator==(const Edge& e) const {
		if (graph_ == e.graph_ ){
			if (std::get<0>(nodes_ind_) == std::get<0>(e.nodes_ind_))
				if (std::get<1>(nodes_ind_) == std::get<1>(e.nodes_ind_))
					return true;
			if (std::get<1>(nodes_ind_) == std::get<0>(e.nodes_ind_))
				if (std::get<0>(nodes_ind_) == std::get<1>(e.nodes_ind_))
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
		if (index_ < e.index_)		//index-based ordering
			return true;
		if (index_ > e.index_)
			return false;
		if (index_ == e.index_ and std::less<graph_type*>{}(graph_, e.graph_))
			return true;
		return false;
    }
	
	// HW2
	/** Return the Eucledian distance between the two nodes of an edge */
  double  length() const {
	  return norm(node1().position() - node2().position());	  
  }
  
	//HW2
	/** Return the edge's value. */
	edge_value_type& value() {
		return (*graph_).get_edge_val(index_);
	}
    const edge_value_type& value() const {
		return (*graph_).get_edge_val(index_);
	}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
	graph_type* graph_;
	size_type index_;
	//tuple of connected nodes' indices
	std::tuple<size_type, size_type> nodes_ind_;
	/** Private Constructor */
	Edge(const graph_type* graph, size_type index, std::tuple<size_type, size_type> nodes_ind)
		: graph_(const_cast<graph_type*>(graph)),  index_(index), nodes_ind_(nodes_ind){
	}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return internal_edges_.size();
	//return active_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
	assert(i < num_edges());
	Edge result_edge = Edge(this, i, (internal_edges_.find(i))->second);
	assert(result_edge.index_ == i);
	return result_edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
	assert(has_node(a));
	assert(has_node(b));
	if ((neighbors_.find(a.uid_)) == neighbors_.end() or (neighbors_.find(b.uid_)) == neighbors_.end()){
		return false;
	}
	for (auto idx : (neighbors_.find(a.uid_))->second){
		Edge result_edge = edge(idx);
		if (std::get<0>(result_edge.nodes_ind_) == b.uid_ or std::get<1>(result_edge.nodes_ind_) == b.uid_)
			return true;
	}
	return false;
  }
  
  //HW2
  /** Return the edge formed by nodes @a a and @a b.
    * @pre @a a and @b have a valid edge in the graph
	* Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	*/
  Edge edge(const Node& a , const Node& b) const {
	  assert(has_edge(a,b));
	  //check all edges of node a to see if they are edges of b too
	  Edge result_edge;
	  for (auto idx : (neighbors_.find(a.uid_))->second){
		Edge e = edge(idx);
		if (std::get<0>(e.nodes_ind_) == b.uid_ or std::get<1>(e.nodes_ind_) == b.uid_)
			return result_edge = e;
	  }
	  return result_edge;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // HW0: YOUR CODE HERE
	assert(has_node(a));
	assert(has_node(b));
	assert(!(a==b));
	// check if edge exists
	for (auto idx : neighbors_[a.uid_]){
		Edge result_edge = edge(idx);
		if (std::get<0>(result_edge.nodes_ind_) == b.uid_ or std::get<1>(result_edge.nodes_ind_) == b.uid_)
			return result_edge;
	}
	//add new edge to the graph
	internal_edges_[next_edge_idx_] = std::make_tuple(a.uid_, b.uid_);
	//HW2
	internal_edge_val_[next_edge_idx_] = val;
	//add new edge to the sets of edges of a and b
	neighbors_[a.uid_].insert(next_edge_idx_);
	neighbors_[b.uid_].insert(next_edge_idx_);
	// update nodes' degrees_
	degrees_[a.uid_] = degrees_[a.uid_] + 1;
	degrees_[b.uid_] = degrees_[b.uid_] + 1;
	++next_edge_idx_;
	//++active_edges_;
	assert(has_edge(a, b));
	return Edge(this, next_edge_idx_-1,   std::make_tuple(a.uid_, b.uid_));
  }
  
  /** Return the value of node with index @a i
  * @pre 0 <= @a i < num_edges()
    *
	* Complexity: O(1).
	*/
	const edge_value_type& get_edge_val(size_type i) const {
		assert(i < num_edges());
		return (internal_edge_val_.find(i))->second;
	}
  /** Return the value of node with index @a i
  * @pre 0 <= @a i < num_edges()
    *
	* Complexity: O(1).
	*/
	edge_value_type& get_edge_val(size_type i)  {
		assert(i < num_edges());
		return (internal_edge_val_.find(i))->second;
	}

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
	// HW2
	while (active_size() != 0){
		remove_node(node(0));
	}
	next_edge_idx_ = 0;
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
    // Supply definitions AND SPECIFICATIONS for: TO DO
	public:
	
	/** @brief Dereference the NodeIterator. 
	 *  @pre graph_ptr_ != nullptr, traversal is not finished
	 *  @return the corresponding Node object.
	 */
    Node operator*() const{
		return Node(graph_ptr_, node_ind_);
	}
	
	/** @brief Increment NodeIterator to traverse to the next Node. 
	 *  @post Either graph_ptr_ ==nullptr or (graph_ptr_ != nullptr and node_ptr_ points to the next Node object).
	 *  @return the updated NodeIterator object.
	 */
    NodeIterator& operator++(){
		if (node_ind_ < num_nodes_ -1){
			node_ind_++;
			Node n = Node(graph_ptr_, node_ind_ );						//HW2
			if (! n.valid()){																																// HW2
				++(*this);	
			}
		}			
		else{
			graph_ptr_ = nullptr;
			node_ind_ = 0;
			num_nodes_ = 0;
		}
		return *this;
	}
	
	/** @brief Tests if two NodeIterator objects are equal.
	 *  @return true if every member of the one iterator is equal to the respective member of the other iterator, else false.
	 */
    bool operator==(const NodeIterator& nditer)  const { 			
		if (graph_ptr_ == nditer.graph_ptr_ and
			node_ind_ == nditer.node_ind_ and num_nodes_ == nditer.num_nodes_)
				return true;
		else
				return false;
	}
	
	/** @brief Tests if two NodeIterator objects are different.
	 *  @return false if every member of the one iterator is equal to the respective member of the other iterator, else true.
	 */
	bool operator!=(const NodeIterator& nditer) const {
		if ( *this == nditer)
			return false;
		else
			return true;
	}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
	Graph* graph_ptr_;
	size_type node_ind_;
	size_type num_nodes_;
	
	//Private constructor that can be accessed by the Graph class.
	NodeIterator(Graph* graph_ptr, size_type node_ind, size_type num_nodes){
		if (graph_ptr != nullptr){
			graph_ptr_ = graph_ptr;
			node_ind_ = node_ind;
			num_nodes_ = num_nodes;
		}
		else{
			graph_ptr_ = nullptr;
			node_ind_ = 0;
			num_nodes_ = 0;
		}
	}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return a NodeIterator object pointing at the first  valid Node object of the graph.*/
  node_iterator node_begin() const{			
		graph_type* graph_ptr = const_cast<graph_type*>(this);
		size_type ind = 0 ;
		Node n = Node(this, ind);																// HW2
		while (! n.valid()){
			++ind;
			n = Node(this, ind);
		}
		return NodeIterator(graph_ptr, ind, size_);
	} 
	
	/** Return a NodeIterator object denoting the end of the Node collection.
	 * @post the iterator returned has all its members set to nullptr and 0, according to their type.
	 */
	node_iterator node_end() const{
		return NodeIterator(nullptr, 0, 0);
	}


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        		= Edge;                     // Element type
    using pointer           			= Edge*;                    // Pointers to elements
    using reference         		= Edge&;                    // Reference to elements
    using difference_type   	= std::ptrdiff_t;           // Signed difference
    using iterator_category 	= std::input_iterator_tag;  // Weak Category, Proxy
	using iterator 					= std::unordered_set< size_type>::iterator;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
	
	/** @brief Dereference the IncidentIterator. 
	 *  @pre graph_ptr_ != nullptr, traversal is not finished
	 *  @return the corresponding Edge object.
	 */
    Edge operator*() const {
		return Edge(graph_ptr_, edge_ind_, std::make_tuple(node1_ind_, node2_ind_));
	}
	
	/** @brief Increment IncidentIterator to traverse to the next Edge. 
	 *  @post Either graph_ptr_ ==nullptr or (graph_ptr_ != nullptr and node_ptr_ points to the next incident Edge object).
	 *  @return the updated IncidentIterator object.
	 */
    IncidentIterator& operator++() {
		Node n1 = Node(graph_ptr_, node1_ind_);
		if (n1.degree() == 0){
			graph_ptr_ = nullptr;
			node1_ind_ = -1;
			node2_ind_ = -1;
			edge_ind_ = -1;
			set_iter_ = (std::unordered_set<size_type>{}).end();
			incid_edges_ = std::unordered_set<size_type>{};
			return *this;
		}
		if (++set_iter_ != incid_edges_.end()) {
			edge_ind_ = *set_iter_;
			if (std::get<0>((*graph_ptr_).get_edge(edge_ind_))   == node1_ind_)
				node2_ind_ = std::get<1>((*graph_ptr_).get_edge(edge_ind_));
			else
				node2_ind_ = std::get<0>((*graph_ptr_).get_edge(edge_ind_));
		}
		else{
			graph_ptr_ = nullptr;
			node1_ind_ = -1;
			node2_ind_ = -1;
			edge_ind_ = -1;
			set_iter_ = (std::unordered_set<size_type>{}).end();
			incid_edges_ = std::unordered_set<size_type>{};
		}
		return *this;
	}
	
	/** @brief Tests if two IncidentIterator objects are equal.
	 *  @return true if every member of the one iterator is equal to the respective member of the other iterator, else false.
	 */
    bool operator==(const IncidentIterator& iter) const {
		if ( graph_ptr_ == iter.graph_ptr_ and node1_ind_ == iter.node1_ind_ and
			  node2_ind_ == iter.node2_ind_ and edge_ind_ == iter.edge_ind_ and
			  set_iter_ == iter.set_iter_ and incid_edges_ == iter.incid_edges_ )
					return true;
		else
					return false;
	}
	
	/** @brief Tests if two IncidentIterator objects are different.
	 *  @return false if every member of the one iterator is equal to the respective member of the other iterator, else true.
	 */
	bool operator!=(const IncidentIterator& iter) const {
		if ( *this == iter)
			return false;
		else
			return true;
	}

   private:
    friend class Graph;
	friend class Node;
    // HW1 #3: YOUR CODE HERE
	Graph* graph_ptr_;
	size_type node1_ind_;					// the node that spawns the iterator
	size_type node2_ind_;
	size_type edge_ind_;
	iterator set_iter_;
	std::unordered_set<size_type> incid_edges_;
	
	//Private constructor that can be accessed by the Node class.
	IncidentIterator(Graph* graph_ptr, size_type node_ind, std::unordered_set<size_type> incid_edges) {
		if (graph_ptr != nullptr) {
			graph_ptr_ = graph_ptr;
			node1_ind_ = node_ind;
			incid_edges_ = incid_edges;
			set_iter_ = incid_edges_.begin();
			edge_ind_ = * set_iter_;
			if (std::get<0>((*graph_ptr_).get_edge(edge_ind_))   == node1_ind_)
				node2_ind_ = std::get<1>((*graph_ptr_).get_edge(edge_ind_));
			else
				node2_ind_ =  std::get<0>((*graph_ptr_).get_edge(edge_ind_));
		}
		else {
			graph_ptr_ = nullptr;
			node1_ind_ = -1;
			node2_ind_ = -1;
			edge_ind_ = -1;
			set_iter_ = (std::unordered_set<size_type>{}).end();
			incid_edges_ = std::unordered_set<size_type>{};
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
    using value_type        			= Edge;                     // Element type
    using pointer           				= Edge*;                    // Pointers to elements
    using reference         			= Edge&;                    // Reference to elements
    using difference_type   		= std::ptrdiff_t;           // Signed difference
    using iterator_category 		= std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
	
	/** @brief Dereference the EdgeIterator. 
	 *  @pre graph_ptr_ != nullptr, traversal is not finished
	 *  @return the corresponding Edge object.
	 */
    Edge operator*() const {
		return *incid_iter_;
	}
	
	/** @brief Increment EdgeIterator to traverse to the next Edge. 
	 *  @post Either graph_ptr_ ==nullptr or (graph_ptr_ != nullptr and node_ptr_ points to the next Edge object).
	 *  @return the updated EdgeIterator object.
	 */
    EdgeIterator& operator++() {
		Node n = *(node_iter_);
		if (n.degree() == 0) {
			if (++(node_iter_) != (*graph_ptr_).node_end()) {
				Node n = *(node_iter_);
				incid_iter_ = n.edge_begin();
			}
			else{
				graph_ptr_ = nullptr;
				node_iter_ = graph_ptr_->node_end();
				Node n;
				incid_iter_ = n.edge_end();
			}
		}
		else if (++incid_iter_ != n.edge_end()) {									
			Edge e = *incid_iter_;
			size_type ind2 = (e.node2()).index();										
			size_type ind1 = (e.node1()).index();																	
			if (ind2 < ind1){
				return ++(*this);
			}
		}
		else if (++(node_iter_) != (*graph_ptr_).node_end()) {
			Node n = *(node_iter_);
			if (n.degree() == 0){
				incid_iter_ = n.edge_end();
				return ++(*this);
			}
			else {
				incid_iter_ = n.edge_begin();
				Edge e = *incid_iter_;
				size_type ind2 = (e.node2()).index();										
				size_type ind1 = (e.node1()).index();																	
				if (ind2 < ind1){
					return ++(*this);
				}
			}
		}
		else {
			graph_ptr_ = nullptr;
			node_iter_ = graph_ptr_->node_end();
			Node n;
			incid_iter_ = n.edge_end();
		}
		return *this;
	}
	
	/** @brief Tests if two EdgeIterator objects are equal.
	 *  @return true if every member of the one iterator is equal to the respective member of the other iterator, else false.
	 */
    bool operator==(const EdgeIterator& iter) const {
		if (iter.graph_ptr_ ==  nullptr and graph_ptr_ == nullptr)
			return true;
		if (graph_ptr_ == iter.graph_ptr_ and node_iter_ == iter.node_iter_ and incid_iter_ == iter.incid_iter_) 
			return true;
		else
			return false;
	}
	
	/** @brief Tests if two EdgeIterator objects are different.
	 *  @return false if every member of the one iterator is equal to the respective member of the other iterator, else true.
	 */
	bool operator!=(const EdgeIterator& iter) const {
		if ( *this == iter)
			return false;
		else
			return true;
	}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
	graph_type* graph_ptr_;
	NodeIterator  node_iter_;
	IncidentIterator incid_iter_;
	
	//Private constructor that can be accessed by the Graph class.
	EdgeIterator (graph_type* graph_ptr, NodeIterator node_iter){
		if (graph_ptr != nullptr) {
			graph_ptr_ = graph_ptr;
			node_iter_ = node_iter;
			Node n = * (node_iter_);
			incid_iter_ = n.edge_begin();
		}
		else {
			graph_ptr_ = nullptr;
			node_iter_ = graph_ptr_->node_end();
			Node n;
			incid_iter_ = n.edge_end();
		}
	}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return an edge_iterator object pointing at the first Edge object of the graph.*/
  edge_iterator edge_begin() const {
	  graph_type* graph_ptr = const_cast<graph_type*>(this);
	  NodeIterator node_iter = node_begin();
	  return EdgeIterator(graph_ptr,  node_iter);
  }  
  
  /** Return an edge_iterator object denoting there are no other edges to visit.
	* @post the iterator returned has all its members set to nullptr.
	*/
  edge_iterator edge_end() const {
	  return EdgeIterator(nullptr, this->node_end());
  }


//
// HW2 - Removal Methods
//

/** @brief Function that removes a given node @a n from the graph
  * @pre @a n is a valid node (belongs to active set of nodes)
   *@post @a n is not a valid node
   *@return the number of nodes removed from graph (1 indicating successful removal of given node)
   *
   * Complexity: O(num_nodes() + num_edges()).
   */
size_type remove_node(const Node& n){
	assert(n.valid());
	size_type idx = n.index();
	size_type uid = i2u_[idx];
	//remove all adj edges
	while (! neighbors_[uid].empty()){
		remove_edge(edge(*neighbors_[uid].begin()));
	}
	size_type swapped_uid = i2u_[i2u_.size() - 1];
	idx_[swapped_uid] = idx;
	std::swap(i2u_[idx], i2u_[i2u_.size() - 1]);
	i2u_.pop_back();
	assert(! n.valid());
	return 1;
}

/** @brief Function that removes a node pointed by the given iterator @a n_it
  * @pre @a *n_it is a valid node (belongs to active set of nodes)
   *@post @a *n_it is not a valid node
   *@return a node iterator to the next unvisited node
   *
   * Complexity: O(num_nodes() + num_edges()).
   */
node_iterator remove_node(node_iterator n_it){
	Node n = * n_it;
	assert(n.valid());
	node_iterator n_it_next;
	remove_node(n);																	// O(num_nodes())
	n_it_next = ++n_it;
	return n_it_next;
}

/** @brief Function that removes the edge formed by nodes @a n1 and @a n2
  * @pre @a n1 and @a n2 are valid nodes (belong to active set of nodes)
   *@post @a n1 and @a n2 don't have an edge together
   *@return 1 if the removal is successful or 0 if the edge doesn't exist
   *
   * Complexity: O(num_nodes() + num_edges()).
   */
size_type remove_edge(const Node& n1, const Node& n2){
	assert(n1.valid());											// O(1)
	assert(n2.valid());
	if (! has_edge(n1,n2))									// O(num_edges())
		return 0;
	Edge e = edge(n1, n2);
	remove_edge(e);
	return 1;
	
}

/** @brief Function that removes edge @a e
  * @pre the nodes that are connected by @a e are valid nodes
   *@post @a e is not an edge of the grpah
   *@return 1 if the removal is successful or 0 if the edge doesn't exist
   *
   * Complexity: O(num_nodes() + num_edges()).
   */
size_type remove_edge(const Edge& e){
	Node n1 = e.node1();
	Node n2 = e.node2();
	size_type idx1 = n1.index();
	size_type idx2 = n2.index();
	size_type uid1 = i2u_[idx1];
	size_type uid2 = i2u_[idx2];
	assert(n1.valid());												// O(1)
	assert(n2.valid());
	if (! has_edge(n1,n2))										// O(num_edges())
		return 0;
	size_type i = e.index();
	
	//update nodes' degrees_
	degrees_[uid1] =degrees_[uid1] - 1;
	degrees_[uid2] =degrees_[uid2] - 1;
	
	auto it_edges1 = neighbors_[uid1].find(i);
	if (it_edges1 != neighbors_[uid1].end()){
		neighbors_[uid1].erase(it_edges1);
	}
	
	auto it_edges2 = neighbors_[uid2].find(i);
	if (it_edges2 != neighbors_[uid2].end()){
		neighbors_[uid2].erase(it_edges2);
	}
	/*
	for (auto it_nodes = neighbors_.begin(); it_nodes != neighbors_.end(); ++it_nodes){
		auto it_edges = it_nodes->second.find(i);
		if (it_edges != it_nodes->second.end()){
			it_nodes->second.erase(it_edges);
		}
	}
	*/
	if (num_edges() > 0){
		Edge final_edge = edge(num_edges()-1);
		size_type n1f_idx = final_edge.node1().index();
		size_type n2f_idx = final_edge.node2().index();
		size_type uidf1 = i2u_[n1f_idx];
		size_type uidf2 = i2u_[n2f_idx];
		auto it_edges_f1 = neighbors_[uidf1].find(num_edges()-1);
		if (it_edges_f1 != neighbors_[uidf1].end()){
			neighbors_[uidf1].erase(it_edges_f1);
		}
		if ( i != num_edges() -1){
			neighbors_[uidf1].insert(i);
		}
		auto it_edges_f2 = neighbors_[uidf2].find(num_edges()-1);
		if (it_edges_f2 != neighbors_[uidf2].end()){
			neighbors_[uidf2].erase(it_edges_f2);
		}
		if ( i != num_edges() -1){
			neighbors_[uidf2].insert(i);
		}
	}
	/*
	for (auto it_nodes = neighbors_.begin(); it_nodes != neighbors_.end(); ++it_nodes){
		auto it_edges = it_nodes->second.find(num_edges()-1);
		if (it_edges != it_nodes->second.end()){
			it_nodes->second.erase(it_edges);
			if (i != num_edges()-1)
				it_nodes->second.insert(i);
		}
	}
	*/
	internal_edge_val_[i] = internal_edge_val_[num_edges()-1];
	internal_edge_val_.erase(num_edges()-1);
	internal_edges_[i] = internal_edges_[num_edges()-1];
	internal_edges_.erase(num_edges()-1);
	-- next_edge_idx_;
	//--active_edges_;
	return 1;
}

/** @brief Function that removes edge pointed by iterator @a e_it
  * @pre the nodes that are connected by @a e are valid nodes
   *@post @a *e_it is not an edge of the grpah
   *@return an edge iterator to the next unvisited edge
   *
   * Complexity: O(num_nodes() + num_edges()).
   */
edge_iterator remove_edge(edge_iterator e_it){
	Edge e = * e_it;
	edge_iterator e_it_next = ++ e_it;
	remove_edge(e);
	return e_it_next;
}


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  //maps node index to point position
  std::unordered_map<size_type, Point> internal_nodes_;
  //maps edge index to pair of nodes' indices
  std::unordered_map<size_type, std::tuple<size_type, size_type>> internal_edges_;
  size_type size_;								//total number of nodes in graph
  size_type next_node_idx_;			//next node index
  size_type next_edge_idx_;			//next node index
  //size_type active_edges_;
  //maps node index to its set of edge indices 
  std::unordered_map<size_type, std::unordered_set<size_type>> neighbors_; 
  //maps node index to its value
  std::unordered_map<size_type, node_value_type> internal_node_val_;
  //maps edge index to its value
  std::unordered_map<size_type, edge_value_type> internal_edge_val_;
  //maps a node's internal uid to its user facing idx_
  std::unordered_map<size_type, size_type> idx_;
  //stores the currently active nodes, indexed by idx_
  std::vector<size_type> i2u_;
  //maps a node's internal uid to its degree
  std::unordered_map<size_type, size_type> degrees_;
};



#endif // CME212_GRAPH_HPP