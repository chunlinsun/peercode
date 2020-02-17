#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */


#include <algorithm>
#include <set>
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
/*  private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.) */

  public:
 
  //using node_value_type = V;

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  
  typedef V node_value_type;
  typedef E edge_value_type;

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

  // new structs
  
  struct EdgeStruct{
    size_type id; //user index
    size_type Node1_uid;
    size_type Node2_uid;
    edge_value_type val;
    EdgeStruct(size_type id_, size_type Node1_uid_, size_type Node2_uid_, edge_value_type val_){
      id = id_;
      Node1_uid = Node1_uid_;
      Node2_uid = Node2_uid_;
      val = val_;
    }
  };
  struct NodeStruct{
    Point p;
    size_type id; // user index
    node_value_type val;
    
    NodeStruct(Point p_, size_type id_, node_value_type val_){
      p = p_;
      id = id_;
      val = val_;
    }
  };  
  /** Construct an empty graph. */

  Graph() {
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
	 

    Node(){
      graph_pointer = nullptr;
      node_index = 0;
    }
	
    Node(graph_type* graph_pointer_, size_type node_index_) {
		
      graph_pointer =graph_pointer_; 
      node_index = node_index_;
	  
    }

    /** Return this node's position as a constant reference. */
    const Point& position() const {
	  
      return graph_pointer->nodes_vector[unique_index()].p;
    }
    
    /** Return the position as a modifiable lvalue. */
    Point& position() {
	  
      return graph_pointer->nodes_vector[unique_index()].p;
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
		
      return size_type(node_index);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
	
    /** Return the number of nodes connected to current node
    */

    size_type degree() const{
    if ( graph_pointer->nodes_connections_map.find(unique_index()) != 
         graph_pointer->nodes_connections_map.end())
      return graph_pointer->nodes_connections_map.at(unique_index()).size();
    return 0;
    }
    
    /** Return value_type as a left value reference
        to be changed.
    */
    
      node_value_type& value(){
        return graph_pointer->nodes_vector[unique_index()].val;
      }
    
    /** Return value_type as a constant left value reference
        to be read only.
    */	
      const node_value_type& value() const{
        return graph_pointer->nodes_vector[unique_index()].val;
      }
     
    /** Return incident iterator. This is used to iterate
    /  the incident nodes. Begin() starts from the first node,
    /  end() is one iteration after the last node. 
    / (end() cannot be dereferenced and is used for exiting loops).
    / @pre degree() > 0
    */
    
    incident_iterator edge_begin() const{
        return IncidentIterator(graph_pointer, node_index, graph_pointer->nodes_connections_map.at(unique_index()).begin());
    }              
    incident_iterator edge_end() const{
        return IncidentIterator(graph_pointer, node_index, graph_pointer->nodes_connections_map.at(unique_index()).end());
    }	

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
		
      bool isequal {false};
      
      if (this->graph_pointer == n.graph_pointer and this->node_index == n.node_index)
      isequal = true;
        return isequal;
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
        
              
      if (graph_pointer == n.graph_pointer)
        if (unique_index() < n.unique_index())
          return true;
        else
          return false;
      else{
        std::set<graph_type*> ordered_set {graph_pointer, n.graph_pointer};
        if (*ordered_set.begin() == graph_pointer)
          return true;
        else
          return false;
      }
        
    }

  
    size_type unique_index() const {
		
      return size_type(graph_pointer->nodes_i2u[node_index]);
    }
    
     private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
	
	graph_type* graph_pointer; 
	size_type node_index;
	
	
	
	
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
	  
    return nodes_i2u.size();
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
    
    size_type node_id = nodes_i2u.size();
    size_type node_uid = nodes_vector.size();
    
    nodes_vector.push_back(NodeStruct(position, node_id, value));

    nodes_i2u.push_back(node_uid);    
    graph_type* graph_ptr = const_cast<graph_type*>(this);
    
    return Node(graph_ptr, node_id);       
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      
    bool doesit {false};
    
    if (n.index() < size())
      doesit = true;
	
    return doesit;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    graph_type* graph_ptr = const_cast<graph_type*>(this);
    node_type node_candidate {Node(graph_ptr,i)};
    if (has_node(node_candidate))
      return node_candidate;
    else
      return Node();
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
	
    Edge(){
      graph_pointer = nullptr;
      index_node1_pointer = nullptr; 
      index_node2_pointer = nullptr;
      // this is current id (not unique id)
      edge_index = 0;
    }
	
    Edge(graph_type* graph_pointer_, 
         size_type* index_node1_pointer_, 
         size_type* index_node2_pointer_, 
         size_type edge_index_) {
      graph_pointer = graph_pointer_;
      index_node1_pointer = index_node1_pointer_; 
      index_node2_pointer = index_node2_pointer_; 
      edge_index = edge_index_;
      
    }


    /** Return a node of this Edge */
    Node node1() const {
		
	  return Node(graph_pointer, *index_node1_pointer);
    }

    /** Return the other node of this Edge */
    Node node2() const {

	  return Node(graph_pointer, *index_node2_pointer);
    }
    
	/** Return value_type as a left value reference
	    to be changed.
	*/
	
    edge_value_type& value(){
      //return graph_pointer->edges_value_vector[graph_pointer->edges_i2u[edge_index]];
      
      return graph_pointer->edges_vector[unique_index()].val;
    }
	
	/** Return value_type as a constant left value reference
	    to be read only.
	*/	
    const edge_value_type& value() const{
      return graph_pointer->edges_vector[unique_index()].val;
    }
 
    
    double length() const{
      return norm(node1().position()-node2().position());
    }
    
     

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
		
	  if (*index_node1_pointer == *e.index_node1_pointer and 
        *index_node2_pointer == *e.index_node2_pointer and
        edge_index == e.edge_index and 
        graph_pointer == e.graph_pointer)
	    return true;
	  else if (*index_node1_pointer == *e.index_node2_pointer and 
             *index_node2_pointer == *e.index_node1_pointer and
             edge_index == e.edge_index and 
             graph_pointer == e.graph_pointer)
	    return true;
    else           // Quiet compiler warning
      return false;
		
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return false;
	  
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      
      if (graph_pointer == e.graph_pointer)
        if (unique_index() < e.unique_index())
          return true;
        else
          return false;
      else{
        std::set<graph_type*> ordered_set {graph_pointer, e.graph_pointer};
        if (*ordered_set.begin() == graph_pointer)
          return true;
        else
          return false;
      }
	  
    }

    
    
    size_type unique_index() const{
      
      return graph_pointer->edges_i2u[edge_index];
    }
    
    private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
	
    size_type* index_node1_pointer;
    size_type* index_node2_pointer;
    size_type edge_index;
    graph_type* graph_pointer;
    
	
	
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
	  
    return edges_i2u.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
   
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
	
    if ( i >= num_edges() ){
      return Edge();
    }
    
    size_type edge_uid = edges_i2u[i];
    
    size_type node1_uid = edges_vector[edge_uid].Node1_uid;
    size_type node2_uid = edges_vector[edge_uid].Node2_uid;
    
    size_type* node1_id_ptr = const_cast<size_type*>(&nodes_vector[node1_uid].id);
    size_type* node2_id_ptr = const_cast<size_type*>(&nodes_vector[node2_uid].id);
    
    graph_type* graph_ptr = const_cast<graph_type*>(this);
    
    return Edge(graph_ptr ,node1_id_ptr ,node2_id_ptr, i);
    
    }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
	
    size_type node1_id = a.index();
    size_type node2_id = b.index();
    
    size_type node1_uid = nodes_i2u[node1_id];
    size_type node2_uid = nodes_i2u[node2_id];
  
    bool doesit {false};
    
      
    if ( nodes_connections_map.find(node1_uid) != nodes_connections_map.end())
      if ( nodes_connections_map.at(node1_uid).find(node2_uid) != nodes_connections_map.at(node1_uid).end() )
        doesit = true;
    
    return doesit;
	
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
    
    size_type node1_id = a.index();
    size_type node2_id = b.index();
    
    size_type node1_uid = nodes_i2u[node1_id];
    size_type node2_uid = nodes_i2u[node2_id];
    
    size_type* node1_id_ptr = &nodes_vector[node1_uid].id;
    size_type* node2_id_ptr = &nodes_vector[node2_uid].id;
    
    graph_type* graph_ptr = this;
       
	// I check if a and b are already connected.
	// has_edge(a,b) and has_edge(b,a) are never both true at the same time
	// because has_edge checks the unique connections and if an unique connection
	// does not exist I add only one in the unique map
    if (has_edge(a, b)){
      
      // I get the unique id from the connection map
      size_type edge_uid = nodes_connections_map[node1_uid][node2_uid];
      size_type edge_id = edges_vector[edge_uid].id;
      
      return Edge(graph_ptr, node1_id_ptr, node2_id_ptr, edge_id);
                         
    }

    // new connection!
    else{
       
      // user index of the new edge
      size_type edge_id = edges_i2u.size();
      // global index of the new edge
      size_type edge_uid = edges_vector.size();

      // add new edge to edge vector
      edges_vector.push_back(EdgeStruct(edge_id,node1_uid,node2_uid,value));

      // new edge index is always the number of total edges befre addition
      edges_i2u.push_back(edge_uid);

      // I add both in the all map
      nodes_connections_map[node1_uid][node2_uid] = edge_uid;
      nodes_connections_map[node2_uid][node1_uid] = edge_uid;
      return Edge(this, node1_id_ptr, node2_id_ptr, edge_id);
    }
	
  }
	

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_i2u = {};
    edges_i2u = {};
    edges_vector = {};
    nodes_vector = {};
    nodes_connections_map = {};
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
		graph_pointer = nullptr;
		current_index = 0;
    }
	

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
	
	/** operator() dereference the iterator and returns current node
	*/
	Node operator*() const{
		return Node(graph_pointer, current_index);
	}

	/** operator++() points at the next node or end, if end is reached
	*/
	
	NodeIterator& operator++(){
		current_index++;
		return *this;
	}

	/** operator==() returns true if two iterators have the same member attributes
	*/
	
	bool operator==(const NodeIterator& node_iterator) const{
		if (this->graph_pointer == node_iterator.graph_pointer and 
			this->current_index == node_iterator.current_index)
			return true;
		else
			return false;
	}

	/** operator!=() is the logical opposit of operator++()
	*/
	
	bool operator!=(const NodeIterator& node_iterator) const{
		if (*this == node_iterator)
			return false;
		else
			return true;
	}
	
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
	
	// thanks to current index I can access the array randomly.
	size_type current_index;
	graph_type* graph_pointer;
	
	
	NodeIterator(graph_type* graph_pointer_, size_type current_index_) {
		graph_pointer = graph_pointer_;
		current_index = current_index_;
		
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
  /** node_begin() is a node_iterator that points
  // to the first node
  */

  node_iterator node_begin() const{
    graph_type* graph_ptr = const_cast<graph_type*>(this);
	  return NodeIterator(graph_ptr,0);
  }

  /** node_begin() is a node_iterator that points
  // to one after the last node (cannot be dereferenced and is used for exiting loops)
  */
  
  node_iterator node_end() const{
    graph_type* graph_ptr = const_cast<graph_type*>(this);
	  return NodeIterator(graph_ptr, nodes_i2u.size());
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
		graph_pointer = nullptr;
		node1_index = 0;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** operator*() dereference the Incident iterator by returning current edge connection.
    // Edge orientation is taken into consideration by using nodes_connections_map.
    // Edge constructor needs pointers to nodes. The nodes Edge points to belong to the graph member 
    // attribute nodes_connections_map which can be access with node indeces and stores all 
    // the connected nodes. 
    */ 
    Edge operator*() const{
      //return Edge(const_cast<const graph_type*>(graph_pointer),const_cast<size_type*>( &((graph_pointer->nodes_connections_map[map_iterator->first].find(node1_index))->first) ), const_cast<size_type*>(&(map_iterator->first)), map_iterator->second);
      
      size_type node1_uid = graph_pointer->nodes_i2u[node1_index];
      size_type* node1_id_ptr = &graph_pointer->nodes_vector[node1_uid].id;
      
      size_type node2_uid = map_iterator->first;
      size_type* node2_id_ptr = &graph_pointer->nodes_vector[node2_uid].id;
      
      size_type edge_uid = map_iterator->second;
      size_type edge_id = graph_pointer->edges_vector[edge_uid].id;
      
      return Edge(graph_pointer, node1_id_ptr, node2_id_ptr, edge_id);    
      
    }

    /** operator++() iterates to the next edge connection
    */
    
    IncidentIterator& operator++(){
      ++map_iterator;
      return *this;
    }

    /** operator ==() checks if two iteratior of this class have the same member attributes
    */
    
    bool operator==(const IncidentIterator& incident_iterator2) const{
      if(this->graph_pointer == incident_iterator2.graph_pointer and
         this->node1_index   == incident_iterator2.node1_index and
         this->map_iterator  == incident_iterator2.map_iterator)
        return true;
      else
        return false;
    }
    /** operator!=() is the logical opposite of operator==()
    */
      
    bool operator!=(const IncidentIterator& incident_iterator) const{
      if (*this == incident_iterator)
        return false;
      else
        return true;
    }
    
     private:
      friend class Graph;
      // HW1 #3: YOUR CODE HERE
    
    graph_type* graph_pointer;
    size_type node1_index;
    std::map<size_type, size_type>::iterator map_iterator;
    
    IncidentIterator(graph_type* graph_pointer_, size_type node1_index_, std::map<size_type, size_type>::iterator map_iterator_){
      graph_pointer = graph_pointer_;
      node1_index = node1_index_;
      map_iterator = map_iterator_;
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
		graph_pointer = nullptr;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
	
	/** Edge operator*() dereference the iterator by returning the current edge
	//  Edge constructor needs pointers to nodes. These pointer are taken from graph member
	//  attribute edges_map which is very easy and quick to access but does not store edge orientation
	//  but it is fine since it is not required here.
	*/
	Edge operator*() const{
    size_type edge_uid = *edges_i2u_iterator;
    size_type edge_id = graph_pointer->edges_vector[edge_uid].id;
    size_type node1_uid = graph_pointer->edges_vector[edge_uid].Node1_uid;
    size_type node2_uid = graph_pointer->edges_vector[edge_uid].Node2_uid;
    
    size_type* node1_id_ptr = &graph_pointer->nodes_vector[node1_uid].id;
    size_type* node2_id_ptr = &graph_pointer->nodes_vector[node2_uid].id;
    
		return Edge(graph_pointer, node1_id_ptr, node2_id_ptr , edge_id );
  }
	
	/** operator++() iterates to the next edge or end if current iterator dereferences to last edge
	*/
	EdgeIterator& operator++(){
		++edges_i2u_iterator;
		return *this;
	}

	/** operator==() returns true if input iterator has the same member attributes as
	// current iterator
	*/
	
	bool operator==(const EdgeIterator& edge_iterator2) const{
		if (this->graph_pointer     ==  edge_iterator2.graph_pointer and
			this->edges_i2u_iterator ==  edge_iterator2.edges_i2u_iterator)
			return true;
		else
			return false;
	}

	/** operator!=() is the logical opposite of operator==()
	*/
	
	bool operator !=(const EdgeIterator& edge_iterator2) const{
		if (*this == edge_iterator2)
			return false;
		else
			return true;
	}
	
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
	
	graph_type* graph_pointer;
  std::vector<size_type>::iterator edges_i2u_iterator;
	
	EdgeIterator(graph_type* graph_pointer_, std::vector<size_type>::iterator edges_i2u_iterator_){
		graph_pointer = graph_pointer_;
    edges_i2u_iterator = edges_i2u_iterator_;
		
    }
	
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  
  /** edge_begin() initializes the edge iterator to point to the 
  // first edge (the first that was added to edge_map)
  */
  edge_iterator edge_begin(){
	  return EdgeIterator(this, edges_i2u.begin());
  }

  /** edge_end() points to one after the last edge. It is used to stop iterations
  // and exit loop.
  */
  
  edge_iterator edge_end(){
	  return EdgeIterator(this, edges_i2u.end());
  }
  
  // Edge Removal
  /**
  * @brief remove_edge removes the edge from the graph
  * @param[in] Node1: one of the two nodes defining the edge
  * @param[in] Node2: the other node defining the edge
  * @param[out] return_value is 0 if the edge does not exist; else,
  *             return_value is the user id of the removed edge, 
  *             (could be 0, user needs to check) 
  *             which is now the user id of the swapped edge.
  * @post (0 <= return_value < number of active edges in the graph)
  * Complexity: O(1)
  */
  
  size_type remove_edge(const Node& Node1, const Node& Node2){
    // I check that the edge exists
    if(edges_i2u.size() == 0 || has_edge(Node1,Node2) == false || has_edge(Node2,Node1) == false){
      return 0;
    }
    // I get the nodes unique indeces
    size_type node1_uid = Node1.unique_index();
    size_type node2_uid = Node2.unique_index();
    
    // I get the edge unique index
    size_type edge_uid = nodes_connections_map.at(node1_uid)[node2_uid];
    // I get the edge current index (the one I want to remove)
    size_type edge_id = edges_vector[edge_uid].id;
    
    // I remove the connection in the map (both entries)
    nodes_connections_map.at(node1_uid).erase(node2_uid);
    nodes_connections_map.at(node2_uid).erase(node1_uid);
    
    // I get the unique index of the element I want to swap
    size_type swap_edge_uid = edges_i2u.back();
    
    // I change the user index in the data structure for bookkeeping
    // now the swapped edge has the id of the one I am removing
    edges_vector[swap_edge_uid].id = edge_id;
    
    // I swap the value (effectively removing the old one)
    edges_i2u[edge_id] = swap_edge_uid;
    
    // I finally pop the last element (which now is a useless copy)
    edges_i2u.pop_back();
    
    return edge_id;
 
  }
  
  /**
  * @brief remove_edge removes the edge from the graph
  * @param[in] e: edge
  * @param[out] return_value is 0 if the edge does not exist; else,
  *             return_value is the user id of the removed edge, 
  *             (could be 0, user needs to check) 
  *             which is now the user id of the swapped edge.
  * @post (0 <= return_value < number of active edges in the graph)
  * Complexity: O(1)
  */
  
  size_type remove_edge ( const Edge & e){
    node_type Node1 = e.node1();
    node_type Node2 = e.node2();
    return remove_edge(Node1,Node2);
    
  }

  /**
  * @brief remove_edge removes the edge from the graph
  * @param[in] e_it: edge iterator
  * @param[out] return_value is 0 if the edge does not exist; else,
  *             return_value is the user id of the removed edge, 
  *             (could be 0, user needs to check) 
  *             which is now the user id of the swapped edge.
  * @pre edge is a dereferencable edge iterator (not end())
  * @post (0 <= return_value < number of active edges in the graph)
  * Complexity: O(1)
  */
  
  edge_iterator remove_edge(edge_iterator e_it){
    edge_type e = *e_it;
    remove_edge(e);
    return (this->edge_begin());
  }

  /**
  * @brief remove_node removes the node from the graph
  * @param[in] n: node to remove
  * @param[out] return_value is 0 if the node does not exist; else,
  *             return_value is the user id of the removed node, 
  *             (could be 0, user needs to check) 
  *             which is now the user id of the swapped node.
  * @post (0 <= return_value < number of active nodes in the graph)
  * Complexity: O(1)
  */  
  
  size_type remove_node ( const Node & n){
    
    // I chack if the node exists
    if (has_node(n) == false) return 0;
    
    // I get the id and uid of the node I am removing
    size_type node_id = n.index();
    size_type node_uid = n.unique_index();
    
    if (nodes_connections_map.find(node_uid) != nodes_connections_map.end())
      for(auto it = n.edge_begin(); it != n.edge_end(); ++it)
        remove_edge(*it);
    
    // I get the uid of the node I need to swap
    size_type swap_node_uid = nodes_i2u.back();
    
    // I change the swap node's user id to the one I am erasing for bookkeeping
    nodes_vector[swap_node_uid].id = node_id;
    
    // I swap the element
    nodes_i2u[node_id] = swap_node_uid;
    
    //Now I can pop out the useless copy
    nodes_i2u.pop_back();
    
    return node_id;
  }
    
  /**
  * @brief remove_node removes the node from the graph
  * @param[in] n_it: node iterator pointing to the node to remove
  * @param[out] return_value is 0 if the node does not exist; else,
  *             return_value is the user id of the removed node, 
  *             (could be 0, user needs to check) 
  *             which is now the user id of the swapped node.
  * @pre n_it must be a dereferencable node iterator (not end())
  * @post (0 <= return_value < number of active nodes in the graph)
  * Complexity: O(1)
  */  
  
  node_iterator remove_node ( node_iterator n_it ){
    remove_node(*n_it);
    return node_begin();
  }
 
 
 private:
 
 

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // vecotr of all nodes
  std::vector<NodeStruct> nodes_vector {};
  
  // vector of unique node indeces
  std::vector<size_type> nodes_i2u {};
  
  // vector of all edges
  std::vector<EdgeStruct> edges_vector {};
  
  // vector of unique edge indeces
  std::vector<size_type> edges_i2u {};
  
  // keys are node unique indeces, and value is a map where the keys are the unique indeces of the nodes connected to the node
  // and value is the edge unique index for that connection
 
  // all the edges connected to the node are saved
  std::map< size_type, std::map<size_type, size_type> > nodes_connections_map {};
  
};

#endif // CME212_GRAPH_HPP
