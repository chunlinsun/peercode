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
class Graph {
/*  private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.) */

  public:
 
  using node_value_type = V;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

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
	
    Node(const graph_type* graph_pointer_, size_type node_index_) {
		
	  graph_pointer = const_cast<graph_type*> (graph_pointer_); 
	  node_index = node_index_;
	  
    }

    /** Return this node's position. */
    const Point& position() const {
	  
      return graph_pointer->nodes_vector[node_index];
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
//--design_1
//--what happens if there are no edges?
//--START
	size_type degree() const{
		return graph_pointer->nodes_connections_map.at(node_index).size();
	}
//--END
	/** Return value_type as a left value reference
	    to be changed.
	*/
	
    node_value_type& value(){
	  return graph_pointer->nodes_value_vector[node_index];
    }
	
	/** Return value_type as a constant left value reference
	    to be read only.
	*/	
    const node_value_type& value() const{
	   return graph_pointer->nodes_value_vector[node_index];
    }
   
	/** Return incident iterator. This is used to iterate
	//  the incident nodes. Begin() starts from the first node,
	//  end() is one iteration after the last node. 
	// (end() cannot be dereferenced and is used for exiting loops).
	*/
	
	incident_iterator edge_begin() const{
			return IncidentIterator(graph_pointer, node_index, graph_pointer->nodes_connections_map.at(node_index).begin());
	}              
	incident_iterator edge_end() const{
			return IncidentIterator(graph_pointer, node_index, graph_pointer->nodes_connections_map.at(node_index).end());
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
	  
	  bool isless {false};
	  
	  if (this->node_index < n.node_index)
		  isless = true; 
	  
      return isless;
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
	  
    return nodes_vector.size();
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
	
	nodes_vector.push_back(position);
	
	//default initialization (will be modified by node class)
	nodes_value_vector.push_back(value);
    //(void) position;      // Quiet compiler warning
    return Node(this, nodes_vector.size()-1);       
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
	  
	bool doesit {false};
	
	if (n.node_index < this->nodes_vector.size())
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

  node_type node_candidate {Node(this,i)};
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
	  edge_index = 0;
	}
	
    Edge(const graph_type* graph_pointer_, size_type* index_node1_pointer_, size_type* index_node2_pointer_, size_type edge_index_) {
	  graph_pointer = const_cast<graph_type*> (graph_pointer_);
	  index_node1_pointer = const_cast<size_type*>(index_node1_pointer_); 
	  index_node2_pointer = const_cast<size_type*>(index_node2_pointer_); 
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

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
		
	  if (*index_node1_pointer == *e.index_node1_pointer and *index_node2_pointer == *e.index_node2_pointer)
	    return true;
	  else if (*index_node1_pointer == *e.index_node2_pointer and *index_node2_pointer == *e.index_node1_pointer)
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
		
		
	  if (edge_index < e.edge_index)
	    return true;
	  else
                 // Quiet compiler warning
        return false;
		
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return false;
	  
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
	  
    return edges_map.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
   
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
//--functionality_1
//--not sure why you're passing refernces to unsigned ints, but you need to cosnt cast them to pass to Edge
//--START	
	if ( edges_map.find(i) != edges_map.end() )
        return Edge(this ,const_cast<size_type*>(&edges_map.at(i)[0]),const_cast<size_type*>( &edges_map.at(i)[1]), i);//--END
	else
    //(void) i;             // Quiet compiler warning
      return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
	
	bool doesit {false};
	
	if ( nodes_connections_map.find(a.node_index) != nodes_connections_map.end())
		if ( nodes_connections_map.at(a.node_index).find(b.node_index) != nodes_connections_map.at(a.node_index).end() )
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE

	// I check if a and b are already connected.
	// has_edge(a,b) and has_edge(b,a) are never both true at the same time
	// because has_edge checks the unique connections and if an unique connection
	// does not exist I add only one in the unique map
	if (has_edge(a, b)){
	  return Edge(this, &edges_map[nodes_connections_map[a.node_index][b.node_index]][0], &edges_map[nodes_connections_map[a.node_index][b.node_index]][1], nodes_connections_map[a.node_index][b.node_index]);
	}

	// new connection!
    else{
	  size_type edge_index = edges_map.size();
	  edges_map[edge_index] = {a.node_index, b.node_index};
	  // I add both in the all map
	  nodes_connections_map[a.node_index][b.node_index] = edge_index;
	  nodes_connections_map[b.node_index][a.node_index] = edge_index;
	  return Edge(this, &edges_map[edge_index][0], &edges_map[edge_index][1], edge_index);
	}
	
  }
	

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
	nodes_vector = {};
    edges_map = {};
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
	const graph_type* graph_pointer;
	
	
	NodeIterator(const graph_type* graph_pointer_, size_type current_index_) {
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
	  return NodeIterator(this,0);
  }

  /** node_begin() is a node_iterator that points
  // to one after the last node (cannot be dereferenced and is used for exiting loops)
  */
  
  node_iterator node_end() const{
	  return NodeIterator(this,this->nodes_vector.size());
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
		return Edge(const_cast<const graph_type*>(graph_pointer),const_cast<size_type*>( &((graph_pointer->nodes_connections_map[map_iterator->first].find(node1_index))->first) ), const_cast<size_type*>(&(map_iterator->first)), map_iterator->second);	
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
	
	IncidentIterator(const graph_type* graph_pointer_, size_type node1_index_, std::map<size_type, size_type>::iterator map_iterator_){
		graph_pointer = const_cast<graph_type*> (graph_pointer_);
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
		return Edge(graph_pointer, &edge_map_iterator->second[0], &edge_map_iterator->second[1], edge_map_iterator->first);
	}
	
	/** operator++() iterates to the next edge or end if current iterator dereferences to last edge
	*/
	EdgeIterator& operator++(){
		++edge_map_iterator;
		return *this;
	}

	/** operator==() returns true if input iterator has the same member attributes as
	// current iterator
	*/
	
	bool operator==(const EdgeIterator& edge_iterator2) const{
		if (this->graph_pointer     ==  edge_iterator2.graph_pointer and
			this->edge_map_iterator ==  edge_iterator2.edge_map_iterator)
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
	
	const graph_type* graph_pointer;
	std::map<size_type, std::vector<size_type>>::iterator edge_map_iterator;
	
	EdgeIterator(const graph_type* graph_pointer_, std::map<size_type, std::vector<size_type>>::iterator edge_map_iterator_){
		graph_pointer = graph_pointer_;
		edge_map_iterator = edge_map_iterator_;
		
    }
//--documentation_0
//--please use doxygen style documentation
//--END	
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  
  /** edge_begin() initializes the edge iterator to point to the 
  // first edge (the first that was added to edge_map)
  */
  edge_iterator edge_begin(){
	  return EdgeIterator(this, edges_map.begin());
  }

  /** edge_end() points to one after the last edge. It is used to stop iterations
  // and exit loop.
  */
  
  edge_iterator edge_end(){
	  return EdgeIterator(this, edges_map.end());
  }
  
 private:
 
 

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // vector of all the nodes
  std::vector<Point> nodes_vector {};
  
  // vector of node_value_type
  std::vector<node_value_type> nodes_value_vector {};
  
  // keys are edge indeces, value is a vector of 2 elements (left and right node indeces) 
  std::map<size_type, std::vector<size_type>> edges_map {};
  
  // keys are node indeces, and value is a map where the keys are the indeces of the nodes connected to the node
  // and value is the edge index for that connection
 
  // all the edges connected to the node are saved
  std::map< size_type, std::map<size_type, size_type> > nodes_connections_map {};
  
};

#endif // CME212_GRAPH_HPP
