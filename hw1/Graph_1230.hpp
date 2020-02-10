#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>
#include <map>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V> //used for node value
class Graph {
 
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /**Node's user specified vale*/
  using node_value_type = V;

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

  	std::vector<Point> node_vec; // makes an empty vector of Points
  	std::vector<std::tuple<int, int>> edge_ends; //makes an empty vector of tuples
 	std::map<int, std::map<int, int>> node_to_edge; //makes an empty map of maps
 	std::vector<node_value_type> value_vec; //makes an empty vector of "value types"


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
   	
    /** Constructs an invalid node**/
    Node(){}
   

    /** Returns this node's position. */
    const Point& position() const {
      return node_graph_ptr-> node_vec[node_idx];
    }

    /** Returns this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_idx;
    }

    /** Returns this node's value. */
    node_value_type& value() {
    	return node_graph_ptr->value_vec.at(node_idx);

    };

    /** Returns this node's value for constant node */
    const node_value_type& value() const {
    	return &(node_graph_ptr->value_vec.at(node_idx));

    };

     /**returns the number of incident edgevalue()s*/
          size_type degree() const{
     	return node_graph_ptr-> node_to_edge.at(node_idx).size();
     }

	//**Start of the incident iterator */
     incident_iterator edge_begin() const {
     	incident_iterator iter;
     	iter.init_node_ptr = const_cast<Node*>(this);
    	iter.init_itr = node_graph_ptr->node_to_edge.at(node_idx).begin();

     	return iter;
     }

     /**End of the incident iterator*/
     incident_iterator edge_end() const {
     	incident_iterator iter;
    	iter.init_itr = node_graph_ptr->node_to_edge.at(node_idx).end();
     	iter.init_node_ptr = const_cast<Node*>(this);
     	return iter;
     }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
    	//tests if node is in the same graph
    	if (n.node_graph_ptr != node_graph_ptr){
    		return false;
    	}
    	//tests if node has same point in graph
    	if (n.node_idx != node_idx){
    		return false;
    	}
    	return true;
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

    //checks if @a n's index is less than or equal 
    if (n.node_idx > node_idx){
    	return true;
    }
    //uses graph pointer to break equality
    if (n.node_idx == node_idx){
    	if (n.node_graph_ptr > node_graph_ptr){
    		return true;
        }
    }
    return false;

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

   	graph_type* node_graph_ptr; //reference to node's graph
    int node_idx; //index in node_vec

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_vec.size();
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
    //adds point to node_vec method in Graph class
  	node_vec.push_back(position);
  	value_vec.push_back(value);


  	//adds a map in the node_to edge map, used later if node connects an edge
	std::map<int, int> m;
  	node_to_edge[node_vec.size()-1]=m;

  	//creates a node to be returned
  	Node n;
  	n.node_idx = node_vec.size()-1;
  	n.node_graph_ptr = const_cast<graph_type*>(this);
    return n;     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
  	//chekcs if node points to this graph
  	if (n.node_graph_ptr == const_cast<graph_type*>(this)){
  		//checks if index in range
  		if (n.node_idx < int(node_vec.size())){
  			return true;
  		}
  		else{
  			return false;
  		}
  	}
  	else{
  		return false;
  	}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
  	//creats a node with this index to be returned
  	Node n;
  	n.node_idx = i;
  	n.node_graph_ptr = const_cast<graph_type*>(this);
    return n;      
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

    /** Constructs an invalid Edge. */
       	Edge(){}

    /** Returns a node of this Edge */
    Node node1() const {
    	Node n;
    	if (main_node == nullptr){
    		n.node_graph_ptr = edge_graph_ptr;
    		n.node_idx = std::get<0>(edge_graph_ptr->edge_ends[edge_idx]);
      	return n;     
      	}
      	else{
      		return *main_node;
      	}
    }

    /** Returns the other node of this Edge */
    Node node2() const {
    	Node n;
    	n.node_graph_ptr = edge_graph_ptr;
    	if (main_node == nullptr){
    		n.node_idx = std::get<1>(edge_graph_ptr->edge_ends[edge_idx]);
  		}
  		else{
  			if ((*main_node).node_idx != std::get<0>(edge_graph_ptr->edge_ends[edge_idx])){
  				n.node_idx = std::get<0>(edge_graph_ptr->edge_ends[edge_idx]);
  			}
  			else{
  				n.node_idx = std::get<1>(edge_graph_ptr->edge_ends[edge_idx]);

  			}
  		}
  		return n;
     }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
    	//comapres @a nodes to first node of edge
    	if (std::get<0>(edge_graph_ptr->edge_ends[e.edge_idx]) != \
    		std::get<0>(edge_graph_ptr->edge_ends[edge_idx])){
    		if (std::get<0>(edge_graph_ptr->edge_ends[e.edge_idx]) != \
    			std::get<1>(edge_graph_ptr->edge_ends[edge_idx])){
    			return false;
    		}
    		else{
    			if (std::get<1>(edge_graph_ptr->edge_ends[e.edge_idx]) != \
    				std::get<0>(edge_graph_ptr->edge_ends[edge_idx])){
    				return false;
    			}
    			else{
    				return true;
    			}
    		}

    	}
    	//compares @a nodes to second node of edge
    	else{
    		if (std::get<1>(edge_graph_ptr->edge_ends[e.edge_idx]) != \
    		 std::get<1>(edge_graph_ptr->edge_ends[edge_idx])){
    			return false;
    		}
    		else{
    			return true;
    		}
    	}
    }
  

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
    	//compares the indices of the edges
      if (edge_idx < e.edge_idx){
      	return true;
      }
      else{
      	return false;
  		}
    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

	graph_type* edge_graph_ptr; //pointer to the graph
  	int edge_idx; //index of edge in edge_ends
 	node_type* main_node; //pointer to node from the incident iterator 
  	
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_ends.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
  	//creates edge to return
  	Edge e;
  	e.edge_idx = i;
  	e.edge_graph_ptr = const_cast<graph_type*>(this);
    return e;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
  		if (node_to_edge.count(a.node_idx !=0)){
  			if (node_to_edge.at(a.node_idx).count(b.node_idx)!=0){
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
  Edge add_edge(const Node& a, const Node& b) {
  	//checks if edge already in graph structures
  	if (has_edge(a,b) == false){
  		//edge doesn't already exist
  		//updates graph structures
  		std::tuple<int, int> tups(a.node_idx, b.node_idx);
    	edge_ends.push_back(tups);
	  	node_to_edge[a.node_idx][b.node_idx] = edge_ends.size()-1;
	  	node_to_edge[b.node_idx][a.node_idx] = edge_ends.size()-1;

    	//creates an edge to return
    	Edge e;
    	e.edge_graph_ptr = const_cast<graph_type*>(this);
		e.edge_idx = edge_ends.size()-1;
		return e;
  	}
  	else{
  		//edge already exists
  		Edge e;
  		e.edge_graph_ptr = const_cast<graph_type*>(this);
  		e.edge_idx = node_to_edge[a.node_idx][b.node_idx];
  		return e;
  		}  	
  	}

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_vec = {}; // Empties the vector of Points
  	edge_ends = {}; //Empties the vector of tuples
 	node_to_edge = {}; //Empties the map of maps
 	value_vec = {}; //Empties vector of "value types"


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
    NodeIterator() {}

    //
    //Operators
    //
 

     /** Defines dereferencing a node iterator.*/
    Node operator*() const{
    	Node node;
    	node.node_idx = node_it_idx;
    	node.node_graph_ptr = node_it_graph_ptr;
    	return node;
    }

    /** Defines incrimnenting a node iterator.*/
    NodeIterator& operator++(){
    	node_it_idx += 1;
    	return *const_cast<NodeIterator*>(this);
    }


     /** Defines equality for node iterators.*/
    bool operator==(const NodeIterator& iter) const{
    	if (node_it_idx == iter.node_it_idx){
    		if (node_it_graph_ptr == iter.node_it_graph_ptr){
    			return true;
    		}
    	}
    return false;

    }

   private:
    friend class Graph;
   int node_it_idx; // index to node in node_vec
   graph_type* node_it_graph_ptr; //pointer to graph 

  };

 
  //**Start of the ndoe iterator */
  node_iterator node_begin() const{
  	node_iterator iter;
  	iter.node_it_idx = 0;
  	iter.node_it_graph_ptr = const_cast<graph_type*>(this);

  	return iter;
  }

  //**End of the incident iterator */
  node_iterator node_end() const{
  	node_iterator iter;
  	iter.node_it_idx = node_vec.size();
  	iter.node_it_graph_ptr = const_cast<graph_type*>(this);
  	return iter;

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

    //
    //Operators
    //

     /** Defines dereferencing a incident iterator.*/
     Edge operator*() const {
     	Edge e;
     	e.edge_graph_ptr = (*init_node_ptr).node_graph_ptr;
     	e.edge_idx = init_itr->second;
     	e.main_node = init_node_ptr;
     	return e;
     }

     /** Defines incrimenting a node iterator.*/
     IncidentIterator& operator++(){
     	++init_itr;
     	return *const_cast<IncidentIterator*>(this);
     }

     /** Defines equality for node iterators.*/
     bool operator==(const IncidentIterator& instance) const {
     	if (init_node_ptr == instance.init_node_ptr){
     		if (init_itr== instance.init_itr){
     			return true;
     		}
     	}
     	return false;
     }

   private:
    friend class Graph;
    std::map<int, int>::iterator init_itr; //key in map in node_to_edge
    node_type* init_node_ptr;// pointer to main node

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

    //
   	//Operators
   	//


    /** Defines dereferencing for an EgeIterator. */
    Edge operator*() const {
    	Edge e;
    	e.edge_idx = edgeit_idx;
    	e.edge_graph_ptr = edgeit_graph_ptr;
    	return e;
    }

/** Defines incrimenting an EgeIterator. */
    EdgeIterator& operator++(){
    	edgeit_idx += 1;
    	return *const_cast<EdgeIterator*>(this);
    }

    /** Defines equality for two EdgeIterators. */
    bool operator==(const EdgeIterator& test) const{
    	if (edgeit_idx == test.edgeit_idx){
    		if(edgeit_graph_ptr == test.edgeit_graph_ptr){
    			return true;
    		}
    	}
    return false;
    }

   private:
    friend class Graph;
    int edgeit_idx;
    graph_type* edgeit_graph_ptr;


  };

//**Start of the edge iterator */
  edge_iterator edge_begin() const {
  	edge_iterator e;
  	e.edgeit_idx = 0;
  	e.edgeit_graph_ptr = const_cast<graph_type*>(this);
  	return e;
  }

//**End of the edge iterator */
  edge_iterator edge_end() const {
  	edge_iterator e;
  	e.edgeit_idx = edge_ends.size();
  	e.edgeit_graph_ptr = const_cast<graph_type*>(this);
  	return e;

  }

private:

 	std::vector<Point> node_vec; //vector of points
 	
 	std::vector<std::tuple<int, int>> edge_ends; //vector of tuples of an edge's nodes

 	// map of nodes to map of nodes that share an edge and the edge index
 	std::map<int, std::map<int, int>> node_to_edge;
 	
 	std::vector<node_value_type> value_vec; // vector of "value types"

};
//--style_1
//--The code compiles with warnings, please make sure your submitted code does not compile with warnings.
//--Otherwise, great job!
//--END
#endif // CME212_GRAPH_HPP
