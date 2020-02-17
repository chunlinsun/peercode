#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include<cmath>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>


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
  using graph_type = Graph<V, E>;

  /**Node's user specified vale*/
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


  struct nodeinfo {
  int idx_;       //< Node index
  Point p_;     //< Node position
  V v_;	         //Node Value
  nodeinfo(int idx, Point p, V v) : idx_(idx), p_(p), v_(v) {}
};

struct edgeinfo {
  int idx_1;       //< Node1 index
  int idx_2;     //< Node2 position
  edge_value_type e_;	         //edge Value
  edgeinfo(int idx1, int idx2, edge_value_type e) : idx_1(idx1), idx_2(idx2), e_(e) {}
};

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {

  	std::vector<nodeinfo> node_vec; // makes an empty vector of node info
  	std::vector<int> i2u_; //empty vector of current node ids
  	std::vector<edgeinfo> edge_ends; //makes an empty vector ofedge info
 	std::map<int, std::map<int, int>> node_to_edge; //makes an empty map of maps


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

    /** Returns this node's position and allows for position to be changed. */
    Point& position (){
    	return node_graph_ptr-> node_vec[node_graph_ptr->i2u_[node_idx]].p_;
    }

    /** Returns this node's position. */
    const Point& position() const {
    	return node_graph_ptr-> node_vec.at(node_graph_ptr->i2u_.at(node_idx)).p_;    	
    }


    /** Returns this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_idx;
    }

    /** Returns this node's value. */
    node_value_type& value() {
    	return node_graph_ptr->node_vec[node_graph_ptr-> i2u_[node_idx]].v_;

    }

    /** Returns this node's value for constant node */
    const node_value_type& value() const {
    	return &(node_graph_ptr->node_vec[node_graph_ptr-> i2u_.at(node_idx)].v_);

    }

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
    int node_idx; //index in i2u_

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
    //adds point to node_vec method in Graph class
  
	nodeinfo info(node_vec.size(), position, value);
  	node_vec.push_back(info);
  	i2u_.push_back(node_vec.size()-1);

  	//adds a map in the node_to edge map, used later if node connects an edge
	std::map<int, int> m;
  	node_to_edge[i2u_.size()-1]=m;

  	//creates a node to be returned
  	Node n;
  	n.node_idx = i2u_.size()-1;
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
  		if (n.node_idx < int(i2u_.size())){
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

    edge_value_type& value() {
    	return edge_graph_ptr->edge_ends[edge_idx].e_;

    }

    /** Returns this edges's value for constant edge */
    const edge_value_type& value() const {
    	return &(edge_graph_ptr->edge_ends.at(edge_idx).e_);
    }


    /** Returns the length of the edge using the position of the nodes */
    double length() const{
    	Point pos1 = edge_graph_ptr-> node_vec.at(edge_graph_ptr->i2u_.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_1)).p_;
    	Point pos2 = edge_graph_ptr-> node_vec.at(edge_graph_ptr->i2u_.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_2)).p_;
    	double length = norm(pos1-pos2);
    	return length;
    }

    /** Returns a node of this Edge */
    Node node1() const {
    	Node n;
    	if (main_node == nullptr){
    		n.node_graph_ptr = edge_graph_ptr;
    		n.node_idx = edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_1).idx_;
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
    		n.node_idx = edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_2).idx_;
  		}
  		else{
  			if ((*main_node).node_idx != edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_1).idx_){
  				n.node_idx = edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_1).idx_;
  			}
  			else{
  				n.node_idx = edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_2).idx_;

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
    	if (edge_graph_ptr == e.edge_graph_ptr){
    	if (edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(e.edge_idx).idx_1).idx_ != \
    		edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_1).idx_){
    		if (edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(e.edge_idx).idx_1).idx_ != \
    			edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_2).idx_){
    			return false;
    		}
    		else{
    			if (edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(e.edge_idx).idx_2).idx_ != \
    				edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_1).idx_){
    				return false;
    			}
    			else{
    				return true;
    			}
    		}

    	}
    	//compares @a nodes to second node of edge
    	else{
    		if (edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(e.edge_idx).idx_2).idx_ != \
    		 edge_graph_ptr->node_vec.at(edge_graph_ptr->edge_ends.at(edge_idx).idx_2).idx_){
    			return false;
    		}
    		else{
    			return true;
    		}
    	}
    }
    else{
    	return false;
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
      	if (edge_idx> e.edge_idx){
      		return false;
  		}
  		else{
  			if (e.edge_graph_ptr > edge_graph_ptr){
  				return true;
  			}
  			else{
  				return false;
  			}
  		}
    }
}

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

	graph_type* edge_graph_ptr; //pointer to the graph
  	int edge_idx; //index of edge in edge_ends
 	node_type* main_node = nullptr; //pointer to node from the incident iterator 
  	
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
  	//checks if edge already in graph structures
  	if (has_edge(a,b) == false){

  		//edge doesn't already exist
  		//updates graph structures
  		edgeinfo info(i2u_[a.node_idx], i2u_[b.node_idx], value);
    	edge_ends.push_back(info);

    	//value_vec_edge.push_back(value);
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
    node_vec = {}; // Empties the vector of node info
    i2u_ = {}; //Empties the vector of unique indexs
  	edge_ends = {}; //Empties the vector of edge info
 	node_to_edge = {}; //Empties the map of maps


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
    	++node_it_idx;
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
   int node_it_idx; // index to node in i2u_
   graph_type* node_it_graph_ptr; //pointer to graph 
  };

 
  //**Start of the node iterator */
  node_iterator node_begin() const{
  	node_iterator iter;
  	iter.node_it_idx = 0;
  	iter.node_it_graph_ptr = const_cast<graph_type*>(this);

  	return iter;
  }

  //**End of the node iterator */
  node_iterator node_end() const{
  	node_iterator iter;
  	iter.node_it_idx = i2u_.size();
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

     /** Defines dereferencing an incident iterator.*/
     Edge operator*() const {
     	Edge e;
     	e.edge_graph_ptr = (*init_node_ptr).node_graph_ptr;
     	e.edge_idx = init_itr->second;
     	e.main_node = init_node_ptr;
     	return e;
     }

     /** Defines incrimenting an incident iterator.*/
     IncidentIterator& operator++(){
     	++init_itr;
     	return *const_cast<IncidentIterator*>(this);
     }

     /** Defines equality for indcident iterators.*/
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
    	++edgeit_idx;
    	return *const_cast<EdgeIterator*>(this);
    }

    /** Defines equality for two EdgeIterators. */
    bool operator==(const EdgeIterator& test) const{
    	if (edgeit_idx== test.edgeit_idx){
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
	
/** 
 * @brief Removes a node and all incident edges from a graph
 * @param[in] n An instance of the node class.
 * @post Number of nodes decreases by 1
 * @post Edges incident to the node are removed
 * @post The incident of the node in node_vec is now invalid
 * @return True if the node is in the graph, false otherwise
 */
  size_type remove_node ( const Node & n){

		//checks that node is in graph
		if (has_node(n) ==true){
			IncidentIterator inc = n.edge_begin();
			while (inc != n.edge_end()){
				remove_edge(*inc);
				++inc;
			}

			int old_idx = i2u_.size()-1;

			//checks if removing the end node		
			if (old_idx!=n.node_idx){
			std::swap(i2u_.at(n.node_idx), i2u_.back());
			i2u_.pop_back();

			int index;
			int edge_i;
			std::map<int, int>::iterator it;

			std::map<int, int>::iterator itr = node_to_edge.at(old_idx).begin();
			while (itr != node_to_edge.at(old_idx).end()){
				index = itr->first;
				edge_i = itr->second;

				it = node_to_edge.at(index).find(old_idx); 
	    		node_to_edge[index].erase(it);


				// node_to_edge[index].erase(old_idx);
				node_to_edge[index][n.node_idx]= edge_i;
				++itr;
			}
			std::map<int, int> can = node_to_edge.at(old_idx);
			node_to_edge[n.node_idx] = can;

			std::map<int, std::map<int, int>>::iterator it1 = node_to_edge.find(old_idx); 
	    	node_to_edge.erase(it1);


			//node_to_edge.erase(old_idx);
			node_vec[i2u_[n.node_idx]].idx_ = n.node_idx;

			return true;
		}
			//removes if already last node
			else{
			i2u_.pop_back();
			std::map<int, std::map<int, int>>::iterator it1 = node_to_edge.find(old_idx); 
	    	node_to_edge.erase(it1);

			//node_to_edge.erase(old_idx);

			return true;

			}
	}
	//if the node isn't in the graph
	else{
		return false;
		}
	}

/** 
 * @brief Removes a node and all incident edges from a graph
 * @param[in] n_it  An instance of the node_iterator class.
 * @pre n_it references a valid node in the graph and is not an end pointer
 * @post number of nodes decreases by 1
 * @post The incident of the ndoe in node_vec is now invalid
 * @post Edges incident to the node are removed
 * @return A node iterator that points to the same location in i2u_,
 * which now references the node that was previously in the back position. 
 * If it previously pointed to the back posision, it now points to i2u_.end()
 */

	node_iterator remove_node ( node_iterator n_it ){

		//removes edges connected to that node
		IncidentIterator inc = (*n_it).edge_begin();
		while (inc != (*n_it).edge_end()){
			remove_edge(*inc);
			++inc;	
		}

		int old_idx = i2u_.size()-1;

		//checks if this already points to last node
		if (old_idx != (*n_it).node_idx){
			std::swap(i2u_.at((*n_it).node_idx), i2u_.back());
			i2u_.pop_back();

			int index;
			int edge_i;
			std::map<int, int>::iterator it;

			std::map<int, int>::iterator itr = node_to_edge.at(old_idx).begin();
			while (itr != node_to_edge.at(old_idx).end()){
				index = itr->first;
				edge_i = itr->second;
				it = node_to_edge.at(index).find(old_idx);
				node_to_edge[index].erase(it);

				node_to_edge[index][(*n_it).node_idx] = edge_i;
				++itr;
			}
			node_to_edge[(*n_it).node_idx] = node_to_edge[old_idx];

			std::map<int, std::map<int, int>>::iterator it1 = node_to_edge.find(old_idx); 
	    	node_to_edge.erase(it1);

			node_vec[i2u_[(*n_it).node_idx]].idx_ = (*n_it).node_idx;
			return n_it;
		}

		//if node is in back position
		else{
		i2u_.pop_back();
		std::map<int, std::map<int, int>>::iterator it1 = node_to_edge.find(old_idx); 
	    node_to_edge.erase(it1);
		return n_it;
		}
	}

/** 
 * @brief Removes an edge from a graph
 * @param[in] a  An instance of the node class, represents one end of an edge
 * @param[in] b  An instance of the node class. representing the other end of an edge
 * @post number of edges decreases by 1
 * @post Edge index in node_to_edge is updated to reflect new order of edge_ends
 * @return True if the edge between the two nodes was in the graph, False otherwise
 */


	size_type remove_edge ( const Node & a, const Node & b){
		//checks if node in graph
		if (has_edge(a,b) == true){
			int e_index = node_to_edge[a.node_idx][b.node_idx];

			std::map<int, int>::iterator it = node_to_edge[a.node_idx].find(b.node_idx); 
    		node_to_edge[a.node_idx].erase(it);

			it = node_to_edge[b.node_idx].find(a.node_idx); 
    		node_to_edge[b.node_idx].erase(it);


			//checks if edge in last position
			if(e_index != int(edge_ends.size()-1)){
				std::swap(edge_ends[e_index], edge_ends.back());
				edge_ends.pop_back();	
				int node1 = node_vec[edge_ends[e_index].idx_1].idx_;
				int node2 = node_vec[edge_ends[e_index].idx_2].idx_;
				node_to_edge[node1][node2] = e_index;
				node_to_edge[node2][node1] = e_index;
			}

			//removes edge in last position
			else{
				edge_ends.pop_back();
			}
			return true;
		}
		//edge isn't in graph
		else{
			return false;}

/** 
 * @brief Removes an edge from a graph
 * @param[in] e  An instance of the edge class.
 * @post Mumber of edges decreases by 1
 * @post Edge index in node_to_edge is updated to reflect new order of edge_ends
 * @return True if the edge was in the graph, False otherwise
 */


	}
	size_type remove_edge ( const Edge & e){
		
		Node a;
		a.node_idx = node_vec.at(edge_ends.at(e.edge_idx).idx_1).idx_;
		a.node_graph_ptr = e.edge_graph_ptr;
		Node b;
		b.node_idx = node_vec.at(edge_ends.at(e.edge_idx).idx_2).idx_;
		b.node_graph_ptr = e.edge_graph_ptr;

		//checks if edge is in the graph
		if (has_edge(a,b) == true){

			int node1 = node_vec.at(edge_ends.at(e.edge_idx).idx_1).idx_;
			int node2 = node_vec.at(edge_ends.at(e.edge_idx).idx_2).idx_;

			std::map<int, int>::iterator it = node_to_edge.at(node1).find(node2); 
    		node_to_edge[node1].erase(it);

			it = node_to_edge[node2].find(node1); 
    		node_to_edge[node2].erase(it);

			//checks if edge in last postion
			if(e.edge_idx != int(edge_ends.size()-1)){
			std::swap(edge_ends.at(e.edge_idx), edge_ends.back());
			edge_ends.pop_back();

			node1 = node_vec.at(edge_ends.at(e.edge_idx).idx_1).idx_;
			node2 = node_vec.at(edge_ends.at(e.edge_idx).idx_2).idx_;

			node_to_edge[node1][node2] = e.edge_idx;
			node_to_edge[node2][node1] = e.edge_idx;
			}
			//edge in last position
			else{
			edge_ends.pop_back();
			}

			return true;
		}
		//edge not in graph
		else{
			return false;}
	}

	/** 
 * @brief Removes an edge from a graph
 * @param[in] e_it  An instance of the edge_iterator class.
 * @pre e_it references a valid edge in the graph and is not an end pointer
 * @post number of nodes decreases by 1
 * @post Edge index in node_to_edge is updated to reflect new order of edge_ends
 * @return A edge iterator that points to the same location in edge_ends,
 * which now references the node that was previously in the back position. 
 * If it previously pointed to the back posision, it now points to edge_ends.end()
 */


	edge_iterator remove_edge ( edge_iterator e_it ){

		int node1 = node_vec.at(edge_ends.at((*e_it).edge_idx).idx_1).idx_;
		int node2 = node_vec.at(edge_ends.at((*e_it).edge_idx).idx_2).idx_;

		std::map<int, int>::iterator it = node_to_edge.at(node1).find(node2); 
    	node_to_edge[node1].erase(it);

		it = node_to_edge.at(node2).find(node1); 
    	node_to_edge[node2].erase(it);

		//checks that edge is not in last position
		if((*e_it).edge_idx != int(edge_ends.size()-1)){
			std::swap(edge_ends.at((*e_it).edgeit_idx), edge_ends.back());
			edge_ends.pop_back();
			node1 = node_vec.at(edge_ends.at((*e_it).edgeit_idx).idx_1).idx_;
			node2 = node_vec.at(edge_ends.at((*e_it).edgeit_idx).idx_2).idx_;
			node_to_edge[node1][node2] = (*e_it).edgeit_idx;
			node_to_edge[node2][node1] = (*e_it).edgeit_idx;
		}

		//edge in last position
		else{
			edge_ends.pop_back();
		}
		return e_it;

		}

private:

 	std::vector<nodeinfo> node_vec; //vector of nodeinfo, uses uid

 	//Stores currently active nodes
 	std::vector<int> i2u_;   // Indexed by node idx
 	
 	std::vector<edgeinfo> edge_ends; //vector of edge info

 	// map of nodes to map of nodes that share an edge and the edge index
 	std::map<int, std::map<int, int>> node_to_edge;

 	
};

#endif // CME212_GRAPH_HPP
