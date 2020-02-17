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
template <typename V,typename E>
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
  using graph_type = Graph<V,E>;
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
  Graph(): node_list(),num_edge_(0),adj_list(),i2u_(){
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
    Node() {
      // HW0: YOUR CODE HERE
      // HW2 updated: notice that idx_ in Node() is the unique identifier, while i in node_info is the current index.
    }

    /** 
     * @brief Return this node's position, which can be modefied. This is for homework2.
     **/
    Point& position() {
    	assert(valid());
      return graph_->node_list[idx_].p;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->node_list[idx_].p;
    }

    /** Return this node's index, a number in the range [0, size of i2u_). */
    //hw2 changed
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->node_list[idx_].i;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

   /** 
   * @brief  Return the value associated with the
             node we are asking, which can be changed.
   **/
    node_value_type& value() {
      return graph_->node_list[idx_].val;
    }

   /** 
   * @brief  Return the value associated with the
             node we are asking, which can't be changed.
   **/    
    const node_value_type& value() const {
      return graph_->node_list[idx_].val;
    }

   /** 
   * @brief  Return the number of edges adjancent to the node.
   **/    
    size_type degree() const{
      return graph_->adj_list[idx_].size();
    }

   /** 
   * @brief  Return incident_iterator corresponding to the first edge
             'explored' that incident to the node.
   **/    
    IncidentIterator edge_begin() const{
      return IncidentIterator(graph_,idx_,0);

    }
  /** 
   * @brief  Return incident_iterator marking the end of visiting 
             incident edges of the node.
            .
   **/   
    IncidentIterator edge_end() const{
      return IncidentIterator(graph_,idx_,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (n.graph_==graph_ && n.idx_==idx_);
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
      return idx_<n.idx_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Here, Node has two private members: the graph it lies in, 
    // as well as its index.
    graph_type* graph_;
    size_type idx_;

    /** 
     * @brief 	Boolean. Determine wwhether a node is  currently valid.
     * @returns True if 0 <= idx_ < node_list.size() && 0 <= i < i2u_.size()
     *				 && i2u_[node_list[idx_].i] == idx_;
     *               i.e it's currently valid and used to be in the graph, also in sync
     **/
    // HW2 added
    bool valid() const {
      return 0<=idx_ && idx_ < graph_->node_list.size() && graph_->node_list[idx_].i >=0 &&
             graph_->node_list[idx_].i < graph_->i2u_.size() && graph_->i2u_[graph_->node_list[idx_].i] == idx_;
    }

    // private constructor of Node

    Node(const graph_type* graph, size_type idx)
    :graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  //HW2 changed
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  // HW2 changed
  Node add_node(const Point& position,const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type unique_id=node_list.size();
    node_list.push_back(node_info(position,num_nodes(),value));// num_nodes() indicating the size of current nodes+1
    //update the i2u_ list
    i2u_.push_back(unique_id);
    //update the adjacency
    std::vector<edge_info> empty;
    adj_list.push_back(empty);

    //notice that, after expanding the node_list, 
    //the length has been added up by 1
    //thus the index of the node should be the current size-1
    return Node(this, unique_id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  // HW2 changed
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (num_nodes()<=n.index()){
      return false;
    }
    //n.index() get the current index of n
    //i2u_[n.index()] get the unique index of n
    node_info node_in_list=node_list[i2u_[n.index()]];

    return node_in_list.p==n.position();
}
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  //HW2 changed
  //Given a 'user known' index i, return the Node with attribute unique index.
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // HW2 CHANGED
    return Node(this, i2u_[i]);
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

    //initialization of your graph.
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, idx_a_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, idx_b_);
    }

    /** Return the length of the edge */
    double length() const{
    	return norm(node1().position()-node2().position());
    }

    /**
    * @brief return the reference of edge value, which can be modified.
    * @pre there is such an edge.
    * @post if edge_1==edge_2, then edge_1.value()==edge_2.value()
	**/
	//hw2 added
	edge_value_type& value(){
		size_type min_idx=std::min(node1().idx_,node2().idx_);
		size_type max_idx=std::max(node1().idx_,node2().idx_);
		for (size_type i=0;i<graph_->adj_list[min_idx].size();++i){
			size_type adj_idx=graph_->adj_list[min_idx][i].idx_other;
			if(adj_idx==max_idx){
				return(graph_->adj_list[min_idx][i].val);
			}
		}
	}

	/**
    * @brief return the reference of edge value, which CAN'T be modified.
    * @pre there is such an edge.
    * @post if edge_1==edge_2, then edge_1.value()==edge_2.value()
	**/
	//hw2 added
	const edge_value_type& value() const{
		return value();
	}

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      if (e.graph_!=graph_){
        return false;
      }

      if (e.idx_a_==idx_a_ and e.idx_b_==idx_b_){
        return true;  
      } 
      else if (e.idx_a_==idx_b_ and e.idx_b_==idx_a_){
        return true;
      }else{
        return false;
      }
      
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    //HW2 changed
    bool operator<(const Edge& e) const {
      size_type min_idx_e=std::min(e.idx_a_,e.idx_b_);
      size_type min_idx_this=std::min(idx_a_,idx_b_);
      size_type max_idx_e=std::max(e.idx_a_,e.idx_b_);
      size_type max_idx_this=std::max(idx_a_,idx_b_);
      if (graph_==e.graph_){
      	if (min_idx_this==min_idx_e){
      		return max_idx_this<max_idx_e;
      	}
      	return min_idx_this<min_idx_e;
      }
      return graph_<e.graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
 	graph_type* graph_;
    size_type idx_a_;
    size_type idx_b_;
    
    //constructor
    Edge(const graph_type* graph,size_type ida,size_type idb)
     : graph_(const_cast<graph_type*>(graph)), idx_a_(ida),idx_b_(idb) {
    }
 };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edge_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
  	size_type count=0;
  	for (auto edge_itr=this->edge_begin(); edge_itr!=this->edge_end();++edge_itr){
  		if (count==i){
  			return *edge_itr;
  		}else{
  			count++;
  		}
  	}
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE    
    auto adj=adj_list[a.idx_];
    for (auto k=adj.begin();k!=adj.end();++k){
    	if ((*k).idx_other==b.idx_){
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
  Edge add_edge(const Node& a, const Node& b,const edge_value_type& value=edge_value_type()) {
    size_type ida=a.idx_;
    size_type idb=b.idx_;

    for(auto inci_itr=a.edge_begin();inci_itr!=a.edge_end();++inci_itr){
      Edge old_inci=*inci_itr;
      // if there is such an edge
      if (old_inci.node2()==b){
        return old_inci;
      }
    }
    // if there isn't
    size_type id_min=std::min(ida,idb);
    size_type id_max=std::max(ida,idb);
    adj_list[id_min].push_back(edge_info(id_max,value));
    adj_list[id_max].push_back(edge_info(id_min,edge_value_type()));
    num_edge_++;
    return Edge(this,ida,idb);  

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_list.clear();
    num_edge_=0;
    adj_list.clear();
    i2u_.clear();
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

    /** 
    * @brief Get access to a node via node_itr_idx_ 

    * @return the node in graph, 'explored' by iterator
    **/
    Node operator*() const{
        return graph_->node(node_itr_idx_);
    }

    /** 
    * @brief  Move the iterator to the next node or if 
              currently at the last node, then move to 1 past 
              the end of the node, which should be node_end(). 

    * @return the node iterator for next node or 1 past the 
              end of the node.
    **/
    NodeIterator& operator++(){
        node_itr_idx_++;
        return *this;
    }

    /** 
    * @brief  Check the equality between two node iterators.

    * @param[in] node_itr_new   The iterator we want to
                                 compare with. 
    * @return Boolean. True if they have the same iterator
                       index in the same graph.

    **/
    bool operator==(const NodeIterator& node_itr_new) const{
        return (graph_==node_itr_new.graph_ && node_itr_idx_==node_itr_new.node_itr_idx_); 
    }


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    // index for node interator
    size_type node_itr_idx_;

    /** 
    * @brief private constructor for NodeIterator.
    * 
    * @param[in] graph            A pointer, pointing to the Graph we are trying to iterate through.
    * @param[in] node_itr_idx     The index for node iterator. Initialized as 0.
    **/
    NodeIterator(const graph_type* graph, const size_type node_itr_idx=0):graph_(const_cast<graph_type*>(graph)),node_itr_idx_(node_itr_idx){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
  /**
  *@ brief Return the iterator for the starting of nodes.

  *@ return An iterator for the starting of nodes.
  **/
  NodeIterator node_begin() const{
    return NodeIterator(this,0);
  }

  /**
  *@ brief Return the iterator for the end of nodes.

  *@ return An iterator for the end of nodes.
  **/
  NodeIterator node_end() const{
    return NodeIterator(this,num_nodes());
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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
    
    /** 
    * @brief Get access to the incident edge via inci_count_. 

    * @return the incident in graph, currently 'explored' by iterator
    **/
    Edge operator*() const{
        size_type node2_idx=graph_->adj_list[node1_idx_][inci_count_].idx_other;
        return Edge(graph_,node1_idx_,node2_idx);
    }

    /** 
    * @brief  Move the iterator to the next incident edge or if 
              currently at the last incident edge of node1, 
              then move to 1 past the end of the adjacency edges,
              which should be equal to the degree of node1 (root node). 

    * @return the incident iterator for next incident edge or 1 past the 
              end of the incident edges of node1.
    **/
    IncidentIterator& operator++(){
        inci_count_++;
        return *this;
    }

    /** 
    * @brief  Check the equality between two incident iterators.

    * @param[in] inci_itr_new   The iterator we want to
                                 compare with. 
    * @return Boolean. True if they have the same pair of iterator
                       indexes in the same graph.

    **/
    bool operator==(const IncidentIterator& inci_itr_new) const{
        return (graph_==inci_itr_new.graph_ && node1_idx_==inci_itr_new.node1_idx_ && inci_count_==inci_itr_new.inci_count_); 
    }
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type node1_idx_;
    size_type inci_count_;
    
    /** 
    * @brief private constructor for IncidentIterator.
    * 
    * @param[in] graph            A pointer pointing to the Graph we are trying to iterate through.
    * @param[in] node1_idx        The index for root node
    * @param[in] inci_count       Number of incident edges iterated. Initialized as 0.
    **/

    IncidentIterator(const graph_type* graph,const size_type node1_idx,const size_type inci_count=0):graph_(const_cast<graph_type*>(graph)),node1_idx_(node1_idx),inci_count_(inci_count){}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */

  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /** 
    * @brief Get access to the current edge via edge_idx_. 

    * @return the edge in graph, currently 'explored' by iterator
    **/
    Edge operator*() const{
        return (*incident_itr_);
    }

    /** 
    * @brief  Move the iterator to the next edge or if 
              currently at the last edge, 
              then move to 1 past the end of the edges,
              which should be equal to total amount of edges.

    * @return the edge iterator for next edge or 1 past the 
              end of the edges in the graph.
    **/
    EdgeIterator& operator++(){
        while (node_itr_!=graph_->node_end()){
        	auto root=(*node_itr_);
        	while (incident_itr_!=root.edge_end()){
        		++incident_itr_;
        		if(incident_itr_==root.edge_end()){
        			break;
        		//This is to avoid double counting. We restrict the edge iterator goes through 
        		// edges that node1.idx_<node2.idx_.
        		}else if (root.idx_<=(*incident_itr_).node2().idx_){
        			return *this;	
        		}
        	}
        	++node_itr_;
        	//re-judge
        	if (node_itr_ != graph_->node_end()) {
        		incident_itr_ = (*node_itr_).edge_begin();
        	}
        }
        return *this;
        
    }

    /** 
    * @brief  Check the equality between two edge iterators.

    * @param[in] edge_itr_new   The iterator we want to
                                 compare with. 
    * @return Boolean. True if they have the same iterator
                       indexes in the same graph.

    **/
    bool operator==(const EdgeIterator& edge_itr_new) const{
        return (graph_==edge_itr_new.graph_ && node_itr_==edge_itr_new.node_itr_);
    }


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    NodeIterator node_itr_;
    IncidentIterator incident_itr_;
    /** 
    * @brief private constructor for EdgeIterator.
    * 
    * @param[in] graph            A pointer pointing to the Graph we are trying to iterate through.
    * @param[in] edge_idx         The index for current iterated edge, initialized as 0.
    **/

    EdgeIterator(const graph_type* graph,NodeIterator node_itr, IncidentIterator incident_itr):
    				graph_(const_cast<graph_type*>(graph)),node_itr_(node_itr),incident_itr_(incident_itr){}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
  *@ brief Return the iterator for the starting of edges.

  *@ return An iterator for the starting of edges.
  **/
  EdgeIterator edge_begin() const{
  	node_iterator node_itr_start=node_begin();
  	incident_iterator incident_itr_start=(*node_itr_start).edge_begin();
    return EdgeIterator(this,node_itr_start,incident_itr_start);
  }

  /**
  *@ brief Return the iterator for the end of edges.

  *@ return An iterator for the end of edges.
  **/

  EdgeIterator edge_end() const{
  	node_iterator node_itr_end=node_end();
  	incident_iterator incident_itr_end=(*node_begin()).edge_begin();
    return EdgeIterator(this,node_itr_end,incident_itr_end);
  }

  //HW2 added
  /**
  * @brief  Given pair of nodes, remove the edge between them from the graph.
  	return -1 if not removed successfully, 1 o.w. Also change the corresponding
  	property

  * 
  * @param[in,out] node1       one end point of the edge
  * @param[in,out] node2       the other end point of the edge
  
  * @return 1 if the edge is sucessfully removed. 0 otherwise.
  * @post 	If successfully removed(return1) then 
  			(1) current num_edges()=old num_edges()-1;
  			(2) has_edge(node1,node2) == false;
  			(3) all other edges remain.
  			(4) Invalidate corresponding iterators.
  * Complexity: O(d), d is the overall maximum degree of nodes in the graph.
  				This is because we only go though the adjacent edges of a node.
  **/
  size_type remove_edge(const Node& node1,const Node& node2){
  	if (!has_edge(node1,node2)){
  		return 0;
  	}
  	size_type degree1=node1.degree();
  	size_type degree2=node2.degree();
  	
  	num_edge_-=1;

  	for (size_type i=0;i<degree1;++i){
  		if (adj_list[node1.idx_][i].idx_other==node2.idx_){
  			adj_list[node1.idx_].erase(adj_list[node1.idx_].begin()+i);
  			break;
  		}
  	}
  	for (size_type i=0;i<degree2;++i){
  		if (adj_list[node2.idx_][i].idx_other==node1.idx_){
  			adj_list[node2.idx_].erase(adj_list[node2.idx_].begin()+i);
  			break;
  		}
  	}

  	assert(node1.degree()==degree1-1);
  	assert(node2.degree()==degree2-1);

  	return 1;
  }

  /**
  * @brief  Given specific edge, remove it.
  		return 0 if not removed successfully, 1 o.w. 
  		Also change the corresponding property.

  * 
  * @param[in,out] edge       The edge we need to remove
  
  * @return 1 if the edge is sucessfully removed. -1 otherwise.
  * @post 	If successfully removed(return1) then 
  			(1) current num_edges()=old num_edges()-1;
  			(2) has_edge(edge.node1(),edge.node2()) == false;
  			(3) all other edges remain.
  			(4) Invalidate corresponding iterators.
  * Complexity: O(d), d is the overall maximum degree of nodes in the graph.
  				This is because we call the remove_edge function on nodes.
  **/
  size_type remove_edge(const Edge& edge){
  	return remove_edge(edge.node1(),edge.node2());
  }

/**
  * @brief  Given specific edge_iterator, remove the edge the iterator pointing to.
  		return 0 if not removed successfully, 1 o.w. 
  		Also change the corresponding property.

  * 
  * @param[in,out] e_it       The edge_iterator pointing to the edge we need to remove
  
  * @return 1 if the edge is sucessfully removed. -1 otherwise.
  * @post 	If successfully removed(return1) then 
  			(1) current num_edges()=old num_edges()-1;
  			(2) has_edge(*e_it) == false;
  			(3) all other edges remain.
  			(4) Invalidate corresponding iterators.
  * Complexity: O(d), d is the overall maximum degree of nodes in the graph.
  				This is because we call the remove_edge function on edge.
  **/
  edge_iterator remove_edge(edge_iterator e_it){
  	Node root=*e_it.node_itr_;
  	incident_iterator inci_itr=e_it.incident_itr_;
  	//remove the whole edge
  	remove_edge(*e_it);
  	//notice that we need to check the boundary case
  	if (inci_itr==root.edge_end()){
  		e_it++;//move further, to the next node.
  	}
  	return e_it;
  }

  /**
  * @brief  Given a specific node, remove it from the graph. 
  			Also change the corresponding property.

  * 
  * @param[in,out] node   		the node need to be removed.
  
  * @return 1 if the node is sucessfully removed. 0 otherwise.
  * @post 	If successfully removed(return1) then 
  			(1) current num_nodes()=old num_nodes()-1;
  			(2) has_node(node) == false;
  			(3) all other nodes remain.
  			(4) Invalidate corresponding iterators.
  *
  * Complexity: O(num_nodes()), because we will adjust remaining nodes.
  **/
  size_type remove_node(const Node& node){
  	if (!has_node(node)){
  		return 0;
  	}
  	size_type old_num_node=num_nodes();
  	size_type idx=node.index();
  	
  	//Adjust the corresponding properties before removing node
  	//Adjust indexes of remaining nodes in the graph.
  	for (auto id=idx+1;id<old_num_node;++id){
  			size_type id_uid=i2u_[id];
  			node_list[id_uid].i=id-1;
  		}
  	
  	//Remove edges incident to node.
  	//First make a copy to store edges need to be removed.
  	std::vector<edge_info> removed_edges=adj_list[node.idx_];
  	size_type node_degree=node.degree();
  	for (size_type i=0;i<node_degree;i++){
  		size_type inci_node=removed_edges[i].idx_other;
  		remove_edge(node,Node(this,inci_node));
  	}
  	//remove the node
  	i2u_.erase(i2u_.begin()+idx);

  	assert(!has_node(node));
  	assert(num_nodes()==old_num_node-1);
  	return 1;	
  }

  /**
  * @brief  Given a specific node_iterator, remove the node
  			pointed by it from the graph. 
  			Also change the corresponding property.

  * 
  * @param[in,out] n_it   		the node iterator, pointing to the 
  *								node that need to be removed.
  
  * @return 1 if the node is sucessfully removed. -1 otherwise.
  * @post 	If successfully removed(return1) then 
  			(1) current num_nodes()=old num_nodes()-1;
  			(2) has_node(*n_it) == false;
  			(3) all other nodes remain.
  			(4) Invalidate corresponding iterators.
  			(5) This iterator pointing to the next node.
  *
  * Complexity: O(num_nodes()), because we call remove_node on node.
  **/
  node_iterator remove_node(node_iterator n_it){
  	remove_node(*n_it);
  	return n_it;
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //define the basic structure of node, which should be private
   //define the information of node structure
  struct node_info;
  // Used to store nodeinfo for any node added, even if later removed.
  std::vector<node_info> node_list;
  //define the information of edge structure
  struct edge_info;
  //store the number of edges
  unsigned num_edge_;
  //store the information of adjancency
  std::vector<std::vector<edge_info>> adj_list;
  //Store the currently "active" set of nodes.
  std::vector<unsigned> i2u_; //Indexed by node idx
  struct node_info {
      Point p;// position of node
      size_type i;// CURRENT index of node
      node_value_type val;// value of node
      node_info(const Point &position, size_type index,node_value_type value):
      p(position), i(index),val(value){}

  };
  struct edge_info {
      size_type idx_other;
      edge_value_type val;
      
      edge_info(size_type idx,edge_value_type value):
      idx_other(idx), val(value){}

  };

};

#endif // CME212_GRAPH_HPP
