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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  //--style_1
  //--Please place these private data members at the end of the class, where the
  //--comment indicates it should go.
  //--START
  //define the information of node structure
  struct node_info;
  //store the list of nodes 
  std::vector<node_info> node_list;

  //define the information of edge structure
  struct edge_info;
  //store the list of edges 
  std::vector<edge_info> edge_list;
  //--END
  //--design_0
  //--There is a better data structure choice for an adjacency list that will
  //--allow you to look up has_edge in average O(1) time. Look into std::map and
  //--std::unordered_map.
  //--START
  //store the list of adjacency
  std::vector<std::vector<unsigned>> adj_list;
  //--END

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;
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
  Graph(): node_list(),edge_list(),adj_list(){
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->node_list[idx_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
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

    // private constructor of Node

    Node(const graph_type* graph, size_type idx)
    :graph_(const_cast<graph_type*>(graph)), idx_(idx) {

    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_list.size();
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
  Node add_node(const Point& position,const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type id=num_nodes();
    node_list.push_back(node_info(position,id,value));
    //update the adjacency
    std::vector<size_type> empty;
    adj_list.push_back(empty);

    //notice that, after expanding the node_list, 
    //the length has been added up by 1
    //thus the index of the node should be the current size-1
    return Node(this, id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //(void) n;            // Quiet compiler warning
    if (num_nodes()<=n.index()){
      return false;
    }
    node_info node_in_list=node_list[n.index()];
    //--functionality_0
    //--Node position is not related to node identity; multiple nodes are allowed
    //--to occupy the same position. So the node position should have no role in
    //--checking whether n is a valid Node of this Graph. (I should have caught
    //--this on the last HW.) Please fix for next HW submission.
    //--START
    return node_in_list.p==n.position();
    //--END
}
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    assert(0 <= i && i < num_nodes());
    return Node(this, i);// Invalid node
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
      return Node(graph_, idx_a_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, idx_b_);      // Invalid Node
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
      //--style_0
      //--More concise way of writing this:
      //-- return (e.idx_a_ == idx_a_ and e.idx_b_ == idx_b_)
      //--     || (e.idx_a_ == idx_b_ and e.idx_b_ == idx_a_);
      //--START
      if (e.idx_a_==idx_a_ and e.idx_b_==idx_b_){
        return true;  
      } 
      else if (e.idx_a_==idx_b_ and e.idx_b_==idx_a_){
        return true;
      }else{
        return false;
      }
      //--END
      
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    //not sure about what does that mean...
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      return std::tie(graph_,idx_a_,idx_b_)<std::tie(e.graph_,e.idx_a_,e.idx_b_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
 	  graph_type* graph_;
    //size_type idx_edge_;
    size_type idx_edge_;
    size_type idx_a_;
    size_type idx_b_;
    
    //constructor
    Edge(const graph_type* graph, size_type eid, size_type ida,size_type idb)
     : graph_(const_cast<graph_type*>(graph)), idx_edge_(eid), idx_a_(ida),idx_b_(idb) {
    }
 };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i and i < num_edges());

    return Edge(this,i,edge_list[i].idx_a,edge_list[i].idx_b);// Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE    
    //--design_0
    //--In add_edge, you check to make sure the edge does not already exist.
    //--Interestingly, in add_edge, your existence check looks faster than it is
    //--here, since you limit your search to neighbors. Can you do that here?
    //--START
    Edge new_edge=Edge(this,0,a.index(),b.index());
    for (unsigned i = 0; i < num_edges(); i=i+1) {
      if (edge(i)==new_edge){
        return true;  
      }
    }
    return false;
    //--END
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
    /*hw1*/
    size_type ida=a.index();
    size_type idb=b.index();
    std::vector<size_type> init_nodes=adj_list[ida];

    for(auto inci_itr=a.edge_begin();inci_itr!=a.edge_end();++inci_itr){
      Edge old_inci=*inci_itr;
      // if there is such an edge
      if (old_inci.node2()==b){
        return old_inci;
      }
    }
    // if there isn't
    adj_list[ida].push_back(idb);
    adj_list[idb].push_back(ida);

    size_type id=num_edges();
    edge_list.push_back(edge_info(id,ida,idb));
    return Edge(this,id,ida,idb);  

  }



  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_list.clear();
    edge_list.clear();
    adj_list.clear();
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
        size_type node2_idx=graph_->adj_list[node1_idx_][inci_count_];
        return Edge(graph_,0,node1_idx_,node2_idx);
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
        size_type node1_idx=graph_->edge_list[edge_idx_].idx_a;
        size_type node2_idx=graph_->edge_list[edge_idx_].idx_b;
        return Edge(graph_,edge_idx_,node1_idx,node2_idx);
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
        edge_idx_++;
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
        return (graph_==edge_itr_new.graph_ && edge_idx_==edge_itr_new.edge_idx_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type edge_idx_;
    /** 
    * @brief private constructor for EdgeIterator.
    * 
    * @param[in] graph            A pointer pointing to the Graph we are trying to iterate through.
    * @param[in] edge_idx         The index for current iterated edge, initialized as 0.
    **/

    EdgeIterator(const graph_type* graph,const size_type edge_idx=0):graph_(const_cast<graph_type*>(graph)),edge_idx_(edge_idx){}

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
    return EdgeIterator(this,0);
  }

  /**
  *@ brief Return the iterator for the end of edges.

  *@ return An iterator for the end of edges.
  **/
  EdgeIterator edge_end() const{
    return EdgeIterator(this,edge_list.size());
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //define the basic structure of node, which should be private

  // in case furthur changes in node information, instead of 
  // adding private members directly in Node(),
  // I choose to use this structure to avoid abundant
  // information in Node()

  struct node_info {
      Point p;
      size_type i;
      node_value_type val;
      // Creat new node
      node_info(const Point &position, size_type index,node_value_type value):
      p(position), i(index),val(value){}

  };
  struct edge_info {
      //the index of the edge, start from 0
      size_type idx_edge;
      size_type idx_a;
      size_type idx_b;
      
      
      // Creat new node
      edge_info(size_type ei, size_type ida,size_type idb):
      idx_edge(ei), idx_a(ida),idx_b(idb){}

  };

};

//--style_0
//--Code formatting could be improved, breaking up long lines, inserting missing
//--blank lines and removing unnecessary blank lines, removing irrelevant
//--comments, making spacing consistent, etc., generally making the code easier
//--to read. This will be a 1-point style deduction on the next HW.
//--END

#endif // CME212_GRAPH_HPP
