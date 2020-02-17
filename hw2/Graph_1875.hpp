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
template <typename V, typename E> 
class Graph {
 private:
  struct Internal_node;
  struct Internal_edge;
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** An optional value to our nodes */
  using node_value_type = V;

  /** An optional value to our edges */
  using edge_value_type = E;

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
  class Node: private totally_ordered<Node>{
   public:
    /** @brief Construct an invalid node.
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
      //Nothing to do, we just create a void node.
    }

    /** @brief Return this node's position. 
    @pre This node is valid and is an element of the graph*/
    Point& position() {
      return graph_->nodes_vector[uid_].pos;
    }

    /** @brief Return this node's position. 
    @pre This node is valid and is an element of the graph*/
    const Point& position() const {
      return graph_->nodes_vector[uid_].pos;
    }

    /** @brief Return this node's value as a reference. 
    @pre This node is valid and is an element of the graph*/
    node_value_type& value(){
       return graph_->nodes_vector[uid_].val ; 
    }

    /** @brief Return this node's value as a constant. 
    @pre This node is valid and is an element of the graph*/
    const node_value_type& value() const {
       return graph_->nodes_vector[uid_].val ; 
    }


    /** @bried Return this node's index, a number in the range [0, graph_size). 
    @pre This node is valid*/
    size_type index() const {
      return uid_;
    }

    /** @brief return the number of incident edges 
    @pre This node is valid*/
    size_type degree() const {
      return graph_-> adjacency_list[uid_].size();
    }

    /** @brief Return the start of the incident iterator */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, uid_, 0);
    }

    /** @brief Return the end of incident iterator */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ((n.uid_ == uid_) and (n.graph_ == graph_)){return true;}
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
      if (n.uid_ < uid_){return true;}
          return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the graph the node lies in
    Graph* graph_;

    // This node's current index number
    // Like in the example proxy class given, our proxy node is defined with a 
    // unique id
    size_type uid_;

    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->nodes_vector.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** @brief Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val 
      = node_value_type()) {
    size_type new_id = nodes_vector.size();
    nodes_vector.push_back(Internal_node(new_id,position,val));
    adjacency_list.push_back(std::vector<std::pair<size_type,size_type>>());
    edge_value_list.push_back(std::vector<edge_value_type>());
    return Node(this, new_id);
  }

  /** @brief Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //We want to check that the node lies in this graph
    return (n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
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
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      //The edge is void so there is nothing to do
    }

    /** @brief Return the value of this Edge 
    	@pre The edge is in the graph and is valid
    */
    edge_value_type& value(){
      std::vector<std::pair<size_type,size_type>> neighbours = graph_ -> adjacency_list[nuid1_];
      //Then we go through our list of neighbours;
      size_type i = 0;
      while(neighbours[i].first != nuid2_) i++;
      return graph_ -> edge_value_list[nuid1_][i];
    }
    
    /** @brief Return the value of this Edge 
    	@pre The edge is in the graph and is valid
    */
    const edge_value_type& value() const{
	  std::vector<std::pair<size_type,size_type>> neighbours = graph_ -> adjacency_list[nuid1_];
      //Then we go through our list of neighbours;
      size_type i = 0;
      while(neighbours[i].first != nuid2_) i++;
      return graph_ -> edge_value_list[nuid1_][i];
    }

    /** @brief Return a node of this Edge */
    Node node1() const {
      return Node(graph_, nuid1_);
    }

    /** @brief Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, nuid2_);
    }

    /** @brief Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * Complexity: O(1)
     */
    bool operator==(const Edge& e) const {
      if ((e.nuid1_ == nuid1_)
        and(e.nuid2_ == nuid2_)
        and(e.graph_ == graph_)){return true;}
      if ((e.nuid2_ == nuid1_)
        and(e.nuid1_ == nuid2_)
        and(e.graph_ == graph_)){return true;}
      return false;
    }

    /** @brief Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * Complexity: O(1)
     */
    bool operator<(const Edge& e) const {
      if (e.graph_ == graph_){
        if (e.nuid1_ == nuid1_){
          return (e.nuid2_ < nuid2_);
        } else {
        return (e.nuid1_ < nuid1_);
      }
      } else {
        return (e.graph_ < graph_);
      }
      return (e.graph_ < graph_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    //Pointer to the overall graph. Once again edge acts as a proxy.
    Graph* graph_ ;

    //The id of the adjacent nodes.
    size_type nuid1_;
    size_type nuid2_;

    // A private constructor for Edge class. 
    Edge(const Graph* edge_graph, size_type n1, size_type n2)
     : graph_(const_cast<Graph*>(edge_graph)), nuid1_(n1), nuid2_(n2){
     }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return this -> edges_vector.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    assert(i<edges_vector.size());
    Internal_edge edge = edges_vector[i];
    return Edge(this, edge.nuid1, edge.nuid2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(max_degree)
   */
  bool has_edge(const Node& a, const Node& b) const {
    // We check that our nodes are in the graph:
    if ((has_node(a)) and (has_node(b)))
    {
      std::vector<std::pair<size_type,size_type>> a_neighbours = adjacency_list[a.uid_];
      //Then we go through our list of neighbours;
      for(size_type i = 0; i < a_neighbours.size(); i++)
      {
          if (a_neighbours[i].first == b.uid_){return true;}
      }
      return false;
    }
    else{
      return false;
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
   * Complexity: O(max_degree)
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val 
      = edge_value_type()) {
    //If the edge already exist we return it
    if (has_edge(a,b)){ //Here we have a complexity of O(max_degree)
      return Edge(this, a.uid_, b.uid_);
    }

    size_type index = edges_vector.size();

    //Otherwise we add it in the adjacency list
    std::pair<size_type, size_type> pairb(b.uid_, index);
    adjacency_list[a.uid_].push_back(pairb);
    edge_value_list[a.uid_].push_back(val);

    std::pair<size_type, size_type> paira(a.uid_, index);
    adjacency_list[b.uid_].push_back(paira);
    edge_value_list[b.uid_].push_back(val);

    //Then we add our edge in our vector of edges
    edges_vector.push_back(Internal_edge(a.uid_, b.uid_));
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_vector.clear() ;
    edges_vector.clear();
    adjacency_list.clear();
    edge_value_list.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator>{
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

    /** @brief Dereference the current node and return it **/
    Node operator*() const{
      return Node(graph_, uid_);
    }

    /** @brief Increment the iterator **/
    NodeIterator& operator++(){
      uid_++;
      return *this;
    }

    /** @brief Test the equality between 2 iterators **/
    bool operator==(const NodeIterator& n) const{
      return ((graph_ == n.graph_) and (uid_ == n.uid_));
    }

   private:
    friend class Graph;
    //The graph the iterator lies in
    Graph* graph_;
    //The id of the current nodes it points to
    size_type uid_;

    /** @brief Basic constructor of the iterator **/
    NodeIterator(const Graph* g, size_type u)
      :graph_(const_cast<Graph*>(g)), uid_(u){}
  };

  /** @brief Return an iterator pointing at the beginning of our node_vector **/
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** @brief Return an iterator pointing at the end of our node_vector **/
  node_iterator node_end() const{
    return NodeIterator(this, nodes_vector.size());
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

    /** @brief Dereference the current node and return it **/
    Edge operator*() const{
      size_type nuid2 = graph_ -> adjacency_list[nuid1_][uid_].first;
      return Edge(graph_, nuid1_, nuid2);
    }

    /** @brief Increment the iterator **/
    IncidentIterator& operator++(){
      uid_++;
      return *this;
    }

    /** @brief Test the equality between 2 iterators **/
    bool operator==(const IncidentIterator& n) const{
      return ((graph_ == n.graph_)
        and (nuid1_ == n.nuid1_)
        and (uid_ == n.uid_));
    }

   private:
    friend class Graph;
    Graph* graph_; //The graph the iterator lies in
    size_type nuid1_; //The node this iterator is based on
    size_type uid_; //The current position in the adjacency list

    /** @brief Basic constructor of the iterator **/
    IncidentIterator(const Graph* g, size_type nuid1, size_type u)
      :graph_(const_cast<Graph*>(g)), nuid1_(nuid1), uid_(u){}
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

    /** @brief Dereference the current node and return it **/
    Edge operator*() const{
      size_type nuid1 = graph_->edges_vector[uid_].nuid1;
      size_type nuid2 = graph_->edges_vector[uid_].nuid2;
      return Edge(graph_, nuid1, nuid2);
    }

    /** @brief Increment the iterator **/
    EdgeIterator& operator++(){
      uid_++;
      return *this;
    }

    /** @brief Test the equality between 2 iterators **/
    bool operator==(const EdgeIterator& n) const{
      return ((graph_ == n.graph_) and (uid_ == n.uid_));
    }

   private:
    friend class Graph;
    Graph* graph_; //The graph the iterator lies in
    size_type uid_; //Our position in the vector of edges

    /** @brief Basic constructor of the iterator **/
    EdgeIterator(const Graph* g, size_type u)
      :graph_(const_cast<Graph*>(g)), uid_(u){}
  };

  /** @brief Return an iterator pointing at the beginning of our edge_vector **/
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** @brief Return an iterator pointing at the end of our edge_vector **/
  edge_iterator edge_end() const{
    return EdgeIterator(this, edges_vector.size());
  }

  /** @brief Remove the edge 
  * @return a boolean stating if the edge was removed
  * @post has_edge(e.node1(),e.node2()) = false;
  * @post if the edge was in the graph, it has been removed
  * 
  * This method invalidates all iterators concerning the incident
  * of the nodes adjacent to this edge
  *
  * Furhtermore, it moves adjacency_list[e.nuid1][adjacency_list[e.nuid1].size()-1]
  * and adjacency_list[e.nuid2][adjacency_list[e.nuid2].size()-1], so anything related to
  * those edges is now invalidate
  *
  * Complexity : O(degree)
  */
  size_type remove_edge(const Edge &e){
    if(!has_edge(e.node1(), e.node2())) return false;

    size_type n1 =e.node1().index();
    size_type n2 =e.node2().index();

    //We clean edges_vector. Complexity : O(degree)

    //Search for this edge in the adjacency list
    for(size_type i=0; i<adjacency_list[n1].size(); i++){
      if (adjacency_list[n1][i].first==n2){
        //Find its index in edge_vector
        size_type index = adjacency_list[n1][i].second;

        //Send it at the back of edge_vector
        edges_vector[index] = edges_vector[edges_vector.size()-1];

        //One edge has been moved, so we need some updates
        //We update the index in the adjacency list
        size_type n1_to_update = edges_vector[index].nuid1;
        size_type n2_to_update = edges_vector[index].nuid2;

        //Search it in the adjacency list of n1_to_update
        for(size_type j=0; j<adjacency_list[n1_to_update].size(); j++){
          if (adjacency_list[n1_to_update][j].first==n2_to_update){
            //Update the index
            adjacency_list[n1_to_update][j].second = index;
          }
        }
        //Search it in the adjacency list of n2_to_update
        for(size_type j=0; j<adjacency_list[n2_to_update].size(); j++){
          if (adjacency_list[n2_to_update][j].first==n1_to_update){
            //Update the index
            adjacency_list[n2_to_update][j].second = index;
          }
        }
        //Finally remove the edge
        edges_vector.pop_back();
      }
    }

    //We clean adjacency_list and edge_value_list. Complexity: O(degree)
    //Search for the edge in the adjacency list of n1
    for(size_type i=0; i<adjacency_list[n1].size(); i++){
      if (adjacency_list[n1][i].first==n2){
        //By swapping 2 edges we invalidate the one at the back
        adjacency_list[n1][i] = adjacency_list[n1][adjacency_list[n1].size()-1];
        adjacency_list[n1].pop_back();
        edge_value_list[n1][i] = edge_value_list[n1][edge_value_list[n1].size()-1];
        edge_value_list[n1].pop_back();
      }
    }

    //Search for the edge in the adjacency list of n1
    for(size_type i=0; i<adjacency_list[n2].size(); i++){
      if (adjacency_list[n2][i].first==n1){
        //By swapping 2 edges we invalidate the one at the back
        adjacency_list[n2][i] = adjacency_list[n2][adjacency_list[n2].size()-1];
        adjacency_list[n2].pop_back();
        edge_value_list[n2][i] = edge_value_list[n2][edge_value_list[n2].size()-1];
        edge_value_list[n2].pop_back();
      }
    }
    return !(has_edge(e.node1(),e.node2()));
  }

  /** @brief Remove the edge 
  * @return a boolean stating if the edge was removed
  * @post has_edge(e.node1(),e.node2()) = false;
  * @post if the edge was in the graph, it has been removed
  * 
  * This method invalidates all iterators concerning the incident
  * of the nodes adjacent to this edge
  *
  * Furhtermore, it moves adjacency_list[n1.uid][adjacency_list[n1.uid].size()-1]
  * and adjacency_list[n2.uid][adjacency_list[n2.uid].size()-1], so anything related to
  * those edges is now invalidate
  *
  * Complexity : O(degree)
  */
  size_type  remove_edge(const  Node &n1, const  Node &n2){
    return remove_edge(Edge{this, n1.index(), n2.index()});
  }

  /** @brief Remove the edge 
  * @return a valid iterator
  * @post has_edge(e.node1(),e.node2()) = false;
  * @post if the edge was in the graph, it has been removed
  * 
  * This method invalidates all iterators concerning the incident
  * of the nodes adjacent to this edge
  *
  * Furhtermore, it moves adjacency_list[(*e_it).nuid1][adjacency_list[(*e_it).nuid1].size()-1]
  * and adjacency_list[(*e_it).nuid2][adjacency_list[(*e_it).nuid2].size()-1], so anything related to
  * those edges is now invalidate
  *
  * Complexity : O(degree)
  */
  edge_iterator  remove_edge(edge_iterator  e_it){
    remove_edge((*e_it));
    return e_it;
  }

  /** @brief Remove the node, and all edges incident to it.
  *   @return a boolean stating if the node was removed 
  *   @post if the node was in the graph it has been removed
  *   @post all edges incident to the node have been removed
  *   Invalidate: The node with index num_nodes-1 now have index 
  *               n.index()
  *               So any edge or iterator related to the removed node
  *               or to the node that has its index changed is now invalid
  *
  *   Complexity: O(degree^2)
  */
  size_type  remove_node(const  Node &n){
    if (!has_node(n)){
      return false;
    }
    size_type uid = n.index();
    size_type size = num_nodes();

    //Clean node_vector
    nodes_vector[uid] = nodes_vector[size-1];
    nodes_vector.pop_back();

    //First we have to erase the incident edges.
    //Complexity O(degree^2)

    size_type deg = n.degree();
    std::vector<size_type> neighbour_vect;

    for(size_type i=0; i < deg; i++){
      neighbour_vect.push_back(adjacency_list[uid][i].first);
    }

    for(auto it = neighbour_vect.begin(); it != neighbour_vect.end(); ++it){
      remove_edge(Edge(this, uid, *it));
    }

    //Now clean this node adjacency list
    adjacency_list[uid] = adjacency_list[size-1];
    adjacency_list.pop_back();
    edge_value_list[uid] = edge_value_list[size-1];
    edge_value_list.pop_back();

    //Now, we have to update all the node of edges concerning the node 
    //that has been moved to the index uid, if it occurs (if uid != size-1)
    if (uid != size-1){
      for(size_type i = 0; i < adjacency_list[uid].size(); i++){
        size_type neighbour = adjacency_list[uid][i].first;
        //Search it in the adjacency list of its neighbour
        for(size_type j = 0; j < adjacency_list[neighbour].size(); j++){
          if(adjacency_list[neighbour][j].first == (size-1)){
            adjacency_list[neighbour][j].first = uid;
            //Update it also in edges_vector
            edges_vector[adjacency_list[neighbour][j].second].nuid1 = neighbour;
            edges_vector[adjacency_list[neighbour][j].second].nuid2 = uid;
          }
        }
      }
      }
    return (!has_node(n));
  }

   /** @brief Remove the node, and all edges incident to it.
  *   @return a boolean stating if the node was removed 
  *   @post if the node was in the graph it has been removed
  *   @post all edges incident to the node have been removed
  *   Invalidate: The node with index num_nodes-1 now have index 
  *               n.index()
  *               So any edge or iterator related to the removed node
  *               or to the node that has its index changed is now invalid
  *
  *   Complexity: O(max_degree^2)
  */
  node_iterator  remove_node(node_iterator  n_it){
    remove_node((*n_it));
    return n_it;
  }

 private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //STL containers. I use vectors to store the data and a mapping to associate ids 
  //with positon in those vectors
  std::vector <Internal_node> nodes_vector;
  std::vector <Internal_edge> edges_vector; 

  //Adjacency list of our graph. Every point is a pair <neighbour, index in edge_vector>
  std::vector <std::vector <std::pair<size_type,size_type>>> adjacency_list;

  //List where we store our edge values.
  std::vector <std::vector <edge_value_type>> edge_value_list;

  // Disable copy and assignment of a graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  /** @brief A struct to represent our node */
  struct Internal_node {
    //The id defining the node
    size_type uid; //The id of the node
    Point pos; //The position of our node
    node_value_type val ; //The value stored in our node

    /** @brief basic constructor */
    Internal_node(size_type uid,Point pos, node_value_type val) 
      :uid(uid),pos(pos), val(val){}

   };

   /** @brief A struct to represent our edges */
   struct Internal_edge {
    size_type nuid1; //The id of an adjacent node
    size_type nuid2; //The id of an adjacent node

    /** @brief basic constructor */
    Internal_edge(size_type nuid1, size_type nuid2) 
      :nuid1(nuid1),nuid2(nuid2) {}

   };
};

#endif // CME212_GRAPH_HPP
