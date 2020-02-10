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
    adjacency_list.push_back(std::vector<size_type>());
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

    /** @brief Return a node of this Edge */
    Node node1() const {
      return Node(graph_, nuid1_);
    }

    /** @brief Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, nuid2_);
    }

//--functionality_1
//--the graph is undirected, so as long as both nodes are the same this should return true
//--START
    /** @brief Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((e.nuid1_ == nuid1_)
        and(e.nuid2_ == nuid2_)
        and(e.graph_ == graph_)){return true;}
      return false;
    }
//--END
    /** @brief Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
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
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this -> edges_vector.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // We check that our nodes are in the graph:
    if ((has_node(a)) and (has_node(b)))
    {
      std::vector<size_type> a_neighbours = adjacency_list[a.uid_];
      //Then we go through our list of neighbours;
      for(size_type i = 0; i < a_neighbours.size(); i++)
      {
          if (a_neighbours[i] == b.uid_){return true;}
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    //If the edge already exist we return it
    if (has_edge(a,b)){
      return Edge(this, a.uid_, b.uid_);
    }
    //Otherwise we add it in the adjacency list
    adjacency_list[a.uid_].push_back(b.uid_);
    adjacency_list[b.uid_].push_back(a.uid_);

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
  }

//--documentation_1
//--no pre or post conditions in the doxygen comments
//--END

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
      size_type nuid2 = graph_ -> adjacency_list[nuid1_][uid_];
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

 private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //STL containers. I use vectors to store the data and a mapping to associate ids 
  //with positon in those vectors
  std::vector <Internal_node> nodes_vector;
  std::vector <Internal_edge> edges_vector; 

  //Adjacency list of our graph.
  std::vector <std::vector <size_type>> adjacency_list;

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
