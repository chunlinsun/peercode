#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
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

  /** Key proxy class for storing the identity (a unique 
  random number that is different for any two nodes) to 
  uniquely identity a node, and a Point class object p
  for storing the location of the node. I use a proxy 
  class because this can be easily extended for more
  complex purposes.*/
  class Node_set{
  	friend class Graph;
    size_type identity;
    Point p;
    Node_set(size_type identity, Point p){
      this->identity = identity;
      this->p = p;
    }
  };

  /**
  Key proxy class for storing the index of the jth neighbor of the
  ith node. I use a proxy class because proxy class can be easily
  extended for more complex purposes. 
  */
  class Edge_set{  
  public:
    Edge_set(){}
    Edge_set(size_type neighbor_index){
      this->neighbor_index = neighbor_index;
    }
  private:
    friend class Graph;
    size_type neighbor_index;
  };

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
  class Node {
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

    Node(size_type identity, const graph_type* Graph_node){
      this->identity = identity;
      this->Graph_node = const_cast<graph_type*>(Graph_node);
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return Graph_node->Node_sets[this->Graph_node->identity_index[identity]].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return Graph_node->identity_index[this->identity];
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return this->Graph_node==n.Graph_node && this->identity==n.identity;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
	 *
     * Since two different nodes have different identities from the way
     * they are constructed and the defintion of identity, comparing their 
     * identities define a global comparison relationship.
     */
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
    	return this->identity<n.identity;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // identity is a unique number that is different for any two nodes
    // when they are added to the graph.
    size_type identity;
    // The graph that the node is on.
    graph_type* Graph_node;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return Node_sets.size();
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    size_type random_number = random();
    while (this->identity_index.count(random_number)){
      random_number = random();
    }
    size_type index = this->size();
    this->identity_index[random_number] = index;
    this->Node_sets.push_back(Node_set(random_number, position));
    this->Edge_sets.push_back(std::vector<Edge_set>(0));
    return Node(random_number,this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   * A node is on the graph if and only if the identity of 
   * the node is a key of the identity_index of the graph.
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this->identity_index.count(n.identity);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(Node_sets[i].identity, this);
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    Edge(const graph_type* Graph_edge, size_type identity_1, size_type identity_2){
      this->identity_1 = identity_1;
      this->identity_2 = identity_2;
      this->Graph_edge = const_cast<graph_type*>(Graph_edge);
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(identity_1, this->Graph_edge);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(identity_2, this->Graph_edge);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->identity_1== e.identity_1) && (this->identity_2== e.identity_2) && (this->Graph_edge== e.Graph_edge);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->identity_1!=e.identity_1){
        return this->identity_1<e.identity_1;
      }
      else {
      	return this->identity_2<e.identity_2;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type identity_1;
    size_type identity_2;
    graph_type* Graph_edge;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type num_edges = 0;
    size_type num_nodes = this->Edge_sets.size();
    for (size_type i = 0; i<num_nodes; i++){
      num_edges+= this->Edge_sets[i].size();
    }
    return num_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
   * First of all, find the first node index j of the ith edge, and then 
   * find which neighbor the ith edge corresponds to the jth node. After
   * find the 
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type index_1 = 0;
    while (i>=Edge_sets[index_1].size()){
      i-= Edge_sets[index_1].size();
      index_1++;
    }
    return Edge(this, Node_sets[index_1].identity, Node_sets[Edge_sets[index_1][i].neighbor_index].identity);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
   * 
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type index_1 = a.index();
    size_type index_2 = b.index();
    if (index_1>index_2){
      size_type temp = index_2;
      index_2 = index_1;
      index_1 = temp;
    }
    for (size_type i= 0; i<Edge_sets[index_1].size(); i++){
    	if (Edge_sets[index_1][i].neighbor_index== index_2){
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
   * When I add edges, I make sure that I did not double count the edges. 
   * In Edge_sets, I make sure that Edge_sets[i], which is the neighbors of 
   * the ith node, the index of all of its neighbors are larger than i. 
   */
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    size_type index_1 = this->identity_index[a.identity];
    size_type index_2 = this->identity_index[b.identity];
    if (index_1>index_2){
      size_type temp = index_2;
      index_2 = index_1;
      index_1 = temp;
    }
    Edge_sets[index_1].push_back(Edge_set(index_2));
    return Edge(this, a.identity, b.identity);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    identity_index.clear();
    Node_sets.clear();
    Edge_sets.clear();
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
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

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
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:
  /** identity_index is the mapping from identity of a node
  to the index of the node on the graph. Since identity is 
  unique for a node, this unordered_map is going to work.
  **/
  std::unordered_map<size_type, size_type> identity_index;
  /** Node_sets are a vector of Node_set, where Node_set is 
  a proxy struct that stores the identity and position of the node.
  Node_sets[i], i is the ith node in terms of index on the graph.**/
  std::vector<Node_set> Node_sets;
  /**
  Edge_sets are a vector of Edge_set, where Edge_set is a proxy class
  that stores the neighbor index of a node. Edge_sets[i][j] is the ith 
  node in terms of index on the graph, and j is the jth neighbor of the 
  ith node. 
  **/
  std::vector<std::vector<Edge_set>> Edge_sets;

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
