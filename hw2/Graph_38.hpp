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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;
  struct internal_edge;

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
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  using edge_value_type = E;

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
  std::vector<internal_node> nodes {}; //default constructor for the std::vector type
  std::vector<internal_edge> edges {};
  std::vector<std::vector<std::pair<size_type, size_type>>> adjacency_matrix {};
  std::vector<size_type> i2u {};
  std::vector<size_type> edges_i2u {};
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
      graph_pointer = nullptr; //does not point to any graph
      unique_id = 0;
    }

    /**
    * Constructor for a node_type belonging to a certain graph with a certain uid
    * (main reason I made this public is so that the class mass_spring::PinConstraint
    * could use it)
    *
    * @param[in] graph  pointer to the graph that the node_type belongs to
    * @param[in] uid  unique_id of the node_type in question
    *
    */
    Node(const graph_type* graph, size_type uid){ //constructor used by the graph 
      graph_pointer = const_cast<graph_type*>(graph);
      unique_id = uid;
      }

    /**
    * specifies whether this node_type is valid meaning it currently belongs to a graph_type
    *
    * @return     True if the node_type is valid false otherwise
    *
    */
    bool valid() const {
      return (graph_pointer != nullptr && unique_id >= 0 && unique_id < graph_pointer->nodes.size() && graph_pointer->nodes[unique_id].index_ < graph_pointer->i2u.size() && graph_pointer->i2u[graph_pointer->nodes[unique_id].index_] == unique_id);
    }

    /** Return this node's position. 
    * @pre this node_type is valid()
    */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(valid());
      return graph_pointer->nodes[unique_id].position_; //position is stored in internal node in graph
    }

    /** Return this node_type's position. 
    * @pre this node_type is valid()
    */
    Point& position() {
      assert(valid());
      return graph_pointer->nodes[unique_id].position_; //position is stored in internal node in graph
    }

    /** Return this node's index, a number in the range [0, graph_size). 
    * @pre this node_type is valid()
    */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(valid());
      return graph_pointer->nodes[unique_id].index_;
    }

    /** return the node_type's unique_id. */
    size_type uid() const {
      return unique_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value();
    // size_type degree();
    // incident_iterator edge_begin();
    // incident_iterator edge_end();

    /**
    * @brief supplies the value associated with this node
    * 
    * @return the value by reference associated with this node
    * @pre this node belongs to a graph and is valid
    */
    node_value_type& value() {
      assert(valid());
      return graph_pointer->nodes[unique_id].value_;
    }

    /**
    * @brief supplies the value associated with this node
    * 
    * @return the value associated with this node as a constant reference
    * @pre this node belongs to a graph and is valid
    */
    const node_value_type& value() const {
      assert(valid());
      return const_cast<const node_value_type>(graph_pointer->nodes[unique_id].value_);
    }
    
    /**
    * @brief supplies the degree of this node (number of edges connected to it)
    * 
    * @return the degree of this node
    * @pre this node belongs to a graph and is valid
    */
    size_type degree() const {
      assert(valid());
      size_type counter = 0;
      for(auto eiter = edge_begin(); eiter != edge_end(); ++eiter) {
        ++counter;
      }
      return counter;
    }

    /**
    * @brief returns an incident_iterator pointing to the first edge coming out of this node_type
    *
    * @return incident_iterator pointing to fist edge of this node_type
    * @post if this node_type is invalid or does not have any edges, result may not be dereferenced
    */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_pointer, unique_id, graph_pointer->adjacency_matrix[unique_id].begin());
    }


    /**
    * @brief returns an incident_iterator pointing to one past the last edge coming out of this node_type
    *
    * @return incident_iterator pointing to one past the last edge of this node_type
    * @post result may not be dereferenced
    */
    incident_iterator edge_end() const {
      return incident_iterator(graph_pointer, unique_id, graph_pointer->adjacency_matrix[unique_id].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (valid() && n.valid() && graph_pointer == n.graph_pointer && unique_id == n.unique_id);
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
      if(graph_pointer == n.graph_pointer)
        return (this->index() < n.index());
      else
        return graph_pointer < n.graph_pointer;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_pointer;
    size_type unique_id;
        
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return i2u.size();
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
    // HW0: YOUR CODE HERE
    internal_node new_node {position, static_cast<size_type>(i2u.size()), value};
    i2u.push_back(nodes.size());
    nodes.push_back(new_node);
    adjacency_matrix.push_back(std::vector<std::pair<size_type, size_type>>());
    return node(new_node.index_);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_pointer == this && n.valid()); 
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0 && i < size());
    assert(nodes[i2u[i]].index_ == i); //graph invariant
    return node_type(this, i2u[i]);        
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
      graph_pointer = nullptr; //does not belong to any graph 
      node1_uid = 0;
      node2_uid = 0;
      unique_id = 0;
    }

     /**
    * specifies whether this edge_type is valid meaning it currently is an edge_type in the graph_type specified by graph_pointer
    *
    * @return     True if the edge_type is valid false otherwise
    *
    */
    bool valid() const {
       return (graph_pointer != nullptr && node_type(graph_pointer, node1_uid).valid() && node_type(graph_pointer, node2_uid).valid() && unique_id >= 0 && unique_id < graph_pointer->edges.size() && graph_pointer->edges[unique_id].index_ < graph_pointer->edges_i2u.size() 
        && graph_pointer->edges_i2u[graph_pointer->edges[unique_id].index_] == unique_id);
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return node_type(graph_pointer, node1_uid);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return node_type(graph_pointer, node2_uid);       
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (valid() && e.valid() && graph_pointer == e.graph_pointer && unique_id == e.unique_id);
    }


    /**
    * @brief supplies the value associated with this edge_type
    * 
    * @return the value by reference associated with this edge_type
    * @pre this edge_type belongs to a graph_type and is valid
    */
    edge_value_type& value(){
      assert(valid());
      return graph_pointer->edges[unique_id].value_;
    }

    /**
    * @brief supplies the value associated with this edge_type
    * 
    * @return the value by reference associated with this edge_type
    * @pre this edge_type belongs to a graph_type and is valid
    */
    const edge_value_type& value() const {
      assert(valid());
      return graph_pointer->edges[unique_id].value_;
    }

    /** 
    * returns the length of the edge_type i.e. the euclidean distance between node1() and node2(). 
    *
    * @return l2 norm of node1().position() and node2().position()
    *
    * @pre the edge_type is valid()
    */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return this edges's index, a number in the range [0, graph_pointer->num_edges()). 
    * @pre the edge_type is valid()
    */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(valid());
      return graph_pointer->edges[unique_id].index_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if(graph_pointer == e.graph_pointer)
        return (index() < e.index()); 
      else
        return graph_pointer < e.graph_pointer;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type node1_uid;
    size_type node2_uid;
    size_type unique_id;
    graph_type* graph_pointer;

    /**
    * Edge private constructor (used mainly by the graph_type, incident_iterator class)
    *
    * @param[in] graph  pointer to the the graph_type
    * @param[in] node1  unique_id of node1
    * @param[in] node2  unique_id of node2
    * @param[in] euid   unique_id of the edge_type in question
    */
    Edge(const graph_type* graph, size_type node1, size_type node2, size_type euid){ //constructor used by the graph
      graph_pointer = const_cast<graph_type*>(graph);
      node1_uid = node1;
      node2_uid = node2;
      unique_id = euid;
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_i2u.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0 && i < num_edges());
    assert(edges[edges_i2u[i]].index_ == i); //graph invariant
    return edge_type(this, edges[edges_i2u[i]].node_1_uid_, edges[edges_i2u[i]].node_2_uid_, edges_i2u[i]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_pointer == b.graph_pointer && a.graph_pointer == this && a.valid() && b.valid()); //verify pre-conditions
    return (std::find_if(a.edge_begin(), a.edge_end(), [&] (edge_type e) {return b.unique_id == e.node2_uid;}) != a.edge_end());
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
  edge_type add_edge(const Node& a, const Node& b,  const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    assert(a.graph_pointer == b.graph_pointer && a.graph_pointer == this && a.valid() && b.valid()); //verify pre-conditions
    edge_type return_edge;
    incident_iterator is_it_the_edge = std::find_if(a.edge_begin(), a.edge_end(), [&] (edge_type e) {return b.unique_id == e.node2_uid;});
    if (is_it_the_edge != a.edge_end()) {
      return_edge = *is_it_the_edge;
    } else {
      adjacency_matrix[a.unique_id].push_back(std::pair<size_type, size_type>(b.unique_id, edges.size()));
      adjacency_matrix[b.unique_id].push_back(std::pair<size_type, size_type>(a.unique_id, edges.size()));
      internal_edge new_edge {a.unique_id, b.unique_id, static_cast<size_type>(edges_i2u.size()), value};
      edges_i2u.push_back(edges.size());
      edges.push_back(new_edge);
      return_edge =  Edge(this, a.unique_id, b.unique_id, new_edge.index_) ;
    }
    return return_edge;
  }

  /**
  * removes a node_type from the graph_type
  * 
  * @param[in] node   node_type to remove
  * @return   1 if remove was succesful in case node_type is valid, otherwise 0
  *
  * @post all graph_type invariants are preserved
  * 
  * Complexity: O(1)
  */
  size_type remove_node(const node_type& node) {
     if(node.valid()) {
      for(auto e = node.edge_begin(); e != node.edge_end(); ++e)
        remove_edge(*e);
      size_type removed_index {node.index()};
      node.graph_pointer->i2u[removed_index] = node.graph_pointer->i2u.back();
      node.graph_pointer->nodes[node.graph_pointer->i2u.back()].index_ = removed_index;
      node.graph_pointer->i2u.pop_back();
      return 1;
    } else {
      return 0;
    }
  }

  /**
  * removes a node_type from the graph_type
  * 
  * @param[in] n_it   node_iterator pointing to node_type to remove
  * @return   node_iterator pointing to the same location but now to the node_type that follows the removed node_type
  *
  * @post naturally all graph_type invariants are preserved
  * 
  * Complexity: O(1)
  */
  node_iterator remove_node(node_iterator n_it) {
    if((*n_it).valid()) {
      for(auto e = (*n_it).edge_begin(); e != (*n_it).edge_end(); ++e)
        remove_edge(*e);
      size_type removed_index  {(*n_it).index()};
      n_it.graph_pointer->i2u[removed_index] = n_it.graph_pointer->i2u.back();
      n_it.graph_pointer->nodes[n_it.graph_pointer->i2u.back()].index_ = removed_index;
      n_it.graph_pointer->i2u.pop_back();
    }
   return n_it;
  }

  /**
  * removes an edge_type from the graph_type
  * 
  * @param[in] edge   edge_type to remove
  * @return   1 if remove was succesful in case edge_type is valid, otherwise 0
  *
  * @post all graph_type invariants are preserved
  * 
  * Complexity: O(1)
  */
  size_type remove_edge(const edge_type& edge) {
    if(edge.valid()) {
      size_type removed_index {edge.index()};
      edge.graph_pointer->edges_i2u[removed_index] = edge.graph_pointer->edges_i2u.back();
      edge.graph_pointer->edges[edge.graph_pointer->edges_i2u.back()].index_ = removed_index;
      edge.graph_pointer->edges_i2u.pop_back();
      return 1;
    } else {
      return 0;
    }
  }


  /**
  * removes an edge_type from the graph_type
  * 
  * @param[in] n1, n2   node_types connecting the edge_type to remove
  * @return   1 if remove was succesful in case edge_type is valid and belongs to graph_type, otherwise 0
  *
  * @post all graph_type invariants are preserved
  * 
  * Complexity: O(1)
  */
  size_type remove_edge(const node_type& n1, const node_type& n2) {
    assert(n1.graph_pointer == n2.graph_pointer && n1.graph_pointer == this); //verify pre-conditions
    incident_iterator is_it_the_edge = std::find_if(n1.edge_begin(), n1.edge_end(), [&] (edge_type e) {return n2.unique_id == e.node2_uid;});
    if (is_it_the_edge != n1.edge_end()) {
      size_type removed_index {(*is_it_the_edge).index()};
      n1.graph_pointer->edges_i2u[removed_index] = n1.graph_pointer->edges_i2u.back();
      n1.graph_pointer->edges[n1.graph_pointer->edges_i2u.back()].index_ = removed_index;
      n1.graph_pointer->edges_i2u.pop_back();
      return 1;
    } else {
      return 0;
    }
  }

  /**
  * removes an edge_type from the graph_type
  * 
  * @param[in] e_it   edge_iterator pointing to the edge_type to remove
  * @return   edge_iterator pointing to the same location but now to the edge_type that follows the removed edge_type
  *
  * @post all graph_type invariants are preserved
  * 
  * Complexity: O(1)
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    if((*e_it).valid()) {
      size_type removed_index {(*e_it).index()};
      e_it.graph_pointer->edges_i2u[removed_index] = e_it.graph_pointer->edges_i2u.back();
      e_it.graph_pointer->edges[e_it.graph_pointer->edges_i2u.back()].index_ = removed_index;
      e_it.graph_pointer->edges_i2u.pop_back();
    }
    return e_it;
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear(); //stl clear method for std::vector class
    edges.clear();
    adjacency_matrix.clear();
    i2u.clear();
    edges_i2u.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator :  private totally_ordered<NodeIterator> {
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
      internal_node_iterator = typename std::vector<size_type>::const_iterator();
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
    * @brief dereferences the NodeIterator to obtain the corresponding node_type
    * 
    * @return node_type behind the NodeIterator
    * @pre this node_iterator points to an actual node_type in a non-empty graph_type
    */
    value_type operator*() const {
      assert(graph_pointer != nullptr && internal_node_iterator!= graph_pointer->i2u.end());
      return graph_pointer->node(graph_pointer->nodes[*internal_node_iterator].index_);
    }

    /**
    * @brief increments node_iterator to point to the next node_type
    * 
    * @return node_iterator reference corresponding to the next node_type
    * @pre this node_iterator points to an actual node_type in a non-empty graph_type
    * @post (*old node_iterator).index() + 1 == (*new node_iterator).index()
    */
    node_iterator& operator++() {
      assert(graph_pointer != nullptr && internal_node_iterator != graph_pointer->i2u.end());
      ++internal_node_iterator;
      return (*this);
    }

    /**
    * @brief checks whether this node_iterator is equal to another node_iterator nodeiter,
    * which is equivalent to checking equality of the underlying node_types
    * 
    * @param[in] nodeiter   other NodeIterator object to compare with this NodeIterator
    * @return true if the NodeIterator's are equal false otherwise
    */
    bool operator==(const NodeIterator& nodeiter) const {
      return (graph_pointer == nodeiter.graph_pointer && internal_node_iterator == nodeiter.internal_node_iterator);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_pointer;
    typename std::vector<size_type>::const_iterator internal_node_iterator; //iterator over unique_id's of the node_types in the i2u vector of the graph_type

    /** 
    * node_iterator private constructor (mainly used by the graph_type)
    *
    * @param[in] graphpointer   pointer to the graph_type
    * @param[in] intni  iterator over the unique_id's of the graph_type's node_types in the i2u vector
    *
    */
    NodeIterator(const graph_type* graphpointer, typename std::vector<size_type>::const_iterator intni) {
      graph_pointer = const_cast<graph_type*>(graphpointer);
      internal_node_iterator = intni;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // NodeIterator node_begin() const
  // NodeIterator node_end() const

  /**
  * @brief returns a node_iterator corresponding to the first node of this graph
  * 
  * @return node_iterator pointing to the first node of this graph, if graph is empty returns a nullptr
  * @post if graph is non-empty result can be derefenced, otherwise the result cannot be dereferenced
  */
  node_iterator node_begin() const {
    return node_iterator(this, i2u.begin());
  } 

  /**
  * @brief returns a node_iterator corresponding to one past the last node of this graph_type
  * 
  * @return node_iterator pointing to one past the last node of this graph_type
  * @post result cannot be dereferenced
  */
  node_iterator node_end() const {
    return node_iterator(this, i2u.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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
      node1_uid = 0;
      node2_uid_iterator = typename std::vector<std::pair<size_type, size_type>>::const_iterator();
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
    * @brief dereferences the incident_iterator to return the corresponding edge_type
    *
    * @result edge_type that points to the corresponding edge
    * @pre the incident_iterator points to an actual edge_type in a non_empty graph_type
    * @post (*incident_iterator).node1() == graph_pointer->node(node1_id)
    */
    value_type operator*() const {
      assert(graph_pointer != nullptr && node2_uid_iterator != graph_pointer->adjacency_matrix[node1_uid].end());
      return edge_type(graph_pointer, node1_uid, (*node2_uid_iterator).first, (*node2_uid_iterator).second);
    }

    /**
    * @brief increments the incident_iterator to point to the next edge incident from the main node 1
    *
    * @result incident_iterator reference that points to the next edge
    * @pre the incident_iterator points to an actual edge_type in a non_empty graph_type
    * @post (*old incident_iterator).node1() == (*new incident_iterator).node1()
    */
    incident_iterator& operator++() {
      assert(graph_pointer != nullptr && node2_uid_iterator != graph_pointer->adjacency_matrix[node1_uid].end());
      ++node2_uid_iterator;
      typename std::vector<std::pair<size_type, size_type>>::const_iterator end_iterator = graph_pointer->adjacency_matrix[node1_uid].end(); //to convert end() from non-const iterator to const
      node2_uid_iterator = std::find_if(node2_uid_iterator, end_iterator, 
        [&] (std::pair<size_type, size_type> pe) {return edge_type(graph_pointer, node1_uid, pe.first, pe.second).valid();});
      return (*this);
    }

     /**
    * @brief checks whether the edges pointed to by the incident_iterator are equal
    *
    * @param[in] incident_iterator iter to which we compare this incident_iterator
    * @result true if equal false otherwise
    */
    bool operator==(const IncidentIterator& iter) const {
      return (graph_pointer == iter.graph_pointer && node1_uid == iter.node1_uid && node2_uid_iterator == iter.node2_uid_iterator);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_pointer;
    size_type node1_uid;
    typename std::vector<std::pair<size_type, size_type>>::const_iterator node2_uid_iterator; //iterator over the possible node_types directly connected to node1_uid

    /**
    * private constructor of incident_iterator (mainly used by node_type)
    *
    * @param[in] graphpointer   pointer to the graph_type
    * @param[in] node1  unique_id of node1 of the edge_type that is pointed to by this incident_iterator, which is the main node_type
    * @param[in] node2iter  iterator over the possible node_types directly connected to _node1_
    *
    */
    IncidentIterator(const graph_type* graphpointer, size_type node1, typename std::vector<std::pair<size_type, size_type>>::const_iterator node2iter){
      graph_pointer = const_cast<graph_type*>(graphpointer);
      node1_uid = node1;
      typename std::vector<std::pair<size_type, size_type>>::const_iterator end_iterator = graph_pointer->adjacency_matrix[node1].end(); //to convert end() from non-const iterator to const
      node2_uid_iterator = std::find_if(node2iter, end_iterator, 
        [&] (std::pair<size_type, size_type> pe) {return edge_type(graph_pointer, node1, pe.first, pe.second).valid();}) ;
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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
      internal_edge_iterator = typename std::vector<size_type>::const_iterator();
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
    * @brief dereferences the edge_iterator to obtain the corresponding edge_type
    * 
    * @return edge_type behind the edge_iterator
    * @pre this edge_iterator points to an actual edge_type in a non-empty graph_type
    */
    value_type operator*() {
      assert(graph_pointer != nullptr && internal_edge_iterator != graph_pointer->edges_i2u.end());
      return graph_pointer->edge(graph_pointer->edges[*internal_edge_iterator].index_);
    }

    /**
    * @brief increments edge_iterator to point to the next edge_type
    * 
    * @return edge_iterator reference corresponding to the next edge_type
    * @pre this edge_iterator points to an actual edge_type in a non-empty graph_type
    */
    edge_iterator& operator++() {
      assert(graph_pointer != nullptr && internal_edge_iterator != graph_pointer->edges_i2u.end());
      ++internal_edge_iterator;
      return (*this);
    }

    /**
    * @brief checks whether this edge_iterator is equal to another edge_iterator nodeiter,
    * which is equivalent to checking equality of the underlying edge_types
    * 
    * @param[in] eiter   other edge_iterator object to compare with this edge_iterator
    * @return true if the edge_iterator's are equal false otherwise
    */
    bool operator==(const edge_iterator& eiter) const {
      //return (**this) == (*eiter); 
      return (graph_pointer == eiter.graph_pointer && internal_edge_iterator == eiter.internal_edge_iterator);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_pointer;
    typename std::vector<size_type>::const_iterator internal_edge_iterator; //iterator over the edge_type unique_id's in the edges_i2u vector


    /**
    * private constructor of edge_iterator (mainly used by graph_type)
    *
    * @param[in] graphpointer   pointer to the graph_type
    * @param[in] iei  iterator over the edge_type unique_id's in the edges_i2u vector
    *
    */
    EdgeIterator(const graph_type* graphpointer, typename std::vector<size_type>::const_iterator iei) {
      graph_pointer = const_cast<graph_type*>(graphpointer);
      internal_edge_iterator = iei;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
  * @brief returns a edge_iterator corresponding to the first edge of this graph
  * 
  * @return edge_iterator pointing to the first edge of this graph, if graph is empty returns a nullptr
  * @post if graph is non-empty result can be derefenced, otherwise the result cannot be dereferenced
  */
  edge_iterator edge_begin() const {
    return edge_iterator(this, edges_i2u.begin());
  }

  /**
  * @brief returns a edge_iterator corresponding to one past the last edge of this graph_type
  * 
  * @return edge_iterator pointing to one past the last edge of this graph_type
  * @post result cannot be dereferenced
  */
  edge_iterator edge_end() const {
    return edge_iterator(this, edges_i2u.end());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //internal node that actually contains all the node data
  struct internal_node { 
    Point position_;   
    size_type index_;
    node_value_type value_;     
  };

  //internal edge that actually contains all the edge data
  struct internal_edge { 
    size_type node_1_uid_;
    size_type node_2_uid_;
    size_type index_;
    edge_value_type value_;
  };

  std::vector<internal_node> nodes; //vector of internal nodes where all the node_type's info is stored
  std::vector<internal_edge> edges; //vector of internal edges where all the edge_type's info is stored
  std::vector<std::vector<std::pair<size_type, size_type>>> adjacency_matrix;
  std::vector<size_type> i2u; //node_type index to unique_id
  std::vector<size_type> edges_i2u; //edge_type index to unique_id
};

#endif // CME212_GRAPH_HPP
