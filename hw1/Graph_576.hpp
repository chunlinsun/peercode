#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <tuple>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  std::vector<Point> nodes;
  std::vector<V> vals;
  std::map<unsigned, std::vector<unsigned>> edges_map;
  std::vector<std::tuple<unsigned, unsigned>> edges_vec;


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
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
  // CONSTRUCTORS AND DESTRUCTORnode_value_
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    this->nodes = std::vector<Point>();
    this->vals = std::vector<node_value_type>();
    this->edges_map = std::map<size_type, std::vector<size_type>>();
    this->edges_vec = std::vector<std::tuple<size_type, size_type>>();
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
      this->idx = size_type();
      this->graph_id = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return (this->graph_id)->nodes[this->idx];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** @return This node's value, of type node_value_type. */
    node_value_type& value() {
      return const_cast<Graph*>(this->graph_id)->vals[this->idx];  
    }

    /** @return This node's value, of type node_value_type. */
    const node_value_type& value() const {
      return (this->graph_id)->vals.at(this->idx);
    }

    /** @return This node's degree, a number in the range [0, total_number_of_edges_in_graph]. */
    size_type degree() const {
      if ((this->graph_id)->edges_map.find(this->idx) == (this->graph_id)->edges_map.end()) {
        return 0;
      }      
      else {
        return (this->graph_id)->edges_map.at(this->idx).size();
      }
    }

    /** @return The IncidentIterator pointing at the first edge incident to this node. */
    incident_iterator edge_begin() const {
      if ((this->graph_id)->edges_map.find(this->idx) == (this->graph_id)->edges_map.end()) {
        return IncidentIterator(this->graph_id, nullptr, this->idx, 0);
      }
      else {
        return IncidentIterator(this->graph_id, 
          &(const_cast<Graph*>(this->graph_id)->edges_map.at(this->idx)), this->idx, 0);        
      }
    }
    
    /** @return The IncidentIterator pointing at the last edge incident to this node. */
    incident_iterator edge_end() const {
      if ((this->graph_id)->edges_map.find(this->idx) == (this->graph_id)->edges_map.end()) {
        return IncidentIterator(this->graph_id, nullptr, this->idx, 0);
      }      
      else {
        return IncidentIterator(this->graph_id, 
          &(const_cast<Graph*>(this->graph_id)->edges_map.at(this->idx)), this->idx, 
          ((this->graph_id)->edges_map.at(this->idx)).size());        
      }
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      bool same_graph = (this->graph_id == n.graph_id);
      bool same_index = (this->index() == n.index());
      return same_graph && same_index;
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
      bool same_graph = (this->graph_id == n.graph_id);
      bool lesser_index = (this->index() < n.index());
      return same_graph && lesser_index;      
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Node(const size_type& idx, const Graph* graph_id) {
      this->idx = idx;      
      this->graph_id = graph_id;
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
//--design_1
//--Node proxy should not store your value
//--START
    size_type idx;
    const Graph* graph_id;
    node_value_type val;
//--END
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->nodes.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    Node new_node((this->nodes).size(), this);
    (this->nodes).push_back(position);
    (this->vals).push_back(val);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_id;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= Graph::num_nodes()) {
      std::cerr << "In Node::node: Index is out of bounds." << std::endl;
      return Node();
    }
    else {
      Node node_i(i, this);
      return node_i;
    }

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
      this->src_idx = size_type();
      this->dest_idx = size_type();
      this->graph_id = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return this->graph_id->node(this->src_idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return this->graph_id->node(this->dest_idx);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool case_1 = (this->node1() == e.node1() && this->node2() == e.node2());
      bool case_2 = (this->node2() == e.node1() && this->node1() == e.node2());
      return (case_1 || case_2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      bool equal = (*this == e);
      bool same_graph = (this->graph_id == e.graph_id);
      bool lesser_src_idx = (this->src_idx <= e.src_idx);
      return !equal && same_graph && lesser_src_idx;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Edge(const size_type& a_idx, const size_type& b_idx, const Graph* graph_id) {
      this->src_idx = a_idx;
      this->dest_idx = b_idx;
      this->graph_id = graph_id;
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type src_idx;
    size_type dest_idx;
    const Graph* graph_id;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->edges_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= Graph::num_edges()) {
      std::cerr << "In Edge::edge: Index is out of bounds." << std::endl;
      return Edge();
    }
    else {
      Edge edge_i(std::get<0>((this->edges_vec)[i]), std::get<1>((this->edges_vec)[i]), this);
      return edge_i;
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
    size_type a_idx = a.index();
    size_type b_idx = b.index();     
    if (this->edges_map.find(a_idx) == this->edges_map.end()) {
      return false;      
    }
    else {
      if (std::find((this->edges_map).at(a_idx).begin(), (this->edges_map).at(a_idx).end(), b_idx) == (this->edges_map).at(a_idx).end()) {
        return false;
      }
      else {
        return true;
      }
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
    // HW0: YOUR CODE HERE
    size_type a_idx = a.index();
    size_type b_idx = b.index();    
    if (Graph::has_edge(a, b)) {
      Edge curr_edge(a_idx, b_idx, this);
      return curr_edge;
    }
    else {
      this->edges_vec.push_back(std::make_tuple(a_idx, b_idx));
      if (this->edges_map.find(a_idx) == this->edges_map.end()) {
        std::vector<size_type> vect;
        vect.push_back(b_idx);
        this->edges_map[a_idx] = vect;
      }
      else {
        this->edges_map[a_idx].push_back(b_idx);
      }
      if (this->edges_map.find(b_idx) == this->edges_map.end()) {
        std::vector<size_type> vect;
        vect.push_back(a_idx);
        this->edges_map[b_idx] = vect;
      }
      else {
        this->edges_map[b_idx].push_back(a_idx);
      }
      Edge new_edge(a_idx, b_idx, this);
      return new_edge;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->nodes.clear();
    this->edges_vec.clear();
    this->edges_map.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      this->graph_id = nullptr;
      this->i = size_type();
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereference operator.
     * @return The correponding Node object _n_.
     * @post _n_ == Node(_i_, _graph_id_).
     */
    Node operator*() const {
      return Node(this->i, this->graph_id);
    }

    /** Increment operator.
     * @return The NodeIterator _node_iter_ for which _i_ has been incremented 
     * by 1.
     */
    NodeIterator& operator++() {
      this->i++;
      return *this;
    }

    /** Equality operator.
     * @return True if this NodeIterator has the same index and is in the 
     * same graph as _node_iter_.
     *         False otherwise.
     */
    bool operator==(const NodeIterator& node_iter) const {
      bool same_graph = (this->graph_id == node_iter.graph_id);
      bool same_i = (this->i == node_iter.i);
      return same_graph && same_i;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* graph_id;
    size_type i;

    /** Private constructor accessible inside Graph class. 
    * @param[in] _graph_id_ The pointer of the associated graph object.
    * @param[in] _i_ The index of the corresponding node this iterator is
    *               pointed at.
    * @pre 0 <= i <= num_nodes.
    */
    NodeIterator(const Graph* graph_id, size_type i){
      assert(i <= graph_id->num_nodes());
      this->graph_id = graph_id;
      this->i = i;
    }     
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** @return The NodeIterator pointing at the first node of this graph. */  
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** @return The NodeIterator pointing at the last node of this graph. */  
  node_iterator node_end() const {
    return NodeIterator(this, nodes.size());
  }

  //
  // Incident Iterator
  //

  // alias pointer for convenience
  using edges_map_ptr = std::vector<unsigned>*;

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
      this->graph_id = nullptr;
      this->p = nullptr;
      this->src_idx = size_type();
      this->incident_idx = size_type();      
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereference operator.
     * @return The Edge object _e_ linking the Node at index _src_idx_
     * with the Node at index _incident_idx_.
     * @post _e_ == Edge(_i_, incident_node_index, _graph_id_).
     */    
    Edge operator*() const {
      size_type dest_idx = (*(this->p))[this->incident_idx];
      return Edge(this->src_idx, dest_idx, this->graph_id);
    }

    /** Increment operator.
     * @return A IncidentIterator _incident_iter_ for which _incident_idx_
     * has been incremented by 1.
     */
    IncidentIterator& operator++() {
      this->incident_idx++;
      return *this;
    }

    /** Equality operator.
     * @return True if this NodeIterator has the same source and incident 
     * node indices, the same pointer to the vector of incident edges and 
     * is in the same graph as _node_iter_.
     *         False otherwise.
     */
    bool operator==(const IncidentIterator& incident_iter) const {
      bool same_graph = (this->graph_id == incident_iter.graph_id);
      bool same_p = (this->p == incident_iter.p);
      bool same_src = (this->src_idx == incident_iter.src_idx);
      bool same_incident = (this->incident_idx == incident_iter.incident_idx);
      return same_graph && same_p && same_src && same_incident;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph_id;
    edges_map_ptr p;
    size_type src_idx;
    size_type incident_idx;

    /** Private constructor accessible inside Graph class. 
    * @param[in] _graph_id_ The pointer of the associated graph object.
    * @param[in] _p_ The pointer of the vector of indices of nodes incident to 
    *               the node of index _src_idx_.
    * @param[in] _src_idx_ The index of the node for which we want to iterate over
    *                     its incident nodes.
    * @param[in] _incident_idx_ The index of the node linked to the node of index 
    *                         _src_idx_.
    * @pre 0 <= _src_idx_ < num_nodes.
    * @pre 0 <= incident_idx < number of nodes incident to the node of index 
    *      _src_idx_.
    */    
    IncidentIterator(const Graph* graph_id, edges_map_ptr p, size_type src_idx, size_type incident_idx) {
      assert(src_idx <= graph_id->num_nodes());
      assert(incident_idx <= (*p).size());
      this->graph_id = graph_id;
      this->p = p;
      this->src_idx = src_idx;
      this->incident_idx = incident_idx;
    }
  };

  //
  // Edge Iterator
  //

  // alias pointer for convenience
  using edges_vec_ptr = const std::vector<std::tuple<unsigned, unsigned>>*; 

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
      this->graph_id = nullptr;
      this->p = nullptr;
      this->i = size_type();
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference operator.
     * @return The Edge object _e_ such that Graph::edge(_i_) == _e_.
     */      
    Edge operator*() const {
      size_type src_idx = std::get<0>((*(this->p)).at(this->i));
      size_type dest_idx = std::get<1>((*(this->p)).at(this->i));
      return Edge(src_idx, dest_idx, this->graph_id);
    }

    /** Increment operator.
     * @return A EdgeIterator _edge_iter_ for which the index _i_
     * has been incremented by 1.
     */
    EdgeIterator& operator++() {
      this->i++;
      return *this;
    }

    /** Equality operator.
     * @return True if this EdgeIterator has the same index, 
     * the same pointer to the vector of edges and is in the 
     * same graph as _node_iter_.
     *         False otherwise.
     */
    bool operator==(const EdgeIterator& edge_iter) const {
      bool same_graph = (this->graph_id == edge_iter.graph_id);
      bool same_p = (this->p == edge_iter.p);
      bool same_i = (this->i == edge_iter.i);
      return same_graph && same_p && same_i;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_id;
    edges_vec_ptr p;
    size_type i;

//--documentation_-1
//--good doxygen style docs
//--END

    /** Private constructor accessible inside Graph class. 
    * @param[in] _graph_id_ The pointer of the associated graph object.
    * @param[in] _p_ The pointer of the vector of indices of edges.
    * @param[in] _i_ The index of the corresponding edge.
    * @pre 0 <= _i_ < num_edges.
    */        
    EdgeIterator(const Graph* graph_id, edges_vec_ptr p, size_type i) {
      assert(i <= graph_id->num_edges());
      this->graph_id = graph_id;
      this->p = p;
      this->i = i;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** @return The EdgeIterator pointing at the first edge of this graph. */    
  edge_iterator edge_begin() const {
    if ((this->edges_vec).size() == 0) {
      return EdgeIterator(this, nullptr, 0);
    }    
    return EdgeIterator(this, &(this->edges_vec), 0);
  }

  /** @return The EdgeIterator pointing at the last edge of this graph. */  
  edge_iterator edge_end() const {
    if ((this->edges_vec).size() == 0) {
      return EdgeIterator(this, nullptr, 0);
    }
    return EdgeIterator(this, &(this->edges_vec), (this->edges_vec).size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
