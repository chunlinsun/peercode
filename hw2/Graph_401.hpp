#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <unordered_map>
#include <tuple>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** Custom helper functions for debugging
 * Printing maps and vectors.
 */
template<typename T>
void print_vec(std::vector<T> const &input)
{
  for (unsigned i = 0; i < input.size(); i++) {
    std::cout << input.at(i) << ' ';
  }
}

template<typename K, typename V>
void print_map(std::unordered_map<K,V> const &m)
{
    for (auto const& pair: m) {
      std::cout << "{" << pair.first << ": ";
      print_vec<K>(pair.second);
      std::cout << "}\n";
    }
}

template<typename K, typename V>
void print_keys(std::unordered_map<K, V> const &m)
{
  for (auto const& pair: m) {
    std::cout << pair.first << ' ';
  }  
}

template<typename K, typename V, typename K_, typename V_>
void print_map_val(std::unordered_map<K,V> const &m)
{
    for (auto const& pair: m) {
      std::cout << "{" << pair.first << ": ";
      print_keys<K_, V_>(pair.second);
      std::cout << "}\n";
    }
}


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V, typename E>
// template<typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  typedef std::tuple<unsigned, E> tuple_edgeIdx_val;
  typedef std::unordered_map<unsigned, tuple_edgeIdx_val> map_val;  
  std::vector<Point> nodes;
  std::unordered_map<unsigned, V> node_vals;
  std::unordered_map<unsigned, map_val> edges_map_val;
  std::vector<std::tuple<unsigned, unsigned>> edges_vec;


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  typedef V node_value_type;
  typedef E edge_value_type;

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
    this->node_vals = std::unordered_map<unsigned, node_value_type>();
    this->edges_map_val = std::unordered_map<unsigned, map_val>() ;
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
    Point& position() {
      return const_cast<Graph*>(this->graph_id)->nodes[this->idx];
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
      return const_cast<Graph*>(this->graph_id)->node_vals[this->idx];  
    }

    /** @return This node's value, of type node_value_type. */
    const node_value_type& value() const {
      return (this->graph_id)->node_vals.at(this->idx);
    }

    /** @return This node's degree, a number in the range [0, total_number_of_edges_in_graph]. */
    size_type degree() const {
      if ((this->graph_id)->edges_map_val.count(this->idx) == 0) {
        return 0;
      }      
      else {
        return (this->graph_id)->edges_map_val.at(this->idx).size();
      }
    }

    /** @return The IncidentIterator pointing at the first edge incident to this node. */
    incident_iterator edge_begin() const {
      if ((this->graph_id)->edges_map_val.count(this->idx) == 0) {
        map_val empty_map_val;
        return IncidentIterator(this->graph_id, empty_map_val.end(), this->idx);
      }
      else {
        return IncidentIterator(this->graph_id, 
          (const_cast<Graph*>(this->graph_id)->edges_map_val.at(this->idx)).begin(), this->idx);        
      }      
    }
    
    /** @return The IncidentIterator pointing at the last edge incident to this node. */
    incident_iterator edge_end() const {
      if ((this->graph_id)->edges_map_val.count(this->idx) == 0) {
        map_val empty_map_val;
        return IncidentIterator(this->graph_id, empty_map_val.end(), this->idx);
      }
      else {
        return IncidentIterator(this->graph_id, 
          (const_cast<Graph*>(this->graph_id)->edges_map_val.at(this->idx)).end(), this->idx);        
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
      bool lesser_graph_ptr = (this->graph_id < n.graph_id);
      bool same_graph_cond = (same_graph && lesser_index);
      bool diff_graph_cond = (!same_graph && lesser_graph_ptr);
      return (same_graph_cond || diff_graph_cond);      
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
    size_type idx;
    const Graph* graph_id;
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
    (this->node_vals)[(this->nodes).size()-1] = val;
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
    if (i >= num_nodes()) {
      std::cerr << "In Node::node: Index is out of bounds." << std::endl;
      return Node();
    }
    else {
      Node node_i(i, this);
      return node_i;
    }
  }

  /** Remove a node from the graph, along with all the edges incident to it,
   * and return the number of removed nodes.
   * @param[in] @a n The node to be removed.
   * @pre @a n need not be a node in the graph.
   * @return 1 if the node exists in the graph, 0 otherwise.
   * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
   *       Else,                  new num_nodes() == old num_nodes().      
   * @post has_node(@a n) == false
   *
   * Can invalidate node and edge indices -- in other words, old node(@a i) and 
   * old edge(@a j) might not equal new node(@a i) and new edge(@a j), resp. 
   *
   * Complexity: O(num_edges() + num_nodes()).
   */  
  size_type remove_node(const Node& n){
    if (has_node(n)) {
      size_type i = n.index();
      auto it = n.edge_begin();
      while (it != n.edge_end()) {
        remove_edge(*it);
        it = n.edge_begin();
      }    
      std::swap((this->nodes)[i], (this->nodes).back());
      size_type last_idx = (this->nodes).size() - 1;
      (this->nodes).pop_back();
      (this->node_vals)[i] = (this->node_vals)[last_idx];
      (this->node_vals).erase(last_idx);
      if ((this->edges_map_val).count(last_idx) != 0) {
        (this->edges_map_val)[i] = (this->edges_map_val)[last_idx];
        (this->edges_map_val).erase(last_idx);    
        for (auto it = (this->edges_map_val)[i].begin(); it != (this->edges_map_val)[i].end(); it++)  {
          size_type k = it->first;
          (this->edges_map_val)[k][i] = (this->edges_map_val)[k][last_idx];
          (this->edges_map_val)[k].erase(last_idx);
          size_type edge_idx = std::get<0>((this->edges_map_val)[k][i]);
          (this->edges_vec)[edge_idx] = std::make_tuple(k, i);
        }       
      }
      return 1;
    }
    return 0;
  }

  /** Wrapper function around Graph::remove_node(const Node&).
   * @param[in] @a n_it The NodeIterator corresponding to the node to be removed.
   * @pre @a *n_it need not be a node in the graph.
   * @return If the node exists in the graph,
   *            the NodeIterator associated with the node that was removed,
   *         Else,
   *            the next NodeIterator.
   * @post If old has_node(@a *n_it), new num_nodes() == old num_nodes() - 1.
   *       Else,                      new num_nodes() == old num_nodes().      
   * @post has_node(@a *n_it) == false
   *
   * Can invalidate node and edge indices -- in other words, old node(@a i) and 
   * old edge(@a j) might not equal new node(@a i) and new edge(@a j), resp. 
   *
   * Complexity: O(num_edges() + num_nodes()).
   */  
  node_iterator remove_node(node_iterator n_it){
    size_type num_removed_nodes = remove_node(*n_it);
    if (num_removed_nodes) {
      return NodeIterator(this, n_it.i);      
    }
    return NodeIterator(this, n_it.i + 1);
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

    /** @return The length of the edge, which is computed as the Euclidean distance 
     * between the Point of the Node of index _src_idx_, and the Point of the Node
     * of index _dest_idx_.
     */
    double length() const {
      Point point_1 = node1().position();
      Point point_2 = node2().position();
      return std::sqrt(std::pow((point_1.x - point_2.x), 2) + 
        std::pow((point_1.y - point_2.y), 2) + std::pow((point_1.z - point_2.z), 2));
    }

    /** @return This edge's value, of type edge_value_type. 
    * Note: updating the value of an edge(src_idx, dest_idx) supposes doing so for
    * edge(dest_idx, src_idx) as well.
    */
    edge_value_type& value() {
      return std::get<1>(const_cast<Graph*>(this->graph_id)->edges_map_val[this->src_idx][this->dest_idx]);  
    }

    /** @return This edge's value, of type edge_value_type. */
    const edge_value_type& value() const {
      return std::get<1>((this->graph_id)->edges_map_val.at(this->idx).at(this->dest_idx));
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
      bool lesser_graph_ptr = (this->graph_id <= e.graph_id);
      bool same_graph_cond = (!equal && same_graph && lesser_src_idx);
      bool diff_graph_cond = (!same_graph && lesser_graph_ptr);
      return (same_graph_cond || diff_graph_cond);
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
    if (i >= num_edges()) {
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
    if (this->edges_map_val.count(a_idx) == 0) {
      return false;      
    }
    else {
      if (this->edges_map_val.at(a_idx).count(b_idx) == 0) {
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
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // HW0: YOUR CODE HERE
    size_type a_idx = a.index();
    size_type b_idx = b.index();    
    if (has_edge(a, b)) {
      Edge curr_edge(a_idx, b_idx, this);
      return curr_edge;
    }
    else {
      this->edges_vec.push_back(std::make_tuple(a_idx, b_idx));
      size_type edge_idx = (this->edges_vec).size() - 1;
      if (this->edges_map_val.count(a_idx) == 0) {
        tuple_edgeIdx_val edgeIdx_val = std::make_tuple(edge_idx, val);
        map_val map;
        map[b_idx] = edgeIdx_val;
        this->edges_map_val[a_idx] = map;
      }
      else {
        this->edges_map_val[a_idx][b_idx] = std::make_tuple(edge_idx, val);
      }
      if (this->edges_map_val.count(b_idx) == 0) {
        tuple_edgeIdx_val edgeIdx_val = std::make_tuple(edge_idx, val);
        map_val map;
        map[a_idx] = edgeIdx_val;
        this->edges_map_val[b_idx] = map;              
      }
      else {
        this->edges_map_val[b_idx][a_idx] = std::make_tuple(edge_idx, val);
      }
      Edge new_edge(a_idx, b_idx, this);
      return new_edge;
    }
  }

  /** Remove an edge from the graph, and return the number of removed edges.
   * @param[in] @a n_1 One of the nodes in the edge to be removed.
   * @param[in] @a n_2 One of the nodes in the edge to be removed.
   * @pre @a n_1 and @a n_2 need not be distinct nodes, nor nodes from the 
   * same graph. 
   * @return 1 if the edge exists in the graph, 0 otherwise.
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   *
   * Complexity: O(num_edges()).
   */
  size_type remove_edge(const Node& n_1, const Node& n_2) {
    if (has_edge(n_1, n_2)) {
      size_type i = n_1.index();
      size_type j = n_2.index();
      size_type edge_idx = std::get<0>(this->edges_map_val[i][j]);
      size_type i_new = std::get<0>((this->edges_vec).back());
      size_type j_new = std::get<1>((this->edges_vec).back());
      std::swap((this->edges_vec)[edge_idx], (this->edges_vec).back());
      (this->edges_vec).pop_back();
      std::get<0>((this->edges_map_val)[i_new][j_new]) = edge_idx;
      std::get<0>((this->edges_map_val)[j_new][i_new]) = edge_idx;
      ((this->edges_map_val)[i]).erase(j);
      ((this->edges_map_val)[j]).erase(i);
      if ((this->edges_map_val)[i].empty()) {
        (this->edges_map_val).erase(i);
      }
      if ((this->edges_map_val)[j].empty()) {
        (this->edges_map_val).erase(j);
      }            
      return 1;      
    }
    return 0;
  }

  /** Wrapper method around Graph::remove_edge(const Node&, const Node&).
   * @param[in] @a e Edge to be removed.  
   * @pre @a e need not exist in the graph.
   * @return 1 if the edge exists in the graph, 0 otherwise.
   * @post has_edge(@a e.node1(), @a e.node1()) == false
   * @post If old has_edge(@a e.node1(), @a e.node1()), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   *
   * Complexity: O(num_edges()).
   */
  size_type remove_edge(const Edge& e) {
    Node n_1 = e.node1();
    Node n_2 = e.node2();
    return remove_edge(n_1, n_2);
  }

  /** Wrapper method around Graph::remove_edge(const Edge&).
   * @param[in] @a e_it EdgeIterator corresponding to the edge to be removed.
   * @pre @a *e_it need not exist in the graph.
   * @return If the edge exists in the graph,
   *            the EdgeIterator associated with the edge that was removed,
   *         Else,
   *            the next EdgeIterator.
   * @post has_edge(@a (*e_it).node1(), @a (*e_it).node2()) == false
   * @post If old has_edge(@a (*e_it).node1(), @a (*e_it).node2()), 
   *            new num_edges() == old num_edges() - 1.
   *       Else,      
   *            new num_edges() == old num_edges().
   *
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   *
   * Complexity: O(num_edges()).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type num_removed_edges = remove_edge(*e_it);
    if (num_removed_edges) {
      return EdgeIterator(this, &(this->edges_vec), e_it.i);  
    }
    return EdgeIterator(this, &(this->edges_vec), e_it.i + 1);  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->nodes.clear();
    this->node_vals.clear();
    this->edges_vec.clear();
    this->edges_map_val.clear();
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
      map_val empty_map_val;
      this->graph_id = nullptr;
      this->map_val_it = empty_map_val.end();
      this->src_idx = size_type();
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereference operator.
     * @return The Edge object _e_ linking the Node at index _src_idx_
     * with the Node at index the first element of the iterator of
     * the map _map_val_it_, whose keys are the indices of the indicent 
     * nodes to the one whose index is _src_idx_.
     * @post _e_ == Edge(_i_, incident_node_index, _graph_id_).
     */    
    Edge operator*() const {
      size_type dest_idx = (this->map_val_it)->first;
      return Edge(this->src_idx, dest_idx, this->graph_id);
    }

    /** Increment operator.
     * @return A IncidentIterator _incident_iter_ for which _map_val_it_
     * has been incremented by 1.
     */
    IncidentIterator& operator++() {
      this->map_val_it++;
      return *this;
    }

    /** Equality operator.
     * @return True if this NodeIterator has the same source,
     * the same iterator to the map whose keys are the indices of
     * the incident nodes and is in the same graph as _node_iter_.
     *         False otherwise.
     */
    bool operator==(const IncidentIterator& incident_iter) const {
      bool same_graph = (this->graph_id == incident_iter.graph_id);
      bool same_it = (this->map_val_it == incident_iter.map_val_it);
      bool same_src = (this->src_idx == incident_iter.src_idx);
      return same_graph && same_it && same_src;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph_id;
    typename map_val::iterator map_val_it;
    size_type src_idx;

    /** Private constructor accessible inside Graph class. 
    * @param[in] _graph_id_ The pointer of the associated graph object.
    * @param[in] _map_val_it_ The iterator of the map whose keys are indices of nodes incident to 
    *               the node of index _src_idx_.
    * @param[in] _src_idx_ The index of the node for which we want to iterate over
    *                     its incident nodes.
    * @pre 0 <= _src_idx_ < num_nodes.
    */    
    IncidentIterator(const Graph* graph_id, typename map_val::iterator map_val_it, size_type src_idx) {
      assert(src_idx <= graph_id->num_nodes());
      this->graph_id = graph_id;
      this->map_val_it = map_val_it;
      this->src_idx = src_idx;
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
