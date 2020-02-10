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
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node. */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return g->Nodes.at(idx);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx;
    }
    
    /** Set a Node's value.
     * @param[in] val of the type V, the type of values in the graph.
     * @post The Node's value has been updated to val.
     *
     * Complexity: O(1) amortized operations.
     */
    node_value_type& value_set(V val){
      g->Values.at(idx) = val;
      return g->Values.at(idx);
    }

    /** Return the value of a node.
     * Complexity: O(1) amortized operations.
     */
    const node_value_type& value() const {
      return g->Values.at(idx);
    }
    
    /** Return the degree of a node
     */
    size_type degree() const {
      //--functionality_1
      //--Code does not compile because Edgemap can only be accessed via the graph pointer.
      //--This fixes the compilation error. However, the code still segfaults when the degree 
      //--method is called. Most likely because you are not adding the node index to the edgemap
      //--when adding nodes. So nodes with zero degree will cause this to segfault.   
      //--START
      return g->EdgeMap.at(idx).size();
      //--END
    }

    /** Return an IncidentIterator that points to the beginning of the node's
    list of adjacent edges. */
    incident_iterator edge_begin() const {
      return IncidentIterator(g, idx, 0);
    }

    /** Return an IncidentIterator that points to the end of the node's
    list of adjacent edges. */
    incident_iterator edge_end() const {
      return IncidentIterator(g, idx, g->EdgeMap.at(idx).size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (g == n.g && idx == n.idx);
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
      return (idx < n.idx);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /**Node attributes*/
    Graph* g;
    size_type idx;
    
    /**Valid Node constructor*/
    Node(const Graph* g, size_type idx) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->idx = idx;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return Nodes.size();
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
    Nodes.push_back(position);
    Values.push_back(val);
    return Node(this, num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (num_nodes() > n.idx);
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (backwards == false) {
        return Node(g, std::get<0>(g->EdgeVector.at(idx)));
      }
      else {
        return Node(g, std::get<1>(g->EdgeVector.at(idx)));
      }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (backwards == false) {
        return Node(g, std::get<1>(g->EdgeVector.at(idx)));
      }
      else {
        return Node(g, std::get<0>(g->EdgeVector.at(idx)));
      }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (g == e.g && idx == e.idx);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (idx < e.idx);
    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    /** Edge attributes */
    Graph* g;
    size_type idx;
    bool backwards;

    /** Valid Edge constructor */
    Edge(const Graph* g, size_type idx, bool backwards) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->idx = idx;
      this->backwards = backwards;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return EdgeVector.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i, false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (a == b){
      return false; //Cannot have self-loops
    }
    else {
      //Check whether a is in the list of nodes adjacent to b.
      if (EdgeMap.count(b.idx) == 0)
        return false;
      else {
        for (unsigned int i = 0; i < EdgeMap.at(b.idx).size(); i++) {
          if (std::get<0>(EdgeMap.at(b.idx).at(i)) == a.idx)
            return true;
        }
      }
    }
    return false; //This means we did not find it.
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
    assert((a == b) == false); //Cannot create self-loops.
    //If we already have the edge, return that edge.
    if (has_edge(a, b)){
      for (unsigned int i = 0; i < EdgeMap.at(a.idx).size(); i++) {
        if (std::get<0>(EdgeMap.at(a.idx).at(i)) == b.idx) {
          return Edge(this, std::get<1>(EdgeMap.at(a.idx).at(i)), (b < a));
          }
        }
      }
    //Add a new edge to the EdgeVector.
    if (a < b) {
      std::tuple <size_type, size_type> tupV = std::make_tuple(a.idx, b.idx);
      EdgeVector.push_back(tupV);
    }
    else {
      std::tuple <size_type, size_type> tupV = std::make_tuple(b.idx, a.idx);
      EdgeVector.push_back(tupV);
    }

    //Add a new edge to the EdgeMap at index a.
    if (EdgeMap.count(a.idx) == 0) {
      std::vector <std::tuple <size_type, size_type>> vec;
      EdgeMap.insert(std::pair<size_type, std::vector<std::tuple<size_type, size_type>>>(a.idx, vec));
    }
    std::tuple <size_type, size_type> tupM = std::make_tuple(b.idx, num_edges() - 1);
    EdgeMap.at(a.idx).push_back(tupM);
    
    //Add a new edge to the EdgeMap at index b.
    if (EdgeMap.count(b.idx) == 0) {
      std::vector <std::tuple <size_type, size_type>> vec;
      EdgeMap.insert(std::pair<size_type, std::vector<std::tuple<size_type, size_type>>>(b.idx, vec));
    }
    tupM = std::make_tuple(a.idx, num_edges() - 1);
    EdgeMap.at(b.idx).push_back(tupM);

    return Edge(this, num_edges() - 1, (b < a)); //Return the edge we added.
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    std::vector<Point> n;
    Nodes = n;
    std::map<size_type, std::vector<std::tuple <size_type, size_type>>> em;
    EdgeMap = em;
    std::vector<std::tuple<size_type, size_type>> ev;
    EdgeVector = ev;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Return the node that a NodeIterator is pointing to.
     * @pre idx is a valid index in the Nodes list in Graph g.
     * @post Graph g has not been changed.
     *
     * Complexity: O(1) amortized operations.
     */  
    Node operator*() const {
      return Node(g, idx);
    }

    /** Increment the NodeIterator and return it.
     * @pre idx is a valid index in the Nodes list in Graph g.
     *
     * Complexity: O(1) amortized operations.
     */
    NodeIterator& operator++() {
      idx++;
      return *this;
    }

    /* Return true if this NodeIterator is equal to NodeIterator other. 
     */
    bool operator==(const NodeIterator& other) const {
      return (g == other.g && idx == other.idx);
    }

   private:
    friend class Graph;
    Graph* g;
    size_type idx;

    /* Valid NodeIterator Constructor
     */
    NodeIterator(const Graph* g, size_type idx) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->idx = idx;
    }
  };

  /** Return an NodeIterator that points to the first Node.
   */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }


  /** Return an NodeIterator that points to the last Node.
   */
  NodeIterator node_end() const {
    return NodeIterator(this, num_nodes());
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Return the Edge that the IncidentIterator is pointing to.
     */
    Edge operator*() const {
      size_type node2_idx = std::get<0>(g->EdgeMap.at(node_idx).at(edge_num));
      size_type edge_idx = std::get<1>(g->EdgeMap.at(node_idx).at(edge_num));
      return Edge(g, edge_idx, (node2_idx < node_idx));
    }

    /** Increment the IncidentIterator and return it.
     */
    IncidentIterator& operator++() {
      edge_num++;
      return *this;
    }

    /** Return true if this IncidentIterator is equal to IncidentIterator other.
     */
    bool operator==(const IncidentIterator& other) const {
      return (g == other.g && node_idx == other.node_idx && edge_num == other.edge_num);
    }

   private:
    friend class Graph;
    Graph* g;
    size_type node_idx;
    size_type edge_num;

    /* Valid IncidentIterator Constructor
     */
    IncidentIterator(const Graph* g, size_type node_idx, size_type edge_num) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->node_idx = node_idx;
      this->edge_num = edge_num;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }
    
    /** Return the Edge that the EdgeIterator is pointing to.
     */
    Edge operator*() const {
      return Edge(g, idx, false);
    }

    /** Increment the EdgeIterator and return it.
     */
    EdgeIterator& operator++() {
      idx++;
      return *this;
    }

    /** Return true if the EdgeIterator is equal to EdgeIterator other.
     */
    bool operator==(const EdgeIterator& other) const {
      return (g == other.g && idx == other.idx);
    }

   private:
    friend class Graph;
    Graph* g;
    size_type idx;

    /* Valid EdgeIterator Constructor
     */
    EdgeIterator(const Graph* g, size_type idx) {
      Graph* h = const_cast<Graph *>(g);
      this->g = h;
      this->idx = idx;
    }
  };
  
  /* Return an EdgeIterator that points to the first Edge in the graph.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  
  /* Return an EdgeIterator that points to the last Edge in the graph.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:
  //--style_0
  //--Please don't capitalize names for member variables.
  //--START
  std::vector<Point> Nodes;
  std::vector<node_value_type> Values;
  std::map<size_type, std::vector<std::tuple <size_type, size_type>>> EdgeMap;
  std::vector<std::tuple<size_type, size_type>> EdgeVector;
  //--END
};

#endif // CME212_GRAPH_HPP
