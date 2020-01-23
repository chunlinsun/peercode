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
class Graph {
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
  class Node {
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
    const Graph* g;
    size_type idx;
    
    /**Valid Node constructor*/
    Node(const Graph* g, size_type idx) {
      this->g = g;
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
  Node add_node(const Point& position) { //shoukld call the valid constyructor
    Nodes.push_back(position);
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
  class Edge {
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
    const Graph* g;
    size_type idx;
    bool backwards;

    /** Valid Edge constructor */
    Edge(const Graph* g, size_type idx, bool backwards) {
      this->g = g;
      this->idx = idx;
      this->backwards = backwards;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
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
    else if (a < b) {
      //Check whether b is in the list of nodes adjacent to a. 
      if (EdgeMap.count(a.idx) == 0)
        return false;
      for (unsigned int i = 0; i < EdgeMap.at(a.idx).size(); i++) {
        if (std::get<0>(EdgeMap.at(a.idx).at(i)) == b.idx)
          return true;
      }
    }
    else {
      //Check whether a is in the list of nodes adjacent to b.
      if (EdgeMap.count(b.idx) == 0)
        return false;
      for (unsigned int i = 0; i < EdgeMap.at(b.idx).size(); i++) {
        if (std::get<0>(EdgeMap.at(b.idx).at(i)) == a.idx)
          return true;
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
      if (a < b) {
        for (unsigned int i = 0; i < EdgeMap.at(a.idx).size(); i++) {
          if (std::get<0>(EdgeMap.at(a.idx).at(i)) == b.idx) {
            return Edge(this, std::get<1>(EdgeMap.at(a.idx).at(i)), false);
          }
        }
      }
      else {
        for (unsigned int i = 0; i < EdgeMap.at(b.idx).size(); i++) {
          if (std::get<0>(EdgeMap.at(b.idx).at(i)) == a.idx) {
            return Edge(this, std::get<1>(EdgeMap.at(b.idx).at(i)), true);
          }
        }
      }
    }
    //Add a new edge to the EdgeMap and EdgeVector.
    if (a < b) {
      std::tuple <size_type, size_type> tupV = std::make_tuple(a.idx, b.idx);
      EdgeVector.push_back(tupV);
      if (EdgeMap.count(a.idx) == 0) {
        std::vector <std::tuple <size_type, size_type>> vec;
        EdgeMap.insert(std::pair<size_type, std::vector<std::tuple<size_type, size_type>>>(a.idx, vec));
      }
      std::tuple <size_type, size_type> tupM = std::make_tuple(b.idx, num_edges() - 1);
      EdgeMap.at(a.idx).push_back(tupM);
    }
    else {
      std::tuple <size_type, size_type> tupV = std::make_tuple(b.idx, a.idx);
      EdgeVector.push_back(tupV);
      if (EdgeMap.count(b.idx) == 0) {
        std::vector <std::tuple <size_type, size_type>> vec;
        EdgeMap.insert(std::pair<size_type, std::vector<std::tuple<size_type, size_type>>>(b.idx, vec));
      }
      std::tuple <size_type, size_type> tupM = std::make_tuple(a.idx, num_edges() - 1);
      EdgeMap.at(b.idx).push_back(tupM);
    }
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

  std::vector<Point> Nodes;
  std::map<size_type, std::vector<std::tuple <size_type, size_type>>> EdgeMap;
  std::vector<std::tuple<size_type, size_type>> EdgeVector;

};

#endif // CME212_GRAPH_HPP
