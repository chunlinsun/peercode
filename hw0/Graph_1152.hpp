#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <functional>

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
  
  /** Invalid node index value */
  static size_type constexpr invalid_index = size_type(-1);

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  {
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
    Node()
    : m_graph(nullptr),
      m_index(invalid_index)
    {
    }

    /** Return this node's position. */
    const Point & position() const
    {
      assert(is_valid());
      return m_graph->m_points[m_index];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return m_index;
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
    bool operator==(const Node & n) const
    {
      return m_graph == n.m_graph && m_index == n.m_index;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node & n) const
    {
      return std::less<const Graph *>{}(m_graph, n.m_graph) || 
             (m_graph == n.m_graph && m_index < n.m_index);
    }

   private:
    
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    /** Check node for validity.
     * @return @p true if node is valid (i.e. represents an actual graph node)
     * or @p false otherwise.
     */ 
    bool is_valid() const
    {
      return m_graph != nullptr &&
             m_index != invalid_index;
    }
    
    /** Private constructor used by Graph to create Node objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] index index of the node in owning Graph
     */
    Node(const Graph * const graph, const size_type index)
      : m_graph(graph),
        m_index(index)
    {
      assert(is_valid()); // sanity check, this should never occur
    }
    
    /** Pointer to the containing Graph */
    const Graph * m_graph;
    
    /** Index of the node in the containing Graph */
    size_type m_index;
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return static_cast<size_type>(m_points.size());
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point & position)
  {
    m_points.push_back(position);
    return Node(this, num_nodes()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node & n) const
  {
    return n.m_graph == this && n.index() < num_nodes();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    assert(i < num_nodes());
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
    Edge()
    : m_graph(nullptr),
      m_index(invalid_index),
      m_node1(invalid_index),
      m_node2(invalid_index)
    {
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      assert(is_valid());
      return Node(m_graph, m_node1);
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      assert(is_valid());
      return Node(m_graph, m_node2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge & e) const
    {
      // Graph class maintains the invariant m_node1 <= m_node2.
      // Therefore we don't have to check the other combination.
      return m_graph == e.m_graph &&
             m_node1 == e.m_node1 &&
             m_node2 == e.m_node2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge & e) const
    {
      return std::less<const Graph *>{}(m_graph, e.m_graph) ||
             (m_graph == e.m_graph && (m_node1 < e.m_node1 || 
                                      (m_node1 == e.m_node1 && m_node2 < e.m_node2)));
    }

   private:
    
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    /** Check if the edge is valid.
     * @return @p true if edge is valid (i.e. represents an actual graph edge)
     * or @p false otherwise.
     */ 
    bool is_valid() const
    {
      return m_graph != nullptr && 
             m_index != invalid_index &&
             m_node1 != invalid_index && 
             m_node2 != invalid_index &&
             m_node1 <= m_node2;
    }
    
    /** Private constructor used by Graph to create Edge objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] node1 index of first node in owning Graph
     * @param[in] node2 index of second node in owning Graph
     */
    Edge(const Graph * const graph, 
         const size_type index,
         const size_type node1, 
         const size_type node2)
    : m_graph(graph),
      m_index(index),
      m_node1(node1),
      m_node2(node2)
    {
      assert(is_valid()); // sanity check, this should never occur
    }
    
    /** Pointer to the containing Graph */
    const Graph * m_graph;
    
    // Some of the members below are not mandatory.
    // We could remove m_index, but assuming in the future there
    // will be data in the graph associated with edges, it will
    // allow for faster access. Or we could remove node indices
    // m_node1 and m_node2, but they save us an extra trip to
    // graph storage, and we still fit into the 32-byte limit.
    
    /** Index of the edge */
    size_type m_index;
    
    /** Index of first node */
    size_type m_node1;
    
    /** Index of second node */
    size_type m_node2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return static_cast<size_type>(m_edges.size());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    assert(i < num_edges());
    return Edge(this, i, m_edges[i].first, m_edges[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node & a, const Node & b) const
  {
    return m_edge_lookup.count(get_node_indices(a, b)) > 0;
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
  Edge add_edge(const Node & a, const Node & b)
  {
    const auto nodes = get_node_indices(a, b);
    const size_type index = num_edges();
    
    const bool added = m_edge_lookup.emplace(nodes, index).second;
    if (added)
    {
      m_edges.push_back(nodes);
    }
    
    return Edge(this, index, nodes.first, nodes.second);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    m_edge_lookup.clear();
    m_edges.clear();
    m_points.clear();
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
  
  // Enforce the size constraint
  static_assert(sizeof(Node) <= 16, "Node class exceeds size limit");
  static_assert(sizeof(Edge) <= 32, "Edge class exceeds size limit");
  
  /** Alias for stored edge element */
  using index_pair = std::pair<size_type, size_type>;
              
  /** Get node indices in correct order for an edge
   * @pre @a a and @a b are valid nodes of this graph
   * @param[in] a first node
   * @param[in] b second node
   * @return pair of node indices in correct order
   */
  index_pair get_node_indices(const Node & a, const Node & b) const
  {
    assert(a.is_valid() && a.m_graph == this);
    assert(b.is_valid() && b.m_graph == this);
    return { std::min(a.index(), b.index()),
             std::max(a.index(), b.index()) };
  }
  
  /** Storage for node points */
  std::vector<Point> m_points;
  
  /** Storage for edges, smaller node number first */
  std::vector<index_pair> m_edges;
  
  /** Fast edge lookup by index nodes */
  std::map<index_pair, size_type> m_edge_lookup;

};

#endif // CME212_GRAPH_HPP
