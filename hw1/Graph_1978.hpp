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
template <typename NVT>
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
  
  /** Aliases for node index */
  using node_idx_type = size_type;
  
  /** Aliases for edge index */
  using edge_idx_type = size_type;
  
  /** Invalid node index value */
  static node_idx_type constexpr invalid_node_index = node_idx_type(-1);
  
  /** Invalid node index value */
  static edge_idx_type constexpr invalid_edge_index = edge_idx_type(-1);
  
  /** Alias for value type stored at nodes **/
  using node_value_type = NVT;

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
    Node()
    : m_graph(nullptr),
      m_index(invalid_node_index)
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
    
    /** Accessor for the value associated with the node.
     * @return reference to the value
     * @pre node is valid, i.e. not default-constructed
     **/
    node_value_type & value()
    {
      assert(is_valid());
      return m_graph->m_values[m_index];
    }
    
    /** Accessor for the value associated with the node.
     * @return const reference to the value
     * @pre node is valid, i.e. not default-constructed
     **/
    const node_value_type & value() const
    {
      assert(is_valid());
      return m_graph->m_values[m_index];
    }
    
    /** 
     * Return number of undirected edges incident on the node 
     * @return number of edges starting at this node
     * @pre node is valid, i.e. not default-constructed
     **/
    size_type degree() const
    {
      assert(is_valid());
      return m_graph->m_inc_lookup[m_index].size();
    }
    
    /** Begin iterator over node incidence.
     * @return iterator to the beginning of node's edge list
     * @pre the node is valid, i.e. not default-constructed
     */
    incident_iterator edge_begin() const
    {
      assert(is_valid());
      return incident_iterator(m_graph, m_index, m_graph->m_inc_lookup[m_index].begin());
    }
    
    /** End iterator over node incidence.
     * @return iterator to the end of node's edge list
     * @pre the node is valid, i.e. not default-constructed
     */
    incident_iterator edge_end() const
    {
      assert(is_valid());
      return incident_iterator(m_graph, m_index, m_graph->m_inc_lookup[m_index].end());
    }

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
      return std::less<Graph *>{}(m_graph, n.m_graph) || 
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
             m_index != invalid_node_index;
    }
    
    /** Private constructor used by Graph to create Node objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] index index of the node in owning Graph
     */
    Node(Graph * const graph, const node_idx_type index)
      : m_graph(graph),
        m_index(index)
    {
      assert(is_valid()); // sanity check, this should never occur
    }
    
    /** Pointer to the containing Graph */
    Graph * m_graph;
    
    /** Index of the node in the containing Graph */
    node_idx_type m_index;
    
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
  Node add_node(const Point & position, node_value_type value = node_value_type{})
  {
    m_points.push_back(position);
    m_values.emplace_back(std::move(value)); // could be large, move if possible
    m_inc_lookup.push_back({});
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
    // HACK - UB if called on a truly const object
    return Node(const_cast<Graph *>(this), i);
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
    Edge()
    : m_graph(nullptr),
      m_node1(invalid_node_index),
      m_node2(invalid_node_index)
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
      return m_graph == e.m_graph &&
             (m_node1 == e.m_node1 && m_node2 == e.m_node2 ||
              m_node1 == e.m_node2 && m_node2 == e.m_node1) ;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge & e) const
    {
      auto const nl = std::minmax(m_node1, m_node2);
      auto const nr = std::minmax(e.m_node1, e.m_node2);
      return std::less<Graph *>{}(m_graph, e.m_graph) ||
             (m_graph == e.m_graph && (nl.first < nr.first || 
                                      (nl.first == nr.first && nl.second < nr.second)));
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
             m_node1 != invalid_node_index && 
             m_node2 != invalid_node_index;
    }
    
    /** Private constructor used by Graph to create Edge objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] node1 index of first node in owning Graph
     * @param[in] node2 index of second node in owning Graph
     */
    Edge(Graph * const graph,
         const size_type node1, 
         const size_type node2)
    : m_graph(graph),
      m_node1(node1),
      m_node2(node2)
    {
      assert(is_valid()); // sanity check, this should never occur
    }
    
    /** Pointer to the containing Graph */
    Graph * m_graph;
    
    /** Index of first node */
    node_idx_type m_node1;
    
    /** Index of second node */
    node_idx_type m_node2;
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
    // HACK - UB if called on a truly const object
    return Edge(const_cast<Graph *>(this), m_edges[i].first, m_edges[i].second);
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
    assert(has_node(a) && has_node(b));
  
    const auto nodes = get_node_indices(a, b);
    const edge_idx_type eidx = num_edges();
    
    const bool added = m_edge_lookup.emplace(nodes, eidx).second;
    if (added)
    {
      // This pair of nodes not seen before; add it everywhere
      m_edges.push_back(nodes);
      m_inc_lookup[nodes.first].push_back(eidx);
      m_inc_lookup[nodes.second].push_back(eidx);
    }
    
    return Edge(this, nodes.first, nodes.second);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    m_points.clear();
    m_values.clear();
    m_edges.clear();
    m_edge_lookup.clear();
    m_inc_lookup.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    : m_graph(nullptr),
      m_index(invalid_node_index)
    {
    }
    
    /** Dereference operator. 
     * @return a proxy object representing pointed-to node
     * @pre iterator is valid (not default-constructed or past-the-end)
     **/
    value_type operator*() const
    {
      assert(is_valid());
      return m_graph->node(m_index);
    }
    
    /** Prefix increment operator.
     * @return reference to this iterator
     * @pre iterator is valid (not default-constructed or past-the-end)
     * @post iterator is pointing to the next node in graph's node collection,
     *       or past-the-end if previously pointed-to node was last in graph
     **/
    NodeIterator & operator++()
    {
      assert(is_valid());
      ++m_index;
      return *this;
    }
    
    /** Equality comparison operator. 
     * @return @p true if iterators are identical, i.e. pointing to
     *         the same graph object and same node in that graph,
     *         or both are default-constructed, or @p false otherwise
     **/
    bool operator==(const NodeIterator & other) const
    {
      return m_graph == other.m_graph && m_index == other.m_index;
    }

   private:
   
    friend class Graph;
    
    /** Check whether the iterator is valid (not default-constructed) **/
    bool is_valid() const
    {
      return m_graph != nullptr && m_index != invalid_node_index;
    }
    
    /** Construct a valid NodeIterator. Used privately by Graph. */
    NodeIterator(Graph * const graph, const node_idx_type index)
    : m_graph(graph),
      m_index(index)
    {
      assert(is_valid());
    }
    
    /** Pointer to the graph **/
    Graph * m_graph;
    
    /** Node index in the graph **/
    node_idx_type m_index;
    
  };

  /** Node begin iterator
   * @return an iterator to the beginning of node collection 
   **/
  node_iterator node_begin() const
  {
    // HACK - UB if called on a truly const object
    return node_iterator(const_cast<Graph *>(this), 0);
  }
  
  /** Node end iterator
   * @return an iterator past-the-end of node collection
   **/
  node_iterator node_end() const
  {
    // HACK - UB if called on a truly const object
    return node_iterator(const_cast<Graph *>(this), num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator()
    : m_graph(nullptr),
      m_src_node(invalid_node_index),
      m_iter{}
    {
    }

    /** Dereference operator.
     * @return a proxy object representing an edge in graph
     * @pre iterator is valid (not default-constructed or past-the-end)
     */
    Edge operator*() const
    {
      assert(is_valid());
      const auto & nodes = m_graph->m_edges[*m_iter];
      const node_idx_type dst_node = (m_src_node == nodes.first) ? nodes.second : nodes.first;
      return Edge(m_graph, m_src_node, dst_node);
    }
    
    /** Prefix increment operator.
     * @return reference to this iterator
     * @pre iterator is valid (not default-constructed or past-the-end)
     * @post iterator is pointing to the next incident edge of the node,
     *       or past-the-end if previous edge was last in node's list
     */
    IncidentIterator & operator++()
    {
      assert(is_valid());
      ++m_iter;
      return *this;
    }
    
    /** Equality comparison operator. 
     * @return @p true if iterators are identical, i.e. pointing to
     *         the same graph object and same incident edge of the 
     *         same node in that graph, or both are default-constructed,
     *         or @p false otherwise
     **/
    bool operator==(const IncidentIterator & other) const
    {
      return m_graph    == other.m_graph && 
             m_src_node == other.m_src_node &&
             m_iter     == other.m_iter;
    }

   private:
   
    friend class Graph;
    
    /** Alias for graph internal data structure iterator held here **/
    using edge_list_iterator = std::vector<size_type>::const_iterator;
    
    /** Check whether the iterator is valid (not default-constructed) **/
    bool is_valid() const
    {
      return m_graph != nullptr && m_src_node != invalid_node_index;
    }
    
    /** Private constructor used by  **/
    IncidentIterator(Graph * const graph, 
                     const node_idx_type src_node,
                     const edge_list_iterator iter)
    : m_graph(graph),
      m_src_node(src_node),
      m_iter(iter)
    {
      assert(is_valid());
    }
    
    /** Pointer to the graph object **/
    Graph * m_graph;
    
    /** Index of the source node (whose incidence is being traversed) **/
    node_idx_type m_src_node;
    
    /** Iterator into the list of edges for a node **/
    edge_list_iterator m_iter;
    
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() 
    : m_graph(nullptr),
      m_index(invalid_edge_index)
    {
    }
    
    /** Dereference operator
     * @return a proxy Edge object representing the pointed-to graph edge
     * @pre the iterator is valid (not past-the-end or default-constructed)
     */
    Edge operator*() const
    {
      assert(is_valid());
      return m_graph->edge(m_index);
    }
    
    /** Prefix increment operator
     * @return reference to this iterator
     * @pre the iterator is valid (not past-the-end or default-constructed)
     * @post iterator is pointing to the next edge in the edge collection,
     *       or past-the-end if currently pointed-to edge is the last one
     */
    EdgeIterator & operator++()
    {
      ++m_index;
      return *this;
    }
    
    /** Equality comparison operator
     * @return @p true if iterators are identical, i.e. pointing to
     *         the same graph object and same edge in that graph,
     *         or both are default-constructed, or @p false otherwise
     */
    bool operator==(const EdgeIterator & other) const
    {
      return m_graph == other.m_graph && m_index == other.m_index;
    }

   private:
   
    friend class Graph;
    
    bool is_valid() const
    {
      return m_graph != nullptr && m_index != invalid_edge_index;
    }
    
    EdgeIterator(Graph * graph, edge_idx_type index)
    : m_graph(graph),
      m_index(index)
    {
      assert(is_valid());
    }
    
    /** Pointer to the graph object **/
    Graph * m_graph;
    
    /** Index of node in the graph **/
    edge_idx_type m_index;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Edge begin iterator
   * @return an iterator to the beginning of edge collection 
   **/
  edge_iterator edge_begin() const
  {
    // HACK - UB if called on a truly const object
    return edge_iterator(const_cast<Graph *>(this), 0);
  }
  
  /** Edge end iterator
   * @return an iterator past-the-end of edge collection 
   **/
  edge_iterator edge_end() const
  {
    // HACK - UB if called on a truly const object
    return edge_iterator(const_cast<Graph *>(this), num_edges());
  }

 private:
  
  // Enforce the size constraint
  static_assert(sizeof(Node) <= 16, "Node class exceeds size limit");
  static_assert(sizeof(Edge) <= 32, "Edge class exceeds size limit");
              
  using node_idx_pair = std::pair<node_idx_type, node_idx_type>;
              
  /** Get node indices in correct order for an edge
   * @pre @a a and @a b are valid nodes of this graph
   * @param[in] a first node
   * @param[in] b second node
   * @return pair of node indices in correct order
   */ 
  node_idx_pair get_node_indices(const Node & a, const Node & b) const
  {
    assert(a.is_valid() && a.m_graph == this);
    assert(b.is_valid() && b.m_graph == this);
    return std::minmax(a.index(), b.index());
  }
  
  /** Storage for node points */
  std::vector<Point> m_points;
  
  /** Storage for node values **/
  std::vector<node_value_type> m_values;
  
  /** Storage for edges, smaller node number first */
  std::vector<node_idx_pair> m_edges;
  
  /** Fast edge lookup by two node indices */
  std::map<node_idx_pair, edge_idx_type> m_edge_lookup;
  
  /** Fast incidence lookup by node **/
  std::vector<std::vector<edge_idx_type>> m_inc_lookup;

};

//--functionality_0
//--Good job!
//--END
#endif // CME212_GRAPH_HPP
