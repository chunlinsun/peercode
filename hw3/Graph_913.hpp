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
template <typename NVT, typename EVT = char>
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
  
  /** Aliases for node index */
  using node_idx_type = size_type;
  
  /** Aliases for edge index */
  using edge_idx_type = size_type;
  
  /** Alias for value type stored at nodes **/
  using node_value_type = NVT;
  
  /** Alias for value type stored at edges **/
  using edge_value_type = EVT;

 private:
 
  /** Alias for node uid type **/
  using node_uid_type = size_type;
  
  /** Alias for edge uid type **/
  using edge_uid_type = size_type;
 
  /** Invalid node uid value */
  static node_uid_type constexpr invalid_node_uid = node_uid_type(-1);
  
  /** Invalid node uid value */
  static edge_uid_type constexpr invalid_edge_uid = edge_uid_type(-1);

 public:

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
      m_uid(invalid_node_uid)
    {
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    node_idx_type index() const
    {
      assert(m_uid < m_graph->m_node_indices.size());
      return m_graph->m_node_indices[m_uid];
    }

    /** Return const reference to this node's position. */
    const Point & position() const
    {
      assert(is_valid());
      return m_graph->m_node_positions[index()];
    }
    
    /** Return reference to this node's position. */
    Point & position()
    {
      assert(is_valid());
      return m_graph->m_node_positions[index()];
    }
    
    /** Accessor for the value associated with the node.
     * @return reference to the value
     * @pre node is valid, i.e. not default-constructed
     **/
    node_value_type & value()
    {
      assert(is_valid());
      return m_graph->m_node_values[index()];
    }
    
    /** Accessor for the value associated with the node.
     * @return const reference to the value
     * @pre node is valid, i.e. not default-constructed
     **/
    const node_value_type & value() const
    {
      assert(is_valid());
      return m_graph->m_node_values[index()];
    }
    
    /** 
     * Return number of undirected edges incident on the node 
     * @return number of edges starting at this node
     * @pre node is valid, i.e. not default-constructed
     **/
    size_type degree() const
    {
      assert(is_valid());
      return m_graph->m_incidence[index()].size();
    }
    
    /** Begin iterator over node incidence.
     * @return iterator to the beginning of node's edge list
     * @pre the node is valid, i.e. not default-constructed
     */
    incident_iterator edge_begin() const
    {
      assert(is_valid());
      return incident_iterator(m_graph, m_uid, m_graph->m_incidence[index()].begin());
    }
    
    /** End iterator over node incidence.
     * @return iterator to the end of node's edge list
     * @pre the node is valid, i.e. not default-constructed
     */
    incident_iterator edge_end() const
    {
      assert(is_valid());
      return incident_iterator(m_graph, m_uid, m_graph->m_incidence[index()].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node & n) const
    {
      return m_graph == n.m_graph && m_uid == n.m_uid;
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
             (m_graph == n.m_graph && m_uid < n.m_uid);
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
             m_uid != invalid_node_uid &&
             m_uid < m_graph->m_node_indices.size() &&
             index() < m_graph->num_nodes() &&
             m_graph->m_node_uids[index()] == m_uid;
    }
    
    /** Convenience accessor for uid **/
    node_uid_type uid() const
    {
      return m_uid;
    }
    
    /** Private constructor used by Graph to create Node objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] index index of the node in owning Graph
     */
    Node(Graph * const graph, const node_uid_type uid)
      : m_graph(graph),
        m_uid(uid)
    {
      assert(is_valid()); // sanity check
    }
    
    /** Pointer to the containing Graph **/
    Graph * m_graph;
    
    /** Unique id of node in the containing Graph **/
    node_uid_type m_uid;
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return num_nodes();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return m_node_uids.size();
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
    const node_uid_type uid = m_node_indices.size();
    const node_idx_type idx = m_node_uids.size();
    
    m_node_indices.push_back(idx);
    m_node_uids.push_back(uid);
    m_node_positions.push_back(position);
    m_node_values.emplace_back(std::move(value)); // could be large, move if possible
    m_incidence.push_back({});
    
    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node & n) const
  {
    return n.m_graph == this && n.is_valid();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(node_idx_type i) const
  {
    assert(i < num_nodes());
    // HACK - UB if called on a truly const object
    return Node(const_cast<Graph *>(this), m_node_uids[i]);
  }
  
  /** Remove a node from the graph and all incident edges.
   *
   * @param[in] n the Node to be removed
   * @return the number of nodes removed (always 1)
   * 
   * @pre has_node(n)
   * @post !has_node(n)
   * @post num_nodes() == old num_nodes() - 1
   * @post !has_edge(n, m) for all Node m s.t. old has_edge(n, m),
   * @post num_edges() == old num_edges() - n.degree()
   *
   * Complexity: O(d), where d is degree of node @a n. 
   * O(1) if nodes are assumed to have a fixed upper limit on degree.
   *
   * Invalidates any IncidentIterator created from node @a n or
   * from any node @a m previously adjacent to @a n.
   *
   * Invalidates any NodeIterator pointing to one of:
   * - node being removed
   * - last node in range
   * - past-the-end of the node range
   * 
   * Invalidates any EdgeIterator pointing to one of:
   * - any edge adjacent to node being removed
   * - last n.degree() edges in range
   * - past-the-end of the edge range
   *
   * Invalidates any outstanding Edge objects representing edges adjacent
   * to node @a n, and any outstanding Node objects representing @a n.
   * Does not invalidate any other Node or Edge objects.
   * However, may change node/edge indices and order of traversal.
   */
  size_type remove_node(const Node & n)
  {
    assert(has_node(n));
    
    // remove edges incident on the node
    for (auto it = n.edge_begin(); it != n.edge_end();)
    {
      remove_edge(*it);
    }
    
    const node_idx_type nidx = n.index();
    
    // swap/pop removed node data
    using std::swap; // ADL lookup in case value type has a custom swap function
    swap(m_node_values[nidx], m_node_values.back()); m_node_values.pop_back();
    std::swap(m_node_positions[nidx], m_node_positions.back()); m_node_positions.pop_back();
    std::swap(m_node_uids[nidx], m_node_uids.back()); m_node_uids.pop_back();
    std::swap(m_incidence[nidx], m_incidence.back()); m_incidence.pop_back();
    m_node_indices[m_node_uids[nidx]] = nidx;
    
    return 1;
  }
  
  /** Remove a node from the graph and all incident edges.
   *
   * @param[in] n_it iterator to the Node to be removed
   * @return iterator pointing to the next node after removed one,
   *         or past-the-end if removed node was last (by index)
   *
   * @pre n_it is a valid node iterator of this graph
   * @pre same as remove_node(*n_it)
   * @post same as remove_node(*n_it)
   */
  node_iterator remove_node(node_iterator n_it)
  {
    remove_node(*n_it);
    return n_it; // pointing to next node now
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
      m_uid(invalid_edge_uid),
      m_node1_uid(invalid_node_uid),
      m_node2_uid(invalid_node_uid)
    {
    }
    
    /** Return index of this edge **/
    edge_idx_type index() const
    {
      assert(m_uid < m_graph->m_edge_indices.size());
      return m_graph->m_edge_indices[m_uid];
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      assert(is_valid());
      return Node(m_graph, m_node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      assert(is_valid());
      return Node(m_graph, m_node2_uid);
    }
    
    /** Returns length of a valid edge **/
    double length() const
    {
      assert(is_valid());
      return norm_2(node1().position() - node2().position());
    }
    
    /** Return const reference to edge's value **/
    const edge_value_type & value() const
    {
      assert(is_valid());
      return m_graph->m_edge_values[index()];
    }
    
    /** Return reference to edge's value **/
    edge_value_type & value()
    {
      assert(is_valid());
      return m_graph->m_edge_values[index()];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge & e) const
    {
      return m_graph == e.m_graph && m_uid == e.m_uid;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge & e) const
    {
      return std::less<Graph *>{}(m_graph, e.m_graph) ||
             (m_graph == e.m_graph && (m_uid < e.m_uid));
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
             m_uid != invalid_edge_uid &&
             m_uid < m_graph->m_edge_indices.size() &&
             index() < m_graph->num_edges() &&
             m_graph->m_edge_uids[index()] == m_uid &&
             m_node1_uid != invalid_node_uid && 
             m_node2_uid != invalid_node_uid;
    }
    
    /** Private constructor used by Graph to create Edge objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] uid uid of the edge
     * @param[in] node1_uid uid of first node
     * @param[in] node2_uid uid of second node
     */
    Edge(Graph * const graph,
         const edge_uid_type uid,
         const node_uid_type node1_uid, 
         const node_uid_type node2_uid)
    : m_graph(graph),
      m_uid(uid),
      m_node1_uid(node1_uid),
      m_node2_uid(node2_uid)
    {
      assert(is_valid()); // sanity check, this should never occur
    }
    
    /** Private constructor used by Graph to create Edge objects.
     * @param[in] graph pointer to the owning Graph
     */
    Edge(Graph * const graph,
         const edge_uid_type uid)
    : Edge(graph, uid, 
           graph->m_edge_nodes[graph->m_edge_indices[uid]].first,
           graph->m_edge_nodes[graph->m_edge_indices[uid]].second)
    {}
    
    /** Pointer to the containing Graph */
    Graph * m_graph;
    
    /** Index of edge **/
    edge_uid_type m_uid;
    
    /** Index of first node */
    node_uid_type m_node1_uid;
    
    /** Index of second node */
    node_uid_type m_node2_uid;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return m_edge_uids.size();
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
    return Edge(const_cast<Graph *>(this), m_edge_uids[i], 
                m_edge_nodes[i].first, m_edge_nodes[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node & a, const Node & b) const
  {
    return edge_lookup(a.uid(), b.uid()).first;
  }
  
  /** Test whether a given edge exists in the graph.
   * This is pretty much the same as the edge is being valid.
   * @param[in] e the edge to test
   */
  bool has_edge(const Edge & e) const
  {
    return e.m_graph == this && e.is_valid();
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
  Edge add_edge(const Node & a, const Node & b, edge_value_type value = edge_value_type{})
  {
    assert(has_node(a) && has_node(b));
    
    const auto lookup = edge_lookup(a.uid(), b.uid());
    const edge_uid_type euid = lookup.first ? *lookup.second : m_edge_indices.size();
    
    if (!lookup.first)
    {
      m_edge_indices.push_back(m_edge_uids.size());
      m_edge_uids.push_back(euid);
      m_edge_nodes.emplace_back(a.uid(), b.uid());
      m_edge_values.push_back(value);
      m_incidence[a.uid()].push_back(euid);
      m_incidence[b.uid()].push_back(euid);
    }
    
    return Edge(this, euid, a.uid(), b.uid());
  }
  
  /** Remove the edge connecting two nodes, if one exists.
   *
   * @param[in] a first node
   * @param[in] b second node
   * @return number of edges deleted (0 or 1)
   *
   * @pre has_node(a) && has_node(b)
   * @post !has_edge(a, b)
   * @post num_edges() == old num_edges() - 1, if old has_edge(a,b)
   * @post num_edges() == old num_edges(), if old !has_edge(a,b)
   *
   * Since edges are unique (no two edges can connect the same pair of nodes), 
   * at most one edge will be deleted.
   *
   * Complexity: O(d), where d = max(a.degree(), b.degree()).
   * O(1) if nodes are assumed to have a fixed upper limit on degree.
   *
   * If old !has_edge(a, b), the graph remains unchanged and
   * no objects are invalidated. Othwerwise:
   * 
   * Invalidates IncidentIterator objects created from node @a a of @a b.
   *
   * Invalidates existing EdgeIterator that points to one of:
   * - edge being removed
   * - last edge in range
   * - past-the-end of edge range
   *
   * Invalidates any outstanding Edge objects representing removed edge.
   * Does not invalidate any other Node or Edge objects.
   * However, may change edge indices and order of traversal.
   */
  size_type remove_edge(const Node & a, const Node & b)
  {
    assert(a.m_graph == this && a.is_valid());
    assert(b.m_graph == this && b.is_valid());
    
    if (!has_edge(a, b))
    {
      return 0;
    }
    
    // check if edge exists
    const auto it1 = edge_lookup(a.uid(), b.uid()).second;
    const auto it2 = edge_lookup(b.uid(), a.uid()).second;
    assert(*it1 == *it2); // sanity check
    
    const edge_uid_type euid = *it1;
    const edge_idx_type eidx = m_edge_indices[euid];
    
    // disconnect from both nodes
    std::swap(*it1, m_incidence[a.index()].back()); m_incidence[a.index()].pop_back();
    std::swap(*it2, m_incidence[b.index()].back()); m_incidence[b.index()].pop_back();
    //m_incidence[a.index()].erase(it1);
    //m_incidence[b.index()].erase(it2);
    
    // swap/pop removed edge data
    using std::swap; // ADL lookup in case value type has a custom swap function
    swap(m_edge_values[eidx], m_edge_values.back()); m_edge_values.pop_back();
    std::swap(m_edge_nodes[eidx], m_edge_nodes.back()); m_edge_nodes.pop_back();
    std::swap(m_edge_uids[eidx], m_edge_uids.back()); m_edge_uids.pop_back();
    m_edge_indices[m_edge_uids[eidx]] = eidx;
    
    return 1;
  }
  
  /** Remove an edge from the graph.
   * 
   * @param[in] e the edge to remove
   * @return the number of edges removed (always 1)
   * 
   * @pre @a e is a valid edge of this graph: has_edge(e)
   * @post !has_edge(e)
   * @post num_edges() = old num_edges() - 1
   *
   * Complexity and behavior same as remove_edge(e.node1(), e.node2()).
   */
  size_type remove_edge(const Edge & e)
  {
    assert(has_edge(e));
    return remove_edge(e.node1(), e.node2());
  }
  
  /** Remove an edge from the graph.
   * 
   * @param[in] e_it iterator to the edge to remove
   * @return iterator pointing to the next edge after removed one,
   *         or past-the-end if removed edge was last (by index)
   *
   * @pre @a e_it is a valid edge iterator for this graph.
   * @post same as remove_edge(*e_it)
   * 
   * Complexity and behavior same as remove_edge(*e_it).
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    remove_edge(*e_it);
    return e_it; // pointing to next edge now
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    m_node_indices.clear();
    m_node_uids.clear();
    m_node_positions.clear();
    m_node_values.clear();
    m_incidence.clear();
    m_edge_indices.clear();
    m_edge_uids.clear();
    m_edge_nodes.clear();
    m_edge_values.clear();
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
      m_iter{}
    {
    }
    
    /** Dereference operator. 
     * @return a proxy object representing pointed-to node
     * @pre iterator is valid (not default-constructed or past-the-end)
     **/
    value_type operator*() const
    {
      assert(is_valid());
      return Node(m_graph, *m_iter);
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
      ++m_iter;
      return *this;
    }
    
    /** Equality comparison operator. 
     * @return @p true if iterators are identical, i.e. pointing to
     *         the same graph object and same node in that graph,
     *         or both are default-constructed, or @p false otherwise
     **/
    bool operator==(const NodeIterator & other) const
    {
      return m_graph == other.m_graph && m_iter == other.m_iter;
    }

   private:
   
    using iterator = std::vector<node_uid_type>::const_iterator;
   
    friend class Graph;
    
    /** Check whether the iterator is valid (not default-constructed) **/
    bool is_valid() const
    {
      return m_graph != nullptr;
    }
    
    /** Construct a valid NodeIterator. Used privately by Graph. */
    NodeIterator(Graph * const graph, const iterator iter)
    : m_graph(graph),
      m_iter(iter)
    {
      assert(is_valid());
    }
    
    /** Pointer to the graph **/
    Graph * m_graph;
    
    /** Iterator through Node uid vector **/
    iterator m_iter;
    
  };

  /** Node begin iterator
   * @return an iterator to the beginning of node collection 
   **/
  node_iterator node_begin() const
  {
    // HACK - UB if called on a truly const object
    return node_iterator(const_cast<Graph *>(this), m_node_uids.begin());
  }
  
  /** Node end iterator
   * @return an iterator past-the-end of node collection
   **/
  node_iterator node_end() const
  {
    // HACK - UB if called on a truly const object
    return node_iterator(const_cast<Graph *>(this), m_node_uids.end());
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
      m_src_node(invalid_node_uid),
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
      const auto & nodes = m_graph->m_edge_nodes[m_graph->m_edge_indices[*m_iter]];
      const node_uid_type dst_node = (m_src_node == nodes.first) ? nodes.second : nodes.first;
      return Edge(m_graph, *m_iter, m_src_node, dst_node);
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
    using iterator = std::vector<edge_uid_type>::const_iterator;
    
    /** Check whether the iterator is valid (not default-constructed) **/
    bool is_valid() const
    {
      return m_graph != nullptr && m_src_node != invalid_node_uid;
    }
    
    /** Private constructor used by  **/
    IncidentIterator(Graph * const graph, 
                     const node_uid_type src_node,
                     const iterator iter)
    : m_graph(graph),
      m_src_node(src_node),
      m_iter(iter)
    {
      assert(is_valid());
    }
    
    /** Pointer to the graph object **/
    Graph * m_graph;
    
    /** Index of the source node (whose incidence is being traversed) **/
    node_uid_type m_src_node;
    
    /** Iterator into the list of edges for a node **/
    iterator m_iter;
    
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
      m_iter{}
    {
    }
    
    /** Dereference operator
     * @return a proxy Edge object representing the pointed-to graph edge
     * @pre the iterator is valid (not past-the-end or default-constructed)
     */
    Edge operator*() const
    {
      assert(is_valid());
      return Edge(m_graph, *m_iter);
    }
    
    /** Prefix increment operator
     * @return reference to this iterator
     * @pre the iterator is valid (not past-the-end or default-constructed)
     * @post iterator is pointing to the next edge in the edge collection,
     *       or past-the-end if currently pointed-to edge is the last one
     */
    EdgeIterator & operator++()
    {
      ++m_iter;
      return *this;
    }
    
    /** Equality comparison operator
     * @return @p true if iterators are identical, i.e. pointing to
     *         the same graph object and same edge in that graph,
     *         or both are default-constructed, or @p false otherwise
     */
    bool operator==(const EdgeIterator & other) const
    {
      return m_graph == other.m_graph && m_iter == other.m_iter;
    }

   private:
   
    friend class Graph;
    
    using iterator = std::vector<edge_uid_type>::const_iterator;
    
    bool is_valid() const
    {
      return m_graph != nullptr;
    }
    
    EdgeIterator(Graph * const graph, 
                 const iterator iter)
    : m_graph(graph),
      m_iter(iter)
    {
      assert(is_valid());
    }
    
    /** Pointer to the graph object **/
    Graph * m_graph;
    
    /** Index of node in the graph **/
    iterator m_iter;
  };
  
  /** Edge begin iterator
   * @return an iterator to the beginning of edge collection 
   **/
  edge_iterator edge_begin() const
  {
    // HACK - UB if called on a truly const object
    return edge_iterator(const_cast<Graph *>(this), m_edge_uids.begin());
  }
  
  /** Edge end iterator
   * @return an iterator past-the-end of edge collection 
   **/
  edge_iterator edge_end() const
  {
    // HACK - UB if called on a truly const object
    return edge_iterator(const_cast<Graph *>(this), m_edge_uids.end());
  }

 private:
  
  // Enforce the size constraint
  static_assert(sizeof(Node) <= 16, "Node class exceeds size limit");
  static_assert(sizeof(Edge) <= 32, "Edge class exceeds size limit");
              
  /** Find an edge between nodes a and b.
   * @param a uid of first a
   * @param b uid of second
   * @return a pair of: bool indicating success/failure
   *         and iterator to a's indicence list that can
   *         be dereferenced (if success) to get edge uid.
   *         
   * Internal method. Sensitive to order of a and b.
   * Needed since the incidence lookup map was removed.
   * Linear search in a small vector is faster than a map.
   */
  std::pair<bool, std::vector<edge_uid_type>::iterator> 
  edge_lookup(const node_uid_type a, const node_uid_type b)
  {
    auto & edges = m_incidence[m_node_indices[a]]; 
    auto it = std::find_if(edges.begin(), edges.end(),
                           [&](edge_uid_type const & euid)
                           { 
                             auto const & nodes = m_edge_nodes[m_edge_indices[euid]];
                             return (nodes.first == a && nodes.second == b) ||
                                    (nodes.first == b && nodes.second == a);
                           });
    return std::make_pair(it != edges.end(), it);
  }
  
  std::pair<bool, std::vector<edge_uid_type>::const_iterator> 
  edge_lookup(const node_uid_type a, const node_uid_type b) const
  {
    auto const lookup = const_cast<Graph *>(this)->edge_lookup(a, b);
    return std::pair<bool, std::vector<edge_uid_type>::const_iterator>(lookup.first, lookup.second);
  }
  
  /** Maps node ordinal index to its persistent uid **/
  std::vector<node_uid_type> m_node_uids;
  
  /** Maps node persistent uid to its ordinal index **/
  std::vector<node_idx_type> m_node_indices;
  
  /** Storage for node points, addressed via node index **/
  std::vector<Point> m_node_positions;
  
  /** Storage for node values, addressed via node index **/
  std::vector<node_value_type> m_node_values;
  
  /** Storage for node incidence, addressed via node index **/
  std::vector<std::vector<edge_uid_type>> m_incidence;
  
  /** Maps edge ordinal index to its persistent uid **/
  std::vector<edge_uid_type> m_edge_uids;
  
  /** Maps edge persistent uid to its ordinal index **/
  std::vector<edge_idx_type> m_edge_indices;
  
  /** Storage for edge nodes, addressed via edge index **/
  std::vector<std::pair<node_uid_type, node_uid_type>> m_edge_nodes;
  
  /** Storage for edge values, addressed via edge index **/
  std::vector<edge_value_type> m_edge_values;

};

#endif // CME212_GRAPH_HPP
