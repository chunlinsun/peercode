#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph
{
public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
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

private:
  //
  // PRIVATE TYPE DEFINITIONS
  //

  /** Type for node indices */
  using node_index_type = size_type;
  /** Type for edge indices */
  using edge_index_type = size_type;
  /** Type for mapping */
  using node_edge_map = std::unordered_map<node_index_type, edge_index_type>;

  /**
   * @struct Graph::NodeData
   * @brief Struct representing each Node's data & adjacent edges.
   */
  struct NodeData
  {
    /**
     * @brief Constructs a NodeData object
     * 
     * Constructs a NodeData object, which holds both the Node's data
     * as well as the mapping to adjacent node and edge indices.
     * 
     * @param[in] _data_  The data to initialize this Node with
     * @post _NodeData.data == _data_
     */
    NodeData(const node_value_type& data = node_value_type()) : data(data)
    {
    }

    node_value_type data;
    node_edge_map neighbors;
  };


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
  class Node : private totally_ordered<Node>
  {
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
    Node() : graph(nullptr) 
    {
    }

    /** Return this node's position. */
    const Point& position() const
    {
      return graph->points[index()];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return idx;
    }

    /** @return a mutable reference to this node's data. */
    node_value_type& value()
    {
      return const_cast<Graph *>(graph)->node_data_map[index()].data;
    }

    /** @return a const reference to this node's data. */
    const node_value_type& value() const
    {
      return graph->node_data_map.at(index()).data;
    }

    /** @return the number of neighbors for this node. */
    size_type degree() const
    {
      return graph->node_data_map.at(index()).neighbors.size();
    }

    /** @return an iterator to the first incident edge to _src_node_. */
    incident_iterator edge_begin() const
    {
      const auto& itr = graph->node_data_map.at(index()).neighbors.begin();
      return incident_iterator(graph, index(), itr);
    }

    /** 
     * @return an iterator to just past the last incident edge to _src_node_ 
     */
    incident_iterator edge_end() const
    {
      const auto& itr = graph->node_data_map.at(index()).neighbors.end();
      return incident_iterator(graph, index(), itr);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator ==(const Node& n) const
    {
      return (graph == n.graph) && (index() == n.index());
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator <(const Node& n) const
    {
      std::less<const Graph *> graph_less;
      return graph_less(graph, n.graph) || index() < n.index();
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /**
     * @brief Constructs a new Node object
     * @param[in] _graph_   A pointer to the Node's graph
     * @param[in] _idx_ The node's index within the Graph _g_
     * 
     * @pre _graph_ is a valid Graph object (e.g. != nullptr)
     * @pre _idx_ is a valid node index between 0 <= _idx_ < _graph_->size()
     */
    Node(const Graph *graph, node_index_type idx) : graph(graph), idx(idx) 
    { 
    }

    const Graph *graph;
    node_index_type idx;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return points.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] _position_ The new node's position
   * @param[in] _value_ The new node's initialized value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == _value_
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type())
  {
    Node new_node(this, points.size());
    points.push_back(position);
    node_data_map[new_node.index()] = NodeData(value);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const
  {
    return n.graph == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
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
  class Edge : private totally_ordered<Edge>
  {
  public:
    /** Construct an invalid Edge. */
    Edge()
    {
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      return n1;
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator ==(const Edge& e) const
    {
      return (n1 == e.n1 && n2 == e.n2) || (n1 == e.n2 && n2 == e.n1);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator <(const Edge& e) const
    {
      return (n1 < e.n1) || (n1 == e.n1 && n2 < e.n2);
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Edge(const Node& a, const Node& b) : n1(a), n2(b) 
    {
    }

    Node n1;
    Node n2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    return edges[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const
  {
    edge_index_type _unused;
    return has_edge(a, b, _unused);
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
  Edge add_edge(const Node& a, const Node& b)
  {
    edge_index_type existing_edge_idx;

    if (has_edge(a, b, existing_edge_idx))
      return edge(existing_edge_idx);

    const node_index_type a_index = a.index();
    const node_index_type b_index = b.index();
    const edge_index_type new_edge_index = edges.size();
    node_data_map[a_index].neighbors[b_index] = new_edge_index;
    node_data_map[b_index].neighbors[a_index] = new_edge_index;

    Edge new_edge(a, b);
    edges.push_back(new_edge);
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    points.clear();
    edges.clear();
    node_data_map.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
    }

    /** 
     * @return the current Node being iterated upon 
     * 
     * Complexity O(1)
     */
    Node operator *() const
    {
      return graph->node(idx);
    }

    /** 
     * Increments the iterator to next Node.
     * @return the modified iterator pointing to next Node
     */
    NodeIterator& operator ++()
    {
      ++idx;
      return *this;
    }

//--documentation_-1
//--well done with the doxygen style comments
//--END

    /**
     * Returns true if the current iterator is equal to _itr_
     * @param[in] _itr_  The iterator against which to check for equality 
     * @return true if the two NodeIterators are equal
     */
    bool operator ==(const NodeIterator& itr) const
    {
      return (graph == itr.graph) && (idx == itr.idx);
    }

  private:
    friend class Graph;

    /**
     * Constructs a NodeIterator instance.
     * @param[in] _graph_ The graph through which to iterate 
     * @param[in] _idx_   The Node index from which to start iteration
     * 
     * @pre _graph_ is a valid Graph object (e.g. != nullptr)
     */
    NodeIterator(const Graph *graph, size_t idx) : graph(graph), idx(idx)
    {
    }

    const Graph *graph;
    size_type idx;
  };

  /** @return an iterator to the first node in the graph */
  node_iterator node_begin() const
  {
    return node_iterator(this, 0);
  }

  /** @return an iterator to just past the last node in the graph */
  node_iterator node_end() const
  {
    return node_iterator(this, points.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator()
    {
    }

    /**
     * @return the incident Edge object currently being iterated upon 
     * @post Returned Edge should have _src_node_ as n1()
     */
    Edge operator *() const
    {
      return Edge(graph->node(src_node), graph->node(neighbor_map->first));
    }

    /**
     * Moves the iterator to the next incident edge for _src_node_.
     * 
     * @return a modified iterator pointing to the next incident edge
     */
    IncidentIterator& operator ++()
    {
      ++neighbor_map;
      return *this;
    }

    /**
     * @param[in] _itr_   The iterator against which to check for equality
     * @return true if this iterator and _itr_ are equivalent
     */
    bool operator ==(const IncidentIterator& itr) const
    {
      return (graph == itr.graph) && (src_node == itr.src_node) 
        && (neighbor_map == itr.neighbor_map);
    }

  private:
    friend class Graph;

    /**
     * Constructs a new IncidentIterator to the specified graph, source node,
     * and node neighbor map.
     * 
     * @param[in] _graph_       The graph containing the source node
     * @param[in] _src_node_    The source node for incident edge iteration
     * @param[in] _neighbor_map An iterator to the neighbors of _src_node_
     * 
     * @pre _graph_ is a valid Graph object (e.g. != nullptr)
     * @pre  0 <= _src_node_ < _graph_->size()
     */
    IncidentIterator(const Graph *graph, node_index_type src_node, 
                     const node_edge_map::const_iterator& neighbor_map)
      : graph(graph), src_node(src_node), neighbor_map(neighbor_map)
    {
    }

    const Graph *graph;
    node_index_type src_node;
    node_edge_map::const_iterator neighbor_map;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()
    {
    }

    /** @return the current edge that this iterator is pointing to */
    Edge operator *() const
    {
      return *edge_itr;
    }

    /** 
     * Moves the iterator to the next edge in the graph.
     * @return the modified iterator
     */
    EdgeIterator& operator ++()
    {
      ++edge_itr;
      return *this;
    }

    /**
     * @param[in] _itr_ another iterator against which to check for equality
     * @return true if the current iterator is equal to _itr_
     */
    bool operator ==(const EdgeIterator& itr) const
    {
      return edge_itr == itr.edge_itr;
    }


  private:
    friend class Graph;

    /**
     * Constructs an EdgeIterator instance.
     * 
     * @param[in] _edge_itr_  An iterator to the vector Graph::_edges_ 
     */
    EdgeIterator(const typename std::vector<Edge>::const_iterator& edge_itr) 
      : edge_itr(edge_itr)
    {
    }

    typename std::vector<Edge>::const_iterator edge_itr;
  };

  /** @return An iterator to the first edge in the graph */
  edge_iterator edge_begin() const
  {
    return edge_iterator(edges.begin());
  }

  /** @return An iterator past the last edge in the graph */
  edge_iterator edge_end() const
  {
    return edge_iterator(edges.end());
  }

private:

  /** Test whether two nodes are connected by an edge, and return the 
   * edge index if it exists.
   * 
   * @param[out] _out_edge_idx_ index of the found edge
   * @pre @a a and @a b are valid nodes of this graph
   * @post graph->edge(_out_edge_idx_) == Edge(a, b)
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   */
  bool has_edge(const Node& a, const Node& b, edge_index_type& out_edge_idx) const
  {
    const node_index_type a_index = a.index();
    const node_index_type b_index = b.index();

    if (node_data_map.find(a_index) == node_data_map.end())
      return false;

    const node_edge_map& a_neighbors = node_data_map.at(a_index).neighbors;
    if (a_neighbors.find(b_index) == a_neighbors.end())
      return false;

    out_edge_idx = a_neighbors.at(b_index);
    return true;
  }

  std::vector<Point> points;
  std::vector<Edge> edges;
  std::unordered_map<node_index_type, NodeData> node_data_map;
};

#endif // CME212_GRAPH_HPP
