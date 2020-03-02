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
template <typename V, typename E>
class Graph
{
public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Synonym for V */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Synonym for E */
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
   * @brief Represents each Node's data & adjacent edges.
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
    NodeData(const Point& position, 
             const node_value_type& data = node_value_type()) 
      : position(position), data(data)
    {
    }

    Point position;
    node_value_type data;
    node_edge_map neighbors;
  };

  /**
   * @struct EdgeData
   * @brief Represents each Edge's data.
   */
  struct EdgeData
  {
    /**
     * @brief Constructs an EdgeData object
     * 
     * Constructs an EdgeData object, which holds the Edge object itself and
     * the associated custom data, like an edge weight for example.
     * 
     * @param[in] _n1_    The first node of the edge
     * @param[in] _n2_    The second node of the edge
     * @param[in] _data_  The associated data for this edge
     * @pre _n1_ and _n2_ are valid nodes in the same graph.
     * @post edge_data.edge == Edge(_n1_, _n2_)
     * @post edge_data.data == _data_
     */
    EdgeData(const Node& n1, const Node& n2, 
             const edge_value_type& data = edge_value_type()) 
      : edge(n1, n2), data(data)
    {
    } 

    Edge edge;
    edge_value_type data;
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
      return graph->node_data.at(index()).position;
    }

    /** @return mutable reference to this point's position */
    Point& position()
    {
      return const_cast<Graph *>(graph)->node_data.at(index()).position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return idx;
    }

    /** @return a mutable reference to this node's data. */
    node_value_type& value()
    {
      return const_cast<Graph *>(graph)->node_data.at(index()).data;
    }

    /** @return a const reference to this node's data. */
    const node_value_type& value() const
    {
      return graph->node_data.at(index()).data;
    }

    /** @return the number of neighbors for this node. */
    size_type degree() const
    {
      return graph->node_data.at(index()).neighbors.size();
    }

    /** @return an iterator to the first incident edge to _src_node_. */
    incident_iterator edge_begin() const
    {
      const auto& itr = graph->node_data.at(index()).neighbors.begin();
      return incident_iterator(graph, index(), itr);
    }

    /** 
     * @return an iterator to just past the last incident edge to _src_node_ 
     */
    incident_iterator edge_end() const
    {
      const auto& itr = graph->node_data.at(index()).neighbors.end();
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
    return node_data.size();
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
    NodeData new_node_data(position, value);
    Node new_node(this, node_data.size());
    node_data.push_back(new_node_data);
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

    /** 
     * @return Mutable reference to this Edge's value 
     * @pre this Edge is a valid edge
     * 
     * Complexity: O(1)
     */
    edge_value_type& value()
    {
      // first, get the edge index
      edge_index_type edge_index;
      assert(n1.graph->has_edge(n1, n2, edge_index));
      return const_cast<Graph *>(n1.graph)->edge_data.at(edge_index).data;
    }

    /** 
     * @return const reference to this Edge's value 
     * @pre this Edge is a valid edge
     * 
     * Complexity: O(1)
     */
    const edge_value_type& value() const
    {
      // first, get the edge index
      edge_index_type edge_index;
      assert(n1.graph->has_edge(n1, n2, edge_index));
      return n1.graph->edge_data.at(edge_index).data;
    }

    /** @return the distance between the edge's nodes. */
    double length() const
    {
      return norm(n1.position() - n2.position());
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

    /**
     * Constructs an Edge instance.
     * 
     * @param[in] _a_ The first node of this Edge
     * @param[in] _b_ The second node of this Edge
     */
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
    return edge_data.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    return edge_data[i].edge;
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
    const edge_index_type new_edge_index = edge_data.size();
    node_data[a_index].neighbors[b_index] = new_edge_index;
    node_data[b_index].neighbors[a_index] = new_edge_index;

    EdgeData new_edge_data(a, b);
    edge_data.push_back(new_edge_data);
    return new_edge_data.edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   * 
   * Complexity: O(num_nodes + num_edges)
   */
  void clear()
  {
    node_data.clear();
    edge_data.clear();
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
    return node_iterator(this, node_data.size());
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
     * Compares if two iterators are equal: if they point to the same 
     * graph, have the same source node, and point to the same element in the 
     * adjacency matrix. 
     * 
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
      return (*edge_itr).edge;
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
     * Compares if two EdgeIterator objects point to the same vector and
     * element.
     * 
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
    EdgeIterator(const typename std::vector<EdgeData>::const_iterator& edge_itr) 
      : edge_itr(edge_itr)
    {
    }

    typename std::vector<EdgeData>::const_iterator edge_itr;
  };

  /** @return An iterator to the first edge in the graph */
  edge_iterator edge_begin() const
  {
    return edge_iterator(edge_data.begin());
  }

  /** @return An iterator past the last edge in the graph */
  edge_iterator edge_end() const
  {
    return edge_iterator(edge_data.end());
  }

  /**
   * Removes the node _n_ passed in by the user.
   * @param[in] _n_ The node to remove
   * @return the number of removed nodes.
   * @post if _n_ in graph, old_graph.size() - 1 == new_graph.size()
   * @post if _n_ not in graph, no change to the graph and no invalidations.
   * 
   * Can invalidate outstanding Node objects. In other words, 
   * old_node.index() may not equal new_node.index(). Will invalidate 
   * outstanding node_end() iterators. 
   * 
   * Can invalidate outstanding EdgeIterator objects.
   * 
   * Complexity: O(degree)
   */
  size_type remove_node(const Node& n)
  {
    if (n.graph != this)
      return 0;

    // removing all incident edges is O(degree)
    incident_iterator itr = n.edge_begin();
    while (itr != n.edge_end())
    {
      remove_edge(*itr);
      itr = n.edge_begin();
    }

    const node_index_type n_idx = n.index();
    if (n_idx != node_data.size() - 1)
    {
      NodeData& last_node_data = node_data.back();

      // update adjacency matrix before moving: O(degree)
      for (auto& dest_pair : last_node_data.neighbors)
      {
        node_index_type last_node_dest_idx = dest_pair.first;
        edge_index_type last_node_edge_idx = dest_pair.second;

        // this is O(1) due to unordered_map avg constant search, update, erase
        // first, set n_idx to the edge and then wipe out the last idx
        node_data.at(last_node_dest_idx).neighbors[n_idx] = last_node_edge_idx;
        node_data.at(last_node_dest_idx).neighbors.erase(node_data.size() - 1);

        // O(1): vector access/overwrite of size_type
        edge_data.at(last_node_edge_idx).edge.n1.idx = n_idx;
        edge_data.at(last_node_edge_idx).edge.n2.idx = last_node_dest_idx;
      } 

      // O(1): move 
      // once the matrix is set, we can transplant last node into current index
      node_data[n_idx] = std::move(node_data.back());
    }

    node_data.pop_back();
    return 1;
  }

  /**
   * Removes the node referenced by the current iterator.
   * @param[in] _n_it_ The iterator whose referenced element should be removed 
   * @return iterator to another node that has not been visited yet.
   * 
   * Can invalidate outstanding Node objects. In other words, 
   * old_node.index() may not equal new_node.index(). Will invalidate 
   * outstanding node_end() iterators. 
   * 
   * Can invalidate outstanding EdgeIterator objects.
   * 
   * Complexity: O(degree)
   */
  node_iterator remove_node(node_iterator n_it)
  {
    const node_index_type current_node_idx = (*n_it).index();
    remove_node(*n_it);
    return node_iterator(this, current_node_idx);
  }

  /**
   * Removes the edge connecting Nodes @a a and @a b. No particular ordering
   * required.
   * 
   * @param[in] _a_ One node of the edge.
   * @param[in] _b_ The other node of the edge.
   * @return Number of removed edges. Can be zero.
   * @post old_graph.num_edges() - 1 == new_graph.num_edges() if edge removed
   * @post old_graph.num_edges() == new_graph.num_edges() if no edge removed
   * 
   * May invalidate outstanding edge_end() iterators.
   * 
   * Complexity: O(1)
   */
  size_type remove_edge(const Node& a, const Node& b)
  {
    return remove_edge(Edge(a, b));
  }

  /**
   * Removes any edges equal to the given Edge @a e. 
   * 
   * @param[in] _e_ The edge to remove.
   * @return Number of removed edges. Can be zero.
   * @post old_graph.num_edges() - 1 == new_graph.num_edges() if edge removed
   * @post old_graph.num_edges() == new_graph.num_edges() if no edge removed
   * 
   * May invalidate outstanding edge_end() iterators.
   * 
   * Complexity: O(1)
   */
  size_type remove_edge(const Edge& e)
  {
    edge_index_type edge_index;
    if (!has_edge(e.n1, e.n2, edge_index))
      return 0;

    // clear out n1 -> n2 & n2 -> n1
    // O(1) as indices in map are unique (no hash collision)
    const node_index_type n1_idx = e.n1.index();
    const node_index_type n2_idx = e.n2.index();    
    node_data.at(n1_idx).neighbors.erase(n2_idx);
    node_data.at(n2_idx).neighbors.erase(n1_idx);

    // if we're killing the last index, just wipe it with O(1) pop_back
    if (edge_index == edge_data.size() - 1)
    {
       edge_data.pop_back();
       return 1;
    }

    // otherwise, take the last index and swap it in
    // so first, we need to take the last edge's idx and update to new edge
    // average O(1) for unordered_map search/update
    EdgeData& last_edge_data = edge_data.back();
    const node_index_type last_n1_idx = last_edge_data.edge.n1.index();
    const node_index_type last_n2_idx = last_edge_data.edge.n2.index();
    node_data.at(last_n1_idx).neighbors[last_n2_idx] = edge_index;
    node_data.at(last_n2_idx).neighbors[last_n1_idx] = edge_index;

    // now we move the last element into the old edge index
    // O(1) copy + pop_back
    edge_data[edge_index] = std::move(last_edge_data);
    edge_data.pop_back();
    return 1;
  }

  /**
   * Remove the edge referenced by the @a e_it.
   * 
   * @param[in] _e_it_ The iterator whose edge should be removed.
   * @return an iterator to an unvisited edge.
   * 
   * Complexity: O(1)
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    // usually linear time, but edge_itr is vector::iterator is of type
    // LegacyRandomAccessIterator, turning each itr operation O(1)
    edge_index_type current_edge_idx = std::distance(
      edge_begin().edge_itr, 
      e_it.edge_itr
    );

    remove_edge(*e_it);
    return edge_iterator(edge_begin().edge_itr + current_edge_idx);
  }

private:
  /** Test whether two nodes are connected by an edge, and return the 
   * edge index if it exists.
   * 
   * @param[out] _out_edge_idx_ index of the found edge
   * @pre @a a and @a b are valid nodes of this graph
   * @post graph->edge(_out_edge_idx_) == Edge(a, b)
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   * 
   * Complexity: O(1) average
   */
  bool has_edge(const Node& a, const Node& b, edge_index_type& out_edge_idx) const
  {
    const node_index_type a_index = a.index();
    const node_index_type b_index = b.index();

    const node_edge_map& a_neighbors = node_data.at(a_index).neighbors;
    if (a_neighbors.find(b_index) == a_neighbors.end())
      return false;

    out_edge_idx = a_neighbors.at(b_index);
    return true;
  }

  std::vector<NodeData> node_data;
  std::vector<EdgeData> edge_data;
};

#endif // CME212_GRAPH_HPP
