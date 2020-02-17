#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

 /** Function to implement vector swap pop with a predicate */
/**
 * @brief Function to remove an element from a vector using swap pop
 * @param[in,out] The vector to perform swap pop on
 * @param[in] An iterator pointing to the element of @vector to remove
 *
 * @post The element pointed by @it is removed, replaced by the element
 * that was at the end of the vector. The vector's size is reduced by 1
 *
 * Complexity: O(1)
 */
template<typename T>
void swap_pop(std::vector<T>& vector, typename std::vector<T>::iterator it) {
  *it = vector[vector.size() - 1];
  vector.pop_back();
}

/**
 * @brief Function to selectively remove elements from a vector
 * that satisfy a predicate. Removal is done using swap pop
 * @param[in,out] The vector to perform selective removal on.
 * @param[in] The functor used as predicate
 *
 * @post Elements of the vector which satisfy @pred are removed.
 * Vector ordering is not reserved
 *
 * Complexity: O(1)
 */
template<typename T, typename Pred>
void erase_if(std::vector<T>& vector, Pred pred) {
  for(auto it{vector.begin()}; it != vector.end();) {
    if (pred(*it))
      swap_pop(vector, it);
    else {
      ++it;
    }
  }
}

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

  /** Type of the nodes' values */
  using node_value_type = V;

  /** Type of the edge's values */
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    m_nextUid = 0;
    m_nextUidEdge = 0;
  }

  /** Destructor*/
  ~Graph() {
    clear();
  };

 private:
  /**
   * @brief Struct to store edge data
   */
  struct EdgeData : private totally_ordered<EdgeData> {
    size_type index;
    edge_value_type value;
    std::pair<size_type, size_type> nodes;

    EdgeData() : index(-1), value(-1), nodes() {}
    EdgeData(size_type i_index, edge_value_type i_value, std::pair<size_type, size_type> i_nodes)
      : index{i_index}, value{i_value}, nodes{i_nodes} {};
  };

  /**
   * @brief Struct to store node data
   */
  struct NodeData : private totally_ordered<NodeData> {
    size_type index;
    Point point;
    node_value_type value;
    std::unordered_map<size_type, size_type> connectingNodes; // {node_uid, edge_uid}

//    NodeData() : uid(-1), index(-1), point(), value(-1), connectingNodes() {}
    NodeData(size_type i_index, const Point& i_point, const node_value_type& i_value)
      : index{i_index}, point{i_point}, value{i_value}, connectingNodes() {};
  };

 public:

  //
  // NODES
  //
  /** Struct for error handling */
  struct InvalidNode {};

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
      m_graph = nullptr;
      m_uid = size_type(-1);
    }

    /** Return this node's position. 
     * @pre: This node is valid and belongs to a graph
     * Complexity: O(1)
     */
    Point& position() {
      if (!valid()) throw InvalidNode{};
      return m_graph->m_uidToNodeData[m_uid].point;
    }

    const Point& position() const {
      // HW0: YOUR CODE HERE
      if (!valid()) throw InvalidNode{};
      return m_graph->m_uidToNodeData[m_uid].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). 
     * Complexity: O(1)
     */
    size_type index() const {
      // HW0: YOUR CODE HERE
      if (!valid()) throw InvalidNode{};
      return m_graph->m_uidToNodeData[m_uid].index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /**
     * @brief Get the node's coressponding value
     * @return Node's value
     *
     * Complexity: O(1)
     */
    node_value_type& value() {
      if (!valid()) throw InvalidNode{};
      return m_graph->m_uidToNodeData[m_uid].value;
    }

    const node_value_type& value() const {
      if (!valid()) throw InvalidNode{};
      return m_graph->m_uidToNodeData[m_uid].value;
    }

    /**
     * @brief Get number of incident nodes from this node (i.e. degree)
     * @return This node's degree
     *
     * Complexity: O(1)
     */
    size_type degree() const {
      if (!valid()) throw InvalidNode{};
      return m_graph->m_uidToNodeData[m_uid].connectingNodes.size();
    }

    /**
     * @brief Get incident iterator to first connection
     * @return Incident iterator to an edge that is this node's first connection
     *
     * Complexity: O(1)
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(m_graph, m_uid, m_graph->m_uidToNodeData[m_uid].connectingNodes.begin());
    }

    /**
     * @brief Get incident iterator to one past last connection
     * @return Incident iterator to an invalid edge that is one past this node's last connection
     *
     * Complexity: O(1)
     */
    incident_iterator edge_end() const {
      return IncidentIterator(m_graph, m_uid, m_graph->m_uidToNodeData[m_uid].connectingNodes.end());
    }

    /**
     * @brief Test whether this node and @a n are equal. Equal nodes have the same graph and the same index.
     * @return true if this node and node a are equal, fase otherwise
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (!valid() || !n.valid()) throw InvalidNode{};
      if (m_graph == n.m_graph && m_uid == n.m_uid)
        return true;
      else
        return false;
    }

    /**
     * @brief Test whether this node is less than @a n in a global order.
     * @return true if this node is less than another node n, false otherwise
     * This ordering function is useful for STL containers such as
     * std::map<>. It does not have any geometric meaning.
     *
     * This node ordering relation obeys trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     *
     * Complexity: O(1)
     */
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (!valid() || !n.valid()) throw InvalidNode{};
      if (m_graph == n.m_graph)
        return m_uid < n.m_uid;
      else {
        std::less<Graph*> less_than_graph;
        return less_than_graph(m_graph, n.m_graph);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    /** Private constructor */
    Node(const Graph* pGraph, size_type uid)
      : m_graph(const_cast<Graph*>(pGraph)), m_uid(uid) {}
    /**
     * @brief Check if the node is valid (satisfy representation invariance
     * 1. UID is in range: 0 <= m_uid <= m_uidToNodeData.size()
     * 2. Index stored in NodeData is in range: 0 <= m_uidToNodeData[m_uid].index < m_nodeIndexToUid.size()
     * 3. Index is synced: m_nodeIndexToUid[m_uidToNodeData[m_uid].index] == m_uid
     *
     * Complexity: O(1)
     */
    bool valid() const {
      return m_uid>=0 && m_uid < m_graph->m_uidToNodeData.size()  // uid in range
          && m_graph->m_uidToNodeData[m_uid].index>=0 && m_graph->m_uidToNodeData[m_uid].index < m_graph->m_nodeIndexToUid.size() // index in range
          && m_graph->m_nodeIndexToUid[m_graph->m_uidToNodeData[m_uid].index] == m_uid;  // uid in sync
    }
    Graph* m_graph;
    size_type m_uid;
  };

  /**
   * @brief Get the total number of nodes in the graph.
   * @return Total number of nodes in the graph
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return m_nodeIndexToUid.size();
  }

  /**
   * @brief Synonym for size().
   */
  size_type num_nodes() const {
    return size();
  }

  /**
   * @brief Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1)
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type currUid = m_nextUid++;
    // Store the node
    m_uidToNodeData.push_back(NodeData(m_nodeIndexToUid.size(), position, value));
    m_nodeIndexToUid.push_back(currUid);
    // Make a node to return
    return Node(this, currUid);

  }

  /**
   * @brief Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.m_graph == this && n.valid();
  }

  /**
   * @brief Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < num_nodes());
    return Node(this, m_nodeIndexToUid[i]);
  }

  //
  // EDGES
  //
 public:
  struct InvalidEdge {}; // For exception handling

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   *
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /**
     * @brief Construct an invalid Edge.
     * @return An invalid Edge
     *
     * Complexity: O(1)
     */
    Edge() {
      // HW0: YOUR CODE HERE
      m_graph = nullptr;
      m_uidEdge = size_type(-1);
    }

    /**
     * @brief Return the first node of this Edge
     * Edges are unordered ih Graph, so this first node is only the first one used to construct this instance of edge
     * @return The first Node of the Edge
     *
     * Comlexity: O(1)
     */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if (!valid()) throw InvalidEdge{};
      if (m_firstIsHome)
        return Node(m_graph, m_graph->m_uidToEdgeData[m_uidEdge].nodes.first);
      else
        return Node(m_graph, m_graph->m_uidToEdgeData[m_uidEdge].nodes.second);
    }

    /**
     * @brief Return the second node of this Edge
     * @return The second Node of the Edge
     *
     * Complexity: O(1)
     */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (!valid()) throw InvalidEdge{};
      if (m_firstIsHome)
        return Node(m_graph, m_graph->m_uidToEdgeData[m_uidEdge].nodes.second);
      else
        return Node(m_graph, m_graph->m_uidToEdgeData[m_uidEdge].nodes.first);
    }

    /**
     * @brief Get this edge's value
     * @return Reference to this edge's value
     *
     * Complexity: O(1)
     */
    edge_value_type& value() {
      if (!valid()) throw InvalidEdge{};
      return m_graph->m_uidToEdgeData[m_uidEdge].value;
    }

    const edge_value_type& value() const {
      if (!valid()) throw InvalidEdge{};
      return m_graph->m_uidToEdgeData[m_uidEdge].value;
    }

    /**
     * @brief Get this edge's length
     * @return The length of this edge
     *
     * Complexity: O(1)
     */
    double length() const {
      if (!valid()) throw InvalidEdge{};
      return norm(node1().position() - node2().position());
    }

    /**
     * @brief Test whether this edge and @a e are equal.
     * Equal edges represent the same undirected edge between two nodes.
     * @return true if this edge and @e are equal, false otherwise
     */
    bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE
      if (!valid() || !e.valid()) throw InvalidEdge{};
      if (m_graph == e.m_graph && m_uidEdge == e.m_uidEdge)
        return true;
      else
        return false;
    }

    /**
     * @brief Test whether this edge is less than @a e in a global order.
     * @return true of this edge is less than @e, false otherwise.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE
      if (!valid() || !e.valid()) throw InvalidEdge{};
      if (m_graph == e.m_graph)
        return m_uidEdge < e.m_uidEdge;
      else {
        std::less<Graph*> less_than_graph;
        return less_than_graph(m_graph, e.m_graph);
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    /* Private constructor */
    Edge(const Graph* pGraph, size_type uid, bool firstIsHome = true)
      : m_graph(const_cast<Graph*>(pGraph)), m_uidEdge(uid) , m_firstIsHome(firstIsHome) {}

    /**
     * @brief Check if the edge is valid (satisfy representation invariance
     * 1. UID is in range: 0 <= m_uidEdge <= m_uidToEdgeData.size()
     * 2. Index stored in NodeData is in range: 0 <= m_uidToEdgeData[m_uidEdge].index < m_edgeIndexToUid.size()
     * 3. Index is synced: m_edgeIndexToUid[m_uidToEdgeData[m_uid].index] == m_uidEdge
     *
     * Complexity: O(1)
     */
    bool valid() const {
      return m_uidEdge>=0 && m_uidEdge < m_graph->m_uidToEdgeData.size()  // uid in range
          && m_graph->m_uidToEdgeData[m_uidEdge].index>=0 && m_graph->m_uidToEdgeData[m_uidEdge].index < m_graph->m_edgeIndexToUid.size() // index in range
          && m_graph->m_edgeIndexToUid[m_graph->m_uidToEdgeData[m_uidEdge].index] == m_uidEdge;  // uid in sync
    }
    Graph* m_graph;
    size_type m_uidEdge;
    bool m_firstIsHome;
  };

  /**
   * @brief Return the total number of edges in the graph.
   * @return The number of unique edges in the graph
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return m_edgeIndexToUid.size();
  }

  /**
   * @brief Return the edge with index @a i.
   * @return The Edge of the graph with index @i
   *
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < num_edges());
    return Edge(this, m_edgeIndexToUid[i]);
  }

  /**
   * @brief Test whether two nodes are connected by an edge.
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * @pre @a a and @a b are valid nodes of this graph
   * If either @a or @b is invalid, throw InvalidNode
   * If they belong to different graphs, return false
   *
   * Complexity: Average O(1), worst O(a.degree())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (!a.valid() || !b.valid()) throw InvalidNode{};
    if (!(a.m_graph == b.m_graph)) return false;
    else
      return get_edge_uid(a, b) != size_type(-1);
  }

  /**
   * @brief Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: Average O(1), worst O(a.degree() + b.degree())
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if (!a.valid() || !b.valid()) throw InvalidNode{};
    size_type uidEdge {get_edge_uid(a, b)};
    if (uidEdge != size_type(-1))
      return Edge(this, uidEdge);
    else {
      size_type currUidEdge {m_nextUidEdge++};
      // Add connections
      m_uidToEdgeData.push_back(EdgeData(m_edgeIndexToUid.size(), value, {a.m_uid, b.m_uid}));
      m_edgeIndexToUid.push_back(currUidEdge);
      m_uidToNodeData[a.m_uid].connectingNodes.insert({b.m_uid, currUidEdge});
      m_uidToNodeData[b.m_uid].connectingNodes.insert({a.m_uid, currUidEdge});
      return Edge(this, currUidEdge);
    }
  }

  /**
   * @brief Remove all nodes and edges from this graph.
   *
   * @post num_nodes() == 0 && num_edges() == 0
   * @post Invalidates all outstanding Node and Edge objects.
   * @code
   * Node n = graph.node(i); // 0 <= i <= graph.num_nodes()
   * Edge e = graph.edge(i); // 0 <= i <= graph.num_edges()
   * n.valid(); // True
   * e.valid(); // True
   * graph.clear();
   * n.valid(); // False
   * e.valid(); // False
   * @code
   */
  void clear() {
    // HW0: YOUR CODE HERE
    NodeIterator it {node_begin()};
    while (m_nodeIndexToUid.size() != 0)
      remove_node(it);
    m_nextUidEdge = 0;
    m_nextUid = 0;
    m_uidToNodeData.clear();
    m_uidToEdgeData.clear();
    m_nodeIndexToUid.clear();
    m_edgeIndexToUid.clear();
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
//    using iterator_category = std::forward_iterator_tag;  // For compiling on OSX

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
     * @brief Overload deference operator
     * @return The node this iterator is pointing to
     *
     * Complexity: O(1)
     */
    Node operator*() const {
      return Node(m_graph, *m_nodeIndexToUidIt);
    }

    /**
     * @brief Overload ++ operator
     * @return Increment this iterator to so that it point to the next node in the graph. Return this (incremented) iterator.
     *
     * @post Iterator now points to node with the next index
     * @code
     * Node n1 {*(It++)};
     * Node n2 {*It};
     * n1.index() == n2.index() - 1; // return true
     * @code
     *
     * Complexity: O(1)
     */
    NodeIterator& operator++() {
      ++m_nodeIndexToUidIt;
      return *this;
    }

    /**
     * @brief Overload == operator
     * @return true if the node this iterator is pointing to is the same as the node @nodeIt is pointing to, false otherwise.
     *
     * Complexity: O(1)
     */
    bool operator==(const NodeIterator nodeIt) const {
      return (m_graph == nodeIt.m_graph && m_nodeIndexToUidIt == nodeIt.m_nodeIndexToUidIt);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    /**
     * @brief Private constructor
     * @param[in] graph The graph this iterator belongs to
     * @param[in] Iterator to graph's m_nodeIndexToUid, pointing to the UID for this Node Iterator
     */
    NodeIterator(const graph_type* graph, std::vector<size_type>::const_iterator nodeIndexToUidIt)
      : m_graph{const_cast<Graph*>(graph)}, m_nodeIndexToUidIt{nodeIndexToUidIt} {}
    graph_type* m_graph;
    std::vector<size_type>::const_iterator m_nodeIndexToUidIt;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * @brief Get iterator to first node in the graph.
   * @return Iterator pointing to the first node in the graph.
   *
   * Complexity: O(1)
   */
  // TODO: can I omit const?
  node_iterator node_begin() const {
    return node_iterator(this, m_nodeIndexToUid.begin());
  }

  /**
   * @brief Get iterator to one beyond last node.
   * @return Iterator pointing to an invalid node that is one beyond the last node in the graph.
   *
   * Complexity: O(1)
   */
  node_iterator node_end() const {
    return node_iterator(this, m_nodeIndexToUid.end());
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
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
     * @brief Overload deference operator
     * @return The connecting edge from the this iterator's home node it is pointing to.
     *
     * Complexity: O(1)
     */
    Edge operator*() const {
      if (m_uid == (m_graph->m_uidToEdgeData[m_connectionIt->second]).nodes.first)
        return Edge(m_graph, m_connectionIt->second);
      else
        return Edge(m_graph, m_connectionIt->second, false);
    }

    /**
     * @brief Overload ++ operator
     * @return Increment this iterator so that it points to the next connecting edge from this iterator's home node.
     * Return this (incremented) iterator.
     *
     * @post Iterator now points to the next connecting edge
     * @code
     * Edge e1 {*(It++)};
     * if (It != e1.node1().end()) {
     *    Edge e2 = *It;
     *    e1.node1() == e2.node(1) && e1.node2 != e2.node2(); // return true
     * }
     * @code
     *
     * Complexity: O(1)
     */
    IncidentIterator& operator++() {
      ++m_connectionIt;
      return *this;
    }

    /**
     * @brief Overload == operator
     * @return true if this iterator points to the same connecting edge from the same home node as @incidentIt does, false otherwise.
     *
     * Complexity: O(1)
     */
    bool operator==(const IncidentIterator& incidentIt)  const {
      return (m_graph == incidentIt.m_graph && m_uid == incidentIt.m_uid && m_connectionIt == incidentIt.m_connectionIt);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    /**
     * @brief Private constructor
     * @param[in] pGraph The graph this iterator belongs to
     * @param[in] uid The UID of the home node
     * @param[in] it Iterator of the map in the graph storing edges starting from home node
     *
     * Complexity O(1)
     */
    IncidentIterator(const graph_type* pGraph, size_type uid, std::unordered_map<size_type, size_type>::iterator it) :
      m_graph{const_cast<graph_type*>(pGraph)}, m_uid{uid}, m_connectionIt{it} {}
    graph_type* m_graph;
    size_type m_uid; // UID of starting node
    std::unordered_map<size_type, size_type>::iterator m_connectionIt; // Iterator of map in starting node's data storing connecting-node uid -> edge uid
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
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
     * @brief Overload dereference operator
     * @return The Edge in the graph this iterator is pointing to.
     *
     * Complexity: O(1)
     */
    Edge operator*() const {
      return Edge(m_graph, *m_edgeIndexToUidIt);
    }

    /**
     * @brief Overload ++ operator
     * @return Increment this iterator so that it now points to the next edge in the graph. Return this (incremented) iterator.
     *
     * @post Iterator now points to the next edge
     * @code
     * Edge e1 {*(It++)};
     * Edge e2 {*It};
     * e2 != e1; // return true
     * @code
     *
     * Complexity: O(1)
     */
    EdgeIterator& operator++() {
      ++m_edgeIndexToUidIt;
      return *this;
    }

    /**
     * @brief Overload == operator
     * @return true if this iterator is pointing to the same edge of the same graph as @edgeIt does, false otherwise.
     *
     * Complexity: O(1)
     */
    bool operator==(const EdgeIterator& edgeIt) const {
      return (m_graph == edgeIt.m_graph && m_edgeIndexToUidIt == edgeIt.m_edgeIndexToUidIt);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    /**
     * @brief Private constructor
     * @param[in] pGraph The graph this iterator belongs to
     * @param[in] m_it Iterator to an element in pGraph->m_edgeIndexToUid
     *
     * Complexity: O(1)
     */
    EdgeIterator(const Graph* pGraph, typename std::vector<size_type>::const_iterator it) :
        m_graph{const_cast<Graph*>(pGraph)}, m_edgeIndexToUidIt{it} {}
    Graph* m_graph;
    std::vector<size_type>::const_iterator m_edgeIndexToUidIt;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * @brief Get iterator to first edge of the graph
   * @return Iterator pointing to the first edge of the graph
   *
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, m_edgeIndexToUid.begin());
  }

  /**
   * @brief Get iterator to one beyond last edge
   * @return Iterator pointing to an invalid edge that is one beyond the last edge of the graph
   *
   * Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, m_edgeIndexToUid.end());
  }

  /** Remove methods */
  /**
   * @brief Function to remove an edge
   * @param[in] The edge to remove
   * @return Number of removed edge (1) if successful
   *
   * @pre @e is a valid edge
   * @post Edge @e is removed from the graph and is now invalid
   * Edge indices are not conserved by this operation
   * @code
   * n1 = e.node1();
   * n2 = e.node2();
   * int deg1 = n1.degree();
   * int deg2 = n2.degree();
   * int nE = graph.num_edges();
   * Edge e1 = graph.edge(i); // 0 <= i < graph.num_edges && i != e.index()
   * graph.remove_edge(e);
   * deg1 == n1.degree() + 1 && deg2 == n2.degree() + 1 && nE = graph.num_edges() + 1; // True
   * e1.index() == i; // Not guaranteed
   * @code
   *
   * Complexity: Average O(1), worst O(n1.degree() + n2.degree()), n1 and n2 are nodes of the edge
   */
  size_type remove_edge(const Edge& e) {
    if(!e.valid()) throw InvalidEdge();
    remove_edge_uid(e.m_uidEdge);
    return 1;
  }

  /**
   * @brief Function to remove an edge
   * @param[in] Iterator to the edge to be removed
   * @return Iterator pointing to an unvisited edge (forward iterator)
   *
   * @pre @e_it is points to a valid edge
   * @post Edge at @e_it is removed from the graph and is now invalid
   * Edge indices are not conserved by this operation
   * @code
   * Edge e = *e_it;
   * n1 = e.node1();
   * n2 = e.node2();
   * int deg1 = n1.degree();
   * int deg2 = n2.degree();
   * int nE = graph.num_edges();
   * Edge e1 = graph.edge(i); // 0 <= i < graph.num_edges && i != e.index()
   * graph.remove_edge(e_it);
   * deg1 == n1.degree() + 1 && deg2 == n2.degree() + 1 && nE = graph.num_edges() + 1; // True
   * e1.index() == i; // Not guaranteed
   * @code
   *
   * Complexity: Average O(1), worst O(n1.degree() + n2.degree()), n1 and n2 are nodes of the edge
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (!(*e_it).valid()) throw InvalidEdge();
    remove_edge_uid((*e_it).m_uidEdge);
    return e_it;
  }

  /**
   * @brief Function to remove an edge
   * @param[in] The first node of the edge
   * @param[in] The second ndoe of the edge
   * @return Number of removed edge (1) if successful
   *
   * @pre @a and @b are valid nodes
   * @post Edge that contains @a and @b is removed from the graph and is now invalid
   * Edge indices are not conserved by this operation
   * @code
   * Edge e; // A valid edge
   * n1 = e.node1();
   * n2 = e.node2();
   * int deg1 = n1.degree();
   * int deg2 = n2.degree();
   * int nE = graph.num_edges();
   * Edge e1 = graph.edge(i); // 0 <= i < graph.num_edges && i != e.index()
   * graph.remove_edge(n1, n2);
   * deg1 == n1.degree() + 1 && deg2 == n2.degree() + 1 && nE = graph.num_edges() + 1; // True
   * e1.index() == i; // Not guaranteed
   * @code
   *
   * Complexity: Average O(1), worst O(n1.degree() + n2.degree()), n1 and n2 are nodes of the edge
   */
  size_type remove_edge(const Node& a, const Node & b) {
    if (!a.valid() || !b.valid()) throw InvalidNode();
    size_type edgeUid {get_edge_uid(a, b)};
    if (edgeUid != size_type(-1)) {
      remove_edge(Edge(this, edgeUid));
      return 1;
    } else
      return 0;
  }

  /**
   * @brief Function to remove a node
   * @param[in] The node to remove
   * @return Number of removed node (1) if successful
   *
   * @pre @n is a valid node
   * @post Node @n and assiciating edges is removed from the graph and is now invalid
   * Node and edge indices are not conserved by this operation
   * @code
   * Edge e; // A valid edge
   * n1 = e.node1();
   * n2 = e.node2();
   * int deg1 = n1.degree();
   * int deg2 = n2.degree();
   * int nI = n2.index();
   * int nE = graph.num_edges();
   * Edge e1 = graph.edge(i); // 0 <= i < graph.num_edges && i != e.index()
   * graph.remove_node(n1);
   * deg2 == n2.degree() + 1 && nE = graph.num_edges() + deg1 && !e.valid(); // True
   * e1.index() == i; // Not guaranteed
   * n2.index() == i; // Not guaranteed
   * @code
   *
   * Complexity: Average O(1), worst O(degree()^2), still much less than O(num_nodes + num_edges) for sparse graphs
   */
  size_type remove_node(const Node& n) {
    if (!n.valid()) throw InvalidNode();
    remove_node_uid(n.m_uid);
    return 1;
  }

  /**
   * @brief Function to remove a node
   * @param[in] The iterator to the node to remove
   * @return The iterator pointing to an unvisited node (forward iteratot)
   *
   * @pre @n_it points to a valid node
   * @post Node referenced by n_it and assiciating edges is removed from the graph and is now invalid
   * Node and edge indices are not conserved by this operation
   * @code
   * Edge e; // A valid edge
   * n1 = e.node1();
   * NodeIterator n1_it; // An iterator pointing to node Sn1
   * n2 = e.node2();
   * int deg1 = n1.degree();
   * int deg2 = n2.degree();
   * int nI = n2.index();
   * int nE = graph.num_edges();
   * Edge e1 = graph.edge(i); // 0 <= i < graph.num_edges && i != e.index()
   * graph.remove_node(n1_it);
   * deg2 == n2.degree() + 1 && nE = graph.num_edges() + deg1 && !e.valid() // True
   * e1.index() == i; // Not guaranteed
   * n2.index() == i; // Not guaranteed
   * @code
   *
   * Complexity: Average O(1), worst O(degree()^2), still much less than O(num_nodes + num_edges) for sparse graphs
   */
  node_iterator remove_node(node_iterator n_it) {
    if (!(*n_it).valid()) throw InvalidNode();
    remove_node_uid((*n_it).m_uid);
    return n_it;
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  /**
   * @brief Get the edge uid from two nodes if it exists. Otherwise return a nullptr
   * @param[in] a One node of the edge
   * @param[in] b The other node of the edge
   * @return The edge uid that connects nodes @a and @b. If the edge does not exist, return nullptr
   *
   * @pre @a a and @a b are distinct valid nodes of this graph
   * Complexity: Average O(1), worst O(a.degree())
   */
  size_type get_edge_uid(const Node& a, const Node& b) const {
    if (!a.valid() || !b.valid()) throw InvalidNode();
    if (a.m_graph != b.m_graph) return size_type(-1);
    auto it = m_uidToNodeData[a.m_uid].connectingNodes.find(b.m_uid);
    if (it != m_uidToNodeData[a.m_uid].connectingNodes.end())
      return it->second;
    return size_type(-1);
  }

  /**
   * @brief Remove a valid edge
   * @param[in] The UID of the edge to be removed
   *
   * @pre The UID is of a valid edge
   * @post The edge is removed. See complete postcondition at specific public remove_edge functions
   *
   * Complexity: Average O(1), worst O(n1.degree() + n2.degree()), n1 and n2 are nodes of the edge
   */
  void remove_edge_uid(size_type uid) {
    // Remove from each node's connectingNodes
    auto nodes = m_uidToEdgeData[uid].nodes;
    m_uidToNodeData[nodes.first].connectingNodes.erase(nodes.second);
    m_uidToNodeData[nodes.second].connectingNodes.erase(nodes.first);
    // Remove from m_edgeIndexToUid, update edge index at this place, and invalidate e
    size_type index {m_uidToEdgeData[uid].index};
    swap_pop(m_edgeIndexToUid, m_edgeIndexToUid.begin() + index);
    m_uidToEdgeData[m_edgeIndexToUid[index]].index = index;
    m_uidToEdgeData[uid].index = size_type(-1);
  }

  /**
   * @brief Remove a valid node
   * @param[in] The UID of the node to be removed
   *
   * @pre The UID is of a valid node
   * @post The node and all associated edge is removed. See complete postcondition at specific public remove_node funcitons
   *
   * Complexity: Average O(1), worst O(degree()^2), still much less than O(num_nodes + num_edges) for sparse graphs
   */
  void remove_node_uid(size_type uid) {
    // Remove node connections
    while (!m_uidToNodeData[uid].connectingNodes.empty()) {
      remove_edge_uid(m_uidToNodeData[uid].connectingNodes.begin()->second);
    }
    // Remove from m_nodeIndexToUid, update node index at this place, and invalidate n
    size_type index ={m_uidToNodeData[uid].index};
    swap_pop(m_nodeIndexToUid, m_nodeIndexToUid.begin() + index);
    m_uidToNodeData[m_nodeIndexToUid[index]].index = index;
    m_uidToNodeData[uid].index = size_type(-1);
  }


  size_type m_nextUid;
  std::vector<NodeData> m_uidToNodeData;
  std::vector<size_type> m_nodeIndexToUid;

  size_type m_nextUidEdge;
  std::vector<EdgeData> m_uidToEdgeData;
  std::vector<size_type> m_edgeIndexToUid;
};

#endif // CME212_GRAPH_HPP
