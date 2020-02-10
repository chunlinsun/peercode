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


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    m_numNodes = 0;
    m_nextUid = 0;
    m_numEdges = 0;
    m_nextUidEdge = 0;
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
      m_graph = nullptr;
      m_uid = size_type(-1);
    }

    /** Return this node's position. 
     * @pre: This node is valid and belongs to a graph
     * Complexity: O(1)
     */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return m_graph->m_points[m_uid];
    }

    /** Return this node's index, a number in the range [0, graph_size). 
     * Complexity: O(1)
     */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return m_graph->m_indices[m_uid];
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
      return m_graph->m_uidToValue[m_uid];
    }

    const node_value_type& value() const {
      return m_graph->m_uidToValue[m_uid];
    }

    /**
     * @brief Get number of incident nodes from this node (i.e. degree)
     * @return This node's degree
     *
     * Complexity: O(1)
     */
    size_type degree() const {
      return m_graph->m_fromNode[m_uid].size();
    }

    /**
     * @brief Get incident iterator to first connection
     * @return Incident iterator to an edge that is this node's first connection
     *
     * Complexity: O(1)
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(m_graph, m_uid, m_graph->m_fromNode[m_uid].begin());
    }

    /**
     * @brief Get incident iterator to one past last connection
     * @return Incident iterator to an invalid edge that is one past this node's last connection
     *
     * Complexity: O(1)
     */
    incident_iterator edge_end() const {
      return IncidentIterator(m_graph, m_uid, m_graph->m_fromNode[m_uid].end());
    }

    /**
     * @brief Test whether this node and @a n are equal. Equal nodes have the same graph and the same index.
     * @return true if this node and node a are equal, fase otherwise
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (m_graph == n.m_graph && m_uid == n.m_uid) {
        return true;
      } else {
        return false;
      }
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
      if (m_uid < n.m_uid) {
        return true;
      } else if (m_uid == n.m_uid) {
        // If two nodes have same uid, compare their graphs
        std::less<Graph*> less_than_graph;
        return less_than_graph(m_graph, n.m_graph);
      } else {
        return false;
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
    return m_numNodes;
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
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type value = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type currUid {m_nextUid++};
    ++m_numNodes;
    // Store the position
    m_points.insert(m_points.end(), {currUid, position});
    // Store the value
    m_uidToValue.insert(m_uidToValue.end(), {currUid, value});
    // Assign index to uid
    m_indexToUid.push_back(currUid);
    m_indices.insert(m_indices.end(), {currUid, m_numNodes - 1});
    // Make a node to return
    return Node(this, currUid);  
    // TODO: node_value_type
  }

  /**
   * @brief Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.m_graph == this) {
      return true;
    } else {
      return false;
    }
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
    return Node(this, m_indexToUid[i]);
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
    /**
     * @brief Construct an invalid Edge.
     * @return An invalid Edge
     *
     * Complexity: O(1)
     */
    Edge() {
      // HW0: YOUR CODE HERE
      m_graph = nullptr;
      m_uid = 0;
      m_uidNode1 = 0;
      m_uidNode2 = 0;
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
      return Node(m_graph, m_uidNode1);
    }

    /**
     * @brief Return the second node of this Edge
     * @return The second Node of the Edge
     *
     * Complexity: O(1)
     */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(m_graph, m_uidNode2);      // Invalid Node
    }

    /**
     * @brief Test whether this edge and @a e are equal.
     * Equal edges represent the same undirected edge between two nodes.
     * @return true if this edge and @e are equal, false otherwise
     */
    bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE
      if (m_graph == e.m_graph 
        && ((m_uidNode1 == e.m_uidNode1 && m_uidNode2 == e.m_uidNode2)
          || (m_uidNode1 == e.m_uidNode2 && m_uidNode2 == e.m_uidNode1))) {
        return true;
      } else {
        return false;
      }
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
      if (m_uid < e.m_uid) {
        return true;
      } else if (m_uid == e.m_uid) {
        std::less<Graph*> less_than_graph;
        return less_than_graph(m_graph, e.m_graph);
      } else {
        return false;
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
    Edge(const Graph* pGraph, size_type uid, size_type uidNode1, size_type uidNode2) 
      : m_graph(const_cast<Graph*>(pGraph)), m_uid(uid), m_uidNode1(uidNode1), m_uidNode2(uidNode2) {}

    Graph* m_graph;
    size_type m_uid;
    size_type m_uidNode1;
    size_type m_uidNode2;
  };

  /**
   * @brief Return the total number of edges in the graph.
   * @return The number of unique edges in the graph
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return m_numEdges;
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
    assert(0 <= i && i < m_numEdges);
    size_type uidEdge {m_indexToUidEdge[i]};
    size_type uidNode1 {m_connections.at(uidEdge).first};
    size_type uidNode2 {m_connections.at(uidEdge).second};

    return Edge(this, uidEdge, uidNode1, uidNode2);
  }

  /**
   * @brief Test whether two nodes are connected by an edge.
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * @pre @a a and @a b are valid nodes of this graph
   *
   * Complexity: O(num_nodes())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.m_graph == this && b.m_graph == this);
    if (get_edge_uid(a, b) != size_type(-1)) {
      return true;
    } else {
      return false;
    }
    
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Complexity: O(1)
   */
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    assert(!(a == b) && a.m_graph == this && b.m_graph == this);
    size_type uid_find {get_edge_uid(a, b)};
    if (uid_find != size_type(-1)) {
      return Edge(this, uid_find, a.m_uid, b.m_uid);
    } else {
      size_type currUidEdge {m_nextUidEdge++};
      ++m_numEdges;
      // Assign index to uidEdge
      m_indexToUidEdge.push_back(currUidEdge);
      // Add connection
      m_connections.insert(m_connections.end(), {currUidEdge, {a.m_uid, b.m_uid}});
      m_fromNode[a.m_uid][b.m_uid] = currUidEdge;
      m_fromNode[b.m_uid][a.m_uid] = currUidEdge;
      // Make a new edge to return
      return Edge(this, currUidEdge, a.m_uid, b.m_uid);
    }
  }

  /**
   * @brief Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    m_nextUid = 0;
    m_numNodes = 0;
    m_nextUidEdge = 0;
    m_numEdges = 0;
    m_indexToUid.clear();
    m_indexToUidEdge.clear();
    m_points.clear();
    m_connections.clear();
    m_fromNode.clear();
    m_indices.clear();
    m_uidToValue.clear();
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
      return Node(m_graph, m_uid);
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
      size_type nextIndex = m_graph->m_indices[m_uid] + 1;
      if (nextIndex == m_graph->size()) {
        // Current iterator is at the last node
        m_uid = size_type(-1);
      } else {
        m_uid = m_graph->m_indexToUid[nextIndex];
      }
      return *this;
    }

    /**
     * @brief Overload == operator
     * @return true if the node this iterator is pointing to is the same as the node @nodeIt is pointing to, false otherwise.
     *
     * Complexity: O(1)
     */
    bool operator==(const NodeIterator nodeIt) const {
      return (m_graph == nodeIt.m_graph && m_uid == nodeIt.m_uid);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    /**
     * @brief Private constructor
     * @param[in] graph The graph this iterator belongs to
     * @param[in] uid The node UID this iterator points to
     */
    NodeIterator(const graph_type* graph, size_type uid) : m_graph{const_cast<Graph*>(graph)}, m_uid{uid} {}
    graph_type* m_graph;
    size_type m_uid;
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
  node_iterator node_begin() const {
    return node_iterator(this, m_indexToUid[0]);
  }

  /**
   * @brief Get iterator to one beyond last node.
   * @return Iterator pointing to an invalid node that is one beyond the last node in the graph.
   *
   * Complexity: O(1)
   */
  node_iterator node_end() const {
    return node_iterator(this, size_type(-1));
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
      return Edge(m_graph, m_it->second, m_uid, m_it->first);
    }


//--documentation_-1
//--well done on documentation
//--END
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
      ++m_it;
      return *this;
    }

    /**
     * @brief Overload == operator
     * @return true if this iterator points to the same connecting edge from the same home node as @incidentIt does, false otherwise.
     *
     * Complexity: O(1)
     */
    bool operator==(const IncidentIterator& incidentIt)  const {
      return (m_graph == incidentIt.m_graph && m_uid == incidentIt.m_uid && m_it == incidentIt.m_it);
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
      m_graph{const_cast<graph_type*>(pGraph)}, m_uid{uid}, m_it{it} {}
    graph_type* m_graph;
    size_type m_uid; // UID of starting node
    std::unordered_map<size_type, size_type>::iterator m_it; // Iterator of map in graph storing edges from starting node
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
      return Edge(m_graph, m_it->first, m_it->second.first, m_it->second.second);
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
      ++m_it;
      return *this;
    }

    /**
     * @brief Overload == operator
     * @return true if this iterator is pointing to the same edge of the same graph as @edgeIt does, false otherwise.
     *
     * Complexity: O(1)
     */
    bool operator==(const EdgeIterator& edgeIt) const {
      return (m_graph == edgeIt.m_graph && m_it == edgeIt.m_it);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    /**
     * @brief Private constructor
     * @param[in] pGraph The graph this iterator belongs to
     * @param[in] m_it Iterator to an element in pGraph->m_connections
     *
     * Complexity: O(1)
     */
    EdgeIterator(const Graph* pGraph, std::unordered_map<size_type, std::pair<size_type, size_type>>::const_iterator it) :
        m_graph{const_cast<Graph*>(pGraph)}, m_it{it} {}
    Graph* m_graph;
    // TODO: is const_iterator fine?
    std::unordered_map<size_type, std::pair<size_type, size_type>>::const_iterator m_it; // Iterator of m_connections in Graph
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
    return EdgeIterator(this, m_connections.begin());
  }

  /**
   * @brief Get iterator to one beyond last edge
   * @return Iterator pointing to an invalid edge that is one beyond the last edge of the graph
   *
   * Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, m_connections.end());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  /**
   * @brief Get an edge uid from two nodes if it exists. Otherwise return an invalid edge uid
   * @param[in] a One node of the edge
   * @param[in] b The other node of the edge
   * @return The edge uid that connects nodes @a and @b. If the edge does not exist, return size_type(-1)
   *
   * @pre @a a and @a b are distinct valid nodes of this graph
   * Complexity: O(1)
   */
  size_type get_edge_uid(const Node& a, const Node& b) const {
    assert(!(a == b) && a.m_graph == this && b.m_graph == this);
    std::unordered_map<size_type, std::unordered_map<size_type, size_type>>::const_iterator it {m_fromNode.find(a.m_uid)};
    if (it != m_fromNode.end()) {
      std::unordered_map<size_type, size_type>::const_iterator it1 {(it->second).find(b.m_uid)};
      if (it1 != it->second.end()) {
        return it1->second;
      }
    }
    return size_type(-1);
  }
  std::unordered_map<size_type, Point> m_points;
  std::vector<size_type> m_indexToUid;
  std::unordered_map<size_type, size_type> m_indices;
  size_type m_nextUid;
  size_type m_numNodes;
  std::unordered_map<size_type, node_value_type> m_uidToValue;

  std::vector<size_type> m_indexToUidEdge;
  size_type m_nextUidEdge;
  size_type m_numEdges;
  std::unordered_map<size_type, std::pair<size_type, size_type>> m_connections; // <uid_edge, <uid_node1, uid_node 2>>
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> m_fromNode; // <uid_node1, <uid_node 2, uid_edge>>
};

#endif // CME212_GRAPH_HPP
