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
    Node() {
      // HW0: YOUR CODE HERE
      m_graph = nullptr;
      m_uid = size_type(-1);
    }

    /** Return this node's position. 
     * @pre: This node is valid and belongs to a graph
     * Complexity: O(log(num_nodes())
     */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return m_graph->m_points[m_uid];
    }

    /** Return this node's index, a number in the range [0, graph_size). 
     * Complexity: O(log(num_nodes())
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

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (m_graph == n.m_graph && m_uid == n.m_uid) {
        return true;
      } else {
        return false;
      }
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

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return m_numNodes;
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    size_type currUid {m_nextUid++};
    ++m_numNodes;
    // Store the position
    m_points.insert(m_points.end(), {currUid, position});
    // Assign index to uid
    m_indexToUid.push_back(currUid);
    m_indices.insert(m_indices.end(), {currUid, m_numNodes - 1});
    // Make a node to return
    return Node(this, currUid);  
  }

  /** Determine if a Node belongs to this Graph
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

  /** Return the node with index @a i.
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      m_graph = nullptr;
      m_uid = 0;
      m_uidNode1 = 0;
      m_uidNode2 = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(m_graph, m_uidNode1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(m_graph, m_uidNode2);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
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

    /** Test whether this edge is less than @a e in a global order.
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

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return m_numEdges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Complexity: O(log(num_edges()))
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < m_numEdges);
    size_type uidEdge {m_indexToUidEdge[i]};
    size_type uidNode1 {m_connections.at(uidEdge).first};
    size_type uidNode2 {m_connections.at(uidEdge).second};

    return Edge(this, uidEdge, uidNode1, uidNode2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
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
   * Complexity: O(log(num_nodes()))
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

  /** Remove all nodes and edges from this graph.
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

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  /* Get an edge uid from two nodes if it exists. Otherwise return an invalid edge uid 
   * @pre @a a and @a b are distinct valid nodes of this graph
   * Complexity: O(log(num_node()))
   */
  size_type get_edge_uid(const Node& a, const Node& b) const {
    assert(!(a == b) && a.m_graph == this && b.m_graph == this);
    std::map<size_type, std::map<size_type, size_type>>::const_iterator it {m_fromNode.find(a.m_uid)};
    if (it != m_fromNode.end()) {
      std::map<size_type, size_type>::const_iterator it1 {(it->second).find(b.m_uid)};
      if (it1 != it->second.end()) {
        return it1->second;
      }
    }
    return size_type(-1);

  }
  std::map<size_type, Point> m_points;
  std::vector<size_type> m_indexToUid;
  std::map<size_type, size_type> m_indices;
  size_type m_nextUid;
  size_type m_numNodes;
  std::vector<size_type> m_indexToUidEdge;
  size_type m_nextUidEdge;
  size_type m_numEdges;
  std::map<size_type, std::pair<size_type, size_type>> m_connections; // <uid_edge, <uid_node1, uid_node 2>>
  std::map<size_type, std::map<size_type, size_type>> m_fromNode; // <uid_node1, <uid_node 2, uid_edge>>
};

#endif // CME212_GRAPH_HPP
