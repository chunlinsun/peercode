#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <tuple>
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
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;

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
    // initialize empty vectors and things
    std::vector<internal_node> v;
    nodes = v;
    next_node_uid = 0;
    std::vector<std::tuple<size_type, size_type>> e;
    edges = e; 
    std::unordered_map<size_type, bool> us;
    uid_status = us;
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
      graph = nullptr;
      uid = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph->get_point_by_uid(uid);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph->get_index_by_uid(uid);
    }

    // for debug / testing
    // size_type get_uid() {
    //   return uid;
    // }

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
      return (uid == n.uid) && (graph == n.graph);
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
      return (uid < n.uid);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // data members
    // Pointer to the Graph
    Graph* graph;
    // uid of the Node (corresponds to internal_node.uid in graph)
    size_type uid;

    // constructor of valid nodes
    Node(Graph* g, size_type id) {
      graph = g;
      uid = id;
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes.size();
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
    // add new node to internal graph
    size_type uid = next_node_uid;
    internal_node inode {position, uid};
    nodes.push_back(inode);
    uid_status[uid] = true;
    next_node_uid = next_node_uid + 1;
    // construct accessor to this node as Node object and return it
    Node n = Node(this, uid);
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph == this) && has_uid(n.uid);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // retrive data for index i node
    // construct Node
    return Node(const_cast<Graph*>(this), nodes[i].uid);
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
      graph = nullptr;
      std::tuple<size_type, size_type> d {0,0};
      uids = d;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph, std::get<0>(uids));     
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph, std::get<1>(uids)); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE 
      // test: same graph, same uids
      bool same_graph = (graph == e.graph);
      bool uid_match = (std::get<0>(uids) == std::get<0>(e.uids)) && (std::get<1>(uids) == std::get<1>(e.uids));
      bool uid_reverse_match = (std::get<0>(uids) == std::get<1>(e.uids)) && (std::get<1>(uids) == std::get<0>(e.uids));
      bool same_uids = uid_match || uid_reverse_match;
      return same_graph && same_uids;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // compare index in graph's vector
      // get index of this edge and index of edge e
      size_type this_e_ind;
      for (size_type i=0; i < graph->num_edges(); ++i){
        bool n1 = std::get<0>(graph->edges[i]) == std::get<0>(uids);
        bool n2 = std::get<1>(graph->edges[i]) == std::get<1>(uids);
        if (n1 && n2) {
          this_e_ind = i;
        }
      }

      size_type other_e_ind;
      for (size_type j=0; j < e.graph->num_edges(); ++j) {
        bool n3 = std::get<0>(e.graph->edges[j]) == std::get<0>(e.uids);
        bool n4 = std::get<1>(e.graph->edges[j]) == std::get<1>(e.uids);
        if (n3 && n4) {
          other_e_ind = j;
        }
      }

      return this_e_ind < other_e_ind;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    
    // data members of Edge
    Graph* graph;
    std::tuple<size_type, size_type> uids;
  
    // constructor of valid Edges
    Edge(Graph* g, std::tuple<size_type, size_type> ids){
      graph = g;
      uids = ids;
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(const_cast<Graph*>(this), edges[i]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // iterate through vector of edges, see if either (a,b)
    // is in the vector, or (b,a) (handled by == for an edge)
    for (size_type i = 0; i < edges.size(); ++i) {
      Edge e_tmp = Edge(const_cast<Graph*>(this), std::make_tuple(a.uid, b.uid));
      Edge existing_e = Edge(const_cast<Graph*>(this), edges[i]);
      if (e_tmp == existing_e) {
        return true;
      }  
    }
    return false;
  }

  /** get edge if present. otherwise return invalid edge.
   */
  Edge get_edge(const Node& a, const Node& b) const {
    for (size_type i = 0; i < edges.size(); ++i) {
      Edge e_tmp = Edge(const_cast<Graph*>(this), std::make_tuple(a.uid, b.uid));
      Edge existing_e = Edge(const_cast<Graph*>(this), edges[i]);
      if (e_tmp == existing_e) {
        return existing_e;
      }  
    }
    return Edge(); // invalid edge
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
    // HW0: YOUR CODE HERE
    Edge e = Edge();
    Edge existing = get_edge(a,b);
    if (existing.graph) { // if not null, return existing edge
      e = existing;
    } else { // if null, add edge
      // push edge to list in graph
      edges.push_back(std::make_tuple(a.uid, b.uid));
      // return edge object
      e = Edge(this, std::make_tuple(a.uid, b.uid));
    }
    return e;  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // clear nodes
    std::vector<internal_node> new_nodes;
    nodes = new_nodes;
    std::unordered_map<size_type, bool> new_uid_status;
    uid_status  = new_uid_status;
    //next_node_uid = 0; < not strictly necessary
    // clear edges
    std::vector<std::tuple<size_type, size_type>> new_edges;
    edges = new_edges;
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
  struct internal_node {
      Point loc; // location of node
      size_type uid; // unique identifier of node
  };
  std::vector<internal_node> nodes;
  size_type next_node_uid;
  std::vector<std::tuple<size_type, size_type>> edges;
  // structure to check which uids are in graph, quickly
  std::unordered_map<size_type, bool> uid_status;

  Point& get_point_by_uid(size_type uid) {
    // iterate through vector named nodes
    // if uid matches, return reference to Point loc
    for (size_type i = 0; i < nodes.size(); ++i) {
      if (uid == nodes[i].uid) {
        return nodes[i].loc;
      }
    }
  }

  size_type get_index_by_uid(size_type uid) {
    // iterate through vector named nodes
    // if uid matches, return index
    for (size_type i = 0; i < nodes.size(); ++i) {
      if (uid == nodes[i].uid) {
        return i;
      }
    }
  }

  // Assumes that the pointer to the graph has been checked to match already
  bool has_uid(size_type node_uid) const {
    // O(num nodes) :
    // for (size_type i = 0; i < nodes.size(); ++i) {
    //     if (node_uid == nodes[i].uid) {
    //       return true;
    //     }
    // }
    // return false;
    // }

    // O(1) with a hashmap:
    if (uid_status.find(node_uid) == uid_status.end()){
      return false;
    } else {
      return true;
    }
  }

// for debug:
// public:
//   bool has_uidd(size_type node_uid) const {
//     return has_uid(node_uid);
//   }

};


#endif // CME212_GRAPH_HPP
