#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <functional>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

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

  struct internal_node;
  struct internal_edge;

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
      Assumed to be hashable since it's used as a key of unordered maps.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {}

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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_->get_internal_node_by_uid(uid_).point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->get_node_index_by_uid(uid_);
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
      return graph_ == n.graph_ && index() == n.index();
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
      size_type this_index = index();
      size_type other_index = n.index();
      // Edge case: If the graphs are different but the indices are the same,
      // we have to break the tie by comparing the graphs pointers.
      if (graph_ != n.graph_ && this_index == other_index) {
        return std::less<const Graph*>()(graph_, n.graph_);
      }
      return this_index < other_index;
    }

    /** Test whether this node is less than or equal to @a n.
     *
     * See operator< and operator== for details.
     */
    bool operator<=(const Node& n) const {
      return (*this < n) || (*this == n);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type uid_;
    
    /** Private constructor intended for use by the Graph class. */
    Node(Graph* graph, size_type uid) : graph_(graph), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_vector_.size();
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
    size_type new_index = num_nodes();
    size_type new_uid = next_node_uid_++;

    internal_nodes_vector_.emplace_back(position);

    Node new_node {this, new_uid};
    nodes_vector_.push_back(new_node);

    // Keep track of the UID's (current) internal index so that we get fast
    // lookups later.
    node_uid_to_vector_index_map_.insert({new_uid, new_index});

    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return nodes_vector_.at(i);
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      const internal_edge& internal_struct = get_internal_edge();
      return (smaller_node_first_) ? internal_struct.smaller_node : internal_struct.larger_node;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      const internal_edge& internal_struct = get_internal_edge();
      return (smaller_node_first_) ? internal_struct.larger_node : internal_struct.smaller_node;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return get_ordered_node_pair() == e.get_ordered_node_pair();
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return get_ordered_node_pair() < e.get_ordered_node_pair();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type uid_;
    // This Edge is really just a view over the actual undirected edge in the graph.
    // The Graph doesn't care what the order is, but the user who created the Edge might
    // expect that the connected nodes are still in the same order, so this Edge is used
    // as a "view" on top of the underlying undirected edge. Internally, the edge's nodes
    // are always stored in order (smaller first), but this Edge has to keep track of
    // what the user's "desired" order was.
    bool smaller_node_first_;

    /** Private constructor intended for use by the Graph class. */
    Edge(Graph* graph, size_type uid, bool smaller_node_first)
      : graph_(graph), uid_(uid), smaller_node_first_(smaller_node_first) {}

    /** Retrieves the internal struct stored in the Graph. */
    const internal_edge& get_internal_edge() const {
      return graph_->get_internal_edge_by_uid(uid_);
    }

    /** Builds an ordered pair of the connected nodes, smaller one first. */
    std::pair<Node, Node> get_ordered_node_pair() const {
      return (smaller_node_first_)
          ? std::pair<Node, Node>(node1(), node2())
          : std::pair<Node, Node>(node2(), node1());
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_vector_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges_vector_.at(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::pair<Node, Node> ordered_pair = build_ordered_node_pair(a, b);
    return has_edge_for_ordered_pair(ordered_pair);
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
    bool smaller_node_first_for_desired_edge = a <= b;
    std::pair<Node, Node> ordered_pair = build_ordered_node_pair(a, b);

    // First, we see whether an equivalent edge already exists.
    if (has_edge_for_ordered_pair(ordered_pair)) {
      size_type existing_index = ordered_node_pair_to_edge_index_map_.at(ordered_pair);
      Edge existing_edge = edges_vector_.at(existing_index);

      // Check if the existing edge is already in the right order, returning it if so.
      if (existing_edge.smaller_node_first_ == smaller_node_first_for_desired_edge) {
        return existing_edge;
      }

      // Otherwise, we have to replace that edge with a new one that points to the same data but has its
      // "view" of the internal edge inverted (i.e., with the nodes in the opposite order). See Edge class.
      // Notably, the existing edge is NOT invalidated, as its UID is not removed from the internal
      // mapping, and the underlying (i.e., internal) edge struct isn't being modified.
      Edge flipped_edge = build_and_track_public_edge(existing_index, smaller_node_first_for_desired_edge);
      edges_vector_[existing_index] = flipped_edge;

      return flipped_edge;
    }

    // If we get here, we're adding a brand new edge.
    size_type new_index = num_edges();

    internal_edges_vector_.emplace_back(ordered_pair.first, ordered_pair.second);

    // We need to keep track of where we're adding this edge for fast lookups later.
    ordered_node_pair_to_edge_index_map_.insert({ordered_pair, new_index});

    Edge new_edge = build_and_track_public_edge(new_index, smaller_node_first_for_desired_edge);
    edges_vector_.push_back(new_edge);

    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_vector_.clear();
    internal_nodes_vector_.clear();
    // Implementation note: Existing nodes will (intentionally) be invalidated because their
    // UIDs will no longer exist in this map, causing subsequent lookups to blow up.
    // (The same holds for edges--see below.)
    node_uid_to_vector_index_map_.clear();

    edges_vector_.clear();
    internal_edges_vector_.clear();
    ordered_node_pair_to_edge_index_map_.clear();
    edge_uid_to_vector_index_map_.clear();
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
  // Internal representation of a node's data.
  // Although this could just use a Point instead of a struct, this abstraction makes
  // adding additional properties to nodes much easier. Style-wise, this is also more
  // consistent with the edge representation below.
  // This intentionally doesn't include the node's index because internally, that's
  // stored in the UID-to-index mapping already--see below.
  struct internal_node {
    Point point;

    internal_node(Point point) : point(point) {}
  };

  // Internal representation of an edge's data.
  // Since all edges are undirected, this underlying representation doesn't care
  // about the order of the nodes, which is a property only of the Edge instances
  // that act as "views" over this struct. See the Edge class for details.
  struct internal_edge {
    Node smaller_node;
    Node larger_node;

    internal_edge(Node smaller_node, Node larger_node)
        : smaller_node(smaller_node), larger_node(larger_node) {}
  };

  size_type next_node_uid_ = 0;
  size_type next_edge_uid_ = 0;

  // The nodes (and edges) are stored in vectors since indexes are (for now) always
  // contiguous and we need O(1) random access.
  std::vector<Node> nodes_vector_;
  std::vector<internal_node> internal_nodes_vector_;
  // Because a node's index may *theoretically* be changed by the Graph at any
  // point (and the existing Node might not be told about it), this map explicitly
  // keeps track of the (current) index for each node's UID.
  std::unordered_map<size_type, size_type> node_uid_to_vector_index_map_;

  // These edge vectors and UID mapping work identically to their node equivalents.
  std::vector<Edge> edges_vector_;
  std::vector<internal_edge> internal_edges_vector_;
  std::unordered_map<size_type, size_type> edge_uid_to_vector_index_map_;
  // This map allows for fast lookups by nodes. The keys are ordered pairs of nodes,
  // smaller first. The value for a pair of nodes (assuming one exists) is the
  // corresponding edge's index in the two vectors. If a valid (i.e., ordered) pair
  // of nodes doesn't have an entry in this map, there's no edge connecting them.
  // Notably, this is actually less efficient (time-wise) than std::unordered_map
  // because the underlying tree causes O(log(n)) lookup time. In practice, the
  // better option would be an unordered map, which entails defining stable hash
  // functions for these node pairs but would give us O(1) lookups at the cost of
  // memory, as we don't really care about the ordering.
  // WARNING: This (and the rest of the class) generally assumes that two distinct
  // Nodes will never be equal (e.g., two nodes with the same index and graph but
  // different positions).
  std::map<std::pair<Node, Node>, size_type> ordered_node_pair_to_edge_index_map_;

  /** Readable convenience method for looking up node indexes. */
  size_type get_node_index_by_uid(size_type uid) const {
    return node_uid_to_vector_index_map_.at(uid);
  }

  /** Readable convenience method for accessing internal nodes. */
  const internal_node& get_internal_node_by_uid(size_type uid) const {
    return internal_nodes_vector_.at(get_node_index_by_uid(uid));
  }

  /** Readable convenience method for looking up edge indexes. */
  size_type get_edge_index_by_uid(size_type uid) const {
    return edge_uid_to_vector_index_map_.at(uid);
  }

  /** Readable convenience method for accessing internal edges. */
  const internal_edge& get_internal_edge_by_uid(size_type uid) const {
    return internal_edges_vector_.at(get_edge_index_by_uid(uid));
  }

  /** Returns whether an edge exists between an ordered pair of nodes (smaller first). */
  bool has_edge_for_ordered_pair(std::pair<Node, Node> ordered_node_pair) const {
    return ordered_node_pair_to_edge_index_map_.count(ordered_node_pair) == 1;
  }

  /** Helper method for creating a new public Edge, tracked by UID.
   *
   * This keeps the internal UID-to-index mapping consistent but does not store anything
   * in the internal vectors (neither the new public edge nor a new internal edge).
   */
  Edge build_and_track_public_edge(size_type index, bool smaller_node_first) {
    size_type new_uid = next_edge_uid_++;
    edge_uid_to_vector_index_map_.insert({new_uid, index});

    return Edge(this, new_uid, smaller_node_first);
  }

  /** Given two nodes, constructs an ordered pair where the smaller node comes first. */
  static std::pair<Node, Node> build_ordered_node_pair(const Node& a, const Node& b) {
    return (a <= b) ? std::pair<Node, Node>(a, b) : std::pair<Node, Node>(b, a);
  }
};

#endif // CME212_GRAPH_HPP
