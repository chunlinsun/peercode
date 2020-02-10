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
 * @tparam V The type of data a node can hold.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:

  struct internal_node;
  struct internal_edge;

 public:

//--documentation_-1
//--well done again, good doxygen style in documentation and concise code
//--END

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Synonym for the node payload type. */
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_->get_internal_node_by_uid(uid_).point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->get_node_index_by_uid(uid_);
    }

    /**
     * @brief Returns a mutable reference to this node's payload.
     *
     * @return A mutable reference to this node's payload.
     */
    node_value_type& value() {
      return graph_->get_internal_node_by_uid(uid_).value;
    }

    /**
     * @brief Returns a read-only reference to this node's payload.
     *
     * @return A read-only reference to this node's payload.
     */
    const node_value_type& value() const {
      // Note that get_internal_node_by_uid is overloaded, so this is not the same as above.
      return graph_->get_internal_node_by_uid(uid_).value;
    }

    /**
     * @brief Returns the degree of this node.
     *
     * @return The number of nodes directly connected to this one by an edge.
     */
    size_type degree() const {
      return graph_->insert_or_get_adjacency_map_for_node(*this).size();
    }

    /** 
     * @brief Returns an "starting" iterator for this node's incident edges.
     *
     * @return "starting" iterator for all of this node's incident edges. This will be equal to
     *     edge_end() if there are no incident edges.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, graph_->insert_or_get_adjacency_map_for_node(*this).begin());
    }

    /** 
     * @brief Returns an "end" iterator for this node's incident edges.
     *
     * @return "end" iterator for all of this node's incident edges.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, graph_->insert_or_get_adjacency_map_for_node(*this).end());
    }

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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    size_type new_index = num_nodes();
    size_type new_uid = next_node_uid_++;

    internal_nodes_vector_.emplace_back(position, value);

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
  class Edge : private totally_ordered<Edge> {
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

    Edge get_inverted_copy() const {
      return Edge(graph_, uid_, !smaller_node_first_);
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
    return node_pair_to_edge_index_map_.count(a) == 1
        && node_pair_to_edge_index_map_.at(a).count(b) == 1;
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

    // First, we see whether an equivalent edge already exists.
    if (has_edge(a, b)) {
      size_type existing_index = node_pair_to_edge_index_map_.at(a).at(b);
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

    const Node& smaller_node = (smaller_node_first_for_desired_edge) ? a : b;
    const Node& larger_node = (smaller_node_first_for_desired_edge) ? b : a;
    internal_edges_vector_.emplace_back(smaller_node, larger_node);

    // Since the graph is undirected, we can safely keep track of the edge in both adjacency maps.
    insert_or_get_adjacency_map_for_node(a).insert({b, new_index});
    insert_or_get_adjacency_map_for_node(b).insert({a, new_index});

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
    node_pair_to_edge_index_map_.clear();
    edge_uid_to_vector_index_map_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. An input iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /**
     * @brief Dereferences this iterator.
     *
     * @return A copy of the node this iterator is pointing to.
     *
     * @pre this is not equal to an "end" iterator (i.e., this != the associated graph's node_end())
     *
     * NOTE: This can actually safely return a const reference to the node, but I'm intentionally
     * not changing the signature given in the starter code.
     */
    Node operator*() const {
      return *vector_iterator_;
    }

    /**
     * @brief Advances this iterator by one step.
     *
     * @return A reference to this iterator, for chaining.
     *
     * @pre this is not equal to an "end" iterator (i.e., this != the associated graph's node_end())
     */
    NodeIterator& operator++() {
      ++vector_iterator_;
      return *this;
    }

    /**
     * @brief Checks whether this iterator is equal to _other_.
     *
     * @param[in] other The other iterator to check equality against.
     * @return true if and only if this iterator is equal to _other_.
     */
    bool operator==(const NodeIterator& other) const {
      return vector_iterator_ == other.vector_iterator_;
    }

   private:
    friend class Graph;
    
    typename std::vector<Node>::const_iterator vector_iterator_;

    NodeIterator(typename std::vector<Node>::const_iterator vector_iterator)
        : vector_iterator_(vector_iterator) {}
  };

  /** 
   * @brief Returns an "starting" iterator for this graph's ndes.
   *
   * @return "starting" iterator for all of this graph's nodes. This will be equal to node_end() if
   *     the graph is empty.
   */
  node_iterator node_begin() const {
    return NodeIterator(nodes_vector_.begin());
  }

  /** 
   * @brief Returns an "end" iterator for this graph's ndes.
   *
   * @return "end" iterator for all of this graph's nodes.
   */
  node_iterator node_end() const {
    return NodeIterator(nodes_vector_.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. An input iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /**
     * @brief Dereferences this iterator.
     *
     * @return A copy of the edge this iterator is pointing to, where node1() is guaranteed to be
     *     the node that spawned this iterator, and node2() is the neighboring node.
     *
     * @pre this is not equal to an "end" iterator (i.e., this != the associated node's edge_end())
     */
    Edge operator*() const {
      Edge edge = graph_->edges_vector_.at(map_iterator_->second);
      return (edge.node1() == map_iterator_->first) ? edge.get_inverted_copy() : edge;
    }

    /**
     * @brief Advances this iterator by one step.
     *
     * @return A reference to this iterator, for chaining.
     *
     * @pre this is not equal to an "end" iterator (i.e., this != the associated node's edge_end())
     */
    IncidentIterator& operator++() {
      ++map_iterator_;
      return *this;
    }

    /**
     * @brief Checks whether this iterator is equal to _other_.
     *
     * @param[in] other The other iterator to check equality against.
     * @return true if and only if this iterator is equal to _other_.
     */
    bool operator==(const IncidentIterator& other) const {
      return graph_ == other.graph_ && map_iterator_ == other.map_iterator_;
    }

   private:
    friend class Graph;
    
    Graph* graph_;
    typename std::map<Node, size_type>::const_iterator map_iterator_;

    IncidentIterator(Graph* graph, typename std::map<Node, size_type>::const_iterator map_iterator)
        : graph_(graph), map_iterator_(map_iterator) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. An input iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /**
     * @brief Dereferences this iterator.
     *
     * @return A copy of the edge this iterator is pointing to.
     *
     * @pre this is not equal to an "end" iterator (i.e., this != the associated graph's edge_end())
     *
     * NOTE: This can actually safely return a const reference to the edge, but I'm intentionally
     * not changing the signature given in the starter code.
     */
    Edge operator*() const {
      return *vector_iterator_;
    }

    /**
     * @brief Advances this iterator by one step.
     *
     * @return A reference to this iterator, for chaining.
     *
     * @pre this is not equal to an "end" iterator (i.e., this != the associated graph's edge_end())
     */
    EdgeIterator& operator++() {
      ++vector_iterator_;
      return *this;
    }

    /**
     * @brief Checks whether this iterator is equal to _other_.
     *
     * @param[in] other The other iterator to check equality against.
     * @return true if and only if this iterator is equal to _other_.
     */
    bool operator==(const EdgeIterator& other) const {
      return vector_iterator_ == other.vector_iterator_;
    }

   private:
    friend class Graph;
    
    typename std::vector<Edge>::const_iterator vector_iterator_;

    EdgeIterator(typename std::vector<Edge>::const_iterator vector_iterator)
        : vector_iterator_(vector_iterator) {}
  };

  /** 
   * @brief Returns a starting iterator for this graph's edges.
   *
   * @return Starting iterator for all of this graph's edges. This will be equal to edge_end() if
   *     the graph has no edges.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(edges_vector_.begin());
  }

  /** 
   * @brief Returns an "end" iterator for this graph's edges.
   *
   * @return Ending iterator for all of this graph's edges.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(edges_vector_.end());
  }

 private:
  // Internal representation of a node's data.
  // Although this could just use a Point instead of a struct, this abstraction makes
  // adding additional properties to nodes much easier. Style-wise, this is also more
  // consistent with the edge representation below.
  // This intentionally doesn't include the node's index because internally, that's
  // stored in the UID-to-index mapping already--see below.
  struct internal_node {
    Point point;
    node_value_type value;

    internal_node(const Point& point, const node_value_type& value) : point(point), value(value) {}
  };

  // Internal representation of an edge's data.
  // Since all edges are undirected, this underlying representation doesn't care
  // about the order of the nodes, which is a property only of the Edge instances
  // that act as "views" over this struct. See the Edge class for details.
  struct internal_edge {
    Node smaller_node;
    Node larger_node;

    internal_edge(const Node& smaller_node, const Node& larger_node)
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
  // This map allows for fast lookups by nodes. Each node is mapped to an "adjacency map," which in
  // turn maps adjacent nodes to (public) edge indexes corresponding to the above vectors.
  // Abstractly, this map represents the connectedness of the graph as follows:
  // - The adjacency map for node A doesn't exist or is empty if and only if node A has no incident
  //   edges.
  // - The adjacency map for node A contains an entry for node B if and only if A and B are
  //   connected by an edge.
  // Notably, this is less efficient (time-wise) than using std::unordered_map since the underlying
  // tree causes O(log(n)) lookup time. In practice, the better option would be an unordered map
  // which entails defining stable hash functions nodes but would give us O(1) lookups. (Another
  // option would be using the node UIDs as keys instead.)
  std::map<Node, std::map<Node, size_type>> node_pair_to_edge_index_map_;

  /** Readable convenience method for looking up node indexes. */
  size_type get_node_index_by_uid(size_type uid) const {
    return node_uid_to_vector_index_map_.at(uid);
  }

  /** Readable convenience method for accessing internal nodes. */
  internal_node& get_internal_node_by_uid(size_type uid) {
    return internal_nodes_vector_.at(get_node_index_by_uid(uid));
  }

  /** Read-only (i.e., const) version of the above method. */
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

  /** Returns the adjacency map for a node.
   *
   * If the adjacency map doesn't exist yet (which can happen if the node has no incident edges),
   * this creates an empty map first.
   */
  std::map<Node, size_type>& insert_or_get_adjacency_map_for_node(const Node& n) {
    if (node_pair_to_edge_index_map_.count(n) == 0) {
      node_pair_to_edge_index_map_.insert({n, std::map<Node, size_type>()});
    }
    return node_pair_to_edge_index_map_.at(n);
  }
};

#endif // CME212_GRAPH_HPP
