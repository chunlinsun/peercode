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

  using node_value_type = V;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    // check the size of Node and Edge as required by HW0
    assert(sizeof(Node) <= 16);
    assert(sizeof(Edge) <= 32);
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
  class Node: private totally_ordered<Node> {
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_list[id_].second;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    /**
    * @return the value associated with the Node.
    */
    node_value_type& value() {
      return graph_->nodes_list[id_].first;
    }
    /**
    * @return the const value associated with this Node.
    */
    const node_value_type& value() const {
      return graph_->nodes_list[id_].first;
    }

    //--documentation_0
    //--This postcondition really isn't necessary, since degree() is a const
    //--method. (Nothing's being modified. The assertion isn't necessary either.)
    //--START
    /**
    * @return the degree/number of incident edges associated with this Node.
    * @post the degree should be < num_nodes();
    */
    size_type degree() const {
      assert(graph_->adj_list.at(id_).size() < num_nodes());
      return graph_->adj_list.at(id_).size();
    }
    //--END
    /**
    * @return the iterator that points to the 1st incident edge.
    */
    IncidentIterator edge_begin() const {
      return IncidentIterator(graph_, id_, graph_->adj_list.at(id_).begin());
    }
    /**
    * @return the iterator that points to the last incident edge.
    */
    IncidentIterator edge_end() const {
      //--functionality_1
      //--You cannot increment the end iterator! This will likely cause a
      //--segfault in your IncidentIterator.
      //--START
      return IncidentIterator(graph_, id_, ++(graph_->adj_list.at(id_).end()));
      //--END
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_ == n.graph_) && (id_ == n.id_));
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
        return (id_ < n.id_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(const Graph* graph, size_type id)
      : graph_(const_cast<Graph*>(graph)), id_(id) {
    }
    Graph* graph_;
    size_type id_;
    node_value_type node_value_;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_list.size();
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
    nodes_list.push_back({nodes_list.size(),position});
    return Node(this, num_nodes()-1); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ((this == n.graph_) && (nodes_list[n.id_].first == n.id_));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    // HW0: YOUR CODE HERE
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_id_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_id_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((graph_ == e.graph_) && (node1_id_ == e.node1_id_) && (node2_id_ == e.node2_id_)) || \
             ((graph_ == e.graph_) && (node2_id_ == e.node1_id_) && (node1_id_ == e.node2_id_));

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning. 
     */
    bool operator<(const Edge& e) const {
      //--style_0
      //--No need for a ternary operator here. This is equivalent to
      //--  return id_ < e.id_;
      //--START
      return (id_ < e.id_) ? true : false;
      //--END
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* graph, size_type id, size_type node1_id, size_type node2_id)
      : graph_(const_cast<Graph*>(graph)), id_(id), node1_id_(node1_id), node2_id_(node2_id) {
    }
    Graph* graph_;
    size_type id_;
    size_type node1_id_;
    size_type node2_id_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    // HW0: YOUR CODE HERE
    return Edge(this, i, edges_list[i].first, edges_list[i].second); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    assert(a.graph_ == b.graph_);

    return (adj_list[a.index()].count(b.index()) == 1);
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
    assert(a.graph_ == b.graph_);
    assert(a.id_ != b.id_);

    size_type i;
    if (has_edge(a, b)) {
      i = adj_list.at(a.index()).at(b.index()); //return the edge index if found
    }
    else {
      i = edges_list.size(); // get the new edge index
      edges_list.push_back({a.index(), b.index()}); // insert the new edge
      adj_list[b.index()][a.index()] = i; // store new edge index at adj_list[a.index][b.index]
      adj_list[a.index()][b.index()] = i; // store new edge index at adj_list[b.index][a.index] (redundancy for performance purposes)
    }
    return Edge(this, i, edges_list[i].first, edges_list[i].second);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_list.clear();
    edges_list.clear();
    adj_list.clear();
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
    /**
    * @return the dereferenced node.
    */
    Node operator*() const {
      return Node(graph_, (*adj_list_iter_).first);
    }
    /**
    * @return the next iterator in the sequence.
    */
    NodeIterator& operator++() {
      ++adj_list_iter_;
      return *this;
    }
    /**
    * @return true if _iter_ equals to this iterator.
    *
    * Equal iterators have the same members, otherwise, they are unequal.
    */
    bool operator==(const NodeIterator& iter) const {
      return ((graph_ == iter.graph_) && (adj_list_iter_ == iter.adj_list_iter_));
    }
    /**
    * @return true if _iter_ does not equal to this iterator.
    *
    * Equal iterators have the same members, otherwise, they are unequal.
    */
    bool operator!=(const NodeIterator& iter) const {
      return ((graph_ != iter.graph_) || (adj_list_iter_ != iter.adj_list_iter_));
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    //--design_0
    //--Why not iterate over the nodes_list instead of the adjacency list? This seems like a case of
    //--passing in too much information.
    //--START
    std::map<size_type, std::map<size_type, size_type>>::const_iterator adj_list_iter_;
    //--END
    NodeIterator(const Graph* graph, std::map<size_type, std::map<size_type, size_type>>::const_iterator adj_list_iter)
      : graph_(const_cast<Graph*>(graph)), adj_list_iter_(adj_list_iter) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /**
  * @return NodeIterator that points to the beginning of the sequence.
  */
  NodeIterator node_begin() const {
    return NodeIterator(this, adj_list.begin());
  }

  /**
  * @return NodeIterator that points to the end of the sequence.
  */
  NodeIterator node_end() const {
    return NodeIterator(this, adj_list.end());
  }
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

    /**
    * @return dereferenced edge.
    */
    Edge operator*() const {
      size_type node2_id_ = adj_list_iter2_->first;
      return Edge(graph_, graph_->adj_list.at(node1_id_).at(node2_id_), node1_id_, node2_id_);
    }

    /**
    * @return the next iterator in the sequence.
    */
    IncidentIterator& operator++() {
      ++adj_list_iter2_;
      return *this;
    }

    /**
    * @return true if _iter_ equals to this iterator.
    *
    * Equal iterators have the same members, otherwise, they are unequal.
    */
    bool operator==(const IncidentIterator& iter) const {
      return ((graph_ == iter.graph_) && (node1_id_ == iter.node1_id_) && (adj_list_iter2_->first == iter.adj_list_iter2_->first));
    }

    /**
    * @return true if _iter_ does not equal to this iterator.
    *
    * Equal iterators have the same members, otherwise, they are unequal.
    */
    bool operator!=(const IncidentIterator& iter) const {
      return ((graph_ != iter.graph_) || (node1_id_ != iter.node1_id_) || (adj_list_iter2_->first != iter.adj_list_iter2_->first));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node1_id_;
    std::map<size_type, size_type>::const_iterator adj_list_iter2_;
    IncidentIterator(const Graph* graph, size_type node1_id, std::map<size_type, size_type>::const_iterator adj_list_iter2)
      : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), adj_list_iter2_(adj_list_iter2) {
    }
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

    /**
    * @return dereferenced edge.
    */
    Edge operator*() const {
      //--design_0
      //--Generally speaking, you shouldn't do math with iterators like this. Since
      //--your edge list is a vector, you might consider just keeping track of
      //--the index as a private member rather than a const_iterator.
      //--START
      return Edge(graph_, edges_list_iter_ - graph_->edges_list.begin(), edges_list_iter_->first, edges_list_iter_->second);
      //--END
    }

    /**
    * @return the next iterator in the sequence.
    */
    EdgeIterator& operator++() {
      ++edges_list_iter_;
      return *this;
    }

    /**
    * @return true if _iter_ equals to this iterator.
    *
    * Equal iterators have the same members, otherwise, they are unequal.
    */
    bool operator==(const EdgeIterator& iter) const {
      return ((graph_ == iter.graph_) && (edges_list_iter_ == iter.edges_list_iter_));
    }

    /**
    * @return true if _iter_ does not to this iterator.
    *
    * Equal iterators have the same members, otherwise, they are unequal.
    */
    bool operator!=(const EdgeIterator& iter) const {
      return ((graph_ != iter.graph_) || (edges_list_iter_ != iter.edges_list_iter_));
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    std::vector<std::pair<size_type, size_type>>::const_iterator edges_list_iter_;
    EdgeIterator(const Graph* graph, std::vector<std::pair<size_type, size_type>>::const_iterator edges_list_iter)
      : graph_(const_cast<Graph*>(graph)), edges_list_iter_(edges_list_iter) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
  * @return EdgeIterator that points to the beginning of the sequence.
  */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, edges_list.begin());
  }

  /**
  * @return EdgeIterator that points to the beginning of the sequence.
  */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, edges_list.end());
  }
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  /* nodes_list: vector that store the nodes information:
  *  nodes_list[node_index].first == node.value()
  *  nodes_list[node_index].second == node.position()
  **/
  //--design_0
  //--I'd really recommend designing a struct instead of a std::pair to hold the
  //--node internals (value and position), so that you can do things like
  //--  Point pos = nodes_list[i].position;
  //--  node_value_type val = nodes_list[i].value;
  //--instead of having to remember which one is first or second.
  //--START
  std::vector<std::pair<node_value_type, Point>> nodes_list;
  //--END

  /* edges_list: vector that store the edges information:
  *  edges_list[edge_index].first == node1().index()
  *  edges_list[edge_index].second == node2().index()
  **/
  std::vector<std::pair<size_type, size_type>> edges_list;

  /* adj_list: vector that store the adjacency of the list:
  *  adj_list[node1][node2] == edge.index()
  **/
  //--design_0
  //--Good! Look up std::unordered_map, too.
  //--START
  std::map<size_type, std::map<size_type, size_type>> adj_list;
  //--END
};


#endif // CME212_GRAPH_HPP
