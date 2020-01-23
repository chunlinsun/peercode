#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>
#include <unordered_map>
#include <iostream>
#include <functional>

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
    // HW0: Constructor, vectors are default initialized to size 0
    num_nodes_ = internal_nodes_.size();
    num_edges_ = internal_edges_.size();
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
      // HW0: Construct invalid node with public constructor
      // Do nothing on purpose, this node is invalid
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: return the Point located at this Node's index within
      // this Node's associated Graph
      return graph_->internal_nodes_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: return the value stored as index_
      return index_;
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
      // HW0: if both do not point to same Graph or have different 
      // indices, return false, otherwise return true
      if ((this->graph_ != n.graph_) or (this->index_ != n.index_))
        return false;
      else
        return true;
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
      // HW0: If both Nodes point to same graph, then just compare index
      if (this->graph_ == n.graph_){
        if (this->index_ < n.index_)
          return true;
        else
          return false;
      }

      // HW0: If Nodes point to different graphs, use ordering on pointers
      // to Graph objects
      std::less<graph_type*> check;
      if (check(this->graph_,n.graph_))
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // HW0: Initialize member variables: pointer to Graph and index 
    // telling us where in Graph internal node data to get Point
    graph_type* graph_;
    size_type index_;

    // HW0: Private constructor for valid Nodes, given Graph and index
    Node(const graph_type* graph, size_type ind)
      : graph_(const_cast<graph_type*>(graph)), index_(ind) {}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: return value stored in num_nodes member variable
    return num_nodes_;
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
    // HW0: add the point to the end of the internal_nodes_ vector
    // then increment num_nodes_, allocate space in adjacency for edges
    // incident to this new node, and call the private constructor of Node to return
    internal_nodes_.push_back(position);
    num_nodes_++;
    internal_adj_.push_back({{}});
    return Node(this, num_nodes_);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: check if n points to the current Graph
    // and ensure the index value is valid
    if ((this == n.graph_) and (n.index_ < this->num_nodes_) and (n.index_ >= 0))
      return true;
    else
      return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: Construct valid Node using this and input index, 
    // ensuring input indes is valid for this graph
    assert(i < this->num_nodes_);
    assert(i >= 0);
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: Public constructor for invalid edge
      // Do nothing on purpose, this edge is invalid
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: return a Node created with the current Graph and index
      // given by the first entry in this Edge's internal tuple
      return Node(graph_, std::get<0>(graph_->internal_edges_[index_])); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: return a Node created with the current Graph and index
      // given by the second entry in this Edge's internal tuple
      return Node(graph_, std::get<1>(graph_->internal_edges_[index_])); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // HW0: check if both point to same graph and have same internal index
      if ((this->graph_ != e.graph_) or (this->index_ != e.index_))
        return false;
      else
        return true;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: If the two have same graph, simply compare the index_ values, 
      // which refer to the order of the edges within the internal structure
      if (this->graph_ == e.graph_){
        if (this->index_ < e.index_)
          return true;
        else
          return false;
      }

      // HW0: Otherwise, use ordering on Graph pointers
      std::less<graph_type*> check;
      if (check(this->graph_,e.graph_))
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // HW0: Initialize member variables: pointer to Graph and index
    // which tells us where to look in internal data to get the endpoints
    graph_type* graph_;
    size_type index_;

    // HW0: Private constructor for valid Edges, given Graph and index
    Edge(const graph_type* graph, size_type ind)
      : graph_(const_cast<graph_type*>(graph)), index_(ind) {}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: return the value stored in internal num_edges_ variable
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: create valid Edge using the current Graph and i,
    // ensuring index is valid
    assert(i < this->num_edges_);
    assert(i >= 0);
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: Get the nodal indices, search through adjacency structure
    // to see if these nodes are connected

    size_type inda = a.index_;
    size_type indb = b.index_;

    auto search = internal_adj_[inda].find(indb);

    if (search != internal_adj_[inda].end())
      return true;
    else
      return false;
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

    // HW0: Check if there is an edge by calling the previous method,
    // if it does exist then get the index associated with this edge and return
    // a valid edge

    bool exists = has_edge(a, b);
    size_type inda = a.index_;
    size_type indb = b.index_;
    if (exists) {
      size_type edge_ind = internal_adj_[inda][indb];
      return Edge(this,edge_ind);
    }

    // Otherwise, create this edge in edge list, modify the nodal
    // adjacencies, increment # edges, and return valid Edge
    else {
      std::tuple <size_type, size_type> add;
      add = std::make_tuple(inda,indb);
      num_edges_++;
      internal_edges_.push_back(add);
      internal_adj_[inda][indb] = num_edges_;
      internal_adj_[indb][inda] = num_edges_;
      return Edge(this,num_edges_);
    }
    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: clear internal vectors and set counters to zero
    internal_nodes_.clear();
    num_nodes_ = 0;
    internal_edges_.clear();
    num_edges_ = 0;
    internal_adj_.clear();

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

  // HW0: Initialize (private) member variables

  // Vector of Point objects defining nodes, index of node is location in vector
  std::vector<Point> internal_nodes_; 
  size_type num_nodes_;

  // Vector of tuples of nodal indices defining edges, index of edge is location in vector
  std::vector< std::tuple<size_type, size_type> > internal_edges_;
  size_type num_edges_;

  // Vector of maps, each element of vector is a map where the keys are the indices
  // of adjacent nodes, and the values are the indices of that edge in internal_edges_
  std::vector< std::unordered_map<size_type,size_type> > internal_adj_;

};

#endif // CME212_GRAPH_HPP
