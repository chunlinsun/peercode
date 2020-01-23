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

  // Declarations of important internal types you need
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
    Node() {
      Gr_ = nullptr; // empty graph
      idx_ = size_type(-1); // invalid index
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(Gr_ != nullptr and Gr_->has_node(*this)); // check if valid graph and node
      assert(idx_ >= 0 and idx_ < Gr_->num_nodes()); // check if valid index
      return Gr_->node_lst[idx_].pos_; // return the position
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(Gr_ != nullptr); // check if valid graph
      // note that idx_ (index directly associated with the node) is different from index_ (index from the graph)
      return Gr_->node_lst[idx_].index_; // return the index from the graph (so keeping up to date)
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
      return Gr_ == n.Gr_ and idx_ == n.idx_;
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
      // First compare the graph pointer
      if (Gr_ < n.Gr_) return true;
      if (Gr_ > n.Gr_) return false;
      // if Gr_ == n.Gr, then compare the index
      return idx_ < n.idx_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* Gr_; // pointer to the graph Gr
    size_type idx_; // node index
    // Valid node constructor (only be constructed by a graph)
    Node(const graph_type* Gr, size_type idx){
      Gr_ = const_cast<graph_type*> (Gr);
      idx_ = idx;
    }
  };

  /** Struct for node detailed info (position and index)
  */
  struct Node_Detail{
    Point pos_;
    size_type index_;
    Node_Detail(const Point& pos, size_type index){
      pos_ = pos;
      index_ = index;        
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return index_lst.size(); // index_lst[index_] = idx_; index_ is contiguous
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
    // new index_ (from graph) and idx_ (associated with node)
    size_type new_idx_, new_index_;
    new_idx_ = node_lst.size();
    new_index_ = this->num_nodes();
    // new node
    Node_Detail new_node_detail = Node_Detail(position, new_index_);
    index_lst.push_back(new_idx_);
    node_lst.push_back(new_node_detail);
    adj_lsts.push_back(std::vector<Edge_Detail>()); // initialize empty adjacency list for the node
    return Node(this, new_idx_);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // check if the graph is the same, the index is valid and belongs to the graph
    return n.Gr_ == this and n.idx_ != size_type(-1) and n.index() < this->num_nodes();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < this->num_nodes());
    return Node(this, index_lst[i]);        
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
      Gr_ = nullptr;
      // node index of two end points
      n1_index_ = size_type(-1);
      n2_index_ = size_type(-1);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Gr_->node(n1_index_);     
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Gr_->node(n2_index_);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // when storing the edges, we will make sure that n1_index_ > n2_index_
      return Gr_ == e.Gr_ and n1_index_ == e.n1_index_ and n2_index_ == e.n2_index_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // again, when storing the edges, we will make sure that n1_index_ > n2_index_
      // First compare the graph pointer
      if (Gr_ < e.Gr_) return true;
      if (Gr_ > e.Gr_) return false;
      // if Gr_ == n.Gr, then compare n1_index_
      if (n1_index_ < e.n1_index_) return true;
      if (n1_index_ > e.n1_index_) return false;
      // if Gr_ == n.Gr, n1_index_ == n.n1_index_, then compare n2_index_
      return n2_index_ < e.n2_index_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type*  Gr_;
    size_type n1_index_;
    size_type n2_index_;

    // Valid edge constructor (only be constructed by a graph)
    Edge(const graph_type* Gr, size_type n1_index, size_type n2_index){
      Gr_ = const_cast<graph_type*> (Gr);
      n1_index_ = n1_index;
      n2_index_ = n2_index;
    }

  };

  /** Struct for edge detailed info (n1_index_ and n2_index_)
  */
  struct Edge_Detail{
    size_type n1_index_;
    size_type n2_index_;
    Edge_Detail(size_type n1_index, size_type n2_index){
      n1_index_ = n1_index;       
      n2_index_ = n2_index;     
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // Directly find the number of edges as the vector size (O(1) time)
    return edge_lst.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Directly access the edge i from the edge list/vector (O(1) time)
    assert(i < this->num_edges());
    return Edge(this, edge_lst[i].n1_index_, edge_lst[i].n2_index_);       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) and has_node(b));
    /** Loop over the adjacency list of a or b (depending on which one is shorter).
    The complexity is O(num_nodes()), or indeed O(maximum node degree). 
    It can be further improved to O(log maximum node degree) if using binary search. 
    Similar improvement for add_edge.
    */
    const std::vector<Edge_Detail>& adj_lst_a = adj_lsts[a.index()];
    const std::vector<Edge_Detail>& adj_lst_b = adj_lsts[b.index()];
    const std::vector<Edge_Detail>* adj_lst_search = &adj_lst_a;
    if (adj_lst_a.size() > adj_lst_b.size()){
      adj_lst_search = &adj_lst_b;
    }
    const std::vector<Edge_Detail>& adj_lst_search1 = *adj_lst_search;
    for (size_type i = 0; i < adj_lst_search1.size(); i++) {
      Edge_Detail Ei = adj_lst_search1[i]; 
      size_type index1_ = Ei.n1_index_;
      size_type index2_ = Ei.n2_index_;
      // remember to check two directions
      if ((index1_ == a.index() and index2_ == b.index()) \
        or (index1_ == b.index() and index2_ == a.index())) {
          return true;
      }
    }
    return false; // edge not found
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
    // check if the nodes are valid
    size_type a_index = a.index();
    size_type b_index = b.index();
    assert(has_node(a) and has_node(b) and a_index != b_index);
    // Ensure that n1_index > n2_index
    size_type n1_index, n2_index;
    if (a_index > b_index){
      n1_index = a_index;
      n2_index = b_index;
    }
    else {
      n1_index = b_index;
      n2_index = a_index;
    }
    // update the related vectors
    if (!has_edge(a, b)) {
      Edge_Detail new_edge_detail = Edge_Detail(n1_index, n2_index);
      edge_lst.push_back(new_edge_detail);
      adj_lsts[a.index()].push_back(new_edge_detail);
      adj_lsts[b.index()].push_back(new_edge_detail);
    }

    return Edge(this, n1_index, n2_index); // return the current edge    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // node vectors clear
    node_lst.clear();
    index_lst.clear();
    // edge vectors clear
    edge_lst.clear();
    adj_lsts.clear();
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
  // Graph class's internals:
  //   helper functions, data members, and so forth.
  // vectors for the nodes
  std::vector<Node_Detail> node_lst; // list of nodes, each entry storing the pos_ and index_
  std::vector<size_type> index_lst; // list of index_ (from graph), each entry storing the idx_ (for the node)
  // vectors for the edges
  std::vector<Edge_Detail> edge_lst; // list of edges, each entry storing n1_index_ and n2_index_
  std::vector<std::vector<Edge_Detail>> adj_lsts; // list of adjacency lists for each node, each entry is a list of edges connected to the node
};

#endif // CME212_GRAPH_HPP
