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

  // JX: Predeclare the internal struct, which is used later in Helper method
  //     of node and edge class. 
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
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): nodes_(), edges_(), adjacency_() {
    // HW0: YOUR CODE HERE
    // JX: use the constructor of node,edge,and ajacency list here
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
      // JX: not sure what to do here????????????
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // return Point();

      // JX: return the position by internal_node
       //    position: internal attribute of node in graph
      return graph_->nodes_[index_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // return size_type(-1);

      // JX: return node's unique index number 
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
      // HW0: YOUR CODE HERE
      // (void) n;          // Quiet compiler warning
      // return false;

      // JX: Equal nodes have the same graph and the same index.
      //     thus we need to retrun true when both conditoins are satisfied
      if (this->graph_ == n.graph_ && this->index_ == n.index_){
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
      // (void) n;           // Quiet compiler warning
      // return false;

      // JX: there are two possible sitation where this node is less than n:
      //     1) graph is less than n's graph in order
      //     2) in same graph, but the index is less than n's index in order
      if (this->graph_ < n.graph_ ){
        return true;
      } else  if (this->graph_ == n.graph_ && this->index_ < n.index_) {
        return true;
      } else{
        return false;
      }

      // Version 2:
      // return index_ < n.index_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // JX: sizeof(Graph::Node) ≤ 16 bytes
    Graph* graph_;       // pointer back to the Graph container
    size_type index_;    // node's unique index number


    // JX: /** Private Constructor */
    Node(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // return 0;

    // JX: #nodes = size of nodes_, vector of nodes in graph
    return nodes_.size();
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
    // (void) position;      // Quiet compiler warning
    // return Node();        // Invalid node

    // JX: 
    //    create the index for new node
    size_type new_nidx = num_nodes();
    //    create the new node with given position : 
    //    Version 2: 
              // internal_node new_node;
              // new_node.position = position;
    internal_node new_node {position};
    //    push the new node to nodes_
    nodes_.push_back(new_node);
    //    add the node in adjacency list
    adjacency_.push_back(std::vector<size_type> ());
    //    return added node; 
    return Node(this, new_nidx);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // (void) n;            // Quiet compiler warning
    // return false;

    // JX:, this is a graph class object, check if n is a Node of this Graph
    if (n.graph_ == this){
      return true;
    } else{
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
    // (void) i;             // Quiet compiler warning
    // return Node();        // Invalid node

    // JX:
    //    Return the node with index i in graph      
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // return Node();      // Invalid Node

      // JX: 
      return  Node(graph_, n1idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // return Node();      // Invalid Node

      // JX: 
      return  Node(graph_, n2idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // (void) e;           // Quiet compiler warning
      // return false;

      // DO we need the edge id to the same??????????
      // JX: since equal edges represent the same "undirected" edge
      //     we have to consider two situations where (a',b')==(b,a)==(a,b)
      if (this->graph_ == e.graph_ && this->n1idx_ == e.n1idx_ && this->n2idx_ == e.n2idx_){
        return true;
      } else if (this->graph_ == e.graph_ && this->n1idx_ == e.n2idx_ && this->n2idx_ == e.n1idx_){
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
      // (void) e;           // Quiet compiler warning
      // return false;

      // JX:
      // Version 1:
      if (this->graph_ < e.graph_ ){
        return true;
      // this->eidx_ & eidx_ should both work
      } else  if (this->graph_ == e.graph_ && this->eidx_ < e.eidx_) {
        return true;
      } else{
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

    // JX: sizeof(Graph::Edge) ≤ 32 bytes
    Graph* graph_;       // pointer back to the Graph container 
    size_type eidx_;     // edge's unique index number
    size_type n1idx_;    // index of the first vertex of edge
    size_type n2idx_;    // index of the second vertex of edge

    // JX: /** Private Constructor */
    Edge(const Graph* graph, size_type eidx, size_type n1idx, size_type n2idx)
        : graph_(const_cast<Graph*>(graph)), eidx_(eidx), n1idx_(n1idx), n2idx_(n2idx) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // return 0;

    // JW: number of edge is size of edge vector, edges_
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // (void) i;             // Quiet compiler warning
    // return Edge();        // Invalid Edge

    // JX:
    // Version 1:
    return Edge(this, i, edges_[i].idx1, edges_[i].idx2);

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // (void) a; (void) b;   // Quiet compiler warning
    // return false;

    // JX: recursion, find if b is in the adj_list of a, and also the other side
    //     for my method, I have check both (a,b), (b,a)
    //     cpmplexity: O(num_nodes)

    size_type a_idx = a.index();
    size_type b_idx = b.index();
    std::vector<size_type> a_adj = adjacency_[a_idx];
    std::vector<size_type> b_adj = adjacency_[b_idx];

    // with a, check if b in adjacency_ of a
    for (auto it = a_adj.begin(); it != a_adj.end(); ++it){
      if (*it == b_idx){
        return true;
      }
    } 
    // with b, check if a in adjacency_ of b
    for (auto it = b_adj.begin(); it != b_adj.end(); ++it){
      if (*it == a_idx){
        return true;
      }
    } 
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
    // HW0: YOUR CODE HERE
    // (void) a, (void) b;   // Quiet compiler warning
    // return Edge();        // Invalid Edge

    // JX: the main part is to update adj_list and update number of edge
    //     complexity: O(num_edges)
    size_type a_idx = a.index();
    size_type b_idx = b.index();

    if (has_edge(a,b)){
      // JX: if the edge exist, find the edge id
      // JX: existe_idx is used here to fix a error happen when only 
      //     return in the two if statements below
      size_type existe_idx;
      for(unsigned int i = 0; i != edges_.size(); i++) {
        if (edges_[i].idx1 == a_idx && edges_[i].idx2 == b_idx){
          existe_idx = i;
        }
        if (edges_[i].idx1 == b_idx && edges_[i].idx2 == a_idx){
          // return Edge(this, i, a_idx, b_idx);
          existe_idx = i;
        }
      }
      return Edge(this, existe_idx, a_idx, b_idx);
    }else {
      //  JX: if the edge not exist, add the edge and update adjacency list

      // create index for new edge
      size_type newe_idx = num_edges();
      // add the new edge to edge list
          // Version 2:
          // internal_edge new_edge ;
          // new_edge.idx1 = a_idx;
          // new_edge.idx2 = b_idx;
      internal_edge new_edge {a_idx, b_idx};
      edges_.push_back(new_edge);
      // add the new adjacency relation to adjacenct list
      adjacency_[a_idx].push_back(b_idx);
      adjacency_[b_idx].push_back(a_idx);
      // return the new edge added
      return Edge(this, newe_idx, a_idx, b_idx);
    }

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    nodes_.clear();
    edges_.clear();
    adjacency_.clear();

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

  // JX:  Internal type for graph node and edge
  struct internal_node {
    Point position; // postion of a node
  };

  struct internal_edge{
    size_type idx1; // index of first vertex of edge
    size_type idx2; // index of second vertex of edge
  };

  // JX: data member
  std::vector<internal_node> nodes_;  // vector of nodes 
  std::vector<internal_edge> edges_;  // vector of edges
  // JX: [[1,5,9],[0,5,7],..], index represent the idx of node1, vector inside
  //                           represent the index of adjacent nodes of node1
  std::vector<std::vector<size_type>> adjacency_; 

};

#endif // CME212_GRAPH_HPP
