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
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  struct node_obj;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  using node_type_value = V;
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
  Graph(): nodes_(), adjacent_(), n_edges_(0) {}

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
  class Node:private totally_ordered<Node> {
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
      return graph_->nodes_[n_id_].pos_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return n_id_;
    }

    /** @brief Returns the value of a node as a mutable reference. */
    node_type_value& value() {
      return graph_->nodes_[n_id_].val_;
    }
    /** @brief Returns the value of a node as a non-mutable reference. */
    const node_type_value& value() const {
      return graph_->nodes_[n_id_].val_;
    }

    /** @brief Returns the number of incident edges to the node. */
    size_type degree() const{
      return graph_->adjacent_[n_id_].size();
    }

    /** @brief Returns the iterator corresponding to the first visited edge. */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, n_id_, 0);
    }

    /** @brief Returns the iterator corresponding to the last visited edge. */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, n_id_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (n.n_id_ == n_id_);
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
      return ((graph_ == n.graph_) && (n_id_ < n.n_id_)) || (graph_ < n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;
    size_type n_id_;

    Node(const Graph* graph, size_type n_id): graph_(const_cast<Graph*>(graph)), n_id_(n_id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& position,
                const node_type_value& val = node_type_value()) { 
    size_type new_id = num_nodes(); // This will increment the size automatically when we push back the node.
    node_obj new_node = node_obj(new_id, position, val);

    // Add the new node along with its empty vector of neighbors
    nodes_.push_back(new_node);
    adjacent_.push_back(std::vector<size_type> ());

    return Node(this, new_id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_ == this) && (n.n_id_ < n.graph_->size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, n1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, n2_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // The undirected nature of the graph requires us to check both nodes
      return (e.graph_ == graph_) && ((e.n1_id_ == n1_id_ && e.n2_id_ == n2_id_) || (e.n1_id_ == n2_id_ && e.n2_id_ == n1_id_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (e.graph_ == graph_){
       return (e.n1_id_ == n1_id_ && e.n2_id_ > n2_id_) || (e.n1_id_ > n1_id_);
      }
      return e.graph_ > graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type n1_id_;
    size_type n2_id_;
    
//--functionality_1
//--same issue as in HW0, use the pointer you feed in ("graph") to initalize the member variable graph_ 
//--START
    // Constructor
    Edge(const Graph* graph, size_type n1_id, size_type n2_id): 
     graph_(const_cast<Graph*> (graph)), n1_id_(n1_id), n2_id_(n2_id) {}
  };
//--END

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return n_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    auto edge_iter = edge_begin();
    for (; i != 0; --i){
      ++ edge_iter;
    }
    return *edge_iter;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for(size_type i=0; i < adjacent_[a.n_id_].size(); i++){
      if(adjacent_[a.n_id_][i] == b.n_id_)  // Ordering b -> a vs. a -> b does not matter: Graph Undirected.
       return true;
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
    if (has_edge(a, b)) {
     return Edge(this, a.n_id_, b.n_id_);
    }

    adjacent_[a.n_id_].push_back(b.n_id_);
    adjacent_[b.n_id_].push_back(a.n_id_);

    n_edges_ ++;

    return Edge(this, a.n_id_, b.n_id_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
   nodes_.clear();
   adjacent_.clear();
   n_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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

    /** @brief Dereference Operator to access the node with its @a index_. */
    Node operator*() const {
      return graph_->node(index_);
    }

    /** @brief Returns the iterator pointing to the next node or @a node_end(). */
    NodeIterator& operator++(){
      index_++;
      return *this;
    }

    /** @brief Check if the NodeIterator @a n_iter is equal to this NodeIterator.
    *
    * @param @a n_iter    NodeIterator to be checked
    * @return             Boolean true if the iterators are equal
    */
    bool operator==(const NodeIterator& n_iter) const {
      return (n_iter.graph_ == graph_) && (n_iter.index_ == index_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type index_;
    NodeIterator(const Graph* graph, const size_type index = 0):
       graph_(const_cast<Graph*>(graph)), index_(index) {}
  };

  /** @brief Returns the iterator pointing at the first node to in the iterator. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** @brief Returns the iterator pointing at the last node to be iterated through. */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    /** @brief Returns the edge currently pointed at by the iterator. */
    Edge operator*() const{
      return Edge(graph_, root_, graph_->adjacent_[root_][traversed_]);
    }
   
    /** @brief Returns the iterator pointing to the next edge or @a edge_end().
     *
     * @pre 0 <= @a num_edges_traversed_ < node.degree()
     *
     * @post @a num_edges_traversed_ = @ num_edges_traversed + 1
     */
    IncidentIterator& operator++(){
      ++traversed_;
      return *this;
    }
   
    /** @brief Verify if an iterator is equal to this incident iterator.
     *
     * @param[in] @a incident  IncidentIterator to be checked.
     * @return                 boolean true if @a incident is equal to this iterator.
     */
    bool operator==(const IncidentIterator& incident) const{
      return (incident.graph_ == graph_) && (incident.root_ == root_) && (incident.traversed_ == traversed_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type root_; // ID of the node for which we are traversing edges.
    size_type traversed_; // Current Number of Traversed Edges.

    IncidentIterator(const Graph* graph, const size_type root, const size_type traversed=0):
      graph_(const_cast<Graph*> (graph)), root_(root), traversed_(traversed) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /** @brief Returns the edge currently iterated through in the EdgeIterator. */
    Edge operator*() const{
      return Edge(graph_, current_, graph_->adjacent_[current_][visited_]);
    }

    /** @brief Returns an iterator pointing to the next edge in or @a this.edge_end().
     *
     * @pre 0 <= @a current_ < @a num_nodes()
     * @pre 0 <= @a visited_ < @a current_.degree()
     * @pre the IDs of the nodes are subsequent and continous. 
     *
     * @post visited_ = visited_ + 1 or visit_ = 0
     * @post 0 <= current_ <= @a num_nodes()
     */
    EdgeIterator& operator++(){
      visited_++;
      while(current_ < graph_->num_nodes()) {
        while(visited_ < graph_->node(current_).degree()){
          if (current_ < graph_->adjacent_[current_][visited_])
           return *this;
          ++visited_;
        }
        visited_ = 0;
        ++current_;
      }
      visited_ = 0;
      return *this;
    }
   
    /** @brief check if an edge iterator input is equal to this edge iterator.
     *
     * @param[in] @a edge_iter  EdgeIterator input for check
     * @return                  Boolean true if iterators are equal
     */
    bool operator==(const EdgeIterator& edge_iter) const{
      return (edge_iter.graph_ == graph_) && (edge_iter.current_ == current_) && (edge_iter.visited_ == visited_);
    }

   private:
    friend class Graph;

    Graph* graph_;
    size_type current_;
    size_type visited_;

    EdgeIterator(const Graph* graph, const size_type current = 0, const size_type visited = 0):
      graph_(const_cast<Graph*>(graph)), current_(current), visited_(visited) {}

  };
 
  /** @brief Returns an iterator indicating where to start visiting edges. */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0, 0);
  }

  /** @brief returns an iterator indicating where to stop visiting edges. */
  edge_iterator edge_end() const {
    return edge_iterator(this, num_nodes(), 0);
  }


 private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct node_obj{
   Point pos_;
   size_type n_id_;
   node_type_value val_;

   // Constructor:
   node_obj(const size_type n_id, const Point& position, node_type_value val)
    : pos_(position), n_id_(n_id), val_(val) {}
  };

  std::vector<node_obj> nodes_;
  std::vector<std::vector<size_type>> adjacent_;
  size_type n_edges_;

};

#endif // CME212_GRAPH_HPP
