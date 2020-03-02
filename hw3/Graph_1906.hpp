#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
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
  using node_value_type = V;
  using edge_value_type = E;

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
  Graph() : size_(0), numedges_(0) {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() = default;
  
  // 
  // nodeinfo struct
  //
  struct nodeinfo
  {
  Point p_;            // node position
  node_value_type v_;  // node value
  size_type idx_;      // node index
  
  nodeinfo(Point p, node_value_type v, size_type idx) : p_(p), v_(v), idx_(idx) { }

  };

  // 
  // edgeinfo struct
  //
  struct edgeinfo
  {
    size_type node1uid_;  // node 1 uid
    size_type node2uid_;  // node 2 uid
    size_type index_;     // edge index

    edgeinfo(size_type node1uid, size_type node2uid, size_type index) : 
    node1uid_(node1uid), node2uid_(node2uid), index_(index) { }
  };

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

    Node() : graph_(nullptr), index_(100) {
      // HW0: YOUR CODE HERE
      // index is set to an arbitrarily large number
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // return the point
      return (graph_->nodes_).at(graph_->i2u_[index_]).p_;
      // return graph_->allpoints_.at(index_);
    }

    /** Return this node's modifiable position */
    Point& position() {
      // return a reference to the point
      // return graph_->allpoints_.at(index_);
      return (graph_->nodes_).at(graph_->i2u_[index_]).p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(index_ < graph_->size_);
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return a modifiable reference to this node's value, a variable whose 
    type is specified by the Graph class declaration */
    node_value_type& value() {return (graph_->nodes_).at(graph_->i2u_[index_]).v_;}

    /** Return a non-modifiable reference to this node's value, const version of
    the value() function */
    const node_value_type& value() const {return (graph_->nodes_).at(graph_->i2u_[index_]).v_;}

    /** Return this node's degree, a measure of how many nodes are incident to it */
    size_type degree() const {
      // pull up the adjacency map for this node
      // find out if node even has an edge, if not 
      // degree = 0
      auto search = graph_->adjacency_.find(graph_->i2u_[index_]);
      if (search == graph_->adjacency_.end()) return 0;

      auto a_map = graph_->adjacency_.at(graph_->i2u_[index_]);
      return a_map.size();
    }
    
    /** Return an iterator to the first edge this node is incident to */
    incident_iterator edge_begin() const {
      // should be an iterator over the adjacency map for this node
      return incident_iterator(0, index_, graph_);
    }

    /** Return an iterator to one past the last edge this node is incident to */
    incident_iterator edge_end() const {
      return incident_iterator(this->degree(), index_, graph_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // Check the Graph pointer and the index
      if (this->graph_ != n.graph_)  return false;
      if (this->index_ != n.index_) return false; 
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
      // HW0: YOUR CODE HERE
      if (index_ == n.index_ && this->graph_ != n.graph_) {
        return std::less<Graph*>{}(this->graph_, n.graph_);
      }
      return (this->index_ < n.index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Graph will contain an index and a pointer to a Graph
    graph_type* graph_;
    size_type index_;

    // Private Constructor
    Node(const graph_type* graph, unsigned index)
        : graph_(const_cast<graph_type*>(graph)), 
	  index_(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
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
                const node_value_type& vl = node_value_type()) {
    // HW0: YOUR CODE HERE
    // allocate memory for position and add point to vector
    Point newpoint = Point(position.x, position.y, position.z);
    nodes_.push_back(nodeinfo(newpoint, vl, size_));
    i2u_.push_back(nodes_.size() - 1);
    Node n = Node(this, size_);
    size_++;
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this != n.graph_) return false;
    // making sure the index is valid
    assert(n.index_ < size_);
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= size_) {
    std::cout << "Invalid index" << std::endl;
    return Node();
    }
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() : graph_(nullptr) {
      // set edgeinfo_ to default value , index to arbitrary number, and graph pointer to null
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, node2_);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Check that they have the same graph pointer
      if (this->graph_ != e.graph_) return false;

      // Check that they have the same nodes
      Node nodea = this->node1();
      Node nodeb = this->node2();
      Node nodec = e.node1();
      Node noded = e.node2();
      
      if ((nodea == nodec) && (nodeb == noded)) return true;
      if ((nodea == noded) && (nodeb == nodec)) return true;
      
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->index_ == e.index_ && this->graph_ != e.graph_) {
        return std::less<Graph*>{}(this->graph_, e.graph_);
      }
      return (this->index_ < e.index_);
    }
    
    edge_value_type& value() {
      size_type n1 = graph_->i2u_[node1_];
      size_type n2 = graph_->i2u_[node2_];
      if (n1 < n2) return graph_->restlength_[n1][n2];
      else return graph_->restlength_[n2][n1];
    }

    const edge_value_type& value() const {
    size_type n1 = (graph_->i2u_).at(node1_);
    size_type n2 = (graph_->i2u_).at(node2_);
    if (n1 < n2) return graph_->restlength_.at(n1).at(n2);
    else return graph_->restlength_.at(n2).at(n1);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // pointer to the parent Graph
    graph_type* graph_;
    // index
    size_type index_;
    // node indices, user facing
    size_type node1_;
    size_type node2_;

    Edge(const graph_type* graph, size_type index, const Node& a, const Node& b) 
        : graph_(const_cast<graph_type*>(graph)), index_(index),
	  node1_(a.index()), node2_(b.index()) {
	}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return numedges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= numedges_) {
      std::cout << "Invalid edge index" << std::endl;
      return Edge();
    }
    edgeinfo e = alledges_.at(ei2u_[i]);
    assert(i == e.index_);
    size_type node1i = nodes_.at(e.node1uid_).idx_;
    size_type node2i = nodes_.at(e.node2uid_).idx_;
    return Edge(this, i, Node(this, node1i), Node(this, node2i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Test validity of nodes
    // if (!(this->has_node(a) && this->has_node(b))) return false;
    // if (numedges_ == 0) return false;

    // find out if node even has an edge
    auto search = adjacency_.find(i2u_[a.index()]);
    if (search == adjacency_.end()) return false;

    // find out if node is attached to b
    auto a_map   = adjacency_.at(i2u_[a.index()]);
    auto search2 = a_map.find(i2u_[b.index()]);
    if (search2 == a_map.end()) return false;
    
    // otherwise we have found our edge
    return true;
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
  Edge add_edge(const Node& a, const Node& b, 
                const edge_value_type& vl = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // a and b are valid nodes

    if (!(this->has_node(a) && this->has_node(b))) {
      std::cout << "Invalid nodes" << std::endl;
      return Edge();}


    // case where edge doesn't exist yet
    if (!(this->has_edge(a, b))) {
      Edge e = Edge(this, numedges_, a, b);
      size_type n1 = i2u_[a.index()];
      size_type n2 = i2u_[b.index()];
      alledges_.push_back(edgeinfo(n1, n2, numedges_));
      ei2u_.push_back(alledges_.size() - 1);
      // update adjacency map in both directions
      // adjacency map has euids
      adjacency_[n1][n2] = (alledges_.size() - 1);
      adjacency_[n2][n1] = (alledges_.size() - 1);
      // update the edge values
      if (n1 < n2) restlength_[n1][n2] = vl;
      else restlength_[n2][n1] = vl;
      // find the rest length
      Point acc  = Point(0,0,0);
      acc += a.position();
      acc -= b.position();
      edge_value_type length = norm(acc);
      // store the length, just once this time
      if (n1 < n2) {
      restlength_[n1][n2] = length;}
      else {
      restlength_[n2][n1] = length;}
      numedges_++;
      return e;
    }

    // case where edge exists, return
    size_type edge_uid   = adjacency_.at(i2u_[a.index()]).at(i2u_[b.index()]);
    size_type edge_index = alledges_.at(edge_uid).index_;

    return Edge(this, edge_index, a, b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // empty containers and set counts to zero
    size_ = 0;
    numedges_ = 0;
    restlength_.clear();
    nodes_.clear();
    i2u_.clear();
    alledges_.clear();
    adjacency_.clear();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the underlying node from the iterator object */ 
    Node operator*() const {
      // create a node object
      return Node(graph_, index_);
    }
    
    /** Step once through the iterator object */
    NodeIterator& operator++() {
      // increment the index
      index_++;
      return *this;
    }
    
    /** Check for equality between two iterators
    @pre _iter2_ is a valid Node Iterator */
    bool operator==(const NodeIterator& iter2) const {
      return (iter2.index_ == index_ && iter2.graph_ == graph_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    // pointer to the graph
    graph_type* graph_;
    // index we are on
    size_type index_;

    // Constructor
    NodeIterator(const graph_type* graph, size_type index) 
                : graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return an iterator to the first node of the graph */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  
  /** Return an iterator to one past the last node of the graph */
  node_iterator node_end() const  {
    return NodeIterator(this, size_);
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the underlying edge referenced by the IncidentIterator object */
    Edge operator*() const {
      // Create the first node
      Node nodea         = Node(graph_, node1_);
      // std::cout << nodea.degree() << std::endl;
      // find the edge and node2 indices
      size_type no2id    = mapitr_->first;
      size_type node2    = graph_->nodes_.at(no2id).idx_;
      size_type edge_uid = mapitr_->second;
      size_type edge_ind = graph_->alledges_.at(edge_uid).index_;
      Node nodeb         = Node(graph_, node2);
      return Edge(graph_, edge_ind, nodea, nodeb);
    }

    /** Increment the iterator in place */
    IncidentIterator& operator++() {
      // step through the map
      mapitr_++;
      counter_++;
      return *this;
    }
    
    /** Check equality between this iterator and a user-supplied one
    @ pre _it2_ is a valid IncidentIterator object */
    bool operator==(const IncidentIterator& it2) const {
      return (counter_ == it2.counter_ && node1_ == it2.node1_ && graph_ == it2.graph_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // counter through keys in the map
    size_type counter_;
    // index of first node
    size_type node1_;
    // pointer to the graph
    graph_type* graph_;
    // Private constructor
    // an iterator over a map
    std::map<unsigned, unsigned>::iterator mapitr_;

    IncidentIterator(size_type counter, size_type node1, const graph_type* graph) : 
                     counter_(counter), node1_(node1), graph_(const_cast<graph_type*>(graph)){ 
                     // Iterator over map only gets initialized if the node has edges
		     Node nodea = Node(graph_, node1_);
		     size_type nodeid = graph_->i2u_[node1_];
		     if (nodea.degree() != 0) {
		       mapitr_ = (graph_->adjacency_.at(nodeid)).begin();
		       // std::cout << mapitr_->second << std::endl;
		       }
		     }
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
    /** Return the underlying edge of the iterator */
    Edge operator*() const {
      // return the edge from the current index
      return graph_->edge(index_);
    }

    /** In-place increment the iterator */
    EdgeIterator& operator++() {
      index_++;
      return *this;
    }
    
    /** Check equality between two iterators
    @pre _e2_ is a valid EdgeIterator object */
    bool operator==(const EdgeIterator& e2) const {
      return (graph_ == e2.graph_ && index_ == e2.index_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;
    
    // Constructor
    EdgeIterator(const graph_type* graph, size_type index) :
                 graph_(const_cast<graph_type*>(graph)), index_(index) { }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return an iterator to the first edge of the graph */
  edge_iterator edge_begin() const {
    // First node, // First edge
    return EdgeIterator(this, 0);
    // return EdgeIterator(node_begin(), (*node_begin()).edge_begin(), size_);
  }

  /** Return an iterator to one past the last edge of the graph */
  edge_iterator edge_end() const {
    // Last node, // last edge
    return EdgeIterator(this, numedges_);
    // return EdgeIterator(node_end(), (*node_end()).edge_end(), size_);
  }

  /**
   * Remove a specified node from the graph
   * 
   * @param[in] n node to delete from the graph
   * @return    The number of nodes removed from the graph
   *
   * @pre  (this->has_node(_n_))
   *
   * @post (!this->has_node(_n_))
   *       for all edges e incident to _n_, (!this->has_edge(e))
   *       new (this->size()) == old (this->size()) - 1
   *       node is invalidated if node.index() >= new this->size()
   *       NodeIterator(node.index()) is invalidated if node.index() == _n_.index()
   *                                               || node.index() >= new this->size()
   *
   * Complexity O(d) armotized where d is the degree of _n_.
   *        
   */
  size_type remove_node(const Node& n) {
    // remove node only if it exists
    if (!this->has_node(n)) return 0;
    // remove all incident edges
    std::vector<Node> v;
    auto dic = adjacency_[i2u_[n.index()]];

    for (auto it = dic.begin(); it != dic.end(); it++) {
      unsigned n_index = nodes_[it->first].idx_;
      v.push_back(Node(this, n_index));
    }
    
    for (unsigned i = 0; i < v.size(); i++) {
      remove_edge(n, v[i]);
    }

    // index of node to be removed
    size_type ind = n.index();
    // change index of node before swapping
    nodes_[i2u_[size_ - 1]].idx_ = ind;
    std::swap(i2u_[ind], i2u_[size_ - 1]);
    i2u_.pop_back();
    size_ -= 1;
    return 1;
  }

  /**
   * Remove the node referenced by an iterator from the graph
   * 
   * @param in  n_it iterator to node to delete from the graph
   * @return    Fresh iterator pointing to the new node at (*n_it).index()
   *
   * @pre  (this->has_node(*_n_it_))
   *
   * @post (!this->has_node(*_n_it_))
   *       for all edges e incident to (*_n_it_), (!this->has_edge(e))
   *       new (this->size()) == old (this->size()) - 1
   *       node is invalidated if node.index() >= new this->size()
   *       NodeIterator(node.index()) is invalidated if node.index() == (*_n_it_).index()
   *                                               || node.index() >= new this->size()
   *
   * Complexity O(d) armotized where d is the degree of _n_.
   *        
   */
  node_iterator remove_node(NodeIterator n_it) {
    Node n = *n_it;
    remove_node(n);
    return n_it;
  }

  /**
   * Remove the edge connecting two nodes from the graph
   *
   * @param [in] a first node to the edge to be removed from the graph
   * @param [in] b second node to the edge to be removed from the grah
   * @return     number of edges removed from the graph
   *
   * @pre (this->has_node(_a_))
   * @pre (this->has_node(_b_))
   * @pre (this->has_edge(_a_, _b_))
   *
   * @post (!this->has_edge(_a_, _b_))
   *       new this->num_edges = old this->num_edges - 1
   *       EdgeIterator eit is invalidated if (*eit).index == Edge(_a_, _b_, this).index
   *                                         || (*eit).index >= new this->num_edges()
   *
   * Complexity O(1) armotized operations.
   *       
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // Remove edge only if it exists
    if (!this->has_edge(a, b)) return 0;
    Edge e = this->add_edge(a, b);
    size_type rval = remove_edge(e);
    return rval;
  }

  /**
   * Remove the edge provided from the graph
   *
   * @param [in] e edge to be removed from the graph
   * @return     number of edges removed from the graph
   *
   * @pre (this->has_edge(_e_.node1(), _e_.node2()))
   *
   * @post (!this->has_edge(_e_.node1(), _e_.node2()))
   *       new this->num_edges = old this->num_edges - 1
   *       EdgeIterator eit is invalidated if (*eit).index == _e_.index
   *                                         || (*eit).index >= new this->num_edges()
   *
   * Complexity O(1) armotized operations.
   *       
   */
  size_type remove_edge(const Edge& e) {
    if (!this->has_edge(e.node1(), e.node2())) return 0;
    // Find the indices of the nodes attached
    size_type idx1 = e.node1().index();
    size_type idx2 = e.node2().index();
    size_type uid1 = i2u_[idx1];
    size_type uid2 = i2u_[idx2];
    // Remove entries from both the adjacency and restlength maps
    adjacency_[uid1].erase(uid2);
    if (adjacency_[uid1].size() == 0) adjacency_.erase(uid1);
    adjacency_[uid2].erase(uid1);
    if (adjacency_[uid2].size() == 0) adjacency_.erase(uid2);
    if (uid1 < uid2) {
      restlength_[uid1].erase(uid2);
      if (restlength_[uid1].size() == 0) restlength_.erase(uid1);
    }
    else if (uid2 < uid1) {
      restlength_[uid2].erase(uid1);
      if (restlength_[uid2].size() == 0) restlength_.erase(uid2);
    }
    else {
      std::cout << "Restlength error" << std::endl;
    }

    // deal with the mapping vector
    size_type old_index = e.index_;
    // change index of edge before swapping
    alledges_[ei2u_[numedges_ - 1]].index_ = old_index;
    std::swap(ei2u_[old_index], ei2u_[numedges_ - 1]);
    ei2u_.pop_back();
    numedges_ -= 1;
    return 1;
  }

 /**
   * Remove the edge referenced by the provided iterator from the graph
   *
   * @param [in] e_it edge iterator to edge to be removed from the graph
   * @return     fresh iterator to new edge at the index referenced by _e_it_
   *
   * @pre (this->has_edge(*_e_it_))
   *
   * @post (!this->has_edge(*_e_it_))
   *       new this->num_edges = old this->num_edges - 1
   *       EdgeIterator et is invalidated if (*et).index == (*_e_it_).index
   *                                         || (*et).index >= new this->num_edges()
   *
   * Complexity O(1) armotized operations.
   *       
   */
  EdgeIterator remove_edge(EdgeIterator e_it) {
    // Dereference iterator and then remove the edge
    Edge e = (*e_it);
    remove_edge(e);
    return e_it;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // a vector containing Point objects
  // std::vector<Point> allpoints_;
  // the size of the graph (number of nodes)
  size_type size_;
  // a vector to hold all the edges ever added
  std::vector<edgeinfo> alledges_;
  // number of Edges
  size_type numedges_;
  // adjacency map for all nodes connected by edges
  std::map<size_type, std::map<size_type, size_type>> adjacency_;
  // hold all values
  // std::vector<node_value_type> allvalues_;
  // adjacency map for all rest lengths
  std::map<size_type, std::map<size_type, edge_value_type>> restlength_;
  // nodeinfo for every node added, even if later removed
  std::vector<nodeinfo> nodes_;
  // store the currently active set of nodes
  std::vector<size_type> i2u_;
  // store the currently active set of edges
  std::vector<size_type> ei2u_;
};

#endif // CME212_GRAPH_HPP
