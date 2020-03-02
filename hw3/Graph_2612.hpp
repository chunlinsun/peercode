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
template <typename V, typename E>
class Graph {

  private:
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

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
  /** Synonym for templated node value V . */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Synonym for templated edge value E . */
  using edge_value_type = E;

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
  Graph()
    : node_size_(0), edge_size_(0) {
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
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(valid());
      return graph_->node_elements_[uid_]->position;
    }

    /** Return a reference of this node's position. */
    Point& position() {
      assert(valid());
      return graph_->node_elements_[uid_]->position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(valid());
      return graph_->node_elements_[uid_]->idx_;
    }

    /** Return this node's value. */
    node_value_type& value() {
      assert(valid());
      return graph_->node_elements_[uid_]->value;
    };

    /** Return this node's value as const. */
    const node_value_type& value() const {
      assert(valid());
      return graph_->node_elements_[uid_]->value;
    }

    /* Return the number of incident edges of this node. */
    size_type degree() const {
      assert(valid());
      return graph_->nodes_degrees_[uid_];
    }

    /* Returns an iterator to the beginning of incident edges/ */
    incident_iterator edge_begin() const {
      assert(valid());
      return IncidentIterator(graph_, uid_, 0);
    }

    /* Returns an iterator to the end of incident edges. */
    incident_iterator edge_end() const {
      assert(valid());
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      assert(valid() && n.valid());
      return n.graph_ == graph_ && n.uid_ == uid_;
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
      assert(valid() && n.valid());
      if(n.graph_ == graph_ && n.uid_ == uid_){
        return false;
      }
      if(graph_ < n.graph_ || (graph_ == n.graph_ && uid_ < n.uid_)){
        return true;
      }
      return false;
    }

   private:
     // returns if this node is valid with representation invariance
     bool valid() const {
       return uid_ >= 0 && uid_<graph_->node_elements_.size()
              && graph_->node_elements_[uid_]->idx_ < graph_->i2u_node_.size()
              && graph_->i2u_node_[graph_->node_elements_[uid_]->idx_] == uid_;
     }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type uid_;    // unique identifier for node information
    Node(const graph_type* graph, size_type uid):
                graph_(const_cast<Graph*>(graph)), uid_(uid){}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node, also support node value.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    internal_node* new_node = new internal_node();
    new_node->position = position;
    new_node->idx_ = node_size_;
    new_node->value = val;
    node_elements_.push_back(new_node);
    node_size_ ++;

    std::vector<size_type> incident_e;
    incident_edges_.push_back(incident_e);
    nodes_degrees_.push_back(0);
    i2u_node_.push_back(node_elements_.size() - 1);

    return Node(this, node_size_ - 1);        // return new node
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
    assert(0 <= i and i < size());   //
    return Node(this, i2u_node_[i]);    // return node with index i
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: at most O(num_nodes()).
   */
  size_type remove_node(const Node& n) {
    (void) n;
    return 0;
  }

  /** Remove the node pointed to by @a n_it.
   *
   * Complexity: at most O(num_nodes()).
   */
  node_iterator remove_node(node_iterator n_it) {
    (void) n_it;
    return node_begin();
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, uid2_);
    }

    /** Return a reference of this edge's value. */
    edge_value_type& value() {
      return graph_->edge_elements_[uid_]->value;
    };

    /** Return this edge's value as const. */
    const edge_value_type& value() const {
      return graph_->edge_elements_[uid_]->value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return e.graph_ == graph_ && e.uid_ == uid_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(this->node1() == e.node1() && this->node2() == e.node2() ){
        return false;
      }
      if(this->node1() < e.node1() || (this->node1() == e.node1() && this->node2() < e.node2())){
        return true;
      }
      return false;
    }

   private:

    friend class Graph;

    graph_type* graph_;
    size_type uid_;    // unique id for this edge
    size_type uid1_;   // unique id for node 1
    size_type uid2_;   // unique id for node 2

    // construct a valid edge
    Edge(const graph_type* graph, size_type uid, size_type uid1, size_type uid2):
                graph_(const_cast<Graph*>(graph)), uid_(uid), uid1_(uid1), uid2_(uid2){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >=0 && i < num_edges());
    size_type uniq = i2u_edge_[i];
    size_type uniq_n1 = i2u_node_[edge_elements_[uniq]->idx1_];
    size_type uniq_n2 = i2u_node_[edge_elements_[uniq]->idx2_];
    return Edge(this, uniq, uniq_n1, uniq_n2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    int res = find_edge(a, b);
    return !(res == -1);
  }


  /** Find the edge between two nodes a and b, or returns -1 if the edge does not exist
   * @pre @a a and @a b are valid nodes of this graph
   * @return the unique id of the edge if the edge exists, -1 otherwise
   *
   * Complexity: No more than O(num_edges())
   */
  int find_edge(const Node& a, const Node& b) const {
    if(!has_node(a) || !has_node(b) || a == b || !a.valid() || !b.valid()){
      return -1;
    }
    size_type uniq_a = i2u_node_[a.index()];
    std::vector<size_type> a_incident = incident_edges_[uniq_a];   // vector of incident edges
    for(size_type i = 0; i < a_incident.size(); i++){
      // check incident edge edge_elements_[a_incident[i]]
      if(edge_elements_[a_incident[i]]->idx2_ == b.index() || edge_elements_[a_incident[i]]->idx1_ == b.index()) {
        return a_incident[i];
      }
    }
    return -1;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {


    int res = find_edge(a, b);

    if (res != -1) {
      // the edge exists already
      return Edge(this, res, a.uid_, b.uid_);
    }

    // not found, add edge
    internal_edge* new_edge = new internal_edge();
    new_edge->idx1_ = a.index();
    new_edge->idx2_ = b.index();
    new_edge->idx_ = edge_size_;
    new_edge->value = val;
    edge_elements_.push_back(new_edge);
    edge_size_ ++;

    // update incident edges
    size_type uniq_e = edge_elements_.size() - 1;
    incident_edges_[a.uid_].push_back(uniq_e);
    incident_edges_[b.uid_].push_back(uniq_e);
    nodes_degrees_[a.uid_]++;
    nodes_degrees_[b.uid_]++;
    i2u_edge_.push_back(edge_elements_.size() - 1);
    return Edge(this, uniq_e, a.uid_, b.uid_);
  }


  /** Remove an edge from the graph
   * @pre @a a and @a b are valid nodes of this graph
   * @pre @a a and @a b have a valid edge of this graph
   * @return 1 if remove edge successully, 0 if remove edge failed
   *
   * @post new num_edges() == old num_edges() - 1;
   * Complexity: at most O(num nodes() + num edges())
   */
  size_type remove_edge(const Node& a, const Node& b) {
    int uniq_e = find_edge(a, b);

    if(uniq_e != -1){
      size_type index_e = edge_elements_[uniq_e]->idx_;
      // remove edge from adj_ mapping
      size_type index_inc = 0;
      // remove from a
      std::vector<size_type> a_incident = incident_edges_[a.uid_];
      for(auto it = a_incident.begin(); it != a_incident.end(); ++it){
        if(*it == (unsigned)uniq_e){
          index_inc = it - a_incident.begin();
          break;
        }
      }
      incident_edges_[a.uid_].erase(index_inc + incident_edges_[a.uid_].begin());
      std::vector<size_type> b_incident = incident_edges_[b.uid_];
      for(auto it = b_incident.begin(); it != b_incident.end(); ++it){
        if(*it == (unsigned)uniq_e){
          index_inc = it - b_incident.begin();
          break;
        }
      }
      incident_edges_[b.uid_].erase(index_inc + incident_edges_[b.uid_].begin());
      // manually change index of following edges
      for(auto it = edge_elements_.begin(); it != edge_elements_.end(); ++it){
          if((*it)->idx_ > (unsigned)index_e){(*it)->idx_ -= 1;}
      }
      i2u_edge_.erase(i2u_edge_.begin() + index_e);  // remove from edge mapping
      edge_size_--;
      return 1;
    }
    // did not find edge, return failure
    return 0;

  }

  /** Remove an edge @a e from the graph
   * @pre @a e is a valid edge of this graph
   * @return 1 if remove edge successully, 0 if remove edge failed
   *
   * @post new num_edges() == old num_edges() - 1;
   * Complexity: at most O(num nodes() + num edges())
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /**
   *
   * Complexity: at most O(num nodes() + num edges())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge((*e_it).node1(), (*e_it).node2());
    return edge_begin();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_size_ = 0;
    edge_size_ = 0;
    size_type all_nodes = node_elements_.size();
    for(size_type i = 0; i < all_nodes; i++){
      delete node_elements_[i];
    }
    node_elements_.clear();

    size_type all_edges = edge_elements_.size();
    for(size_type i = 0; i < all_edges; i++){
      delete edge_elements_[i];
    }
    edge_elements_.clear();

    incident_edges_.clear();
    nodes_degrees_.clear();
    i2u_edge_.clear();
    i2u_node_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : totally_ordered<NodeIterator>{
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

    /** deference a node iterator
     * @return a node the iterator is pointing to
     */
    Node operator*() const {
      // return node object with uniq id
      return Node(graph_, graph_->i2u_node_[ni_id_]);
    }

    /** forward increment of the NodeIterator
     * @return the next NodeIterator
     */
    NodeIterator& operator++() {
      ni_id_++;
      return *this;
    }

    /** test if this iterator is equal to @a ni.
     *
     * Equal NodeIterator points to the same node.
     */
    bool operator==(const NodeIterator& ni) const {
      return *(*this) == *ni;
    }

   private:
    friend class Graph;

    graph_type* graph_;
    size_type ni_id_;    // index of node its pointing to

    /* NodeIterator constructor */
    NodeIterator(const graph_type* g, size_type i) : graph_(const_cast<graph_type*>(g)), ni_id_(i) {
    }

  };

  /* return the first node iterator. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /* return the iterator pointing to one past the last node.  */
  node_iterator node_end() const {
    return NodeIterator(this, node_size_);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator :totally_ordered<IncidentIterator> {
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

    /* Deference an IncidentIterator.
     * @return the Edge the iterator is pointing to.
     */
    Edge operator*() const {
      // Edge should know its own node1 and node2 （which may not be represented in graph data）
      // should not modify graph data. Consider the case where two ii are instantiated.
      size_type id;
      auto e = graph_->edge_elements_[graph_->incident_edges_[node_uid_][incident_index_]];   // an undirected edge in graph
      if (graph_->i2u_node_[e->idx1_] == node_uid_) {
        id = graph_->i2u_node_[e->idx2_];
      }else{
        id = graph_->i2u_node_[e->idx1_];
      }
      return Edge(graph_, graph_->incident_edges_[node_uid_][incident_index_],
                    node_uid_, id);
    }

    /* Forward increment the incident iterator. */
    IncidentIterator& operator++() {
      incident_index_++;
      return *this;
    }

    /* Check if the given incident iterator is equal to this incident iterator.
     * Equal incident iterator belongs to the same graph, the same node and points to the same edge.
     */
    bool operator==(const IncidentIterator& ii) const {
      return (graph_ == ii.graph_) && (node_uid_ == ii.node_uid_) && (*(*this) == *ii);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type node_uid_;    // unique id for the parent node
    size_type incident_index_;   // index among all incident edges from the same parent

    IncidentIterator(const graph_type* g, size_type ni, size_type ii):
        graph_(const_cast<graph_type*>(g)), node_uid_(ni), incident_index_(ii){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /* Deference edge iterator. Returns the edge its pointing to. */
    Edge operator*() const {
      size_type uniq_e = graph_->i2u_edge_[idx_];
      return Edge(graph_,
                  graph_->i2u_edge_[idx_],
                  graph_->i2u_node_[graph_->edge_elements_[uniq_e]->idx1_],
                  graph_->i2u_node_[graph_->edge_elements_[uniq_e]->idx2_]);
    }

    /* Forward increment the iterator. */
    EdgeIterator& operator++() {
      idx_++;
      return *this;
    }

    /* Check if the iterator is equal to a given edge iterator.
     * Equal edge iterators point to the same edge.
     */
    bool operator==(const EdgeIterator& ei) const {
      return *(*this) == *ei;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type idx_;      // index of the edge
    EdgeIterator(const graph_type* g, size_type i) : graph_(const_cast<graph_type*>(g)), idx_(i) {}
  };

  /* Returns an iterator to the beginning. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /* Returns an iterator to the end. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_size_);
  }

 private:

   /*
   *  Node informations
   */
  struct internal_node {
    Point position;   // The text held by an element
    Point::size_type idx_;      // the user facing index of a node
    node_value_type value;
  };

  std::vector<internal_node*> node_elements_;
  size_type node_size_;
  std::vector<size_type> i2u_node_;    // indexed by node index, ith element is the unique id for node (i)


  /*
  *  Edge informations
  */
  struct internal_edge {        // undirected edges
    Point::size_type idx1_;      // The index for node 1
    Point::size_type idx2_;      // The index for node 2
    size_type idx_;            // the user facing index of an edge
    edge_value_type value;
  };
  std::vector<internal_edge*> edge_elements_;     // a vector to store the undirected edges
  size_type edge_size_;
  std::vector<size_type> i2u_edge_;    // indexed by edge index, ith element is the unique id for edge (i)

  // incident_edges_.size() == node_elements_.size()
  // incident_edges_[i] is a list of uniq id of edges incident to uniq id node i
  std::vector<std::vector<size_type>> incident_edges_;

  // nodes_degrees_.size() == node_elements_.size()
  std::vector<size_type> nodes_degrees_;

};

#endif // CME212_GRAPH_HPP
