#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * Note: discussed graph design (node, edges) with Lauren Pendo
 */

#include <algorithm>
#include <functional>
#include <map>
#include <set>
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

template  <typename V, typename E>
class Graph {
 private:

  // create an internal node with value, position, and index
  struct internal_node;
  // create an internal edge with value, node 1 and node 2, and index
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Node value type */
  using  node_value_type = V;

  /** Edge value type */
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
  Graph() 
    : nodes_(), edges_(), node_map_(), i2u_nodes_(), i2u_edges_() {
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
     Point& position() {
      assert(valid());
      return graph_->nodes_.at(uid_).point;
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(valid());
      return graph_->nodes_.at(uid_).point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(valid());
      return graph_->nodes_.at(uid_).idx;
    }

    /* Return the value of the graph; this value can be changed here */
    node_value_type& value() {
      assert(valid());
      return graph_->nodes_.at(uid_).value;
    }

    /* Return a read-only value of the node */
    const node_value_type& value() const {
      assert(valid());
      return graph_->nodes_.at(uid_).value;
    }

    /* Return the number of nodes that the current node is connected to */
    size_type degree() const {
      assert(valid());
      return graph_->node_map_.at(uid_).size();
    }
    
    /* Return the beginning of the incident interator */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, graph_->node_map_.at(uid_).begin(), uid_);
    }

    /* Return the end of the incident iterator */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, graph_->node_map_.at(uid_).end(), uid_);
    }
  
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return uid_ == n.uid_ && graph_ == n.graph_;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * compare graphs using std::less template
     * return some comparison if graphs are not equal, otherwise return
     * the less than indices
     */
    bool operator<(const Node& n) const {
      if (graph_ != n.graph_) return std::less<Graph*>()(graph_, n.graph_);
      else return uid_ < n.uid_; 
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    // Node members a pointer to the actual graph and the index of the node
    Graph* graph_;
    size_type uid_;

    /** Test whether this node is valid */
    bool valid() const
    {
      // uid in range, idx in range, uid in sync
      return uid_>= 0 && uid_< graph_->nodes_.size() 
        && graph_->nodes_.at(uid_).idx < graph_->i2u_nodes_.size() 
        && graph_->i2u_nodes_.at(graph_->nodes_.at(uid_).idx) == uid_; 
    }

    // create a valid constructor of the node
    Node(const Graph* graph, size_type uid) 
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
      }

    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_nodes_.size();
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
  Node add_node(const Point& point, 
    const node_value_type& value = node_value_type()) {

    // get the size of the overall nodes_, which gives it a unique id
    // get the size of the overall graph, which is the current index
    size_type uid = nodes_.size();
    size_type idx = num_nodes();

    // add the point position, value, and index to nodes_
    internal_node nd(point, value, idx);
    nodes_.push_back(nd);

    // add the uid to i2u vector
    i2u_nodes_.push_back(uid);

    // insert the node into the incident iterator
    std::vector<std::vector<size_type>> vec;
    node_map_.insert(std::pair<size_type, 
      std::vector<std::vector<size_type>>>(uid, vec));
    
    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // get the index of the node
    size_type idx = nodes_.at(n.uid_).idx;

    // check if the node is in the graph, has the right size
    // and its uid matches the index
    return this == n.graph_ 
      && idx < i2u_nodes_.size() 
      && n.uid_ == i2u_nodes_.at(idx);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // return node with pointer to a graph and index
    // get the uid from the i2u vector using the index
    assert(i < size());
    return Node(this, i2u_nodes_.at(i));
  }

  /** Remove node @n from the graph
   * @return 1 if node is successfully removed
   */
  size_type remove_node(const Node& n) {
    // if the node is invalid, return 0 (node was not ermoved)
    if (!has_node(n)) return 0;

    // edge removal for all nodes incident to it
    auto it = n.edge_begin();
    while (it != n.edge_end()) {
      remove_edge((*it));
    }

    // erase the instance of uid in the map
    node_map_.erase(n.uid_);

    // pop and swap on the i2u nodes
    size_type cur_idx = nodes_.at(n.uid_).idx;
    size_type last_uid = i2u_nodes_.back(); // pop
    i2u_nodes_.at(cur_idx) = last_uid; // set i2u[idx] as the last uid (swap)
    nodes_.at(last_uid).idx = cur_idx; // change the idx of the last uid (swap)
    i2u_nodes_.pop_back(); // pop

    return 1;
  }

  /** Remove node @n from the graph
   * @return the current iterator
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      assert(valid());
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(valid());
      return Node(graph_, node2_uid_);
    }

    /* Return the value of the edge; this value can be changed here */
    edge_value_type& value() {
      assert(valid());
      return graph_->edges_.at(edge_uid_).value;
    }

    /* Return a read-only value of the edge */
    const edge_value_type& value() const {
      assert(valid());
      return graph_->edges.at(edge_uid_).value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(valid());
      return edge_uid_ == e.edge_uid_ and graph_ == e.graph_;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * compare graphs using std::less template
     */
    bool operator<(const Edge& e) const {
      assert(valid());
      if (graph_ != e.graph_) return std::less<Graph*>()(graph_, e.graph_);
      else return edge_uid_ < e.edge_uid_;

      return false;
    }

    
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // graph member types like graph pointer and indices to nodes and edge uid
    Graph* graph_;
    size_type edge_uid_;
    size_type node1_uid_;
    size_type node2_uid_;

    // Edge constructor
    Edge(const Graph* graph, 
      size_type edge_uid_, 
      size_type node1_uid, 
      size_type node2_uid) 
      : graph_(const_cast<Graph*>(graph)), 
      edge_uid_(edge_uid_), 
      node1_uid_(node1_uid), 
      node2_uid_(node2_uid) {
    }

    /** Test whether this edge is valid */
    bool valid() const
    {
      // uid in range, idx in range, uid in sync.
      return edge_uid_>= 0 && edge_uid_< graph_->edges_.size() 
        && graph_->edges_.at(edge_uid_).idx < graph_->i2u_edges_.size() 
        && graph_->i2u_edges_.at(graph_->edges_.at(edge_uid_).idx) == edge_uid_; 
  }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return i2u_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // check if index exists 
    assert(i < num_edges());

    // get uid of the edge
    size_type uid = i2u_edges_.at(i);
    // translate the set type into a vector for easy indexing
    std::vector<size_type> 
      e(edges_.at(uid).nodes.begin(), edges_.at(uid).nodes.end());
    return Edge(this, i, e[0], e[1]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check if a vector exists at all in this map index
    if (!has_node(a) || !has_node(b)) return false;

    size_type uid_a = i2u_nodes_.at(a.index());
    size_type uid_b = i2u_nodes_.at(b.index());

    // check if uid a is inside the node map
    if (node_map_.find(uid_a) == node_map_.end()) return false;

    // if node b exists in the vector of node a, return true
    // adjacency list stores vectors <edge e, node b>, so index into 1 to
    // check if b exists
    std::vector<std::vector<size_type>> np = node_map_.at(uid_a);

    for (std::size_t i = 0; i < np.size(); i++) {
      if (np.at(i).at(1) == uid_b)  {
        return true;
      }
    }

    // otherwise return false
    return false;
  }


  /** Add an edge to the node adjacency list
   * @pre @a a and @a b are distinct valid uid of this graph
   * @pre @uid is the current index of the of the edge
   * @post add a node to the node adjacency list by pushing a vector of 
   * <edge e, node b> into the vector at node a
   */
  void add_edge_helper(const size_type uid, 
    const size_type uid_a, 
    const size_type uid_b) {
    // push a vector of <edge e, node b> to node a
    std::vector<size_type> node_edge {uid, uid_b};
    node_map_.at(uid_a).push_back(node_edge);
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
  Edge add_edge(const Node& a, 
    const Node& b, 
    const edge_value_type& value = edge_value_type()) {

    size_type uid_a = i2u_nodes_.at(a.index());
    size_type uid_b = i2u_nodes_.at(b.index());

    // if there exists an edge, return an edge
    // iterate through node_map_ at a to get indices of already existing edges
    if (has_edge(a, b)) {
      size_type uid_edge = 0;
      std::vector<std::vector<size_type>> np = node_map_.at(uid_a);

      for (std::size_t i = 0; i < np.size(); i++) {
        if (np.at(i).at(1) == uid_b) uid_edge = np.at(i).at(0);
      }
      return Edge(this, uid_edge, uid_a, uid_b);
    }
    // create a new edge
    else {
      // get the index of the edge
      size_type uid = edges_.size();
      size_type idx = num_edges();

      // add edges to the adjacency list
      // this repetition is to allow for easy look up of nodes from 
      // a -> b and b -> a
      add_edge_helper(uid, uid_a, uid_b);
      add_edge_helper(uid, uid_b, uid_a);

      // push back a and b nodes as a set
      std::set<size_type> nodes = {uid_a, uid_b};
      internal_edge e(nodes, value, idx);
      edges_.push_back(e);

      // add uid to the edge
      i2u_edges_.push_back(uid);

      return Edge(this, uid, uid_a, uid_b);
    }
  }


  /* Remove all adjacent nodes from the Graph given node 1 and 2 */
  void adjacency_removal(const Node& n1, const Node& n2) {
    auto it = node_map_.at(n1.uid_).begin();
    while (it != node_map_.at(n1.uid_).end()) {
      if ((*it).at(1) == n2.uid_) {
        it = node_map_.at(n1.uid_).erase(it);
        break;
      }
      else {
        ++it;
      }
    }
  }

  /* Swap the edge with the last value, update the index */
  void swap_pop_edge(const size_type idx) {
    // set value of current edges to the last value
    size_type last_uid = i2u_edges_.back();
    i2u_edges_.at(idx) = last_uid;
    edges_.at(last_uid).idx = idx;
    // set the idx value of the struct to the current index
    i2u_edges_.pop_back();
  }


  /* Remove edge from the graph
   * @ pre e edge to remove
    * if nodes or edge do not exist, return 0
    Remove from the adjacency list
   */
  size_type remove_edge(const Edge& e) {
    Node a = Node(this, e.node1_uid_);
    Node b = Node(this, e.node2_uid_);

    if (!has_node(a) or !has_node(b)) return 0;
    if (!has_edge(a, b)) return 0;

    size_type uid = e.edge_uid_;
    size_type cur_idx = edges_.at(uid).idx;


    // remove from the adjacency list
    Node n1 = e.node1();
    Node n2 = e.node2();
    adjacency_removal(n1, n2);
    adjacency_removal(n2, n1);

    swap_pop_edge(cur_idx);

    return 1;
  }

  /* Remove edge from the graph
   * given two nodes
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // exit if edge from node 1 to node 2 does not exist
    if (!has_node(n1) || !has_node(n2) || !has_edge(n1, n2)) return 0;

    // get uid of the edge
    size_type uid;
    std::vector<std::vector<size_type>> np = node_map_.at(n1.uid_);

    for (std::size_t i = 0; i < np.size(); i++) {
      if (np.at(i).at(1) == n2.uid_) uid = np.at(i).at(0);
    }

    return remove_edge(Edge(this, uid, n1.uid_, n2.uid_));
  }

  /* Remove edge from the graph
   * given edge iterator
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = (*e_it);
    remove_edge(e);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   * Use stl functions to clear everything in graph.
   */
  void clear() {
    i2u_nodes_.clear();
    i2u_edges_.clear();
    node_map_.clear();
    nodes_.clear();
    edges_.clear();
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
    
   /** Increase node index and return a reference to the current object */
    NodeIterator& operator++() {
      ++idx_;
      return (*this);
    }

    /* Check whether the current iterator is equal to another iterator */ 
    bool operator==(const NodeIterator& node_iter) const {
      return idx_ == node_iter.idx_ and graph_ == node_iter.graph_;
    }

    /* Check whether the current iterator is not equal to another iterator */
    bool operator!=(const NodeIterator& node_iter) const {
      return !(*this == node_iter);
    }

    /* Return a node pointing to the graph and the current index */
    Node operator*() const {
      return Node(graph_, graph_->i2u_nodes_.at(idx_));
    }
  
   private:
    friend class Graph;
    Graph* graph_;
    size_type idx_;

    // constructor takes in graph and index
    NodeIterator(const Graph* graph, size_type idx) : 
      graph_(const_cast<Graph*>(graph)), idx_(idx) {
      }
  };

  /* Return a node iterator object at index 0 to begin the iterator */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /* Return a node iterator object at index node size to end the iterator */
  NodeIterator node_end() const {
    return NodeIterator(this, i2u_nodes_.size());
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
    
    /* Increment the incident iterator by incrementing the vector iterator*/
    IncidentIterator& operator++() {
      ++incident_iter_;
      return (*this);
    }

    /* check if the current iterator is equal to another iterator */
    bool operator==(const IncidentIterator& incident_iter) const {
      return incident_iter_ == incident_iter.incident_iter_;
    }

    /* check if the current iterator is not equal to another iterator */
    bool operator!=(const IncidentIterator&incident_iter) const {
      return incident_iter_ != incident_iter.incident_iter_;
    }

    /* return the edge object, where incident[0] is the index of the edge,
    * idx_ is the index of node a, and incident[1] is node b */
    Edge operator*() const {
      std::vector<size_type> incident = *incident_iter_;
      return Edge(graph_, incident[0], idx_, incident[1]);
    }
  
   private:
    friend class Graph;
    
    // track graph, incident iterator over a vector, and index of the node
    Graph* graph_;
    std::vector<std::vector<size_type>>::iterator incident_iter_;
    size_type idx_;

    // constructor for incident iterator
    IncidentIterator(
      const Graph* graph, 
      std::vector<std::vector<size_type>>::iterator incident_iter,
      size_type idx
    ) : 
      graph_(const_cast<Graph*>(graph)), 
        incident_iter_(incident_iter), 
        idx_(idx) {
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

    /* Increment the vector index */
    EdgeIterator& operator++() {
      ++idx_;
      return (*this);
    }

    /* check if the current edge iterator is equal to another edge iterator */
    bool operator==(const EdgeIterator& edge_iter) const {
      return idx_ == edge_iter.idx_ and graph_ == edge_iter.graph_;
    }

    /*check if the current edge iterator is not equal to another edge iterator*/
    bool operator!=(const EdgeIterator& edge_iter) const {
      return !(*this == edge_iter);
    }

    /* Return the current edge index the iterator is on */
    Edge operator*() const {
      // convert set into vector for indexing edge 0 and edge 1
      size_type uid = graph_->i2u_edges_.at(idx_);
      std::vector<size_type> 
        e(graph_->edges_.at(uid).nodes.begin(), graph_->edges_.at(uid).nodes.end());
      return Edge(graph_, uid, e[0], e[1]);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type idx_;

    // constructor for edge iterator
    EdgeIterator(
      const Graph* graph, 
      size_type idx
    ) : 
      graph_(const_cast<Graph*>(graph)), idx_(idx) {
      }
  };

  /* Return an edge iterator object at index 0 to begin the iterator */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /* Return an edge iterator object at index edge size to begin the iterator */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, i2u_edges_.size());
  }

 private:

  struct internal_node {
    Point point;
    node_value_type value;
    size_type idx;

    internal_node(Point point, node_value_type value, size_type idx): 
      point(point), value(value), idx(idx) {}
  };

  struct internal_edge {
    std::set<size_type> nodes;
    edge_value_type value;
    size_type idx;

    internal_edge(std::set<size_type> nodes, edge_value_type value, size_type idx) : 
      nodes(nodes), value(value), idx(idx) {}
  };

  // get a vector of points that are nodes
  std::vector<internal_node> nodes_;
    // get a mapping of node sets to its indices
  std::vector<internal_edge> edges_;


  // create a map between a node and a vector of <node 2, edge> pairs that
  // connect between a and b
  std::map<size_type, std::vector<std::vector<size_type>>> node_map_;

  // node size
  // size_type node_size_;

  // a vector of indices to uids for nodes
  // this is used for node deletion
  std::vector<size_type> i2u_nodes_;
  // a vector of indices to uids for edges
  // this is used for edge deletion
  std::vector<size_type> i2u_edges_;

  // edge size
  // size_type edge_size_;

};

#endif // CME212_GRAPH_HPP
