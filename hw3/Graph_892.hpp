#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <algorithm>    // std::swap
#include <vector>
#include <cassert>
#include <unordered_map>


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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type for class template */
  using node_value_type = V;
  using edge_value_type = E;

  /** Type of this graph. */
  using graph_type = Graph<V,E>;  //modified

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


 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;
  struct internal_edge;

  std::vector<internal_node> nodes_;  // node proxy
  std::vector<internal_edge> edges_;  // edge proxy
 
  // Store the currently "active" set of nodes and edges.
  std::vector<size_type> i2u_;   // Indexed by node idx
  std::vector<size_type> i2e_;   // Indexed by edge idx



 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() :  nodes_(),edges_(), i2u_(), i2e_() { //, neighbors_() {
    // HW0: YOUR CODE HERE
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
      assert(valid());
      return graph_->nodes_[nid_].pos;  
    }

    /** Return this node's modifiable position. */
    Point& position() {
      assert(valid());
      return graph_->nodes_[nid_].pos;  
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(valid());
      return graph_->nodes_[nid_].idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for: 

    /**
     * @brief Access the user-specified value stored supported by Node
     *
     * @return Reference to the user-specified value
     */
    node_value_type& value() {
      assert(valid());
      return graph_->nodes_[nid_].val;
    }

    /**
     * @brief Access the user-specified value stored supported by Node
     *
     * @return Const reference to the user-specified value
     * 
     * Promising not to modify the value
     */
    const node_value_type& value() const {
      assert(valid());
      return graph_->nodes_[nid_].val;
    }

    /**
     * @brief Get the degree of this Node
     *
     * @return Degree of this Node
     */
    size_type degree() const{
      assert(valid());
      return graph_->nodes_[nid_].adjkey_.size();
    }

    /**
     * @brief Obtain start of the iterator over all edges incident to this Node
     
     * @return Start of the incident iterator
     *
     * Complexity: O(1).
     */
    incident_iterator edge_begin() const{
      assert(valid());
      return IncidentIterator(graph_, nid_, 0);
    }

    /**
     * @brief Obtain end of the iterator over all edges incident to this Node
     *
     * @return End of the incident iterator
     *
     * Complexity: O(1). 
     * Iterating over all the edges incident to a Node n has maximum complexity O(n.degree())
     */
    incident_iterator edge_end() const{
      assert(valid());
      return IncidentIterator(graph_, nid_, graph_->nodes_[nid_].adjkey_.size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      bool eq_id = nid_==n.nid_;
      bool eq_graph = graph_==n.graph_;
      return (eq_id && eq_graph);
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
      return ((graph_==n.graph_) && (nid_<n.nid_)) || (graph_<n.graph_);
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE - modified
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph 
    Graph* graph_;
    // Node index
    size_type nid_; 
    /** Private Constrctor */
    Node(const graph_type* graph, size_type nid)
     : graph_(const_cast<graph_type*>(graph)), nid_(nid) {  
    }


    /** Check if the node is valid.
     * 
     * @return true if the unique ID of node is in range, 
     *                 the user-seen ID is in range, and
     *                 the unique ID is in sync;
     *         false otherwise.
     *
     * Complexity: O(1) 
     */
    bool valid() const {
      return nid_>= 0 && nid_< graph_->nodes_.size()         // uid in range.
             && graph_->nodes_[nid_].idx_ < graph_->i2u_.size()  // idx in range.
             && graph_->i2u_[graph_->nodes_[nid_].idx_] == nid_; // uid in sync.
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
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

    size_type new_id = num_nodes();  //user-seen id
    size_type new_nid = nodes_.size();  //unique id

    // create an internal node
    internal_node new_node(position, new_id, value, std::unordered_map<size_type, size_type>(), std::vector<size_type>());
    nodes_.push_back(new_node);  // O(1)
    i2u_.push_back(new_nid);  // O(1)

    return Node(this,new_nid);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE - modified
    //check if node index within range 
    if (num_nodes()<=n.index()) {  
      return false;
    }
    internal_node node_found = nodes_[i2u_[n.index()]];
    return (node_found.pos==n.position()) && (node_found.idx_==n.index());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {

    assert(i<num_nodes()); // index in range

    return Node(this,i2u_[i]);        
  }


  /** Remove a node (and all indicent edges) from the graph, returning the number of nodes removed.
   * @param[in] n             Reference of the Node to remove
   * @return    num_removed   Number of nodes removed.
   * 
   * @post if result==1, new num_nodes()==old num_nodes()-1, new num_edges()==old num_edges()-n.degree(),
   *                     Node _n_ can no longer be seen by the user, and
   *                     old _n_.index() was now the index of the previous last node in graph
   * @post if result==0, no nodes removed, the graph stays the same
   *
   * Complexity: O(num_nodes()) average (for sparse graph)
   */
  size_type remove_node ( const Node & n) {

    unsigned num_removed= 0;

    if (has_node(n)) {

      size_type n_idx = n.index();
      // remove incident edges, take O(num_nodes()) on average
      while (nodes_[i2u_[n_idx]].adjkey_.size()>0) {

          size_type neighid = nodes_[i2u_[n_idx]].adjkey_[0];

          size_type eid = nodes_[i2u_[n_idx]].adj_[nodes_[i2u_[n_idx]].adjkey_[0]];
          remove_edge(Edge(this,eid,i2u_[n_idx],neighid)); //each edge removal takes O(|E|)
      }
      
      if (i2u_.size()>1){ // more than one node left
        size_type original_last_id = i2u_.back();
        size_type original_last_index = nodes_[original_last_id].idx_;
               
        // update all its adjacency edges
        for (auto it=node(original_last_index).edge_begin(); it!=node(original_last_index).edge_end();++it){
          if (edges_[i2e_[(*it).index()]].nid1==original_last_index){
            edges_[i2e_[(*it).index()]].nid1 = n_idx;
          }else{
            edges_[i2e_[(*it).index()]].nid2 = n_idx;
          }
        }

        // update the internally-stored index of the previous last node in graph
        nodes_[original_last_id].idx_ = n_idx; 

        // swap the node-to-remove and the previous last node in index list
        i2u_[n_idx]=original_last_id; 
      } 
      // then pop_back: take O(1)
      i2u_.pop_back();
      
      // increased the number of node removed
      ++num_removed;
    }

    return num_removed;

  }

  /** Remove a node (and all edges incident to it) from the graph, 
   * returning the NodeIterator of original address of the removed node
   * @param[in] n      NodeIterator of the node to remove
   * @return    n_it   NodeIterator of original address of the removed node
   *
   * Complexity: O(num_nodes()) average
   */
  node_iterator remove_node ( node_iterator n_it ){
    // call remove_node(const Node &), take O(num_nodes()) average
    node_iterator origin = n_it;
    remove_node(*n_it);  
    return origin;
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
      return Node(graph_,nid1_); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(valid());
      return Node(graph_,nid2_); 
    }

    /**
     * @brief Get the length of Edge
     *
     * @return Length of the edge
     *
     * @post 0<=result, result=norm(distance between two endpoints (Node) of this Edge)
     *
     * The complexity of is O(1)
     */
    double length() const {
      assert(valid());
      return norm(node1().position()-node2().position());
    }

    /**
     * @brief Access the user-specified value stored supported by Edge
     *
     * @return Reference to the user-specified value
     */
    edge_value_type& value() {
      assert(valid());
      return graph_->edges_[eid_].val;
    }

    /**
     * @brief Access the user-specified value stored supported by Edge
     *
     * @return Const reference to the user-specified value
     * 
     * Promising not to modify the value
     */
    const edge_value_type& value() const {
      assert(valid());
      return graph_->edges_[eid_].val;
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     *
     * Complexity: O(1).
     */
    bool operator==(const Edge& e) const {
      bool eq_graph = graph_==e.graph_;
      bool eq_id = eid_==e.eid_;
      return (eq_graph && eq_id);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((graph_==e.graph_) && (eid_ < e.eid_)) || (graph_<e.graph_);
    }

    /** Return this edge's index. */
    size_type index() const {
      return graph_->edges_[eid_].idx_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph 
    Graph* graph_;
    // Edge index
    size_type eid_; 
    // Node indicies
    size_type nid1_;
    size_type nid2_;
    /** Private Constrctor */
    Edge(const graph_type* graph, size_type eid, size_type nid1, size_type nid2)
     : graph_(const_cast<graph_type*>(graph)), eid_(eid), nid1_(nid1), nid2_(nid2) {}

     /** Check if the edge is valid.
     * 
     * @return true if the unique ID of edge is in range, 
     *                 the user-seen ID is in range, and
     *                 the unique ID is in sync;
     *         false otherwise.
     *
     * Complexity: O(1) 
     */
    bool valid() const {
      return eid_>= 0 && eid_< graph_->edges_.size()         // eid in range.
             && graph_->edges_[eid_].idx_ < graph_->i2e_.size()  // idx in range.
             && graph_->i2e_[graph_->edges_[eid_].idx_] == eid_; // eid in sync.
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return i2e_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i<num_edges()); // the index needs to be within range
    return Edge(this, i2e_[i], i2u_[edges_[i2e_[i]].nid1], i2u_[edges_[i2e_[i]].nid2]);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1)
   */
  bool has_edge(const Node& a, const Node& b) const {
    // unordered_map::count is O(1) for simple graph
    return nodes_[i2u_[a.index()]].adj_.count(i2u_[b.index()])>0;  
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
   * Complexity: O(1)
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type e_value={}) {
    if (nodes_[i2u_[a.index()]].adj_.count(i2u_[b.index()])>0){ // unordered_map::count is O(1) for simple graph
      return Edge(this, nodes_[i2u_[a.index()]].adj_[i2u_[b.index()]], i2u_[a.index()], i2u_[b.index()]);
    }else{
      // create edge
      size_type new_id = num_edges();  // user-seen id
      size_type new_eid = edges_.size();  // unique id
      // create new internal edge
      internal_edge new_edge(a.index(), b.index(), new_id, e_value);
      edges_.push_back(new_edge);  // O(1)
      i2e_.push_back(new_eid);  // O(1)
      nodes_[i2u_[a.index()]].adjkey_.push_back(i2u_[b.index()]);
      nodes_[i2u_[a.index()]].adj_[i2u_[b.index()]] = new_eid; 
      nodes_[i2u_[b.index()]].adjkey_.push_back(i2u_[a.index()]);
      nodes_[i2u_[b.index()]].adj_[i2u_[a.index()]] = new_eid; // unordered_map::operator[] is O(1)
      return Edge(this, new_eid, i2u_[a.index()],i2u_[b.index()]);
    }

  };


  // HW2

  /** Remove an edge from the graph, returning the number of edges removed.
   * @param[in] a             Reference of the first Node of the Edge
   * @param[in] b             Reference of the second Node of the Edge
   * @return    num_removed   Number of edges removed.
   * 
   * @post if result==1, new num_edges()==old num_edges()-1,
   *                     Edge(_a_,_b_) can no longer be seen by the user.
   * @post if result==0, no edge removed, the graph stays the same
   *
   * Complexity: O(|E|) average
   */
  size_type remove_edge ( const Node & a, const Node & b) {

    size_type num_removed = 0;

    // unordered_map::count is O(1) for simple graph
    if (nodes_[i2u_[a.index()]].adj_.count(i2u_[b.index()])>0){ 

      // unordered_map::operator[] is O(1)
      size_type eid_to_remove = nodes_[i2u_[a.index()]].adj_[i2u_[b.index()]]; 

      nodes_[i2u_[a.index()]].adj_.erase(i2u_[b.index()]); 
      nodes_[i2u_[b.index()]].adj_.erase(i2u_[a.index()]); // unordered_map::erase is O(1)

      nodes_[i2u_[a.index()]].adjkey_.erase(std::find(nodes_[i2u_[a.index()]].adjkey_.begin(), nodes_[i2u_[a.index()]].adjkey_.end(), i2u_[b.index()]));
      nodes_[i2u_[b.index()]].adjkey_.erase(std::find(nodes_[i2u_[b.index()]].adjkey_.begin(), nodes_[i2u_[b.index()]].adjkey_.end(), i2u_[a.index()]));
      // vector::erase is O(1), so std::find is O(|E|) 

      if (i2e_.size()>1) {
        // get the previous last edge in graph
        size_type original_last_id = i2e_.back();
        // get the user-seen id of the edge-to-remove
        size_type id_to_remove = edges_[eid_to_remove].idx_;
        // update the internally-stored index of the previous last edge in graph
        edges_[original_last_id].idx_ = id_to_remove;
        // swap the edge-to-remove and the previous last edge in index list
        // then pop_back: take O(1)
        i2e_[id_to_remove]=original_last_id;
      }
      i2e_.pop_back();

      // increase the number of edge removed
      ++num_removed;
    }
    return num_removed; 

  }

  /** Remove an edge from the graph, returning the number of edges removed.
   * @param[in] e             Reference of the Edge to remove
   * @return    num_removed   Number of edges removed.
   * 
   * @post if result==1, new num_edges()==old num_edges()-1,
   *                     Edge _e_ can no longer be seen by the user.
   * @post if result==0, no edge removed, the graph stays the same
   *
   * Complexity: O(|E|) average
   */
  size_type remove_edge ( const Edge & e) {
    Node a = e.node1();
    Node b = e.node2();

    size_type num_removed = 0;

    // unordered_map::find is O(1)
    std::unordered_map<size_type,size_type>::const_iterator got = nodes_[i2u_[a.index()]].adj_.find(i2u_[b.index()]);
    if ( got == nodes_[i2u_[a.index()]].adj_.end() ) {
      return 0; //no edge found
    } else {

      size_type eid_to_remove = got->second;
      if (eid_to_remove==i2e_[e.index()]) {

        nodes_[i2u_[a.index()]].adj_.erase(i2u_[b.index()]); 
        nodes_[i2u_[b.index()]].adj_.erase(i2u_[a.index()]); // unordered_map::erase is O(1)

        nodes_[i2u_[a.index()]].adjkey_.erase(std::find(nodes_[i2u_[a.index()]].adjkey_.begin(), nodes_[i2u_[a.index()]].adjkey_.end(), i2u_[b.index()]));
        nodes_[i2u_[b.index()]].adjkey_.erase(std::find(nodes_[i2u_[b.index()]].adjkey_.begin(), nodes_[i2u_[b.index()]].adjkey_.end(), i2u_[a.index()]));
        // vector::erase is O(1), so std::find is O(|E|) 

        if (i2e_.size()>1){
          // get the previous last edge in graph
          size_type original_last_id = i2e_.back();
          // get the user-seen id of the edge-to-remove
          size_type id_to_remove = e.index();
          // update the internally-stored index of the previous last edge in graph
          edges_[original_last_id].idx_ = id_to_remove;
          // swap the edge-to-remove and the previous last edge in index list
          // then pop_back: take O(1)
          i2e_[id_to_remove]=original_last_id;
        }
        i2e_.pop_back();
        // increase the number of edge removed
        ++num_removed;
      }
    }

    return num_removed; 
  }

  /** Remove an edge from the graph, returning the number of edges removed.
   * @param[in] e_it   EdgeIterator of the edge to remove
   * @return    e_it   EdgeIterator pointing to original address of the removed edge
   *
   * Complexity: O(1) average
   */
  edge_iterator remove_edge ( edge_iterator e_it ){
    remove_edge(*e_it);  // call remove_edge(const Edge &), take O(|E|) average
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   * Invalidates all outstanding NodeIndex and EdgeIndex objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    i2u_.clear();
    i2e_.clear();
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

    /**
     * @brief Dereference operator of node iterator
     *
     * @return Node object that the current iterator points to
     *
     * Complexity: O(1)
     */
    Node operator*() const 
    {
      return graphPtr_->node(nodeIterIdx_);
    }

    /**
     * @brief Increment operator of node iterator
     *
     * @param[out] nodeIterIdx_ Index of node the iterator currently points to 
     * @return Reference to the current iterator, with pointer pointing to the next
     *
     * Complexity: O(1)
     */
    NodeIterator& operator++()
    {
      ++nodeIterIdx_;
      return *this;
    }

    /**
     * @brief Equality operator of node iterator

     * @param[in] node_iter Another node's iterator
     * @return True if _this_ and _node_iter_ are the same
     *         (belonging to the same graph and pointing to the same node)
     *
     * Complexity: O(1)
     */
    bool operator==(const NodeIterator& node_iter) const
    {
      return (graphPtr_==node_iter.graphPtr_)&&(nodeIterIdx_==node_iter.nodeIterIdx_);
    }

    /**
     * @brief Inequality operator of node iterator

     * @param[in] node_iter Another node's iterator
     * @return True if _this_ and _node_iter_ are not the same
     *         (not belonging to the same graph or not pointing to the same node)
     *
     * Complexity: O(1)
     */
    bool operator!=(const NodeIterator& node_iter) const
    {
      return (graphPtr_!=node_iter.graphPtr_)||(nodeIterIdx_!=node_iter.nodeIterIdx_);
    }

   private:
    friend class Graph;
    // Pointer to graph
    graph_type* graphPtr_;
    // Node index
    size_type nodeIterIdx_;
    /** Private Constrctor */
    NodeIterator(const graph_type* graph, size_type nid) : 
    graphPtr_(const_cast<graph_type*>(graph)), nodeIterIdx_(nid) {}

  };


  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:


  /**
   * @brief Get the start of the iterator over all nodes
   
   * @return Start of the node iterator
   *
   * Complexity: O(1).
   */
  node_iterator node_begin() const {
    return NodeIterator(this,0);
  }

  /**
   * @brief Get the end of the iterator over all nodes
   
   * @return End of the node iterator
   *
   * Complexity: O(1).
   */
  node_iterator node_end() const {
    return NodeIterator(this,num_nodes());
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

    /**
     * @brief Dereference operator of incident iterator
     *
     * @return Edge object that the current iterator points to
     *
     * Complexity: O(1)
     */
    Edge operator*() const {
      size_type neighbor_id = graphPtr_->nodes_[nid1_].adjkey_[adjkeyIdx_];
      return Edge(graphPtr_, graphPtr_->nodes_[nid1_].adj_[neighbor_id], nid1_, neighbor_id);
    }

    /**
     * @brief Increment operator of incident iterator
     *
     * @param[out] adjIdx_ Index of node in adjacency list that the iterator currently points to 
     * @return Reference to the current incident iterator, with pointer pointing to the next in adjacency list
     *
     * Complexity: O(1)
     */
    IncidentIterator& operator++() {
      ++adjkeyIdx_;
      return *this;
    }

    /**
     * @brief Equality operator of incident iterator

     * @param[in] incidentIter Another incident iterator
     * @return True if _this_ and _incidentIter_ are the same
     *         (belonging to the same graph, incidenting to the same node and pointing to the same neighbor)
     *
     * Complexity: O(1)
     */
    bool operator==(const IncidentIterator& incidentIter) const {
      return ((graphPtr_==incidentIter.graphPtr_) && (nid1_==incidentIter.nid1_) && (adjkeyIdx_==incidentIter.adjkeyIdx_));
    }
    
    /**
     * @brief Inequality operator of incident iterator

     * @param[in] incidentIter Another incident iterator
     * @return True if _this_ and _incidentIter_ are not the same 
     *         (not belonging to the same graph, not incidenting to the same node or not pointing to the same neighbor)
     *
     * Complexity: O(1)
     */
    bool operator!=(const IncidentIterator& incidentIter) const {
      return ((graphPtr_!=incidentIter.graphPtr_) || (nid1_!=incidentIter.nid1_) || (adjkeyIdx_!=incidentIter.adjkeyIdx_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Pointer to graph
    graph_type* graphPtr_;
    // Index of the current node
    size_type nid1_;
    // Iterator of the adjacent node
    size_type adjkeyIdx_;
    /** Private Constructor */
    IncidentIterator(const graph_type* graph, const size_type nid1, size_type adjkeyIdx): 
      graphPtr_(const_cast<graph_type*>(graph)), nid1_(nid1), adjkeyIdx_(adjkeyIdx) {}

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

    /**
     * @brief Dereference operator of edge iterator
     *
     * @return Edge object that the current iterator points to
     *
     * Complexity: O(1)
     */
    Edge operator*() const {
      return Edge(graphPtr_, graphPtr_->i2e_[edgeIterIdx_], 
        graphPtr_->i2u_[graphPtr_->edges_[graphPtr_->i2e_[edgeIterIdx_]].nid1], 
        graphPtr_->i2u_[graphPtr_->edges_[graphPtr_->i2e_[edgeIterIdx_]].nid2]);
    }

    /**
     * @brief Increment operator of incident iterator
     *
     * @param[out] edgeIterIdx_ Index of edge the iterator currently points to
     * @return Reference to the current edge iterator, with pointer pointing to the next 
     *
     * Complexity: O(1)
     */
    EdgeIterator& operator++() {
      ++edgeIterIdx_;
      return *this;
    }

    /**
     * @brief Equality operator of edge iterator

     * @param[in] edgeIter Another edge iterator
     * @return True if _this_ and _edgeIter_ are the same
     *         (belonging to the same graph and pointing to the same edge)
     *
     * Complexity: O(1)
     */
    bool operator==(const EdgeIterator& edgeIter) const {
      return ((graphPtr_==edgeIter.graphPtr_) && (edgeIterIdx_==edgeIter.edgeIterIdx_));
    }

    /**
     * @brief Inequality operator of edge iterator

     * @param[in] edgeIter Another edge iterator
     * @return True if _this_ and _edgeIter_ are not the same
     *         (not belonging to the same graph or not pointing to the same edge)
     *
     * Complexity: O(1)
     */
    bool operator!=(const EdgeIterator& edgeIter) const {
      return ((graphPtr_!=edgeIter.graphPtr_) || (edgeIterIdx_!=edgeIter.edgeIterIdx_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Pointer to graph
    graph_type* graphPtr_;
    // Edge index
    size_type edgeIterIdx_;
    /** Private Constructor */
    EdgeIterator(const graph_type* graph, size_type eid) : 
    graphPtr_(const_cast<graph_type*>(graph)), edgeIterIdx_(eid) {}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /**
   * @brief Get the start of the iterator over all edges
   
   * @return Start of the edge iterator
   *
   * Complexity: O(1).
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this,0);
  }

  /**
   * @brief Get the end of the iterator over all edges
   
   * @return End of the edge iterator
   *
   * Complexity: O(1).
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this,num_edges());
  }



 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    Point pos; // position
    size_type idx_; // user-seen node id
    node_value_type val; // value to be specified
    std::unordered_map<size_type, size_type> adj_;  //map of <node id, corresponding edge id> for each neighboring node
    std::vector<size_type> adjkey_;  //list of node id for each neighboring node, used as indexing into the map above
    /** Constructor */
    internal_node(const Point &position, size_type node_index, node_value_type value, std::unordered_map<size_type, size_type> adj, std::vector<size_type> adjkey):
      pos(position), idx_(node_index), val(value), adj_(adj), adjkey_(adjkey) {}
  };

  struct internal_edge {
    size_type nid1;  // first node's id
    size_type nid2;  // second node's id
    size_type idx_;  // user-seen edge id
    edge_value_type val; // value to be specified
    /** Constructor */
    internal_edge(size_type node_index_1, size_type node_index_2, size_type edge_index, edge_value_type value):
      nid1(node_index_1), nid2(node_index_2), idx_(edge_index), val(value) {}
  };

};

#endif // CME212_GRAPH_HPP
