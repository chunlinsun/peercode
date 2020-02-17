#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <functional>
#include <memory>

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

  /** Type of data carried in node. */
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
  Graph() {
    // HW0: YOUR CODE HERE
    //values default initialized below
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

  /** Remove the node referenced by the iterator
      @param[in] _n_  A valid node object

      @return The number of nodes deleted. This will be 1 if the node is valid and 0 otherwise.

      @post Node object with be removed from internal container. Nodes stored before the node removed
            will keep the same ordering. Nodes stored after the node removed may be reindexed and moved.

            If the node referenced is valid, new num_nodes() == old num_nodes() - 1. new num_edges() = old num_edges() - node_degree()
            Else,                            new num_nodes() == old num_nodes().     new num_edges() = old num_edges()

            Invalidates the node being deleted.
            Invalidates Node iterators pointing to the last element or beyond of the edges container
            Invalidates Incident iterators referencing map objects containing the deleted edges (so 
            invident iterators spawned this node object and node objects that are adjacent to this node). 
            Invalidates the deleted Edge objects. Possibly reindexes Edge objects.
            Invalidates Edge iterators pointing to edges pointing to any of the last node_degree() edges in the graph or beyond.

      Complexity: O(node_degree() * O(remove_edges())) = O(node_degree()). 
   */
  size_type remove_node(const Node& n){
    size_type idx = n.index();
    if (idx < num_nodes()){
      node_iterator new_it = NodeIterator(this, idx);
      remove_node(new_it);
    }
    else{
      std::cout << "\n\n TRYING TO REMOVE AN INVALID NODE \n\n";
      return 0;
    }
    return 1;
  };

  /** Remove the node referenced by the iterator
      @param[in, out] _n_it_  As input, this is an node iterator referencing a valid Node object. 
                              As output, this is an node iterator referecing a valid Node object or
                              the end interator for edges.

      @post Node object with be removed from internal container. Nodes stored before the node removed
            will keep the same ordering. Nodes stored after the node removed may be reindexed and moved.

            If the node referenced is valid, new num_nodes() == old num_nodes() - 1. new num_edges() = old num_edges() - node_degree()
            Else,                            new num_nodes() == old num_nodes().     new num_edges() = old num_edges()

            Invalidates the node being deleted.
            Invalidates Node iterators pointing to the last element or beyond of the edges container
            Invalidates Incident iterators referencing map objects containing the deleted edges (so 
            invident iterators spawned this node object and node objects that are adjacent to this node). 
            Invalidates the deleted Edge objects. Possibly reindexes Edge objects.
            Invalidates Edge iterators pointing to edges pointing to any of the last node_degree() edges in the graph or beyond.

      Complexity: O(node_degree() * O(remove_edges())) = O(node_degree()). 
   */
  node_iterator remove_node(node_iterator n_it){
    size_type idx = (*n_it).index();
    Node n = *n_it;
    while(n.edge_begin() != n.edge_end()){
      Edge e = Edge(this, idx, nodes_[idx]->adj_.begin()->first, nodes_[idx]->adj_.begin()->second);
      remove_edge(e);
    }



    auto last = nodes_.rbegin();
    if (nodes_[idx] != *last){
      size_type temp = nodes_[idx]->uuid_;
      nodes_[idx]->adj_ = (*last)->adj_;
      nodes_[idx]->value = (*last)->value;
      nodes_[idx]->point = (*last)->point;
      nodes_[idx]->uuid_ = (*last)->uuid_;

      (*last)->uuid_ = temp;
      for (auto it = (nodes_[idx]->adj_).begin(); it != (nodes_[idx]->adj_).end(); ++it){
        edges_[it->second]->node1 = idx;
        edges_[it->second]->node2 = it->first;  //edit indices that have lastidx as a node
        nodes_[it->first]->adj_.erase(num_nodes() - 1);
        nodes_[it->first]->adj_[idx] = it->second;
      }
      n_uuidToidx[nodes_[idx]->uuid_] = idx;
      n_uuidToidx[temp] = num_nodes() - 1;
    }

    nodes_.pop_back();

    return n_it;
  };

  /** Remove the edge connecting the nodes passed as parameters
      @param[in] _a_  A valid Node object.
      @param[in] _b_  A valid Node object.

      @return The number of edges deleted. This will be 1 if there is an edge connecting _a_ and _b_ and 0 otherwise.

      @post Edge object with be removed from internal container. Edges stored before the edge removed
            will keep the same ordering. Edges stored after the edge removed may be reindexed and moved.

            If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
            Else,                       new num_edges() == old num_edges().

            Invalidates the edge being deleted.
            Invalidates Edge iterators pointing to the last element or beyond of the edges container
            Invalidates Incident iterators referencing map objects containing the deleted edge (so 
            invident iterators spawned from a and b. 
            Possibly reindexes Edge objects coming after the deleted edge in the container ordering.

      Complexity: O(1)
   */
  size_type remove_edge(const Node& a, const Node& b){
    if (!has_edge(a, b)) {
      std::cout << "TRYING TO REMOVE EDGE INVALID EDGE\n";
      return 0;
    }

    size_type edge_idx = nodes_[a.index()]->adj_.at(b.index()); //get edge index
    if (edge_idx >= num_edges()) {
      std::cout << "INDEX IS MORE THAN NUM OF EDGES (INVALID)" << edge_idx << "\n";
      return 0;
    }
    edge_iterator e_it = EdgeIterator(this, edge_idx);
    remove_edge(e_it);
    return 1;
  };

  /** Remove the edge object passed as parameter
      @param[in] _e_  A valid Edge object.

      @return The number of edges deleted. This will be 1 if there is an edge e in the graph and 0 otherwise.

      @post Edge object with be removed from internal container. Edges stored before the edge removed
            will keep the same ordering. Edges stored after the edge removed may be reindexed and moved.

            If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
            Else,                       new num_edges() == old num_edges().

            Invalidates the edge being deleted.
            Invalidates Edge iterators pointing to the last element or beyond of the edges container
            Invalidates Incident iterators referencing map objects containing the deleted edge (so 
            invident iterators spawned from a and b. 
            Possibly reindexes Edge objects coming after the deleted edge in the container ordering.

      Complexity: O(1)
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  };

  /** Remove the edge referenced by the iterator
      @param[in, out] _e_it_  As input, this is an edge iterator referencing a valid Edge object. 
                              As output, this is an edge iterator referecing a valid Edge object or
                              the end interator for edges.

      @post Edge object with be removed from internal container. Edges stored before the edge removed
            will keep the same ordering. Edges stored after the edge removed may be reindexed and moved.

            If old has_edge(@a (*e_it).node1(), @a (*e_it).node1()), new num_edges() == old num_edges() - 1.
            Else,                                                     new num_edges() == old num_edges().

            Invalidates the edge being deleted.
            Invalidates Edge iterators pointing to the last element or beyond of the edges container
            Invalidates Incident iterators referencing map objects containing the deleted edge (so 
            invident iterators spawned from (*e_it).node1() and (*e_it).node2()). 
            Possibly reindexes Edge objects coming after the deleted edge in the container ordering.

      Complexity: O(1)
   */
  edge_iterator remove_edge(edge_iterator e_it){
    Edge e = *e_it;
    size_type edge_idx = nodes_[e.node1().index()]->adj_.at(e.node2().index()); //get edge index

    nodes_[e.node1().index()]->adj_.erase(e.node2().index());// = 0; //remove adjacencies
    nodes_[e.node2().index()]->adj_.erase(e.node1().index());


    auto last = edges_.rbegin();
    if (edges_[edge_idx] != *last){
      edges_[edge_idx]->node1 = (*last)->node1;
      edges_[edge_idx]->node2 = (*last)->node2;
      edges_[edge_idx]->value = (*last)->value;

      nodes_[(*last)->node1]->adj_[(*last)->node2] = edge_idx;  //relabel swapped index adjacencies
      nodes_[(*last)->node2]->adj_[(*last)->node1] = edge_idx;
    }

    edges_.pop_back();//delete last element
    return e_it;
  };


  class Node : private totally_ordered<Node>{
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
      //params are default initialized to something below
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return *(this->graph_->nodes_[index()]->point);
    }

    /** Return this node's position. */
    Point& position() {
      return *(this->graph_->nodes_[index()]->point);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->n_uuidToidx.at(index_);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    /** Return this node's value */
    node_value_type& value(){
      return this->graph_->nodes_[index()]->value;
    }

    /** Return this node's value */
    const node_value_type& value() const{
      return this->graph_->nodes_[index()]->value;
    }

    /** Return this node's degree */
    size_type degree() const{
      return graph_->nodes_[index()]->adj_.size();
    }

    /** @return the begin iterator for edges incident with this node */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, graph_->nodes_[index()]->adj_.begin(), index());
    }

    /** @return the end iterator for edges incident with this node */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, graph_->nodes_[index()]->adj_.end(), index());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
	  return (n.index() == this->index()) && (n.graph_ == this->graph_);
      //(void) n;          // Quiet compiler warning
      //return false;
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
	  if (n.graph_ != this->graph_) return compare(this->graph_, n.graph_);
	  return  this->index() < n.index();
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects 
    //Point * point = nullptr;
    size_type index_ = 0;
    Graph* graph_ = nullptr;
    Node(Graph * graph){
      graph_ = graph;
      index_ = graph_->size();
    }

    Node(Graph * graph, size_type idx){
      graph_ = graph;
      index_ = idx;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();
    //return 0;
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
  Node add_node(const Point& position,  const node_value_type& v = node_value_type()) {
    // HW0: YOUR CODE HERE
    Node n(this);
     std::shared_ptr<Point> pos = std::make_shared<Point>(position.elem[0], position.elem[1], position.elem[2]);
     std::shared_ptr<internal_node> new_internal_node = std::make_shared<internal_node>(this, pos, size(), v);
  	nodes_.push_back(new_internal_node);
    return n;
    //(void) position;      // Quiet compiler warning
    //return Node();        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
	  return n.graph_ == this;
    //(void) n;            // Quiet compiler warning
    //return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i > size()){
      std::cout << "REQUESTING AND INVALID NODE\n";
    }
    return Node(const_cast<Graph*>(this), nodes_[i]->uuid_);
    //(void) i;             // Quiet compiler warning
    //return Node();        // Invalid node
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
  class Edge :  private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      //values default initialized below
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      size_type n1 = graph_->edges_[index_]->node1;
      return graph_->node(n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type n2 = graph_->edges_[index_]->node2;
      return graph_->node(n2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
		// return (e.index_ == this->index_) && (e.graph_ == this->graph_);
      bool matching_nodes = ((e.node1() == node1()) && (e.node2() == e.node2())) || ((e.node1() == node2()) && (e.node2() == e.node1()));
      return (e.graph_ == this->graph_) && matching_nodes;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
		if (e.graph_ != this->graph_) return compare(this->graph_, e.graph_);
		  return this->index_ < e.index_;
    
    }

    edge_value_type& value(){return graph_->edges_[index_]->value;};
    const edge_value_type& value() const {return graph_->edges_[index_]->value;};

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
	Graph * graph_ = nullptr;
	size_type index_ = 0;

  Edge(Graph * graph, size_type a, size_type b, size_type index) {
    this->graph_ = graph;
    this->index_ = index;
    graph_->edges_[index_]->node1 = a;
    graph_->edges_[index_]->node2 = b;
  }
  
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
	  return edges_.size();
    //return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(const_cast<Graph*>(this), edges_[i]->node1, edges_[i]->node2, i);
	  // return *(edges_[i]->edge);
    //(void) i;             // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    std::unordered_map<size_type, size_type> map = nodes_[a.index()]->adj_;
    std::unordered_map<size_type, size_type>::const_iterator got = map.find(b.index());

    if (got == map.end()) {
      return false;
    }
    return true;
	  
   // //(void) a; (void) b;   // Quiet compiler warning
   // return false;
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

    if (has_edge(a,b)){
      size_type idx = nodes_[a.index()]->adj_.at(b.index());
      return Edge(this, a.index(), b.index(), idx);
    }
    std::shared_ptr<internal_edge> new_internal_edge = std::make_shared<internal_edge>(a.index(), b.index());
    edges_.push_back(new_internal_edge);

    nodes_[a.index()]->adj_[b.index()] = num_edges()- 1;
    nodes_[b.index()]->adj_[a.index()] = num_edges()- 1;

    return Edge(this, a.index(), b.index(), num_edges() - 1);
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Dereference iterator to get node. */
    Node operator*() const{
      return graph_->node(graph_->nodes_[idx_]->index_);
    }

    /** Increment iterator to get next. */
    NodeIterator& operator++(){
      idx_++;
      return *this;
    }

    /** Compare equality of iterators (this and _ni_).
    
        Iterators are equal if they belong to the same graph and point to the same node.
    */
    bool operator==(const NodeIterator& ni) const{
      return ((ni.idx_ == idx_) && (ni.graph_ == graph_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* graph_;
    size_type idx_;
    NodeIterator(const Graph* graph, int idx){
      graph_ = graph;
      idx_ = idx;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return the begin iterator for nodes in the graph. */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** Return the end iterator for nodes in the graph. */
  node_iterator node_end() const{
    return NodeIterator(this, nodes_.size());
  }

  class IncidentIterator : private totally_ordered<IncidentIterator> {
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

    /** Dereference iterator to get edge. 
  
        @return an Edge object e with e.node1() == @a a and e.node2() == @a b

    */
    Edge operator*() const{
      size_type edge_idx = it_->second;
      size_type n1 = graph_->edge(edge_idx).node1().index();
      size_type n2 = graph_->edge(edge_idx).node2().index();

      if (n1 == node_idx_) return Edge(const_cast<Graph*>(graph_), n1, n2, edge_idx);
      return Edge(const_cast<Graph*>(graph_), n2, n1, edge_idx);
    }

    /** Increment iterator to get next. */
    IncidentIterator& operator++(){
      ++it_;
      return *this;
    }

    /** Compare equality of iterators (this and _ii_). 

        Iterators are equal if they have the same map iterator, spawn from the same nodes, and belong to the same graph.

    */
    bool operator==(const IncidentIterator& ii) const {
      return (ii.it_ == it_) && (ii.node_idx_ == node_idx_) && (ii.graph_ == graph_);
    }

   private:
    friend class Graph;
    const Graph * graph_;
    size_type node_idx_;
    std::unordered_map<size_type, size_type>::iterator it_;

    IncidentIterator(const Graph* graph, std::unordered_map<size_type, size_type>::iterator it, size_type node_idx){
      graph_ = graph;
      node_idx_ = node_idx;
      it_=it;
    }
    // HW1 #3: YOUR CODE HERE
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Dereference iterator to get edge. */
    Edge operator*() const{
      return Edge(const_cast<Graph*>(graph_), graph_->edges_[idx_]->node1, graph_->edges_[idx_]->node2, idx_);
    }

    /** Increment iterator to get edge. */
    EdgeIterator& operator++(){
      idx_++;
      return *this;
    }

    /** Compare equality of iterators (this and _ei_). 

        Iterators are equal if they point to the same edge and share the same graph.

    */
    bool operator==(const EdgeIterator& ei) const{
      return (graph_ == ei.graph_) && (idx_ == ei.idx_);
    }

   private:
    friend class Graph;
    const Graph * graph_;
    size_type idx_;
    EdgeIterator(const Graph * graph, size_type idx){
      graph_ = graph;
      idx_ = idx;
    }
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** @return the begin iterator for edges in the graph. */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** @return the end iterator for edges in the graph. */
  edge_iterator edge_end() const{
    return EdgeIterator(this, edges_.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
	 struct internal_node {
    std::shared_ptr<Point> point;
     size_type index_;
     size_type uuid_;
     V value {};
     std::unordered_map<size_type, size_type> adj_ {};
     internal_node(Graph* g,std::shared_ptr<Point> p, size_type idx, V v){
      point = p;
      index_ = idx;
      value = v;
      uuid_ = g->n_uuid;
      g->n_uuid = g->n_uuid + 1;
      g->n_uuidToidx[uuid_] = index_;
     }
  };
	 std::vector<std::shared_ptr<internal_node>> nodes_ {};

	 struct internal_edge {
    size_type node1;
    size_type node2;
    E value {};
     internal_edge(size_type n1, size_type n2){
      node1 = n1;
      node2 = n2;
     }
	 };
	 std::vector<std::shared_ptr<internal_edge>> edges_ {};

   size_type n_uuid = 0;

   std::unordered_map<size_type, size_type> n_uuidToidx {};

   static std::less<Graph*> compare; // stable ordering
};

#endif // CME212_GRAPH_HPP
