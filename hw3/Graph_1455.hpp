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


unsigned int compteur()
{
  static unsigned int count = 0;
  return count++;
}

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


  typedef V node_value_type;
  typedef E edge_value_type;
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


//  Graph(size_type uid)
  Graph()
    : size_(0) {}

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
      // HW0: YOUR CODE HERE
    }

    /**
    *@return The value of the node
    *As it returns a reference, this function allows the user to access the
    *value and modify it.
    */
    node_value_type& value() {
      return graph_->nodes_[index_].value_ ;
    }

    /**
    *@return The value of the node
    *This function return a constant reference, thus only allows to access the
    *value but not to modify it.
    */
    const node_value_type& value() const {
      return graph_->nodes_[index_].value_ ;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return get_position();
    }

    Point& position() {
      // HW0: YOUR CODE HERE
      return get_position();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_ ;
      //return size_type(-1);
    }

    /** Return a pointer to this node's graph*/
    graph_type* graph() const {
     return(this->graph_);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /**
    *@brief Compute the number of incident edges to the node.
    *@return The number of incident edges.
    */
    size_type degree() const{
     return ((graph_->adj_map_).at(index_)).size();
    }

    /**
    *@return An iterator pointing at the first incident edge of the node
    */
    // Start of the incident iterator.
    incident_iterator edge_begin() const {
     return IncidentIterator(this->graph_, index_, 0);
    }

    /**
    *@return An iterator pointing to the end of the list of incident edges of
    *the node
    */
    // End of incident iterator.
    incident_iterator edge_end() const {
      return IncidentIterator(this->graph_, index_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
     // return the number of incident edges.
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((n.graph() == this->graph())&&(n.index()==this->index())) {
        return true;
      }
      (void) n;          // Quiet compiler warning
      return false;
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
     if (n.index()>this->index()) {
       return true;
     }
     (void) n;           // Quiet compiler warning
     return false;
    }

    bool change_index(size_type new_index) {

      if (this->index() != new_index ) {
        assert(!graph_->has_node(graph_->node(new_index))) ;

        internal_node in = graph_->nodes_.at(this->index()) ;
        graph_->nodes_[new_index] = in ;
        std::vector<Edge> inc_edges ;
        graph_->adj_map_.insert(std::make_pair(new_index, inc_edges));

        for (auto it = this->edge_begin(); it != this->edge_end();) {

          Edge e = (*it) ;
          Node node2 = e.node2() ;
          edge_value_type val = e.value() ;
          graph_->remove_edge(*this, node2) ;

          graph_->add_edge(graph_->node(new_index), node2, val) ;
          it = this->edge_begin();
        }

        graph_->nodes_.erase(this->index()) ;
        graph_->adj_map_.erase(this->index()) ;
        index_ = new_index ;
        return 1;
      }
      return 0 ;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    size_type index_ ;
    graph_type* graph_;

    /*Private constructor*/
    Node(const graph_type* graph, size_type index) : index_(index), graph_(const_cast<graph_type*>(graph)){
    }

    Point& get_position() const {
      return graph_->nodes_[index_].position_;
    }

    friend class Graph;
    friend class NodeIterator;
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
   Node add_node(const Point& position, const node_value_type& node_val = node_value_type()) {
     // HW0: YOUR CODE HERE
     size_type ind_node = num_nodes();

     internal_node new_node ;
     new_node.position_ = position ;
     new_node.value_ = node_val ;

     nodes_[ind_node] = new_node;
     ++size_;
     std::vector<Edge> inc_edges ;
     adj_map_.insert(std::make_pair(ind_node, inc_edges));
     return Node(this,size_-1);
   }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
   bool has_node(const Node& n) const {
     // HW0: YOUR CODE HERE
     if (nodes_.count(n.index())>0) {
       return true;
     }
     return false;
   }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
   Node node(size_type i) const {
     // HW0: YOUR CODE HERE
     return Node(this, i);
   }
///////////////////////////////////////////////////////////////////

  /**
  @brief Remove a node from the graph, and all its adjacent edges.
  @param n Node to be removed.
  @return 1 if the node have been removed, ie if it was on the graph, 0 otherwise.

  @post 0 <= i < num_nodes() for each index i of the graph nodes.
  @post new num_nodes() = old num_nodes - 1 if @return = 1;
        new num_nodes() = old num_nodes if @return = 0;

  The node n is removed and the node with the highest index in the graph takes the index of n instead, if n wasn't the node with the highest index.

  If this method returns 1, the following objects are invalidated:
  - Former node with index old num_nodes() - 1
  - Edges incident to former node with index old num_nodes() - 1
  - NodeIterator iterating on this graph nodes.
  - IncidentIterator iterating on old n.
  - IncidentIterator iterating on old node with the highest index.
  - EdgeIterator iterating on this graph edges.

  The complexity of this function is, considering the maximum degree of a node as O(1) :
  O(log(num_node()) + log(num_edges()))
  */

  size_type remove_node(const Node& n){
    if (has_node(n)) {
      size_type ind = n.index() ;
  //    IncidentIterator i_i = n.edge_begin();
      auto it_end = adj_map_[ind].end() ;

      for (auto it = adj_map_[ind].begin(); it != adj_map_[ind].end();) {
        remove_edge((*it)) ;
        it = adj_map_[ind].begin();
        it_end = adj_map_[ind].end() ;
      }
      nodes_.erase(ind) ;
      adj_map_.erase(ind) ;
      --size_ ;
      if (ind != num_nodes()) {
        Node last_node = Node(this, num_nodes());
        last_node.change_index(ind) ;
      }
      return 1;
    }
    return 0 ;
  }

  /**
  @brief Remove a node from the graph, and all its adjacent edges.
  @param n_it NodeIterator pointing at the node to be removed.
  @return NodeIterator pointing at the element following the one that has been removed.

  @pre n_it != this->end_nodes().
  @post 0 <= i < num_nodes() for each index i of the graph nodes.
  @post new num_nodes() = old num_nodes - 1 if @return = 1;

  Complexity : O(log(num_node()) + log(num_edges()))
  See specifications of "size_type remove_node(const Node& n)" for more precisions.
  */
  node_iterator remove_node(node_iterator n_it){
    size_t ind = (*n_it).index() ;
    remove_node(*n_it) ;
    if (ind == this->num_nodes()-1) {
      return this->node_end() ;
    }
    return NodeIterator(this, ind);
  }
///////////////////////////////////////////////////////////////////
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    edge_value_type& value(){
      return graph_->edges_[index()].value_ ;
    }

    const edge_value_type& value() const{
      return graph_->edges_[index()].value_ ;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node2_);      // Invalid Node
    }

    size_type index() const {
      size_type ind_n1 = std::min(node1_, node2_);
      size_type ind_n2 = std::max(node1_, node2_);
      auto r = graph_->edges_lookup_.at(std::make_pair(ind_n1, ind_n2));
      return  r;
    }
    graph_type* graph() const {
     return(this->graph_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
     bool operator==(const Edge& e) const {
       if ((e.node1() == graph_->node(node1_))&&(e.node2()==graph_->node(node2_))){
         return true;
       }
       if ((e.node1() == graph_->node(node2_))&&(e.node2()==graph_->node(node1_))){
         return true;
       }
       (void) e;           // Quiet compiler warning
       return false;
     }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
     bool operator<(const Edge& e) const {
       if (this->graph() != e.graph()) {
         return true ;
       }
       if (this->index() < e.index()) {
         return true;
       }
       (void) e;           // Quiet compiler warning
       return false;
     }

     bool change_index(size_type new_index) {
       if (this->index() != new_index ) {

         internal_edge ie = graph_->edges_.at(this->index()) ;
         size_type ind_n1 = std::min(node1_, node2_);
         size_type ind_n2 = std::max(node1_, node2_);
         graph_->edges_.erase(this->index()) ;

         graph_->edges_lookup_.at(std::make_pair(ind_n1, ind_n2)) = new_index;

         graph_->edges_[new_index] = ie ;

         index_ = new_index ;
         return 1 ;
       }
     return 0 ;
     }

   private:

    Edge(const graph_type* graph, Node node1, Node node2, size_type index)
      : graph_(const_cast<graph_type*>(graph)), node1_(node1.index()), node2_(node2.index()), index_(index){
        index_ = this->index() ;
    }

    Edge(const graph_type* graph, Node node1, Node node2)
      : graph_(const_cast<graph_type*>(graph)), node1_(node1.index()), node2_(node2.index()), index_(graph->num_edges()){

    }

    graph_type* graph_;
    size_type node1_ ;
    size_type node2_ ;
    size_type index_;

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
   size_type num_edges() const {
     // HW0: YOUR CODE HERE
     return edges_.size();
   }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
   Edge edge(size_type i) const {
     // HW0: YOUR CODE HERE
     assert((0 <= i)&&(this->num_edges()>i));
     return Edge(this, (edges_.at(i).n1_), (edges_.at(i).n2_), i);
   }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
   bool has_edge(const Node& a, const Node& b) const {
     // HW0: YOUR CODE HERE

      assert(a.graph() == this);
      assert(b.graph() == this);

      size_type ind_a_ = a.index();
      size_type ind_b_ = b.index();

      size_type ind_a = std::min(ind_a_, ind_b_);
      size_type ind_b = std::max(ind_a_, ind_b_);

      std::pair<size_type, size_type> ind_edge = std::make_pair(ind_a, ind_b);

      return (edges_lookup_.count(ind_edge) > 0);
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
    Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_val = edge_value_type()) {
     // HW0: YOUR CODE HERE
      assert(has_node(a)&&has_node(b));

      size_type ind_a_ = a.index();
      size_type ind_b_ = b.index();
      size_type ind_a = std::min(ind_a_, ind_b_) ;
      size_type ind_b = std::max(ind_a_, ind_b_) ;
      Node a_new = node(ind_a);
      Node b_new = node(ind_b);
      Node* a_new_p;
      Node* b_new_p;

      if (ind_a == ind_a_) {
        a_new_p = const_cast<Node*>(&a);
        b_new_p = const_cast<Node*>(&b);
      }
      else {
        a_new_p = const_cast<Node*>(&b);
        b_new_p = const_cast<Node*>(&a);
      }

      size_type ind_edge = num_edges() ;
      bool added;

      std::pair<size_type, size_type> ind_nodes = std::make_pair(ind_a, ind_b);
      added = edges_lookup_.emplace(ind_nodes, ind_edge).second ;

      if (added) {
        internal_edge new_edge ;
//        new_edge.n1_ = *a_new_p;
        new_edge.n1_ = node(ind_a);
//        new_edge.n2_ = *b_new_p;
        new_edge.n2_ = node(ind_b);
        new_edge.value_ = edge_val;
        edges_[ind_edge] = new_edge;
//        adj_map_[ind_a].push_back(Edge(this, a_new_p, b_new_p, ind_edge));

        adj_map_[ind_a].push_back(Edge(this, a_new, b_new, ind_edge));
//        adj_map_[ind_b].push_back(Edge(this, b_new_p, a_new_p, ind_edge));
        adj_map_[ind_b].push_back(Edge(this, b_new, a_new, ind_edge));
      }
      else {
        ind_edge = edges_lookup_[ind_nodes] ;
      }

      return Edge(this, a_new, b_new, ind_edge) ;
   }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    edges_lookup_.clear();
    edges_.clear();
    nodes_.clear();
    size_ = 0;
  }

/////////////////////////////////////////////////////////////
  /**
  @brief Remove and edge from the graph.
  @param n1 Endpoint node 1 of the edge.
  @param n2 Endpoint node 2 of the edge.
  @return 1 if the edge (n1, n2) have been removed, ie if it was on the graph, 0 otherwise.

  @pre n1 and n2 are in the same graph
  @post 0 <= i < num_edge() for each index i of the graph edges.
  @post new num_edges() = old num_edges - 1 if @return = 1;

  If this method returns 1, the following objects are invalidated:
  - Former edge with index old num_edge() - 1
  - EdgeIterator iterating on this graph edges.
  - IncidentIterator iterating on n1 and n2.
  - IncidentIterator iterating on endpoints on the former edge with the highest index.

  The complexity of this function is, considering the maximum degree of a node as O(1) :
  O(log(num_edges()))
  */

  size_type remove_edge(const Node& n1, const Node& n2) {
    if (has_edge(n1, n2)) {
      size_type ind_n1_ = n1.index() ;
      size_type ind_n2_ = n2.index() ;
      size_type ind_n1 = std::min(ind_n1_, ind_n2_);
      size_type ind_n2 = std::max(ind_n1_, ind_n2_);
      std::pair<size_type, size_type> ind_nodes = std::make_pair(ind_n1, ind_n2);

      size_type ind_e = edges_lookup_[ind_nodes] ;

      bool not_found = true ;
      auto it = n1.edge_begin() ;

      while ((not_found)&&(it != n1.edge_end())) {

        if (*it == edge(ind_e) ) {

          adj_map_[ind_n1_].erase(adj_map_[ind_n1_].begin()+it.curr_edge_);
          not_found = false ;
        }
        ++it ;
      }
      not_found = true  ;

      it = n2.edge_begin() ;
      while ((not_found)&&(it != n2.edge_end())) {

        if (*it == edge(ind_e) ) {
          adj_map_[ind_n2_].erase(adj_map_[ind_n2_].begin()+it.curr_edge_);
          not_found = false ;
        }
        ++it ;
      }
      edges_lookup_.erase(ind_nodes) ;
      edges_.erase(ind_e) ;

      if (ind_e != num_edges()) {

        internal_edge last_iedge = edges_[num_edges()];
        Edge last_edge = Edge(this, (last_iedge.n1_), (last_iedge.n2_));
        last_edge.change_index(ind_e) ;
      }
      return 1 ;
    }

    return 0 ;
  }
  /**
  @brief Remove and edge from the graph.
  @param e Edge to be removed.
  @return 1 if the edge e have been removed, ie if it was on the graph, 0 otherwise.

  @post 0 <= i < num_edge() for each index i of the graph edges.
  @post new num_edges() = old num_edges - 1 if @return = 1;

  Complexity :   O(log(num_edges()))
  See specifications of "size_type remove_edge(const Node& n1, const Node& n2)" for more precisions.
  */
  size_type remove_edge(const Edge& e) {
    Node n1 = e.node1() ;
    Node n2 = e.node2() ;
    return remove_edge(n1, n2);
  }

  /**
  @brief Remove and edge from the graph.
  @param e_it  EdgeIterator pointing at the edge to be removed.
  @return 1 if the edge e have been removed, ie if it was on the graph, 0 otherwise.

  @pre e_it != this->end_edges()
  @post 0 <= i < num_edge() for each index i of the graph edges.
  @post new num_edges() = old num_edges - 1 if @return = 1;

  Complexity : O(log(num_edges()))
  See specifications of "size_type remove_edge(const Node& n1, const Node& n2)" for more precisions.
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *(e_it) ;

    ++e_it ;

    if (e_it == this->edge_end()){
      remove_edge(e) ;
      return this->edge_end() ;
    }

    NodeIterator n_i = e_it.n_i_ ;
    Node n = *n_i ;
    IncidentIterator i_i = e_it.i_i_;
    size_type curr_edge = i_i.curr_edge_;

    remove_edge(e) ;

    if (curr_edge != 0) {
        --curr_edge ;
    }

    NodeIterator new_n_i = NodeIterator(this, n.index()) ;
    IncidentIterator new_i_i = IncidentIterator(this, n.index(), curr_edge) ;

    return EdgeIterator(new_n_i, new_i_i, this) ;
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
    NodeIterator() {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
    *Increments to the next node of the graph.
    *return This iterator
    *@pre the iterator is not already at the end
    */
    NodeIterator& operator++() {
      index_node_++ ;
      return *this;
    }

    /**
    *Dereferences operator
    *return The node to which the iterator points.
    */
    Node operator*() const {
      return Node(graph_, index_node_);
    }

    /**Defines equality between two node iterators*/
    bool operator==(const NodeIterator& ni) const {
      return ((ni.graph_==this->graph_)&&(ni.index_node_ == this->index_node_));
    }

    /**Defines inequality between two node iterators*/
    bool operator!=(const NodeIterator& ni) const {
      return ((ni.graph_!=this->graph_)||(ni.index_node_ != this->index_node_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* graph_ ;
    size_type index_node_;

    NodeIterator(const Graph* graph_, size_type index_node_): graph_{graph_}, index_node_{index_node_} {}

    };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
  *@return An iterator pointing at the first node of the graph
  */
   node_iterator node_begin() const {
     return NodeIterator(this, 0) ;
   }

   /**
   *@return An iterator pointing at the end of the nodes of the graph.
   */
   node_iterator node_end() const {
     return NodeIterator(this, this->size_) ;
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

    /**Dereferences operator*/
    Edge operator*() const {
      Edge edge = (graph_->adj_map_.at(ind_node_)).at(curr_edge_);
      return(edge);
    }

    /**
    *Increments to the next edge of the node this iterator iterates on.
    *return This iterator
    *@pre the iterator is not already at the end
    */
    IncidentIterator& operator++() {
      curr_edge_++;
      return *this ;
    }

    /**Defines equality between two incident iterators*/
    bool operator==(const IncidentIterator& ii) const {
      if ((this->graph_ == ii.graph_)&&(this->ind_node_ == ii.ind_node_)&&(this->curr_edge_==ii.curr_edge_)) {
        return true;
      }
    return false;
    }

    /**Defines inequality between two incident iterators*/
    bool operator!=(const IncidentIterator& ii) const {
      if ((this->graph_ != ii.graph_)||(this->ind_node_ != ii.ind_node_)||(this->curr_edge_!=ii.curr_edge_)) {
        return true;
      }
    return false;
    }

   private:
    friend class Graph;
    friend class Node;
    friend class EdgeIterator;
    // HW1 #3: YOUR CODE HERE

    IncidentIterator(const Graph* graph_, size_type ind_node_, size_type curr_edge_) : graph_{graph_}, ind_node_{ind_node_}, curr_edge_{curr_edge_} {}

    const Graph* graph_;
    size_type ind_node_ ;
    size_type curr_edge_ ; //index of the current edge in the adjacency list;
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

    /**Dereferences operator*/
    Edge operator*() const {
      return *i_i_ ;
    }

    /**Defines equality between two incident iterators*/
    bool operator==(const EdgeIterator& ei) const {
      if ((this->n_i_ == ei.n_i_)&&(this->i_i_ == ei.i_i_)) {
        return true;
      }
      return false;
    }

    bool operator!=(const EdgeIterator& ei) const {
      if ((this->n_i_ != ei.n_i_)||(this->i_i_ != ei.i_i_)) {
        return true;
      }
      return false;
    }

    /**
    *Increments to the next edge of the graph.
    *It iterates over all the incident edges of a node before starting with
    *another node.
    *It iterates over each node only once.
    *return This iterator
    *@pre the iterator is not already at the end
    */
    EdgeIterator& operator++() {
      ++i_i_ ;

      if (n_i_ == graph_->node_end()) {
        return *this;
      }

      while (i_i_ == (*n_i_).edge_end()) {
        ++n_i_;

        if (n_i_ == graph_->node_end()) {
          return *this;
        }
        else {
          i_i_ = (*n_i_).edge_begin();
        }
      }

      size_type edge_ind = (*i_i_).index();

      if (check_map_[edge_ind]) {
        this->operator++() ;
        return *this ;
      }

      check_map_[edge_ind] = true ;
      return *this;
    }

   private:
   // HW1 #5: YOUR CODE HERE

    friend class Graph ;

    NodeIterator n_i_ ;
    IncidentIterator i_i_ ;
    std::map<size_type, bool> check_map_;
    const Graph* graph_ ;

    EdgeIterator(NodeIterator n_i_, IncidentIterator i_i_, const Graph* graph_): n_i_{n_i_}, i_i_{i_i_}, graph_{graph_} {
      for (size_type i = 0; i < graph_->num_edges(); i++) {
        check_map_[i] = false;
      }
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
  *@return An iterator pointing at the first edge of the graph
  */
  // Start of the edge iterator.
  edge_iterator edge_begin() const {
    Node node = Node(this, 0) ;
    return EdgeIterator(this->node_begin(), node.edge_begin(), this);
  }

  /**
  *@return An iterator pointing at the end of the edges of the graph.
  */
  // End of the edge iterator.
  edge_iterator edge_end() const {
    Node node = Node(this, (size_-1)) ;
    return EdgeIterator(this->node_end(), node.edge_end(), this);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node{
    Point position_;
    node_value_type value_;
  };

  struct internal_edge{
    Node n1_ ;
    Node n2_ ;
    edge_value_type value_;
  };

  std::map<size_type, internal_node> nodes_ ;
  std::map<size_type, std::vector<Edge>> adj_map_;
  size_type size_ ;
  std::map<size_type, internal_edge> edges_ ;
  std::map<std::pair<size_type, size_type>, size_type> edges_lookup_;

};

#endif // CME212_GRAPH_HPP
