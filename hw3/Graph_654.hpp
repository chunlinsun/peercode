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
  /** Synonym for Node value. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Synonym for Edge value. */
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
    Node(): gr_(nullptr), idx_(size_type(-1)) {} // empty graph and invalid index

    /** Return this node's position. */
    const Point& position() const {
      assert(gr_ != nullptr and gr_->has_node(*this)); // check if valid graph and node
      assert(idx_ >= 0); // check if valid index
      return gr_->node_lst[idx_].pos_; // return the position
    }

    /** Return this node's position (modifiable). */
    Point& position() {
      assert(gr_ != nullptr and gr_->has_node(*this)); // check if valid graph and node
      assert(idx_ >= 0); // check if valid index
      return gr_->node_lst[idx_].pos_; // return the position
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(gr_ != nullptr); // check if valid graph
      // note that idx_ (index directly associated with the node) is different from index_ (index from the graph)
      return gr_->node_lst[idx_].index_; // return the index from the graph (so keeping up to date)
    }

    /** Return the value of the node (by reference)
    * @return value of node */
    node_value_type& value() {
      assert(gr_ != nullptr and gr_->has_node(*this)); // check if valid graph and node
      assert(idx_ >= 0); // check if valid index
      return gr_->node_lst[idx_].val_;
    }

    /** Return the value of the node (by const reference)
    * @return value of node */
    const node_value_type& value() const {
      assert(gr_ != nullptr and gr_->has_node(*this)); // check if valid graph and node
      assert(idx_ >= 0); // check if valid index
      return gr_->node_lst[idx_].val_;
    }

    /** Return the degree of the node  
    * @return degree of node */
    size_type degree() const {
      return (gr_->adj_lsts[this->index()]).size();
    }

    /** Return the begin iterator of the incident iterator of the node
    * @return the begin iterator of the incident iterator of node */
    incident_iterator edge_begin() const {
      return IncidentIterator(gr_, this->index(), 0);
    }

    /** Return the end iterator of the incident iterator of the node
    * @return the end iterator of the incident iterator of node */
    incident_iterator edge_end() const {
      return IncidentIterator(gr_, this->index(), degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return gr_ == n.gr_ and idx_ == n.idx_;
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
      // use std::less to compare pointers
      if (std::less<graph_type*>{}(gr_, n.gr_)) return true;
      if (std::less<graph_type*>{}(n.gr_, gr_)) return false;
      // if gr_ == n.gr_, then compare the index
      return idx_ < n.idx_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* gr_; // pointer to the graph Gr
    size_type idx_; // node index
    // Valid node constructor (only be constructed by a graph)
    Node(const graph_type* Gr, size_type idx)\
    : gr_(const_cast<graph_type*> (Gr)), idx_(idx) {}
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
  Node add_node(const Point& position, const node_value_type & val = node_value_type()) {
    // new index_ (from graph) and idx_ (associated with node)
    size_type new_idx_, new_index_;
    new_idx_ = node_lst.size();
    new_index_ = this->num_nodes();
    // new node
    Node_Detail new_node_detail = Node_Detail(position, new_index_, val);
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
    return n.gr_ == this and n.idx_ != size_type(-1) and n.index() < this->num_nodes();
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    // graph pointer and node index of two end points
    Edge(): gr_(nullptr), n1_idx_(size_type(-1)), n2_idx_(size_type(-1)) {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(gr_, n1_idx_);     
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(gr_, n2_idx_);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return gr_ == e.gr_ and n1_idx_ == e.n1_idx_ and n2_idx_ == e.n2_idx_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // First compare the graph pointer
      // use std::less
      if (std::less<graph_type*>{}(gr_, e.gr_)) return true;
      if (std::less<graph_type*>{}(e.gr_, gr_)) return false;
      // if gr_ == n.Gr, then compare n1_idx_
      if (n1_idx_ < e.n1_idx_) return true;
      if (n1_idx_ > e.n1_idx_) return false;
      // if gr_ == n.Gr, n1_idx_ == n.n1_idx_, then compare n2_index_
      return n2_idx_ < e.n2_idx_;
    }

    double length() const {
      // length of an edge
      return norm_2(node1().position() - node2().position());
    }

    edge_value_type& value() {
      // find the end node with a shorter adjacency list
      bool tmp = true; // indicator of whether node1 or node2 adjacency list is chosen
      Node n1 = node1();
      Node n2 = node2();
      const Node* node_tmp = &n1;
      if ((gr_->adj_lsts[n1.index()]).size() > (gr_->adj_lsts[n2.index()]).size()){
        node_tmp = &n2;
        tmp = false;
      }
      // loop over the (shorter) adjacency list
      std::vector<Edge_Detail>& adj_lst_tmp = gr_->adj_lsts[(*node_tmp).index()];
      for (auto iter=adj_lst_tmp.begin(); iter!=adj_lst_tmp.end(); ++iter) {
        if (tmp) {
          if ((*iter).n2_idx_ == n2_idx_) return (*iter).e_val_;
        }
        else {
          if ((*iter).n2_idx_ == n1_idx_) return (*iter).e_val_;
        }
      }
      // if edge not found in the adjacency list, which should not happen
      return gr_->e_val_invalid;
    }

    const edge_value_type& value() const {
      // find the end node with a shorter adjacency list
      bool tmp = true; // indicator of whether node1 or node2 adjacency list is chosen
      Node n1 = node1();
      Node n2 = node2();
      const Node* node_tmp = &n1;
      if ((gr_->adj_lsts[n1.index()]).size() > (gr_->adj_lsts[n2.index()]).size()){
        node_tmp = &n2;
        tmp = false;
      }
      // loop over the (shorter) adjacency list
      std::vector<Edge_Detail> adj_lst_tmp = gr_->adj_lsts[(*node_tmp).index()];
      for (auto iter=adj_lst_tmp.begin(); iter!=adj_lst_tmp.end(); ++iter) {
        if (tmp) {
          if ((*iter).n2_idx_ == n2_idx_) return (*iter).e_val_;
        }
        else {
          if ((*iter).n2_idx_ == n1_idx_) return (*iter).e_val_;
        }
      }
      // if edge not found in the adjacency list, which should not happen
      return gr_->e_val_invalid; 
    }



   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type*  gr_;
    size_type n1_idx_;
    size_type n2_idx_;

    // Valid edge constructor (only be constructed by a graph)
    Edge(const graph_type* Gr, size_type n1_idx, size_type n2_idx)\
    : gr_(const_cast<graph_type*> (Gr)), n1_idx_(n1_idx), n2_idx_(n2_idx) {}

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
    return Edge(this, edge_lst[i].n1_idx_, edge_lst[i].n2_idx_);       
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

    /** use incident iterators (also for unittest purpose)
    */
    const Node* node_tmp = &a;
    if ((adj_lsts[a.index()]).size() > (adj_lsts[b.index()]).size()){
      node_tmp = &b;
    }
    // remember to check two directions 
    // (can slightly improve by only checking one direction based on 
    //  whether a or b adj_lst is chosen) 
    for (auto iter=(*node_tmp).edge_begin(); iter!=(*node_tmp).edge_end(); ++iter){
      if (((*iter).node1() == a and (*iter).node2() == b) 
        or ((*iter).node1() == b and (*iter).node2() == a)) {
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
    assert(has_node(a) and has_node(b) and (a_index != b_index));
    // update the related vectors
    if (!has_edge(a, b)) {
      Edge_Detail new_edge_detail_ab = Edge_Detail(a_index, b_index, edge_value_type());
      Edge_Detail new_edge_detail_ba = Edge_Detail(b_index, a_index, edge_value_type());
      edge_lst.push_back(new_edge_detail_ab);
      // for adjacency lists, the edge is included ith the current node as the first node
      adj_lsts[a.index()].push_back(new_edge_detail_ab); 
      adj_lsts[b.index()].push_back(new_edge_detail_ba);
    }

    return Edge(this, a_index, b_index); // return the current edge  
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
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator(): gr_(nullptr), pt_(size_type(-1)) {}

    /** Dereferencing the iterator
    * @return the node that is the dereference of the iterator */
    value_type operator*() const {
      assert((gr_ != nullptr) and (pt_<gr_->num_nodes())); // check input validity
      return Node(gr_, gr_->index_lst[pt_]);
    }

    /** Increment the iterator to the next iterator 
    * @return the next iterator */
    NodeIterator& operator++() {
      assert((gr_ != nullptr) and (pt_<gr_->num_nodes())); // check input validity
      ++pt_;
      return *this;
    }

    /** Check if two iterators are equal 
    * @param[in] const reference of an iterator
    * @return true if the two iterators are equal, otherwise return false */
    bool operator==(const NodeIterator& node_iter1) const {
      return ((gr_ == node_iter1.gr_) and (pt_ == node_iter1.pt_));
    }

   private:
    friend class Graph;
    /** private constructor that is only available through Graph
    * @param[in] gr: pointer to the graph 
    * @param[in] pt: index of the node in the iterator */
    NodeIterator(const graph_type* gr, size_type pt): gr_(const_cast<graph_type*> (gr)), pt_(pt) {}

    graph_type* gr_; // pointer to the graph
    size_type pt_; // index of the node (from the graph, not the unique id) in the iterator
  };

  /** Return the begin iterator of the node iterator of the graph
  * @return the begin iterator of the node iterator of the graph */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return the end iterator of the node iterator of the graph
  * @return the end iterator of the node iterator of the graph */
  node_iterator node_end() const {
    return NodeIterator(this, this->num_nodes());
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(): gr_(nullptr), node_index_(size_type(-1)), pt_(size_type(-1)) {}


    /** Dereferencing the iterator
    * @return the edge that is the dereference of the iterator */
    value_type operator*() const {
      assert((gr_ != nullptr) and (node_index_<gr_->num_nodes()));
      assert(pt_ < (gr_->adj_lsts[node_index_]).size());
      value_type edge_tmp = Edge(gr_,
        gr_->adj_lsts[node_index_][pt_].n1_idx_,
        gr_->adj_lsts[node_index_][pt_].n2_idx_);
      return edge_tmp;
    }

    /** Increment the iterator to the next iterator 
    * @return the next iterator */
    IncidentIterator& operator++() {
      assert((gr_ != nullptr) and (node_index_<gr_->num_nodes()));
      assert(pt_ < (gr_->adj_lsts[node_index_]).size());
      ++pt_;
      return *this;
    }

    /** Check if two iterators are equal 
    * @param[in] const reference of an iterator
    * @return true if the two iterators are equal, otherwise return false */
    bool operator==(const IncidentIterator& inc_iter1) const {   
      return ((gr_ == inc_iter1.gr_) and 
        (node_index_ == inc_iter1.node_index_) and 
        (pt_ == inc_iter1.pt_));
    }

   private:
    friend class Graph;
    /** private constructor that is only available through Graph
    * @param[in] gr: pointer to the graph 
    * @param[in] node_index: node index
    * @param[in] pt: index of the incident edge in the adjacency list */
    IncidentIterator(const graph_type* gr, size_type node_index, size_type pt): \
                     gr_(const_cast<graph_type*> (gr)), node_index_(node_index), pt_(pt) {}

    graph_type* gr_; // pointer to the graph
    size_type node_index_; // node index (from the graph, not the unique id)
    size_type pt_; // index of the incident edge in the adjacency list
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(): gr_(nullptr), pt_(size_type(-1)) {}

    /** Dereferencing the iterator
    * @return the edge that is the dereference of the iterator */
    value_type operator*() const {
      assert((gr_ != nullptr) and (pt_<gr_->num_edges())); // check input validity
      value_type edge_tmp = Edge(gr_, 
      gr_->edge_lst[pt_].n1_idx_,
      gr_->edge_lst[pt_].n2_idx_);
      return edge_tmp;
    }

    /** Increment the iterator to the next iterator 
    * @return the next iterator */
    EdgeIterator& operator++() {
      assert((gr_ != nullptr) and (pt_<gr_->num_edges())); // check input validity
      ++pt_;
      return *this;
    }

    /** Check if two iterators are equal 
    * @param[in] const reference of an iterator
    * @return true if the two iterators are equal, otherwise return false */
    bool operator==(const EdgeIterator& edge_iter1) const {
      return ((gr_ == edge_iter1.gr_) and (pt_ == edge_iter1.pt_));
    }

   private:
    friend class Graph;
    /** private constructor that is only available through Graph
    * @param[in] gr: pointer to the graph 
    * @param[in] pt: index of the edge in the iterator */
    EdgeIterator(const graph_type* gr, size_type pt): gr_(const_cast<graph_type*> (gr)), pt_(pt) {}

    graph_type* gr_; // pointer to the graph
    size_type pt_; // index of the edge
  };

  /** Return the begin iterator of the edge iterator of the graph
  * @return the begin iterator of the edge iterator of the graph */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Return the end iterator of the edge iterator of the graph
  * @return the end iterator of the edge iterator of the graph */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->num_edges());
  }

  /** Remove the node n
  * @param[in] n: node to be removed
  * @return number of nodes in the updated graph (after removal)
  *
  * @post new num_nodes() = old num_nodes() - 1 if the node n exists
  * @post new num_nodes() = old num_nodes() if the node does not exist
  * @post invalidate the node n as well as all edges incident to it 
  *
  * Complexity: O(num_nodes())
  */

  size_type remove_node(const Node&n) {
    // check if the node is valid
    if (n.gr_ != this or n.index()>=num_nodes()) return num_nodes();

    while (n.edge_begin() != n.edge_end()) {
      remove_edge(*(n.edge_begin()));
    }

    // remove the empty adjacency list of node n: O(num_nodes()) time at most
    adj_lsts.erase(adj_lsts.begin() + n.index());


    // remove the node n: O(num_nodes()) time at most
    index_lst.erase(index_lst.begin() + n.index()); 
    for (auto iter=index_lst.begin()+n.index(); iter!=index_lst.end(); ++iter) {
      // remove all nodes with index larger than the current node
      node_lst[*iter].index_ -= 1;
    }

    return num_nodes();
  }

  /** Remove the node n with iterator n_it
  * @param[in] n_it: a valid node iterator
  * @return the next iterator
  *
  * @post the same as size_type remove_node(const Node&);
  *
  * Complexity: O(num_nodes())
  */
  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    remove_node(n);
    return n_it;
  }



  /** Remove the edge with end nodes n1 and n2
  * @param[in] n1: node 1 of the edge
  * @param[in] n2: node 2 of the edge
  * @return number of edges in the graph after removing the edge (n1, n2)
  *
  * @post new num_edges() = old num_edges() - 1 if the edge (n1, n2) exists
  * @post new num_edges() = old num_edges() if the edge does not exist
  * @post invalidate the incident iterator from n1 and n2 by removing the edge from the
  *       corresponding adjacency lists in adj_lsts; also invalidate the edge iterator corresponding to the
  *       removed edge by removing it from the edge list edge_lst;
  *
  * Complexity: O(max degree of the graph + num_edges()) 
  */
  size_type remove_edge(const Node&n1, const Node&n2) {
    // check if edge is valid
    if (!has_edge(n1, n2)) return num_edges();

    // remove the edge from edge_lst: O(num_edges()) search, and O(1) removal
    for (auto iter=edge_lst.begin(); iter!=edge_lst.end(); iter++) {
      if (((*iter).n1_idx_==n1.idx_ and (*iter).n2_idx_==n2.idx_) or
          ((*iter).n1_idx_==n2.idx_ and (*iter).n2_idx_==n1.idx_)) {
        *iter = edge_lst.back();
        edge_lst.pop_back();
        break; 
      }
    }

    // remove the edge from n1's adjacency list: O(max degree of the graph)
    for (auto iter=adj_lsts[n1.index()].begin(); iter!=adj_lsts[n1.index()].end(); iter++) {
      if ((*iter).n2_idx_==n2.idx_) {
        adj_lsts[n1.index()].erase(iter);
        break;
      }
    }

    // remove the edge from n2's adjacency list: O(max degree of the graph)
    for (auto iter=adj_lsts[n2.index()].begin(); iter!=adj_lsts[n2.index()].end(); iter++) {
      if ((*iter).n2_idx_==n1.idx_) {
        adj_lsts[n2.index()].erase(iter);
        break;
      }
    }

    return num_edges();
  }

  /** Remove the edge e
  * @param[in] e: a valid edge e 
  * @return number of edges in the updated graph
  * 
  * @post the same as size_type remove_edge(const Node&, const Node&);
  *
  * Complexity: O(max degree of the graph + num_edges()) 
  */
  size_type remove_edge(const Edge& e) {
    Node n1 = e.node1();
    Node n2 = e.node2();
    return remove_edge(n1, n2);
  }

  /** Remove the edge with edge iterator e_it
  * @param[in] e_it: a valid edge iterator
  * @return the next edge iterator
  * 
  * @post the same as size_type remove_edge(const Node&, const Node&);
  *
  * Complexity: O(max degree of the graph + num_edges()) 
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    remove_edge(e);
    return e_it;
  }

 private:
  // Graph class's internals:
  //   helper functions, data members, and so forth.
  /** Struct for node detailed info (position and index)
  */
  struct Node_Detail{
    Point pos_;
    size_type index_;
    node_value_type val_;
    Node_Detail(const Point& pos, size_type index, node_value_type val): \
                pos_(pos), index_(index), val_(val) {}
  };
  /** Struct for edge detailed info (n1_idx_ and n2_idx_)
  */
  struct Edge_Detail{
    size_type n1_idx_;
    size_type n2_idx_;
    edge_value_type e_val_;
    Edge_Detail(size_type n1_idx, size_type n2_idx, edge_value_type e_val): \
                n1_idx_(n1_idx), n2_idx_(n2_idx), e_val_(e_val) {}
  };
  // vectors for the nodes
  std::vector<Node_Detail> node_lst; // list of nodes, each entry storing the pos_ and index_ and val_
  std::vector<size_type> index_lst; // list of index_ (from graph), each entry storing the idx_ (for the node)
  // vectors for the edges
  std::vector<Edge_Detail> edge_lst; // list of edges, each entry storing n1_idx_ and n2_idx_ and e_val_
  std::vector<std::vector<Edge_Detail>> adj_lsts; // list of adjacency lists for each node, each entry is a list of edges connected to the node
  edge_value_type e_val_invalid{}; // invalid edge value type
};

#endif // CME212_GRAPH_HPP
