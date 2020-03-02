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
template <typename V = int, typename E = int>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node;
  struct internal_edge;
  struct internal_adjac;
  
  // nodes and edges in the graph
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

  // unique ids for nodes and edges
  std::vector<unsigned int> n_id2u;
  std::vector<unsigned int> e_id2u;

  // size of the graph in terms of nodes and edge numbers
  unsigned size_= 0;
  unsigned num_edges_= 0;

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

  // Type of nodes
  typedef V node_value_type;
  typedef E edge_value_type;
  
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
      // HW0: YOUR CODE HERE
      
      return g->nodes_[idx_].position;
      
    }
    /** Return this node's position. */
    Point& position() {
      return g->nodes_[idx_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return g->nodes_[idx_].idx;
    }
     
    /** Return this node's value. */
     node_value_type& value() {
       return g->nodes_[idx_].n_val;
     }

     /** Return this node's value. */
     const node_value_type& value() const {
       return g->nodes_[idx_].n_val;
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
      if((g ==n.g) && (idx_== n.idx_))
           return true;
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
       
      if(idx_< n.idx_)
          return true;
      return false;
    }
 
    // return the number of incident edges.
    size_type degree() const {
       return g->nodes_[idx_].adjac.size();
    }

    // start of the incident iterator.
    incident_iterator edge_begin() const {
       return IncidentIterator(g, this, 0);
    }

    // End of the incident iterator.
    incident_iterator edge_end() const {
       return IncidentIterator(g, this, this->degree());
    }

    private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph*  g;
    size_type idx_;
    
    // constructor
    Node(const graph_type* graph, size_type idx) 
        : g(const_cast<graph_type*>(graph)), idx_(idx){
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
  Node add_node(const Point& position, const node_value_type& n_val = node_value_type()) {
    // HW0: YOUR CODE HERE
    // creating vector for incident edges
    std::vector<internal_adjac> adjac{};

    size_type idx = num_nodes();
    nodes_.push_back(internal_node(position, idx, n_val, adjac));
    n_id2u.push_back(nodes_.size()-1);
    ++size_;
    
    return Node(this,num_nodes()-1);
    
   }
  
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return node(n.index())==n;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < size_) {
        Node nd =  Node(this, n_id2u[i]);
        assert(nd.index() == i);
        return nd;
    }
    return Node();
  }

   /** Remove node _n_ from the graph
    * @para n
    * @ pre 0 <= n.index() <= num_nodes(0
    * @post new num_nodes < old num_nodes()
    *
    *Complexity: O(num_nodes()).
   */
  size_type remove_node (const Node& n) {
     if(has_node(n)) {
         while (n.degree() > 0) {
             remove_edge((*n.edge_begin()));
         }
     n_id2u.erase(n_id2u.begin() + n.index());

     for(size_type i = n.index(); i < n_id2u.size(); ++i){
         nodes_[n_id2u[i]].idx = i;
     }
     --size_;
     return 1;
     }
  return 0;
  }
  
   /** Remove node pointed to by iterator _ n_it_
    * @ para n_it 
    * @ pre 0 <= *it.index() <= num_nodes()
    * @post new num_nodes < old num_nodes()
    *
    *Complexity: O(1).
   */
  node_iterator remove_node (node_iterator n_it) {
      remove_node((*n_it));
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
      // HW0: YOUR CODE HERE
      return Node(g, node1_idx); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g, node2_idx);
    }
    
    /** Return the length of this edge */ 
    double length() const{
        return norm(node1().position() - node2().position());
    }    
      
    /** Return this edge's value */ 
    edge_value_type& value() {
        return g->edges_[edge_idx].e_val;
    }

    /** Return this edge's value */
    const edge_value_type& value() const {
        return g->edges_[edge_idx].e_val;
    }
   
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      
      //HW0: YOUR CODE HERE
      if((node1_idx==e.node1_idx && node2_idx==e.node2_idx) || (node1_idx == 
          e.node2_idx && node2_idx == e.node1_idx))
         return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      
      //HW0: YOUR CODE HERE
     
      return (edge_idx < e.edge_idx);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    
    graph_type* g;
    size_type node1_idx;
    size_type node2_idx;
    size_type edge_idx;
       
    Edge(const graph_type* graph, size_type id_1, size_type id_2, size_type e_id)
       :g(const_cast<graph_type*>(graph)), node1_idx(id_1), node2_idx(id_2),edge_idx(e_id){};    

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_edges_){
        size_type id = e_id2u[i];
        return Edge(this, edges_[id].node1_idx, edges_[id].node2_idx, id);
    }
    return Edge();   
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(has_node(a) && has_node(b));
    
    for(unsigned i = 0; i < nodes_[a.idx_].adjac.size(); ++i) {
        if(nodes_[a.idx_].adjac[i].id2==b.idx_) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& e_val = edge_value_type()) {
    // HW0: YOUR CODE HERE
    assert(has_node(a) && has_node(b));
    e_id2u.push_back(num_edges_-1);
    if(!has_edge(a, b)){
        edges_.push_back(internal_edge(a.idx_, b.idx_, num_edges_, e_val));
        // adding ajacent edges
        nodes_[a.idx_].adjac.push_back(internal_adjac(num_edges_-1, b.idx_));
        nodes_[b.idx_].adjac.push_back(internal_adjac(num_edges_-1,a.idx_));
        ++num_edges_;
        return Edge(this, a.idx_, b.idx_, num_edges_-1) ;
    }else{
        return Edge(this, a.idx_, b.idx_, num_edges_-1);
    }
  }

   /** Remove edge represented by two nodes from the graph
    * @para _n1_, a distinct valid node of the graph
    * @para _n2_, a distinct valid node of the graph
    * @ pre has_edge(_n1_,_ n2_)==True
    * @post new num_edges < old num_edges()
    *
    *Complexity: O(num_nodes() + num_edges()).
    * can invalidate nodes connected to the removed nodes of the graph
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
      if(has_edge(n1, n2)) {
          for(size_type i = 0; i < n1.degree(); ++i) {
             internal_adjac e = nodes_[n1.idx_].adjac[i];
             if(e.id2 == n2.idx_) {
                size_type e_id = edges_[e.edge_idx].edge_idx;

                nodes_[n1.idx_].adjac.erase(nodes_[n1.idx_].adjac.begin()+i);
                for(size_type j =0; j < n2.degree(); ++j) {
                   if(nodes_[n2.idx_].adjac[j].id2 == n1.idx_) {
                      nodes_[n2.idx_].adjac.erase(nodes_[n2.idx_].adjac.begin()+j);
                   }
                }
                e_id2u.erase(e_id2u.begin() + e_id);
                for(size_type k = e_id; k < e_id2u.size(); ++k) {
                  edges_[e_id2u[k]].edge_idx = k;
                }
                -- num_edges_;
                return 1;
             }
          }
      }
      return 0;
  }
    
  /** Remove edge from the graph
    * _e_ a distict valid edge object
    * @para _e_ 
    * 
    * @ pre has_edge(n1, n2)
    * @post new num_edges < old num_edges()
    *
    *Complexity: O(num_edges()).
    * can invalidate all edges with index above removed e.idx_
   */
  size_type remove_edge (const Edge& e) {
      return remove_edge(e.node1(), e.node2());
  }

  /** Remove edge pointed by edge iterator _e_it_ from the graph
    * @pre  _e_it_ is a distict valid edge iterator
    * @para _e_it_ 
    * 
    * @post new num_edges < old num_edges()
    *
    *Complexity: O(1)
    *Can invalidate all outstanding edges with idx above removed edge_idx.
   */
  edge_iterator remove_edge (edge_iterator e_it ) {
      remove_edge((*e_it));
      return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    size_= 0;
    num_edges_ = 0;
    nodes_.clear();
    edges_.clear();
    n_id2u.clear();
    e_id2u.clear();
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

    /** Points to the next node*/
    NodeIterator& operator++() {
           ++cur_nod;
           return *this;
    }

   /** A comparison operator to test if two iterators are equal*/
    bool operator == (const NodeIterator& node_iter) const {
        return (node_iter.g == g && node_iter.cur_nod == cur_nod); 
    }

   /**Return the node pointed to by the iterator*/
    Node operator*() const {
        return g->node(cur_nod);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* g;
    size_type cur_nod;
    
    /** private constructor accessible by the Node class*/
    NodeIterator(const graph_type* graph, unsigned idx_): g(const_cast<graph_type*>(graph)), cur_nod(idx_) {
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  //Start of the node iterator pointing to the start of the nodes_ vector

  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }

  // End of a node iterator pointing to the end of the nodes_ vector
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
     
    /** increment operator for incident edges*/
    IncidentIterator& operator++(){
        ++cur_nod;
        return *this;
    }

   /** Checks whether this iterator and inc_iter are equala nd returnxss a bool*/

    bool operator == (const IncidentIterator&  inc_iter) const {
        return((inc_iter.g == g) && (inc_iter.cur_nod == cur_nod) && (inc_iter.node_==node_));
    }

   // returns the incident edge pointed to*/
    Edge operator*() const {
        size_type n_idu = node_->idx_;
        return Edge(g, n_idu, g->nodes_[n_idu].adjac[cur_nod].id2, g->nodes_[n_idu].adjac[cur_nod].edge_idx);
    }
    
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* g;
    node_type* node_;
    size_type cur_nod;

    /** private constructor forthe incident_iterator*/
    IncidentIterator(const graph_type* graph, const node_type* n, size_type it):
       g(const_cast<graph_type*>(graph)), node_(const_cast<node_type*>(n)), cur_nod(it) {} 
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
 
    /** points to the next edge*/
    EdgeIterator& operator++() {
        ++cur_edge;
        return *this;
    }
    
    /** test if this edge iterator is equal to edge_iter and returns a bool*/
    bool operator == (const EdgeIterator& edge_iter) const {
        return ((edge_iter.g ==g) && (edge_iter.cur_edge == cur_edge));
    }

    /** returns an edge pointed to by the iterator*/
    Edge operator*() const {
        return g->edge(cur_edge);
    }

      
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    graph_type* g;
    size_type cur_edge;

    /** private constructor accessible by the edge class*/
    EdgeIterator(const graph_type* graph, size_type e_it):
       g(const_cast<graph_type*>(graph)), cur_edge(e_it) {}
   
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Start of Edge iterator pointing at the start of the edges_ vector*/
  EdgeIterator edge_begin() {
       return EdgeIterator(this, 0);
  }

  /** End of Edge iterator pointing to the end of the edges_vector*/
  EdgeIterator edge_end() {
       return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
   
  /** internal storage for incident edges*/ 
  struct internal_adjac {
       size_type edge_idx;
       size_type id2;
       internal_adjac(size_type edge_idx_, size_type id2_) : edge_idx(edge_idx_), id2(id2_) {};
  };
    
  /**internal storage for a node*/
  struct internal_node {
        Point position;
        size_type idx;
        node_value_type n_val;
        std::vector<internal_adjac> adjac;
        internal_node(Point position_, size_type idx_, node_value_type val, 
          std::vector<internal_adjac> incid) : position(position_), idx(idx_), n_val(val), adjac(incid) {};
  };

  /** internal storage for edge*/
  struct internal_edge {
        size_type node1_idx;
        size_type node2_idx;
        size_type edge_idx;
        edge_value_type e_val;
        internal_edge(size_type n1, size_type n2, size_type e1, edge_value_type val)
          :node1_idx(n1), node2_idx(n2), edge_idx(e1), e_val(val) {};
  };

};

#endif // CME212_GRAPH_HPP
