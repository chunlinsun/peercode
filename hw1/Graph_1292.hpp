#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>
#include <set>
#include <iterator>

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

 public:


  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;

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

  using P = Point;
  using pair_type = std::pair<size_type,size_type>;
  using ptr_smaller = std::less<const graph_type*>;

  struct NodeData{
    P position_;
    V value_;
    std::vector<size_type> adj_;
  };
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_->elements_[uid_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /**
     * @brief Allows to access and modify the value of a node
     * @return The value of the node by reference
     * @pre This is called from a valid node
     */
    node_value_type& value(){
      return graph_->elements_[uid_].value_;
    }

    /**
     * @brief Allows to access the value of a node without
     * modifying it
     * @return The value of the node as a constant reference
     * @pre This is called from a valid node
     */
    const node_value_type& value() const{
      return graph_->elements_[uid_].value_;
    }


    /**
     * @brief Getter function for the degree of a node
     * @return The degree of the node
     * @pre This is called from a valid node
     */
    size_type degree() const{
      return graph_->elements_[uid_].adj_.size();
    }

    /**
     * @return The begin iterator that allows iteration over the edges incident
     * to the node
     * @pre This is called from a valid node
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,uid_,graph_->elements_[uid_].adj_.cbegin());
    }

    /**
     * @return The end iterator that allows iteration over the edges incident
     * to the node
     * @pre This is called from a valid node
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_,uid_,graph_->elements_[uid_].adj_.cend());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((this->graph_ == n.graph_) && (this->uid_ == n.uid_));
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
      return (this->uid_< n.uid_) || ( (this->uid_ == n.uid_) && (ptr_smaller{}(this->graph_, n.graph_)));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    graph_type* graph_;
    size_type uid_;

    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return elements_.size();
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
    //--style_0
    //--This can actually be written more concisely as
    //-- elements_.push_back({position,value,{}});
    //--which works because the compiler is able to infer all the types. Of
    //--course it's also fine as you wrote it.
    //--START
    elements_.push_back(NodeData{position,value,std::vector<size_type>()});
    //--END
    return node(elements_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return ( (this==n.graph_) && (n.uid_ < this->elements_.size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);
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
      return graph_->node(n1id_);   
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(n2id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->graph_ == e.graph_) 
              && (std::minmax(this->n1id_,this->n2id_) == std::minmax(e.n1id_,e.n2id_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (std::minmax(this->n1id_,this->n2id_) < std::minmax(e.n1id_,e.n2id_))
              || ( (std::minmax(this->n1id_,this->n2id_) == std::minmax(e.n1id_,e.n2id_)) 
                    && (ptr_smaller{}(this->graph_, e.graph_)));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    const graph_type* graph_;
    size_type n1id_;
    size_type n2id_;

    Edge(const graph_type* graph, size_type n1id, size_type n2id)
      : graph_(graph), n1id_(n1id), n2id_(n2id){
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    pair_type c = *(std::next(edges_.cbegin(),i));
    return Edge(this,c.first,c.second); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    pair_type c = std::minmax(a.index(),b.index());
    return edges_.find(c)!=edges_.end();
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
    auto res = edges_.insert(std::minmax(a.index(),b.index()));
    if (res.second) {
      elements_[a.index()].adj_.push_back(b.index());
      elements_[b.index()].adj_.push_back(a.index());

    }
    return Edge(this,a.index(),b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_.clear();
    elements_.clear();

  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>{
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

    //--documentation_0
    //--Critical precondition for both * and ++ is that this is not equal to the
    //--end iterator.
    //--START
    /**
     * @brief Dereference operator for NodeIterator
     * @pre This is called from a valid NodeIterator
     * @return The node where the iterator is pointing
     */
    value_type operator*() const{
      return graph_->node(index_);
    } 


    /**
     * @brief ++ operator for NodeIterator
     * @pre This is called from a valid NodeIterator
     * @return The updated iterator
     */
    NodeIterator& operator++(){
      index_++;
      return *this;
    }
    //--END

    /**
     * @brief Test whether @ it and this iterator are equal
     * Two NodeIterators are equal if they have the same underlying graph
     * and are at the same index
     */
    bool operator==(const NodeIterator& it) const{
      return (this->graph_== it.graph_) && (this->index_ == it.index_);
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type index_;

    NodeIterator(const graph_type* graph, size_type index) :
      graph_(graph), index_(index){}


  };

  /**
   * @brief Provide begin NodeIterator to iterate over the nodes
   * of this graph
   */
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }

  /**
   * @brief Provide end NodeIterator to iterate over the nodes
   * of this graph
   */
  node_iterator node_end() const{
    return NodeIterator(this,this->num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>{
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

    /**
     * @brief Dereference operator for IncidentIterator
     * @pre This is called from a valid IncidentIterator
     * @post node1() of the returned Edge is always node
     * over whose incident Edges we are iterating
     * @return The node where the iterator is pointing
     */
    Edge operator*() const{
      return Edge(graph_, n1id_, *it_);
    }

    /**
     * @brief ++ operator for IncidentIterator
     * @pre This is called from a valid IncidentIterator
     * @return The updated iterator
     */
    IncidentIterator& operator++(){
      it_++;
      return *this;
    }

    /**
     * @brief Test whether @ it2 and this iterator are equal
     * Two IncidentIterator are equal if they have the same underlying graph,
     * the same "origin" node and they are pointing to the same edge
     */
    bool operator==(const IncidentIterator& it2) const{
      return (this->graph_== it2.graph_) && (this->it_ == it2.it_)
             && (this->n1id_ == it2.n1id_);
    }


   private:
    friend class Graph;
    const graph_type* graph_;
    size_type n1id_;
    std::vector<size_type>::const_iterator it_;
    
    IncidentIterator(const graph_type* graph, size_type n1id,
                     std::vector<size_type>::const_iterator it) :
                    graph_(graph),n1id_(n1id),it_(it){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
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


    /**
     * @brief Dereference operator for EdgeIterator
     * @pre This is called from a valid EdgeIterator
     * @return The edge where the iterator is pointing
     */
    Edge operator*() const{
      return Edge(graph_, (*it_).first, (*it_).second);
    }

    /**
     * @brief ++ operator for EdgeIterator
     * @pre This is called from a valid EdgeIterator
     * @return The updated iterator
     */
    EdgeIterator& operator++(){
      it_++;
      return *this;
    }

    /**
     * @brief Test whether @ it2 and this iterator are equal
     * Two EdgeIterator are equal if they have the same underlying graph,
     * and they are pointing to the same edge
     */
    bool operator==(const EdgeIterator& it2) const{
      return (this->graph_ == it2.graph_) && (this->it_ == it2.it_);
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    std::set<pair_type>::const_iterator it_;

    EdgeIterator(const graph_type* graph, std::set<pair_type>::const_iterator it):
                 graph_(graph),it_(it){}
  };

  //--style_0
  //--Remove this stuff
  //--START
  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  //--END

  /**
   * @brief Provide begin EdgeIterator to iterate over the edges
   * of this graph
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, edges_.cbegin());
  }

  /**
   * @brief Provide end EdgeIterator to iterate over the edges
   * of this graph
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, edges_.cend());
  }

 private:

  //--documentation_0
  //--Provide some short comments on what the private data members represent.
  //--END
  std::vector<NodeData> elements_;
  //--design_0
  //--Why use a std::set here rather than a std::vector? Using a set means that
  //--retrieving an Edge by index is slow (O(num_edges)). I think you can find
  //--a way to make has_edge and edge retrieval both O(1), though it's not
  //--required.
  //--START
  std::set<pair_type> edges_;
  //--END
};

//--design_0
//--Great work!
//--END

#endif // CME212_GRAPH_HPP
