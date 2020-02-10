#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * @collab Yue Li
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <set>

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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  // HW0: YOUR CODE HERE
  //Private Graph Constructor
    : nodes_(), adj_list_(), edge_tuples_() {}

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
      //Invalid node constructor takes no arguments
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_.at(idx_).position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:


    /** Return this node's value, an argument of the Graph Template
    * This value can be updated
    * @return node value type of node
    */
    node_value_type& value(){
      return graph_->nodes_.at(idx_).node_value_;
    };

    /** Return this node's value as a const, an argument of the Graph Template
    * This value cannot be updated
    * @return node value type of node
    */
    const node_value_type& value() const{
      return graph_->nodes_.at(idx_).node_value_;
    };

    /** Return this node's degree,
    * The degree of a node is defined as the number of edges incident to a node
    * @return size type indicating degree of node
    */
//--functionality_1
//--out of bounds error on degree test, what if there are no edges?
//--START
    size_type degree() const{
      return graph_->adj_list_.at(idx_).size();
    };
//--END
    /** Return an incident iterator object indicating the beginning of an
    * iterator
    * @pre  edge_begin != edge_end
    * @return and incident iterator object
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, idx_,
        graph_->adj_list_.at(idx_).begin(),
        graph_->adj_list_.at(idx_).end());
    };

    /** Return an incident iterator object indicating the end of an
    * iterator
    * @return and incident iterator object
    */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, idx_,
        graph_->adj_list_.at(idx_).end(),
        graph_->adj_list_.at(idx_).end());
    };

    /** Test whether this node and @a n are equal.
     * @param _n_ node to compare this node with
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      return graph_ == n.graph_ && idx_ == n.idx_;
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
      (void) n;           // Quiet compiler warning
      if (graph_ == n.graph_){ //check if graphs are equal
        return idx_ < n.idx_; // return index ordering
      }
      else {
        //compare pointers to maintain trichotomy
        return std::less<Graph*>()(graph_, n.graph_);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_; //pointer to graph object
    size_type idx_; //index of Point object in data structure


    /**Private Constructor */
    Node(const Graph* grph, size_type idx)
      : graph_(const_cast<Graph*>(grph)), idx_(idx){}


  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return nodes_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
     const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE

    nodes_.emplace_back(internal_node(position, node_value));
    (void) position;      // Quiet compiler warning
    return Node(this, nodes_.size() - 1 );        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    assert(n.idx_ < nodes_.size());
    (void) n;            // Quiet compiler warning
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < nodes_.size());
    (void) i;             // Quiet compiler warning
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      //Empty Constructor
    }

    /** Return a node of this Edge */
    Node node1() const {
      //Get first node from edge object
      return Node(graph_, node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //Get second node from edge object
      return Node(graph_, node2_idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return ((graph_ == e.graph_) &&
              ((node1_idx_ == e.node1_idx_
              && node2_idx_ == e.node2_idx_) ||
              (node1_idx_ == e.node2_idx_ &&
              node2_idx_ == e.node1_idx_)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if (graph_ == e.graph_ ){ //Check if graphs are equal
        return edge_idx_ < e.edge_idx_; //compare indices
      }
      else {
        return std::less<Graph*>()(graph_, e.graph_); //compare pointers
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type edge_idx_;
    size_type node1_idx_;
    size_type node2_idx_;


    /**Private Constructor */
    Edge(const Graph* grph, size_type idx,
      size_type node1_idx, size_type node2_idx)
      : graph_(const_cast<Graph*>(grph)),
        edge_idx_(idx),
        node1_idx_(node1_idx),
        node2_idx_(node2_idx){};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_tuples_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    assert(i < edge_tuples_.size());
    return Edge(this, i, edge_tuples_.at(i).at(0), edge_tuples_.at(i).at(1));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning

    if ( adj_list_.find(a.index()) == adj_list_.end() ) {
      //Did not find node A, return false
      return false;
    } else {
      //Found node A, return bool if finds node B
      return adj_list_.at(a.index()).count(b.index());
    }

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

    if (has_edge(a, b)){
      //If edge exists, call constructor with edge index
      return Edge(this, adj_list_.at(a.index()).find(b.index())->second,
                  a.index(), b.index());
    }

    else {
    //Add Pairing to vector of tuples
    std::vector<size_type> index_tuples = {a.index(), b.index()};
    edge_tuples_.emplace_back(index_tuples);

    //Add to adjacency list
    adj_list_[a.index()][b.index()] = edge_tuples_.size();
    adj_list_[b.index()][a.index()] = edge_tuples_.size();

    (void) a, (void) b;   // Quiet compiler warning
    return Edge(this, edge_tuples_.size() - 1, a.index(), b.index());
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edge_tuples_.clear();
    adj_list_.clear();

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

    /** Dereference Operator
    * @pre iterator != end
    * @return a node object
    */
    Node operator*() const{
      return Node(graph_,curr_idx_);
    }

    /** Increment Operator
    * @pre iterator not at end
    * @return a node iterator
    */
    NodeIterator& operator++(){
        ++curr_idx_;
        return (*this);
    }

    /** Equality Operator
    * @param _node_iter_ other node iterator to compare this with
    * @return boolean indicating equality
    */
    bool operator==(const NodeIterator& node_iter) const{
      return (curr_idx_ == node_iter.curr_idx_
      && graph_ == node_iter.graph_);
    }

    /** Inequality Iterator
    * @param _node_iter_ other node iterator
    * @return boolean indicating inequality
    */
    bool operator!=(const NodeIterator& node_iter) const{
      return !(*this == node_iter);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    Graph* graph_;
    size_type curr_idx_;

    //Node Iterator Constructor
    NodeIterator(const Graph* grph,
      size_type idx)
    				: graph_(const_cast<Graph*>(grph)), curr_idx_(idx){};
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return the beginning of a node iterator
  * @pre nodes_.size() >= 0
  * @return the beginning of a node iterator
  */
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }

  /**Return the end of a node iterator
  * @pre nodes.size() >= 0
  * @return the end of a node iterator
  */
  node_iterator node_end() const{
    return NodeIterator(this, nodes_.size());
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

    /** Dereference Operator for Incident Iterators
    * @pre iterator != end
    * @return an edge object
    */
    Edge operator*() const{
      return Edge(graph_, incident_iter_begin_->second,
                  node1_idx_, incident_iter_begin_->first);
    };

    /** Increment Operator for Incident Iterators
    * @pre iterator != end
    * @return an incident iterator object
    */
    IncidentIterator& operator++(){
      ++incident_iter_begin_;
      return(*this);
    };

    /** Equality Operator for Incident Iterators
    * @param _inc_iter_ other incident iterator
    * @return boolean indicating equality
    */
    bool operator==(const IncidentIterator& inc_iter) const{
      return (graph_ == inc_iter.graph_ &&
              node1_idx_ == inc_iter.node1_idx_ &&
              incident_iter_begin_ == inc_iter.incident_iter_begin_ &&
              incident_iter_end_ == inc_iter.incident_iter_end_);
    };

    /** Inequality Operator for Incident Iterators
    * @param _inc_iter_ other incident iterator
    * @return boolean indicating inequality
    */
    bool operator!=(const IncidentIterator& inc_iter) const{
      return !(*this == inc_iter);
    };

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node1_idx_;
    std::map<size_type, size_type>::iterator incident_iter_begin_;
    std::map<size_type, size_type>::iterator incident_iter_end_;

    //Incident Iterator Constructor
    IncidentIterator( const Graph* grph,
      size_type node_idx,
      std::map<size_type, size_type>::iterator incident_iter_begin,
      std::map<size_type, size_type>::iterator incident_iter_end)
						: graph_(const_cast<Graph*>(grph)),
              node1_idx_(node_idx),
              incident_iter_begin_(incident_iter_begin),
						  incident_iter_end_(incident_iter_end){};
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

    /** Dereference Operator for Edge Iterators
    * @pre iterator not at end
    * @return an edge object
    */
    Edge operator*() const{
      return Edge(graph_, curr_edge_idx_,
                  graph_->edge_tuples_[curr_edge_idx_][0],
                  graph_->edge_tuples_[curr_edge_idx_][1]);
    };


    /** Increment Operator for Edge Iterators
    * @pre iterator != end
    * @return a reference to an edge iterator object
    */
    EdgeIterator& operator++(){
      ++curr_edge_idx_;
      return (*this);
    };


    /** Equality Operator for Edge Iterators
    * @pre edge_iter valid iterator
    * @param _edge_iter_ other edge iterator
    * @return boolean indicating equality between iterators
    */
    bool operator==(const EdgeIterator& edge_iter) const{
      return (graph_ == edge_iter.graph_ &&
              curr_edge_idx_ == edge_iter.curr_edge_idx_);
    };


    /** Inequality Operator for Edge Iterators
    * @param _edge_iter_ other edge iterator
    * @pre _edge_iter_ valid iterator
    * @return boolean indicating inquality between iterators
    */
    bool operator!=(const EdgeIterator& edge_iter) const{
      return !(*this == edge_iter);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type curr_edge_idx_;


    //Edge iterator constructor
    EdgeIterator(const Graph* grph, size_type curr_edge_idx)
						: graph_(const_cast<Graph*>(grph)),
            curr_edge_idx_(curr_edge_idx){};
  };


  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Returns beginning of an edge iterator
  * @pre num _edge_tuples_.size() >= 0 
  * @return edge iterator at beginning of iterator
  */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** Returns beginning of an edge iterator
  * @pre _edge_tuples_.size() >= 0
  * @return edge iterator at end of iterator
  */
  edge_iterator edge_end() const{
    return EdgeIterator(this, edge_tuples_.size());
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.

  //Internal type for node elements
  struct internal_node{
    const Point position_;
    node_value_type node_value_;

  //Internal Node Struct representation
  internal_node(const Point position,
     const node_value_type& node_value = node_value_type())
      : position_(position), node_value_(node_value){};

  };

  std::vector<internal_node> nodes_;

  //map of maps representing node relations in an adjacency list with
  //an additional map to edge index
  std::map<size_type, std::map<size_type, size_type>> adj_list_;

  //vector of node pairs in an edge
  std::vector<std::vector<size_type>> edge_tuples_;

};

#endif // CME212_GRAPH_HPP
