#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>

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

  struct node_element;
  struct edge_element;
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V ;
  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
    node_value_type & value ();
    const node_value_type & value () const ;

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
    : nodes(), edges(), next_node_id(0) {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {
  }
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
    Node()
      : g(nullptr), n(nullptr) {
    }

    /** Return this node's position. */
    const Point& position() const {
      return n->p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return n->id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if(this->g == n.g && this->n->id == n.n->id){
        return true;
      }
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
      if(this->n->id < n.n->id){
        return true;
      }
      return false;
    }

    void update_value(node_value_type v_){
      this->n->v = v_;
      return;
    }

    node_value_type& value(){
      return this->n->v;
    }

    const node_value_type& value() const{
      return this->n->v;
    }

    size_type degree() const {
      return this->n->s.size();
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(this->n,this->n->s.begin(),this->g);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(this->n,this->n->s.end(),this->g);
    }

   private:
    const graph_type* g; // pointer to graph
    node_element* n; // pointer to node element

    Node(const graph_type* g_, node_element* n_)
      : g(g_), n(n_) {
      }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->nodes.size();
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
  Node add_node(const Point& position, const node_value_type & val = node_value_type ()) {
    std::set<node_element*> s;
    node_element* e = new node_element{position,this->next_node_id,s,val}; // allocate on heap

    nodes.push_back(e);
    ++next_node_id;

    // HW0: YOUR CODE HERE
    // (void) position;      // Quiet compiler warning
    return Node(this,e);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n_) const {
    // HW0: YOUR CODE HERE

    // check that id is less than length of graph
    // check that point value is equivalent
    if (n_.n->id < size() && this->nodes[n_.n->id]->p == n_.n->p){
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
    // (void) i;             // Quiet compiler warning
    return Node(this, this->nodes[i]);        // Invalid node
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
    Edge() :
      g(nullptr), e(nullptr) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->g, this->e->n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->g, this->e->n2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //--functionality_0
      //--The other order is not checked.
      //--START
      if(this->g == e.g
        && this->e->n1->id == e.e->n1->id
        && this->e->n2->id == e.e->n2->id) {
        return true;
      }
      return false;
    }
      //--END

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //--functionality_1
      //--No trichotomy upheld.
      //--START
      if(this->g == e.g){
        return true;
      }
      return false;
      //--END
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    const graph_type* g;
    edge_element* e;

    Edge(const graph_type* g_, edge_element* e_)
      : g(g_), e(e_) {
      }

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
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // // HW0: YOUR CODE HERE
    return Edge(this, this->edges[i]);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    typename std::set<node_element*>::iterator it = a.n->s.find(b.n);//.find(b->n);

    if (it != a.n->s.end()) {
      return true;
    }
    else {
      return false;
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

    edge_element* e;

    e = new edge_element{a.n,b.n};

    if (!(has_edge(a,b))){
      edges.push_back(e);

      a.n->s.insert(b.n);
      b.n->s.insert(a.n);
    }

    return Edge(this,e);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  //--functionality_0
  //--Graph not cleared properly.
  //--You should clear the member vectors.
  //--START
  void clear() {
    int n_nodes = size();
    for(int i = 0; i < n_nodes; i++){
      delete nodes[i]; // delete memory allocated on the heap
    }

    int n_edges = num_edges();
    for(int i = 0; i < n_edges; i++){
      delete edges[i]; // delete memory allocated on the heap
    }
  //--END
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
   class NodeIterator : private totally_ordered<NodeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() = default;

    NodeIterator(typename std::vector<node_element*>::const_iterator it_,
      const graph_type* g_)
      : it(it_), g(g_) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return Node(g,(*it));
    }

    NodeIterator& operator++() {
      ++(this->it);
      return *this;
    }

    bool operator==(const NodeIterator& ni) const {
      if (this->it == ni.it) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    typename std::vector<node_element*>::const_iterator it;
    const graph_type* g;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const {
    return NodeIterator(nodes.begin(),this);
  }

  node_iterator node_end() const {
    return NodeIterator(nodes.end(),this);
  }


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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() = default;

    IncidentIterator(node_element* n_,
      typename std::set<node_element*>::const_iterator it_,
      const graph_type* g_)
      : n(n_), it(it_), g(g_) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      edge_element* e = new edge_element{this->n,*(this->it)};
      return Edge(this->g,e);
    }
    IncidentIterator& operator++() {
      ++(this->it);
      return *this;
    }
    bool operator==(const IncidentIterator& ii) const {
      if (this->it == ii.it) {
        return true;
      }
      return false;
    }

   private:

    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    node_element* n;
    typename std::set<node_element*>::iterator it;
    const graph_type* g;
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
    EdgeIterator() = default;

    EdgeIterator(typename std::vector<edge_element*>::const_iterator it_,
      const graph_type* g_)
      : it(it_), g(g_) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return Edge(this->g,*(this->it));
    }
    EdgeIterator& operator++() {
      ++(this->it);
      return *this;
    }
    bool operator==(const EdgeIterator& ii) const {
      if (this->it == ii.it) {
        return true;
      }
      return false;
    }

   private:

    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    typename std::vector<edge_element*>::const_iterator it;
    const graph_type* g;

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this->edges.begin(),this);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this->edges.end(),this);
  }

 private:
   struct node_element {
     Point p;
     size_type id;
     std::set<node_element*> s;
     V v;
   };

   struct edge_element{
     node_element* n1;
     node_element* n2;
   };

   std::vector<node_element*> nodes;
   std::vector<edge_element*> edges;
   size_type next_node_id;
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
