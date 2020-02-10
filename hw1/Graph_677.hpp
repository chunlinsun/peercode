#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

const int PRINT_DEBUG = 0;
#define DPRINT(msg)                               \
  do {                                                  \
    if (PRINT_DEBUG) { \
      std::cerr << __FILE__ << ":" << __LINE__ << ": ";   \
      std::cerr << msg;\
      std::cerr << std::endl;\
    }\
  } while (0)


#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <stdio.h>
#include <utility>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V>
class Graph {
 private:


  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct NodeData;

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

  /** Type of node value */
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    numEdges = 0;
  }

  /** New Destructor **/
  ~Graph() {
    clear();
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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return data_->position; 
    }

    node_value_type& value() { return data_->value; };
    const node_value_type& value() const { return data_->value; };

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      int id = data_->index; 
      return size_type(id);
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
      // HW0: YOUR CODE HERE (DONE)
      bool isequal = n.graph_ == this->graph_ && 
                     n.index() == this->index() &&
                     n.position() == this->position();
      return isequal;
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
      // HW0: YOUR CODE HERE (DONE)
      if (this->graph_ == n.graph_) {
        return this->index() < n.index(); 
      }
      else { // Can't compare nodes of different graphs (return error?)
        std::cerr << "WARNING: trying to compare nodes of different graphs\n";
        return 0;
      }
    }

    /** @brief Get the number of edges to/from the node
     *  @return number of edges to/from the node
     * */
    size_type degree() const { return data_->edges.size(); };

    /** @brief Get the starting iterator for all edges connected to the node
     * */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, data_);
    }

    /** @brief Get the ending iterator for all edges connected to the node
     * */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, data_, data_->edges.end());
    }

    Point::value_type distance(const Node& node) {
      return norm(position() - node.position());
    }

    Point::value_type distance(const Point& p) {
      return norm(position() - p);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    friend class Edge;

    const Graph* graph_; // pointer to parent graph
    NodeData* data_;     // point to data in parent graph

    /** Private Constructor **/
    Node(const Graph* graph, size_type index) : graph_(graph) {
      assert(index < graph->num_nodes());
      assert(graph->nodes_[index]->index == (int) index);
      data_ = graph->nodes_[index];
    } 

    Node(const Graph* graph, NodeData* data) : graph_(graph), data_(data) {}

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
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    
    // Create NodeData
    int index = this->size();
    std::set<NodeData*> edges;
    NodeData* ndata = new NodeData{position, index, value, edges}; 
    nodes_.push_back(ndata);

    // Add entry to edges_ vector
    std::vector<int> v;

    // Create the node
    return Node(this, index);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE (DONE)
    int index = n.index();
    if (node(index) == n)
      return true;
    else
      return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE (DONE)
    Node n = Node(this, i);
    assert(n.index() == i);
    return n;
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE (DONE)
      int idx = n1_->index;
      Node n1 = Node(graph_, idx);
      return n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE (DONE)
      int idx = n2_->index;
      Node n2 = Node(graph_, idx);
      return n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      int e11 = this->n1_->index;
      int e12 = this->n2_->index;

      int e21 = e.n1_->index;
      int e22 = e.n2_->index;

      auto e1 = std::minmax(e11,e12);
      auto e2 = std::minmax(e21,e22);
      if (e1.first == e2.first && e1.second == e2.second)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
//--functionality_1
//--it seems like you're using a method edge_index which doesn't exist?
//--START
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (graph_ != e.graph_) {
        std::cerr << "WARNING: Comparing edges from two different graphs" 
                  << std::endl;
        return false;
      }
      int i1 = graph_->edge_index(*this);
      int i2 = graph_->edge_index(e);
      return i1 < i2; 
    }
//--END

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    friend class Node;

    const Graph* graph_;
    const NodeData* n1_;
    const NodeData* n2_;
    
    /* Private Constructor */
    Edge(const Graph* graph, const Node* n1, const Node* n2) 
        : graph_(graph),  n1_(n1->data_), n2_(n2->data_) {
      assert(n1->graph_ == graph);
      assert(n2->graph_ == graph);
    }

    Edge(const Graph* graph, const NodeData* nd1, const NodeData* nd2)
        : graph_(graph), n1_(nd1), n2_(nd2) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE (DONE)
    // Complexity: O(1)
    return edges_.size(); 
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE (DONE)
    // Complexity: O(1)
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE (DONE)
    std::set<NodeData*> edges_a = a.data_->edges;
    std::set<NodeData*> edges_b = b.data_->edges;

    auto a_to_b = edges_a.find(b.data_);
    auto b_to_a = edges_b.find(a.data_);
    if (a_to_b != edges_a.end() && b_to_a != edges_b.end()) {
      assert((*a_to_b)->index == b.data_->index);
      assert((*b_to_a)->index == a.data_->index);
      return true;
    }
    return false;
  }

  /* Print all the index pairs of all the edges in the graph, in order */
  void print_edges() const {
    std::cout << "\nPrinting " << num_edges() << " edges" << std::endl;
    for (size_type i = 0; i < num_edges(); i++) {
      int i1 = edges_[i].node1().index();
      int i2 = edges_[i].node2().index();
      std::cout << "Edge " << i << ": (" << i1 << ", " << i2 << ")" << std::endl;
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
    assert(a.graph_ == this);
    assert(b.graph_ == this);

    // Create edge 
    Edge e = Edge(this, &a, &b);

    if (has_edge(a, b)) {
      // std::cout << "Edge " << a.index() << "-" << b.index() << " already exists\n";
    } else {
      // std::cout << "Adding edge " << a.index() << "-" << b.index() << std::endl;
      
      // Add to list
      edges_.push_back(e);

      // Add edge to node elements
      a.data_->edges.insert(b.data_);
      b.data_->edges.insert(a.data_);
    }

    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    // Clear all edges
    edges_.clear();

    // Clear all nodes
    for (size_type i = 0; i < size(); i++)
      delete nodes_[i];
    nodes_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. 
   * Simply wraps the iterator for the vector of nodes stored in the Graph
   * */
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    
    /** De-reference the node iterator to get a node 
     * */
    Node operator*() const {
      return Node(graph_, *node_itr);
    }

    /** Move the node iterator forward. 
     * */
    NodeIterator& operator++() {
      ++node_itr;
      return *this;
    }

    // Equality comparison
    bool operator==(const NodeIterator& iter) const {
      return (graph_ == iter.graph_) && (node_itr == iter.node_itr);
    }

    bool operator!=(const NodeIterator& iter) const { return !(*this == iter); }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* graph_;  // pointer to the parent graph

    // Iterator for the vector of node data in the parent graph
    typename std::vector<NodeData*>::const_iterator node_itr; 

    /** @brief Private constructor for NodeIterator 
     * @param[in] graph Pointer to the parent graph
     * @param[in] itr Iterator for the graph's vector of node data
     * */
    NodeIterator(const Graph<V>* graph, 
                 typename std::vector<NodeData*>::const_iterator itr)
        : graph_(graph), node_itr(itr) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
  NodeIterator node_begin() const {
    return NodeIterator(this, nodes_.begin());
  }

  NodeIterator node_end() const {
    return NodeIterator(this, nodes_.end());
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
//--documentation_1
//--no documentation for some methods, and not in doxygen style
//--START
    Edge operator*() const {
      return Edge(graph_, node_data_, *itr);
    }

    IncidentIterator& operator++() {
      ++itr;
      return *this;
    }

    bool operator==(const IncidentIterator& i_itr) const {
      return (graph_ == i_itr.graph_) && (itr == i_itr.itr);
    }

    bool operator!=(const IncidentIterator& i_itr) const { 
      return !(*this == i_itr); 
    }
//--END

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph_;    // Pointer to parent graph
    NodeData* node_data_;  // Pointer to node data

    // Iterator over set node (data) connected to the current node
    typename std::set<NodeData*>::const_iterator itr;

    /** Private constructor for Incident Iterator */
    IncidentIterator(const Graph* graph, NodeData* data) 
        : graph_(graph), node_data_(data) {
      // Initialize with iterator at the beginning of the set
      itr = node_data_->edges.begin();
    }

    IncidentIterator(const Graph* graph, NodeData* data, 
                     typename std::set<NodeData*>::const_iterator itr) 
        : graph_(graph), node_data_(data), itr(itr) {} 

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
    
    Edge operator*() const { return *itr; }

    EdgeIterator& operator++() {
      ++itr;
      return *this;
    }

    bool operator==(const EdgeIterator& e_itr) const {
      return (graph_ == e_itr.graph_) && (itr == e_itr.itr);
    }

    bool operator!=(const EdgeIterator& e_iter) const { 
      return !(*this == e_iter); 
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    typename std::vector<Edge>::const_iterator itr;

    /** Private Constructor */
    EdgeIterator(const Graph* graph, 
                 typename std::vector<Edge>::const_iterator itr) 
      : graph_(graph), itr(itr) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** Internal storage representation of a node
   * Contains the following:
   * @param _position_ location of the node in 3D space
   * @param _index_ unique index identifier
   * @param _value_ node value (template type)
   * @param _edges_ std::set of pointers to others nodes to which the node is
   *    connected by an edge
   * */
  struct NodeData {
    const Point position;
    int index;
    node_value_type value;
    std::set<NodeData*> edges;
  };
  
  int numEdges;


  /** Vector of nodes in the graph. 
   * Nodes are allocated on the heap and stored as pointers
   * Invariant: nodes_[i] == nodes_[i]->index
   * */
  std::vector<NodeData*> nodes_;

  /** Vector of edges in the graph
   * Contains all of the edges between nodes of the graph
   * Invariant: if there exists and edge between node i and j, 
   *    THEN the Edge between i,j is in edges_ 
   *    AND j is in nodes_[i]->edges
   *    AND i is in nodes_[j]->edges
   * */
  std::vector<Edge> edges_;

};

#endif // CME212_GRAPH_HPP
