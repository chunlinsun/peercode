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
#include <string>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

#define DEFAULT_EDGE_LENGTH 0.1 

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V, typename E>
class Graph {
 private:


  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct NodeData;
  struct EdgeData;

  // Exceptions
  struct InvalidEdge { std::string msg; };

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
  typedef V node_value_type;

  /** Type of edge value */
  typedef E edge_value_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    numEdges = 0;
    numNodes = 0;
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
    /** Default constructor
     * Create an invalid node
     * */
    Node() {};

    // 
    // Getters
    //

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return data_->position; 
    }

    /** Return this node's position (modifiable). */
    Point& position() { 
      return data_->position;
    }

    /** Return the node's value */
    node_value_type& value() { return data_->value; };
    const node_value_type& value() const { return data_->value; };

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      int id = data_->index; 
      return size_type(id);
    }

    /** @brief Get the number of edges to/from the node
     *  @return number of edges to/from the node
     * */
    size_type degree() const { return data_->edges.size(); };


    //
    // Comparators
    //

    /** Test whether this node and @a n are equal.
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
      } else { 
        long i1 = long(graph_) + long(this->index());
        long i2 = long(n.graph_) + long(n.index());
        return i1 < i2;
      }
    }

    //
    // Iterators
    //

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

    //
    // Extra Methods
    //

    /** Calculate the distance between the node and another node */
    Point::value_type distance(const Node& node) const {
      return norm(position() - node.position());
    }

    /** Calculate the distance between the node and a point*/
    Point::value_type distance(const Point& p) const {
      return norm(position() - p);
    }

    /** Determine if a node is connected to the current node */
    bool is_neighbor(const Node& n) const {
      return graph_->is_neighbor(data_, n.data_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class Edge;

    const Graph* graph_; // pointer to parent graph
    NodeData* data_;     // pointer to data in parent graph

    /** Private Constructor **/
    Node(const Graph* graph, NodeData* data) : graph_(graph), data_(data) {}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return numNodes;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] _position_ The new node's position
   * @param[in] (optional) _value_ the new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    
    // Create NodeData
    int index = nodes_.size(); 
    std::set<EdgeData*> edges;
    NodeData* ndata = new NodeData{position, index, value, edges}; 
    nodes_.push_back(ndata);

    // Swap with first invalid node
    swap_nodes(index, numNodes);
    ++numNodes;

    // Create the node
    return Node(this, ndata);
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
    Node n = Node(this, nodes_[i]);
    assert(n.index() == i);
    return n;
  }

  /** Remove a node
   * @param[in] _node_ to be removed
   * See remove_node(NodeData*) for more details
   */
  size_type remove_node(const Node& node) {
    return remove_node(node.data_);
  }

  /** Remove a node using a node iterator
   * @param[in] node iterator pointing to node to be deleted 
   * @return node iterator pointing to the NEXT valid node, which is
   *         at the same position in nodes_ as before
   * See remove_node(NodeData*) for more details
   *
   * Example: Removing all nodes
   *   for (auto it = g.node_begin(); it != g.node_end(); g.remove_node(it));
   */
  node_iterator remove_node(node_iterator n_it) {
    node_type node = *n_it; 
    remove_node(node);
    return n_it;
  }

  /** Print all nodes (including removed ones) */
  void print_nodes_debug() const {
    std::cout << "\nPrinting " << num_nodes() << " nodes\n";
    for (size_t i = 0; i < nodes_.size(); i++) {
      NodeData* n = nodes_[i];
      Node node(this, n);
      std::cout << "  Node " << n->index << "(" << node.index() << ") val = "
                             << n->value << " degree = "
                             << node.degree() << " pos = "
                             << node.position();
      if (i >= num_nodes())
        std::cout << " (deleted)";
      std::cout << std::endl;
    }
  }

  /** Print all nodes and their positions */
  void print_positions() const {
    std::cout << "Printing Node Positions\n";
    for (auto it = node_begin(); it != node_end(); ++it) {
      node_type n = *it;
      std::cout << "  Node " << n.index() << " : " << n.position() << std::endl;
    }
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
    Edge() {}

    //
    // Getters
    //

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE (DONE)
      if (flip_)
        return Node(graph_, data_->n2);
      else
        return Node(graph_, data_->n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE (DONE)
      if (flip_)
        return Node(graph_, data_->n1);
      else
        return Node(graph_, data_->n2);
    }

    /** Return the resting length of the edge
     * */
    double& length() { return data_->length; };

    /** Return the distance between node1() and node2() */
    double distance() const {
      return norm(node1().position() - node2().position());
    }

    /** Return the edge value*/
    edge_value_type& value() { return data_->value; }
    const edge_value_type& value() const { return data_->value; }

    //
    // Comparators
    //

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Check if they're from the same graph
      if (graph_ != e.graph_)
        return false;

      int e11 = this->data_->n1->index;
      int e12 = this->data_->n2->index;

      int e21 = e.data_->n1->index;
      int e22 = e.data_->n2->index;

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
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // if (graph_ != e.graph_) {
      //   std::cerr << "WARNING: Comparing edges from two different graphs" 
      //             << std::endl;
      //   return false;
      // }
      int i1 = data_->idx;
      int i2 = e.data_->idx;
      long I1 = long(i1) + long(graph_);
      long I2 = long(i2) + long(e.graph_);
      return I1 < I2; 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class Node;

    const Graph* graph_; // pointer to the parent graph
    EdgeData* data_;     // pointer to the edge data
    bool flip_;          // reverse order of node1, node2 from data_
    
    /* Private Constructor */
    Edge(const Graph* graph, EdgeData* data)
        : graph_(graph), data_(data), flip_(false) {}

    /** Private Constructor 
     * Ensure that _n1_ corresponds to node1(), or throw InvalidEdge
     * */
    Edge(const Graph* graph, NodeData* n1, EdgeData* data)
        : graph_(graph), data_(data) {
      if (n1 == data_->n1)
        flip_ = false;
      else if (n1 == data_->n2)
        flip_ = true;
      else
        throw InvalidEdge{
            "Cannot create edge. Given edge does not match the given data"};
    } 
  };

  /** Return the total number of edges in the graph.
   * Complexity: O(1)
   */
  size_type num_edges() const { return numEdges; }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1) 
   */
  Edge edge(size_type i) const { return Edge(this, edgedata_[i]); }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *         Throw InvalidEdge exception if the edge is not symmetric
   *
   * Complexity: O(log(degree)), where degree is the maximum degree
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (a.is_neighbor(b) && b.is_neighbor(a)) {
      return true;
    } else if (a.is_neighbor(b)) {
      throw InvalidEdge{"Edge mismatch"};
    } else if (b.is_neighbor(a)) {
      throw InvalidEdge{"Edge mismatch"};
    } else {
      return false;
    }
  }

  /* Print all the index pairs of all the edges in the graph, in order */
  void print_edges() const {
    std::cout << "\nPrinting " << num_edges() << " edges" << std::endl;
    for (int i = 0; i < numEdges; i++) {
      EdgeData* e = edgedata_[i];
      std::cout << "  Edge " << e->idx << ": " 
                             << e->n1->index << " - "
                             << e->n2->index << " val = " 
                             << e->value << " idx = " 
                             << e->idx << std::endl;
    }
  }

  /** Print all edges, including removed edges */
  void print_edges_debug() const {
    std::cout << "\nPrinting " << num_edges() << " edges" << std::endl;
    for (size_t i = 0; i < edgedata_.size(); i++) {
      EdgeData* e = edgedata_[i];
      std::cout << "  Edge " << e->idx << ": " 
                             << e->n1->index << " - "
                             << e->n2->index << " val = " 
                             << e->value << " idx = " 
                             << e->idx;
      if (i >= size_t(numEdges))
        std::cout << " (deleted)";
      std::cout << std::endl;
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
   * Complexity: O(1)
   */
  Edge add_edge(const Node& a, const Node& b, 
      const edge_value_type& value = edge_value_type()) {

    // HW0: YOUR CODE HERE
    assert(a.graph_ == this);
    assert(b.graph_ == this);

    // Default length
    double len = a.distance(b);

    // Create EdgeData
    int totalEdges = edgedata_.size();
    EdgeData* edata = new EdgeData{a.data_, b.data_, value, len, totalEdges};

    // Create edge 
    Edge e = Edge(this, edata);

    if (has_edge(a, b)) {
      // std::cout << "Edge " << a.index() << "-" << b.index() << " already exists\n";
    } else {
      // std::cout << "Adding edge " << a.index() << "-" << b.index() << std::endl;
      
      // Add edge to list, swapping with first invalid edge
      edgedata_.push_back(edata);
      swap_edges(totalEdges, numEdges);
      numEdges++;

      // Add edge to node elements
      a.data_->edges.insert(edata);
      b.data_->edges.insert(edata);
    }

    return e;
  }


  /** Remove an edge
   * @param[in] _edge_ to be removed
   * See remove_edge(EdgeData*) for more details
   */
  size_type remove_edge(const Edge& edge) {
    EdgeData* edata = edge.data_;
    return remove_edge(edata);
  }

  /** Remove an edge between two edges
   * @param[in] _edge_ to be removed
   * See remove_edge(EdgeData*) for more details
   * Complexity: O(log(degree))
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // Find edge between _a_ and _b_, according to _a_
    auto it = find_edge(a.data_, b.data_); // O(log(degree))
    if (it == a.data_->edges.end()) {
      return false;
    } else {  // if found, delete
      EdgeData* edata = *it;
      return remove_edge(edata);
    }
  }

  /** Remove an edge using an edge iterator
   * @param[in] edge iterator pointing to edge to be deleted 
   * @return edge iterator pointing to the NEXT valid edge, which is
   *         at the same position in edgedata__ as before
   * See remove_node(EdgeData*) for more details
   *
   * Example: Removing all edges
   *   for (auto it = g.edge_begin(); it != g.edge_end(); g.remove_edge(it));
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    edge_type edge = *e_it;
    remove_edge(edge);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    // Clear all nodes
    for (size_type i = 0; i < nodes_.size(); i++)
      delete nodes_[i];
    nodes_.clear();

    // Clear all edges
    for (size_type i = 0; i < edgedata_.size(); i++)
      delete edgedata_[i];
    edgedata_.clear();

    numEdges = 0;
    numNodes = 0;
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
    const Graph* graph_;  // pointer to the parent graph

    // Iterator for the vector of node data in the parent graph
    typename std::vector<NodeData*>::const_iterator node_itr; 

    /** @brief Private constructor for NodeIterator 
     * @param[in] graph Pointer to the parent graph
     * @param[in] itr Iterator for the graph's vector of node data
     * */
    NodeIterator(const Graph* graph, 
                 typename std::vector<NodeData*>::const_iterator itr)
        : graph_(graph), node_itr(itr) {}
  };

  
  /** Beginning iterator for all nodes in graph */
  NodeIterator node_begin() const {
    return NodeIterator(this, nodes_.begin());
  }

  /** Terminal iterator for all nodes in graph */
  NodeIterator node_end() const {
    return NodeIterator(this, nodes_.begin() + numNodes);
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. 
   * Simply wraps an iterator to the set of EdgeData* stored by each node
   * */
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

    /** De-reference incident iterator to get the Edge */
    Edge operator*() const {
      return Edge(graph_, node_data_, *itr);
    }

    /** Increment incident iterator */
    IncidentIterator& operator++() {
      ++itr;
      return *this;
    }

    // Equality comparisons
    bool operator==(const IncidentIterator& i_itr) const {
      return (graph_ == i_itr.graph_) && (itr == i_itr.itr);
    }

    bool operator!=(const IncidentIterator& i_itr) const { 
      return !(*this == i_itr); 
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph_;    // Pointer to parent graph
    NodeData* node_data_;  // Pointer to node data

    // Iterator over set node (data) connected to the current node
    typename std::set<EdgeData*>::const_iterator itr;

    /** Private constructor for Incident Iterator */
    IncidentIterator(const Graph* graph, NodeData* data) 
        : graph_(graph), node_data_(data) {
      // Initialize with iterator at the beginning of the set
      itr = node_data_->edges.begin();
    }

    /** Generic Private Constructor */
    IncidentIterator(const Graph* graph, NodeData* data, 
                     typename std::set<EdgeData*>::const_iterator itr) 
        : graph_(graph), node_data_(data), itr(itr) {} 

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. 
   * Simply wraps an iterator for edgedata_
   * */
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

    /** De-reference edge iterator */
    Edge operator*() const { 
      return Edge(graph_, *itr); 
    }

    /** Increment edge iterator */
    EdgeIterator& operator++() {
      ++itr;
      return *this;
    }

    // Equality comparisons
    bool operator==(const EdgeIterator& e_itr) const {
      return (graph_ == e_itr.graph_) && (itr == e_itr.itr);
    }

    bool operator!=(const EdgeIterator& e_iter) const { 
      return !(*this == e_iter); 
    }

   private:
    friend class Graph;

    const Graph* graph_;  // pointers to parent graph

    // iterator for edgedata_
    typename std::vector<EdgeData*>::const_iterator itr;  

    /** Private Constructor */
    EdgeIterator(const Graph* graph, 
                 typename std::vector<EdgeData*>::const_iterator itr) 
      : graph_(graph), itr(itr) {}
  };
  
  /** Begining iterator for all edges in the graph */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edgedata_.begin());
  }

  /** Terminal iterator for all edges in the graph */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edgedata_.begin() + numEdges);
  }

  /** Print the current size of the graph */
  void print_size() const {
    std::cout << "Num Nodes: " << num_nodes() << " "
              << "Num Edges: " << num_edges() << std::endl;
  }


 private:

  //
  // Private data members
  //

  /** Internal storage representation of a node
   * Contains the following:
   * @param _position_ location of the node in 3D space
   * @param _index_ unique index identifier
   * @param _value_ node value (template type)
   * @param _edges_ std::set of pointers to others nodes to which the node is
   *    connected by an edge
   * */
  struct NodeData {
    Point position;
    int index;
    node_value_type value;
    std::set<EdgeData*> edges;
  };

  /** Internal storage representation of an edge
   * Contains the following:
   * @param _n1_ pointer to first node data
   * @param _n2_ pointer to second node data
   * @param _value_ edge value
   * @param _length_ nominal length of the edge
   * @param _idx_ index of the edge. Invariant: edgedata_[i]->idx == i
   * */
  struct EdgeData {
    NodeData* n1;
    NodeData* n2;
    edge_value_type value;
    double length;
    int idx;  // Index in edgedata_
  };
  
  int numEdges;  // number of edges
  int numNodes;  // number of nodes

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
  std::vector<EdgeData*> edgedata_;

  //
  // Private methods
  //

  /** Swap the edges at indices _i_ and _j_
   * @param[in] _i_, _j_ are valid indices to edgedata_ (includes removed edges)
   * @post edge(i) will refer to edge(j) prior to the call, and vice versa
   *       indices are update to maintain invariant edgedata_[i]->idx == i 
   * */
  void swap_edges(int i, int j) {
    auto itr_i = edgedata_.begin() + i;
    auto itr_j = edgedata_.begin() + j;
    std::iter_swap(itr_i, itr_j);
    EdgeData* edge_i = edgedata_[i];
    EdgeData* edge_j = edgedata_[j];
    std::swap(edge_i->idx, edge_j->idx);
  }

  /** Swap the nodes at indices _i_ and _j_
   * @param[in] _i_, _j_ are valid indices to nodes_ (includes removed nodes)
   * @post node(i) will refer to node(j) prior to the call, and vice versa
   *       indices are update to maintain invariant node(i).index() == i
   * */
  void swap_nodes(int i, int j) {
    auto itr_i = nodes_.begin() + i;
    auto itr_j = nodes_.begin() + j;
    std::iter_swap(itr_i, itr_j);
    NodeData* node_i = nodes_[i];
    NodeData* node_j = nodes_[j];
    std::swap(node_i->index, node_j->index);
  }

  /** Find if _b_ is in _a_'s adjacency list  
   * @param[in] _a_, _b_ NodeData
   * @return _it_ iterator to the set of edges in _a_.
   *         if there is an edge between _a_ and _b_, *it return the EdgeData*
   *         if there is no edge between _a_ and _b_, it == a->edges.end()
   * */
  typename std::set<EdgeData*>::iterator find_edge(const NodeData* a, 
                                                   const NodeData* b) const {
    auto it = a->edges.begin();
    for (; it != a->edges.end(); ++it) {
      EdgeData* e = *it;  
      if ((e->n1 == b) || (e->n2 == b))
        break;
    }
    return it;
  }

  /** Check if b is in a's adjacency list
   * @return true if there exists an edge between _a_ and _b_, according to _a_
   * */
  bool is_neighbor(const NodeData* a, const NodeData* b) const {
    auto it = find_edge(a, b);
    if (it == a->edges.end())
      return false;
    else
      return true;
  }

  /** Remove edge from graph
   * @param EdgeData* _e_, the edge to be removed
   * @post edgedata_.size() remains the same size
   * @post if _e_ has already been removed (e->idx >= numEdges), return false 
   *         with no other effect
   *       otherwise,
   *         numNodes is decremented by 1
   *         _e_ will be swapped with the last valid edge (at index numEdges-1) 
   *         _e_->idx == numEdges
   *         If i is the original index of _e_, edgedata_[i]->idx = i and
   *           corresponds to a valid edge
   *         All other edges remain unchanged
   * @return 1 if the edge is removed, 0 if no action is taken      
   *
   * Complexity: O(log(degree))
   * */
  size_type remove_edge(EdgeData* e) {
    // Check if the edge needs to be removed
    int idx = e->idx;
    if (idx >= numEdges) {
      std::cout << "Edge already deleted\n";
      // Make sure the edge doesn't exist in the node lists
      bool atob = is_neighbor(e->n1, e->n2);
      bool btoa = is_neighbor(e->n2, e->n1);
      assert(!atob);
      assert(!btoa);
      return 0;
    }

    swap_edges(idx, numEdges-1);  // O(1) complexity
    numEdges--;

    // Remove from adjacency lists
    std::set<EdgeData*>& s1 = e->n1->edges;
    std::set<EdgeData*>& s2 = e->n2->edges;    
    s1.erase(e);  // O(log(degree)) complexity
    s2.erase(e);

    return 1;
  }

  /** Remove a node from the graph
   * @param NodeData* _n_, the node to be removed
   * @post All edges connected to _n_ are removed (see remove_edge)
   *       numNodes is decremented by 1
   *       _n_ is swapped with the last valid node (at index numNodes-1)
   *       _n_->index == numNodes
   *       If i is the original index of _n_, nodes_[i]->index = i
   *       All other nodes remain unchanged
   * 
   * Complexity: O(degree * log(degree)) where degree is the max node degree
   * */
  size_type remove_node(NodeData* n) {
    // Remove adjacent edges
    std::set<EdgeData*> edges = n->edges; // copy so that iterators work
    
    // Delete all adjacent edges ( O(degree) )
    for (auto it = edges.begin(); it != edges.end(); ++it) {
      EdgeData* e = *it;
      remove_edge(e);  // O(log(degree))
    }

    // Swap with the last valid node 
    int idx = n->index;
    swap_nodes(idx, numNodes-1);

    // Increase size of deleted portion
    --numNodes;

    return true;
  }

};


#endif // CME212_GRAPH_HPP
