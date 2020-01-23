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
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
 
 // Intern_Node is defined at the end of this code. 
 struct Intern_Node ;  

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): nodes_(), edge_adj_() ,edges_idx_(), num_edges_(0)
            {
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
  class Node {
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
      return node_graph_->nodes_[node_idx_].position_ ;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // We return the attribute `node_idx_` from the `Node` class.
      return node_idx_ ; 
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
      (void) n;          // Quiet compiler warning
      bool res = false ; 
      // We compare both the pointers value and the index of the nodes.
      if ((n.node_graph_ == node_graph_) && (n.node_idx_ == node_idx_)){res = true ;}
          return res;
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
      // Two cases are possible : if the nodes belong to the same graph,
      // we compare their IDs. 
      if (n.node_graph_ == node_graph_){
          return node_idx_ < n.node_idx_ ; }
      // Otherwise we compare the address of the graphs.
      return node_graph_ < n.node_graph_ ;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
   
    // We declare a pointer to the graph associated with the node.
    Graph* node_graph_ ; 
   
    // Node ID for the node.  
    size_type node_idx_ ; 
   
    // A private constructor for Node class. 
    Node(const Graph* node_graph, size_type node_idx)
     : node_graph_(const_cast<Graph*>(node_graph)), node_idx_(node_idx){
     }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // We return the size of the container `nodes_`. 
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    (void) position;     // Quiet compiler warning
    // The new node ID.
    size_type new_node_id = nodes_.size();
    // We declare an `Intern_Node` => the new node. 
    Intern_Node newNode = Intern_Node(new_node_id, position) ; 
    // We update the nodes container. 
    nodes_.push_back(newNode);
    // We update the adjacency vector. 
    edge_adj_.push_back(std::vector<size_type>());
    // We return the new node. 
    return Node(this, new_node_id);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    // We test for the equality of the pointers value. 
    return n.node_graph_ == this ;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    return Node(this, i);        // Invalid node
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(edge_graph_, node1_idx_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(edge_graph_, node2_idx_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
       (void) e;           // Quiet compiler warning
      // The edge is assumed to be undirected. 
      // We test for the nodes index equality. 
      bool nodes_match = (e.node1_idx_ == node1_idx_ && e.node2_idx_ == node2_idx_)
                      || (e.node1_idx_ == node2_idx_ && e.node2_idx_ == node1_idx_);
      // We test for the pointers value equality. 
      bool graphs_match = e.edge_graph_ == edge_graph_;

      return (nodes_match && graphs_match);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e ; // Quiet compiler warning.
      // HW0: YOUR CODE HERE

      // Do th edges belong to the same graph? 
      if (edge_graph_ ==e.edge_graph_){
      // Test if the nodes are equal. 
      if (e.node1_idx_ == node1_idx_){
      // then we compare the other pair of nodes. 
          return node2_idx_ < e.node2_idx_;
        }
        // Then compare first nodes
        return node1_idx_ < e.node1_idx_;
      }
      // Otherwise, let us compare the pointer value. 
      return edge_graph_ < e.edge_graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
   
    // We define a pointer to the Graph to which the edge belongs.
    Graph* edge_graph_ ; 
   
   // Two nodes are asociated t an edge, we identify them with their ID.
    size_type node1_idx_ ; 
    size_type node2_idx_ ; 
   
    // Constructor for `Edge` class.
    Edge(const Graph* edge_graph, size_type node1_idx, size_type node2_idx)
     : edge_graph_(const_cast<Graph*>(edge_graph)), node1_idx_(node1_idx),
       node2_idx_(node2_idx){
       }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // We return the attribute of Graph `num_edges_`. 
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    assert(i<num_edges()) ; 
    std::pair<size_type,size_type> edge_nodes_pair = edges_idx_[i];
    return Edge(this, edge_nodes_pair.first, edge_nodes_pair.second);   // Invalid Edge
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
    // We check if the nodes are in the graph. 
    assert(has_node(a)) ; 
    assert(has_node(b)) ; 
    // We define the neighbors of node a. 
    std::vector<size_type> neighbors = edge_adj_[a.node_idx_];
    // We go through each neighbor of node a. 
    for(size_type l = 0; l < neighbors.size(); l++){
        if (neighbors[l] == b.node_idx_){return true;}
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    (void) a, (void) b;   // Quiet compiler warning
   
    // a) Check if edge (a,b) already exists. 
    // We stock the neighbor of node a.
    std::vector<size_type> neighbors = edge_adj_[a.node_idx_];
    for(size_type l = 0; l < neighbors.size(); l++){
     if (neighbors[l] == b.node_idx_){
        // We return the existing edge. 
        return Edge(this, a.node_idx_, b.node_idx_);
      }
    }
   
    // b) In this case, we add the new edge to the graph. 
    // We update the container for edge adjencies.
    edge_adj_[a.node_idx_].push_back(b.node_idx_);
    edge_adj_[b.node_idx_].push_back(a.node_idx_);

    // We add a new edge index. The pair is ordered with respect to the index.
    size_type min_node_idx = std::min(a.node_idx_,b.node_idx_);
    size_type max_node_idx = std::max(a.node_idx_,b.node_idx_);
    std::pair<size_type,size_type> new_node_idx = std::make_pair(min_node_idx,max_node_idx);
    edges_idx_.push_back(new_node_idx);

    // Updated num. of edges
    num_edges_ ++;
   
    return Edge(this, a.node_idx_, b.node_idx_); // Invalid Edge     
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // We clear the atributes of the clas `Graph`. 
    nodes_.clear() ; 
    edges_idx_.clear() ;
    num_edges_ = 0 ; 
    edge_adj_.clear();
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

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

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
 
  // A struct that enables the storing of the data related to each node. 
 struct Intern_Node {
  
     // Node Index
     size_type node_idx_ ; 
  
     // Space position 
     Point position_ ; 
 
     // Private Constructor for Intern_Node
     Intern_Node(const size_type node_idx, const Point& position) 
       : node_idx_(node_idx), position_(position) {
       }
  } ; 
  
 // A container from STL that stocks every node in the graph.
 std::vector<Intern_Node> nodes_ ;  
 
 // A container for edges adjacency data. 
 std::vector<std::vector<size_type>> edge_adj_;
 
 // A container for edges index. 
 std::vector<std::pair<size_type, size_type>> edges_idx_ ; 
 
 
 // The total number of edges in the graph.
 size_type num_edges_ ; 
  


};

#endif // CME212_GRAPH_HPP
