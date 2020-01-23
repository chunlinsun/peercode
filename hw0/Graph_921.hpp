#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <array>
#include <functional>
#include <unordered_map>

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
  Graph() 
    : uid_points_(), idx_to_uid_(), uid_to_idx_(), num_nodes_(0), edge_by_node1_(), 
    edge_by_uid_(), uid_to_idx_edge_(), idx_to_uid_edge_(), num_edges_(0), next_uid_(0), next_uid_edge_(0) {
    // HW0: YOUR CODE HERE
    //simply initialize all variables to 0 or empty
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
      //empty constructor
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      
      //positions are stored by node uid in graph class
      return graph_->uid_points_[this->uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE

      //indices are stored in graph container
      assert(graph_->uid_to_idx_.count(uid_));
      return graph_->uid_to_idx_.at(uid_);
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

      //check if uid_'s are equal and in same graph
      return (uid_ == n.uid_ && this->graph_ == n.graph_);
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

      //uid track's global order
      //if same graph, compare graph pointers using std::less 
      if (uid_ == n.uid_ && this->graph_ != n.graph_) {
        return std::less<Graph*>{}(this->graph_, n.graph_);
      } else {
        return uid_ < n.uid_;
      }     
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    //Node maintains a unique id (within this graph) and graph pointer
    Graph* graph_;
    size_type uid_;

    Node(const Graph* graph, size_type uid) 
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE

    //num_nodes tracks size
    return num_nodes_;
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

    //add the point itself and mappings to track the index
    uid_points_.insert(uid_points_.end(), {next_uid_, position});
    uid_to_idx_.insert({next_uid_, num_nodes_});

    //currently only add to the end of the index tracker
    idx_to_uid_.push_back(next_uid_);

    //increment state variables
    num_nodes_++;
    next_uid_++;

    //return Node object with correct uid
    return Node(this, next_uid_ - 1);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    //obtain the index in O(1) amortized time and make sure part of appropriate graph
    assert(uid_to_idx_.count(n.uid_));
    return (uid_to_idx_.at(n.uid_) < num_nodes_ && n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

    //check pre-condition
    assert(0 <= i < num_nodes_);

    //return a node with appropriate uid stored in graph container    
    return Node(this, idx_to_uid_.at(i));        
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
      //empty constructor
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE

      //edges obtained by uid_ and contain node uid_s in size two array
      return Node(graph_, graph_->edge_by_uid_.at(uid_)[0]);    
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, graph_->edge_by_uid_.at(uid_)[1]);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {        
      return (this->node1() == e.node1() && this->node2() == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {      
      if (this->uid_ == e.uid_ && this->graph_ != e.graph_) {
        return std::less<Graph*>{}(this->graph_, e.graph_);
      } else {
        return this->uid_ < e.uid_;
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
    size_type uid_;

    Edge(const Graph* graph, size_type uid) 
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
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

    //pre-condition
    assert(0 <= i < num_edges_);
    //std::cout << i << std::endl;           
    return Edge(this, idx_to_uid_edge_[i]);        
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
    if(edge_by_node1_.count(a.uid_)) {
      return edge_by_node1_.at(a.uid_).count(b.uid_);
    } else {
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
    assert(!(a == b));
    assert(has_node(a) && has_node(b));
    if (has_edge(a, b)) {
      return Edge(this, edge_by_node1_.at(a.uid_).at(b.uid_));
    } else {

      //update all the containers
      //store the edge by node
      edge_by_node1_[a.uid_][b.uid_] = num_edges_;
      edge_by_node1_[b.uid_][a.uid_] = num_edges_;

      //store the mapping of edge index to edge uid
      uid_to_idx_edge_[next_uid_edge_] = num_edges_;
      idx_to_uid_edge_.push_back(next_uid_edge_);

      //store the edge by uid
      std::array<size_type, 2> tmp {{a.uid_, b.uid_}};
      edge_by_uid_.insert(edge_by_uid_.end(), {next_uid_edge_, tmp});
      
      //update state variables
      num_edges_++;
      next_uid_edge_++;

      return Edge(this, next_uid_edge_ - 1);
    }     
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    //clear all variables
    num_nodes_ = 0;
    num_edges_ = 0;
    next_uid_edge_ = 0;
    next_uid_ = 0;
    uid_points_.clear();
    uid_to_idx_.clear();
    idx_to_uid_.clear();
    edge_by_node1_.clear();
    edge_by_uid_.clear();
    idx_to_uid_edge_.clear();
    uid_to_idx_edge_.clear();
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

  //map from node uid_ to point to maintain point data
  std::unordered_map<size_type, Point> uid_points_;

  //node index to node uid for fast node lookup
  //vector implementation to maintain graph size/index ordering
  std::vector<size_type> idx_to_uid_;

  //map node uid to index for fast index lookup
  std::unordered_map<size_type, size_type> uid_to_idx_;

  //state variable for size of graph
  size_type num_nodes_;

  //map from "node in edge" to "other node in edge" to edge
  //for fast lookup of edge by node members
  //two entries for each edge (a->b and b->a) for user flexibility
  std::unordered_map<size_type, std::unordered_map<size_type, size_type> > edge_by_node1_;
  
  //map from edge uid to edge (stored as size 2 array, ex: [a, b])
  std::unordered_map<size_type, std::array<size_type, 2> > edge_by_uid_;

  //map from edge uid to edge index 
  //for fast lookup of edge index
  std::unordered_map<size_type, size_type> uid_to_idx_edge_;

  //edge index to edge uid
  //vector implementation to maintain size/index ordering
  std::vector<size_type> idx_to_uid_edge_;

  //state variables for number of edges and maintaining uid ordering
  size_type num_edges_;
  size_type next_uid_;
  size_type next_uid_edge_;
};

#endif // CME212_GRAPH_HPP
