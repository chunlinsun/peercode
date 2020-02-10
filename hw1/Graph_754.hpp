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
template <typename V>
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

  //template adjustments
  using node_value_type = V;

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

    /** Return the value associated with this node.
     * 
     * @return this nodes value
     *
     * Complexity: O(1) amortized operations.
     */
    node_value_type& value() {
      return graph_->uid_to_val_[this->uid_];
    }

    /** Return the value associated with this node.
     *  Constant overwrite of above function.
     * @return this nodes value
     *
     * Complexity: O(1) amortized operations.
     */
    const node_value_type& value() const {
      return graph_->uid_to_val_[this->uid_];
    }

    /** Return an iterator across all edges incident to this node. 
     *  Returns the beginning value for iterator
     * 
     * @pre 
     * @return incident iterator for this node at the beginning
     * @return if no indcident edges, returns incident iterator s.t. begin() == end()
     *
     * Complexity: O(1) amortized operations.
     */
    IncidentIterator edge_begin() const {

      if (graph_->edge_by_node1_.find(uid_) == graph_->edge_by_node1_.end()) {
        return IncidentIterator(graph_, uid_, graph_->uid_to_idx_edge_.end());
      }
      return IncidentIterator(graph_, uid_, graph_->edge_by_node1_.at(uid_).begin());
    }

    /** Return an iterator across all edges incident to this node. 
     *  Returns the end value for iterator
     * 
     * @pre 
     * @return incident iterator for this node at the end
     * @return if no indcident edges, returns incident iterator s.t. begin() == end()
     *
     * Complexity: O(1) amortized operations.
     */
    IncidentIterator edge_end() const {
      if (graph_->edge_by_node1_.find(uid_) == graph_->edge_by_node1_.end()) {
        return IncidentIterator(graph_, uid_, graph_->uid_to_idx_edge_.end());
      }
      return IncidentIterator(graph_, uid_, graph_->edge_by_node1_.at(uid_).end());
    }

    /** Return the total number of edges incident to this node.
     * 
     * @return number of edges incident on this node
     *
     * Complexity: O(m) amortized operations, m is number of incident edges.
     */
    size_type degree() const {
      size_type count(0);
      for (auto it = this->edge_begin(); it != this->edge_end(); ++it) {
         count++;
       }
       return count;
     }

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
  Node add_node(const Point& position, const node_value_type& v = node_value_type()) {
    // HW0: YOUR CODE HERE

    //add the point itself and mappings to track the index
    uid_points_.insert(uid_points_.end(), {next_uid_, position});
    uid_to_idx_.insert(uid_to_idx_.end(), {next_uid_, num_nodes_});
    uid_to_val_.insert(uid_to_val_.end(), {next_uid_, v});

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
    assert(0 <= i);
    assert(i < num_nodes_);

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
  class Edge : private totally_ordered<Edge> {
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
      return Node(graph_, node1_);    
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    //--functionality_1
    //--The second comparison is supposed to check e.node2().
    //--START
    bool operator==(const Edge& e) const {        
      return ((this->node1() == e.node1() && this->node2() == e.node2()) ||
              (this->node2() == e.node1() && this->node1() == e.node1()));
    }
    //--END

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
    size_type node1_;
    size_type node2_;

    Edge(const Graph* graph, size_type uid, size_type node1, size_type node2) 
      : graph_(const_cast<Graph*>(graph)), uid_(uid), node1_(node1), node2_(node2) {
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
    assert(0 <= i);
    assert(i < num_edges_);
    //std::cout << i << std::endl;  
    size_type uid = idx_to_uid_edge_[i];         
    return Edge(this, uid, edge_by_uid_.at(uid)[0], edge_by_uid_.at(uid)[1]);        
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
      return Edge(this, edge_by_node1_.at(a.uid_).at(b.uid_), a.uid_, b.uid_);
    } else {
      //update all the containers
      //store the edge by node
      edge_by_node1_[a.uid_][b.uid_] = next_uid_edge_;
      edge_by_node1_[b.uid_][a.uid_] = next_uid_edge_;

      //store the mapping of edge index to edge uid
      uid_to_idx_edge_[next_uid_edge_] = num_edges_;
      idx_to_uid_edge_.push_back(next_uid_edge_);

      //store the edge by uid
      std::array<size_type, 2> tmp {{a.uid_, b.uid_}};
      edge_by_uid_.insert(edge_by_uid_.end(), {next_uid_edge_, tmp});
      
      //update state variables
      num_edges_++;
      next_uid_edge_++;

      return Edge(this, next_uid_edge_ - 1, a.uid_, b.uid_);
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

    /** Return the node associated with the current iterator position.
     * 
     * @return Node associated with this iterators uid/index
     *
     * Complexity: O(1) amortized operations.
     */
    Node operator*() const {
      size_type index = graph_->uid_to_idx_[uid_];
      return graph_->node(index);
    }

    /** Increment the iterator position.
     * 
     * @return Node iterator at the next uid in the graph
     *
     * Complexity: O(1) amortized operations.
     */
    NodeIterator& operator++() {
      uid_++;
      return *this;
    }

    /** Equality check with another iterator. True if all iterator
     *  fields are equal (same graph and at same position.)
     * 
     * @return equality boolean
     *
     * Complexity: O(1) amortized operations.
     */
    bool operator==(const NodeIterator& ni) const {
      return (graph_ == ni.graph_ && uid_ == ni.uid_);
    }

    /** Inquality check with another iterator. True if all iterator
     *  fields are equal (same graph and at same position.)
     * 
     * @return Inequality boolean
     *
     * Complexity: O(1) amortized operations.
     */
    bool operator!=(const NodeIterator& ni) const {
      return !(*this == ni);
    }

   private:
    friend class Graph;

    // HW1 #2: YOUR CODE HERE

    /** Private constructor to be used within graph class
     *  
     * @param graph to be associated with
     * @param uid of node to be associated with
     */
    NodeIterator(const Graph* graph, size_type uid) 
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    
    Graph* graph_;
    size_type uid_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return a node iterator at the beginning of this graph.
     * 
     * @return node iterator pointing to first element added
     *
     * Complexity: O(1) amortized operations.
     */
  NodeIterator node_begin() const {
    return NodeIterator(this, size_type(0));
  }

  /** Return a node iterator at the end of this graph.
     * 
     * @return node iterator pointing to one past the 
     *          most recent element added
     *
     * Complexity: O(1) amortized operations.
     */
  NodeIterator node_end() const {
    return NodeIterator(this, next_uid_);
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

    /** Returns the edge associated with the current iterator position.
     * 
     * @pre ensures the edge iterator points to a valid edge
     *        *this condition is helpful in case node associated
     *          with edge iterator has no incident edges
     * @return edge at current position
     *
     * Complexity: O(1) amortized operations.
     */
    Edge operator*() const {
      assert(it_ != graph_->uid_to_idx_edge_.end());
      return Graph::Edge(graph_, it_->second, node_uid_, it_->first);
    }

    /** Returns the an incident iterator at the next position
     *  in this nodes set of edges.
     * 
     * @return iterator at next position
     *
     * Complexity: O(1) amortized operations.
     */
    IncidentIterator& operator++() {
      it_++;
      return *this;
    }

    /** Equality check between two iterators. Equal if at some position and 
     *  for same node and in same graph.
     * 
     * @return true if equal false if not
     *
     * Complexity: O(1) amortized operations.
     */
    bool operator==(const IncidentIterator& it2) const {
      return (graph_ == it2.graph_ && node_uid_ == it2.node_uid_ && it_ == it2.it_);
    }

    /* not equal version of equality check */
    bool operator!=(const IncidentIterator& it2) const {
      return !(*this == it2);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    //state variables
    Graph* graph_;
    size_type node_uid_;
    std::unordered_map<size_type, size_type>::iterator it_;

    /** Private constructor to be used within graph class
     *  
     * @param graph to be associated with
     * @param node_uid to be associated with
     * @param iterator for set of incident edges
     */
    IncidentIterator(const Graph* graph, size_type node_uid, std::unordered_map<size_type, size_type>::iterator it) 
      : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid), it_(it) {
    }
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

    /** Dereferencing operator overload. Current Edge
     *  
     * @return edge at current iterator position
     */
    Edge operator*() const {
      return Graph::Edge(graph_, uid_, graph_->edge_by_uid_.at(uid_)[0], graph_->edge_by_uid_.at(uid_)[1]);
    }

    /** Increment operator overload. Moves iterator to next edge.
     *  
     * @return Edge Iterator at next position.
     */
    EdgeIterator& operator++() {
      uid_++;
      return *this;
    }

    /** Equality operator overloading
     *  
     * @return true if same graph and same edge associated with iterator
     */
    bool operator==(const EdgeIterator& ei) const {
      return (graph_ == ei.graph_ && uid_ == ei.uid_);
    }

    /* Not equal operator overload */
    bool operator!=(const EdgeIterator& ei) const {
      return !(*this == ei);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    /** Private constructor to be used within graph class
     *  
     * @param graph to be associated with
     * @param uid edge to be associated with
     */
    EdgeIterator(const Graph* graph, size_type uid) 
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    
    Graph* graph_;
    size_type uid_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Beginning iterator to move across edges.
     *  
     * @return EdgeIterator at beginning of all edges
     */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, size_type(0));
  }

  /** End iterator for edges
     *  
     * @return Iterator pointing one past last element
     */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, next_uid_edge_);
  }

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

  //map node uid to value for template version
  std::unordered_map<size_type, node_value_type> uid_to_val_;

  //state variable for size of graph
  size_type num_nodes_;

  //map from "node in edge" to "other node in edge" to edge
  //for fast lookup of edge by node members
  //two entries for each edge (a->b and b->a) for user flexibility,
  std::unordered_map<size_type, std::unordered_map<size_type, size_type> > edge_by_node1_;
  
  //map from edge uid to edge (stored as size 2 array, ex: [a, b])
  std::unordered_map<size_type, std::array<size_type, 2> > edge_by_uid_;

  //map from edge uid to edge index 
  //for fast lookup of edge index
  std::unordered_map<size_type, size_type> uid_to_idx_edge_;

  //edge index to edge uid
  //vector implementation to maintain size/index ordering
  //only maps to first edge stored (a->b)
  std::vector<size_type> idx_to_uid_edge_;

  //state variables for number of edges and maintaining uid ordering
  size_type num_edges_;
  size_type next_uid_;
  size_type next_uid_edge_;
};

#endif // CME212_GRAPH_HPP
