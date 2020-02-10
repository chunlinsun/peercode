#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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
  
  std::vector<Point> points;  //stores our vector of points objects
  std::map<unsigned, std::vector<unsigned>> edge_key;   //stores a map with the edge index as the key, lower node index is first one in the tuple
  std::map<unsigned, std::vector<unsigned>> node_idx_key;  //stores a map with the node index as the key, edges are stores both ways (low: high) and (high:low)
  std::vector<V> values;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  using node_value_type = V;
 
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

  class Node: private totally_ordered<Node>{
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

    size_type idx;  //the index of the node
    Graph* g;  //pointer to the Graph (proxy programming)

    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g->points.at(idx);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return size_type(idx);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the value at a certain node.
     * @return A node_value_type reference with the value at a node.
     *
     * @post For all i, i must be a valid value.
     */

    node_value_type& value() {
    	return g->values.at(idx);
    };

    /** Return the const value at a certain node.
     * @return A const node_value_type reference with the value at a node.
     *
     * @post For all i, i must be a valid value and const.
     */

    const node_value_type& value() const {
    	return g->values.at(idx);
    };
    
    /** Return the number of incident edges
     * @return A size_type value that is the number of incident edges from a Node.
     *
     * @post For all i, 0 <= i <= g->size().
     */

    size_type degree() const {
      return g->node_idx_key[idx].size();
    };

    /** Start of the incident iterator.
     * @return An IncidentIterator of the first Edge to be iterated on.
     *
     * @post For all IncidentIterator, the index = 0.
     */

    incident_iterator edge_begin() const {
      return IncidentIterator((*this), 0);
    };
    
    /** End of the incident iterator.
     * @return An IncidentIterator one past the last Edge to be iterated on.
     *
     * @post For all IncidentIterator, the index = this->degree().
     */

    incident_iterator edge_end() const {
      return IncidentIterator((*this), this->degree());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */

    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (n.g == this->g && n.index() == this->index()) {
        return true; //return true if the indices are equal and in the same graph
      }
      else {
        return false; //return false if the indices are not equal or different graphs
      }
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
      if (n.g == this->g) {
        if (this->index() < n.index()) {
          return true; //return true if the index is less than n's index and in the same graph
        }
      }
      else {
        return std::less<Graph* >{}(this->g, n.g); //used for comparing pointers for different graphs
      }
      return false; //return false if same graph but index is greater to or equal
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(size_type gindex, Graph* gpointer){ //private Node constructor
      g = gpointer;
      idx = gindex;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return points.size();
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
    // HW0: YOUR CODE HERE
    Node new_node(num_nodes(), this);
    points.push_back(position);
    values.push_back(value);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.g && this->num_nodes() >= n.idx) {
        return true;
    }
    else {
        return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_nodes()) {
      Node result_node(i, const_cast<Graph*>(this));
      return result_node;
    }
    else {
      return Node();
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
    
    Graph* g;
    size_type node1_idx;
    size_type node2_idx;
    
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(node1_idx, g);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(node2_idx, g);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (e.g == this->g && e.node1_idx == this->node1_idx && e.node2_idx == this->node2_idx) {
        return true; //return true if same graph and both node indices are equal
      }
      else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (e.g == this->g) {
        if (std::min(this->node1_idx, this->node2_idx) < std::min(e.node1_idx, e.node2_idx)) {
          return true; //return true if same graph and the minimum node index of this is less than the minimum node index of e
        }
        else if (std::min(this->node1_idx, this->node2_idx) == std::min(e.node1_idx, e.node2_idx)) {
            return std::max(this->node1_idx, this->node2_idx) < std::max(e.node1_idx, e.node2_idx);
            //if minimum node index is the same and same graph, check for the other node index
        }
      }
      else { 
        return std::less<Graph* >{}(this->g, e.g); //for when we are comparing two different graphs
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(size_type e_n1_idx, size_type e_n2_idx, Graph* epointer) { //private constructor for Edge
        node1_idx = e_n1_idx;
        node2_idx = e_n2_idx;
        g = epointer;
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_key.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_edges()){
      Edge result_edge(edge_key.at(i)[0], edge_key.at(i)[1], const_cast<Graph* >(this));
      return result_edge;
    }
    else {
      return Edge();
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type lower_idx = std::min(a.idx, b.idx); //our edges are added with the lower node index first
    size_type higher_idx = std::max(a.idx, b.idx);
    if (node_idx_key.count(lower_idx)) {
      std::vector<unsigned> vec_at_index = node_idx_key.at(lower_idx);
      if (std::find(vec_at_index.begin(), vec_at_index.end(), higher_idx) != vec_at_index.end()) {
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    size_type lower_idx = std::min(a.idx, b.idx); //our edges are added with the lower node index first
    size_type higher_idx = std::max(a.idx, b.idx);

    if (has_edge(a, b)) { //check to see if current edge already exists
      Edge same_edge(lower_idx, higher_idx, this);
      return same_edge;
    }
    else {
      edge_key.insert(std::pair<size_type, std::vector<size_type>>(num_edges(), {lower_idx, higher_idx}));
      std::vector<unsigned> low_idx_vec_at_index;
      std::vector<unsigned> high_idx_vec_at_index;

      if (node_idx_key.count(lower_idx)) {
        low_idx_vec_at_index = node_idx_key.at(lower_idx);
      }

      low_idx_vec_at_index.push_back(higher_idx);
      node_idx_key[lower_idx] = low_idx_vec_at_index;

      if (node_idx_key.count(higher_idx)) {
        high_idx_vec_at_index = node_idx_key.at(higher_idx);
      }

      high_idx_vec_at_index.push_back(lower_idx);
      node_idx_key[higher_idx] = high_idx_vec_at_index;

      Edge new_edge(a.idx, b.idx, this);
      return new_edge;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    points.clear();
    edge_key.clear();
    node_idx_key.clear();
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

//forward_iterator_tag
    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Deference a NodeIterator.
     * @return The Node that the iterator is pointing to.
     *
     * @post For all returned Node, Node index < size of items in the NodeIterator.
     */
    Node operator*() const {
    	return Node(idx, g);
    }

    /** Increment a NodeIterator.
     * @return The NodeIterator that has an index one larger than before.
     *
     * @post For all returned Node, Node index < size of items in the NodeIterator.
     */
    NodeIterator& operator++() {
    	idx = idx + 1;
    	return (*this);
    }

    /** Tests whether this NodeIterator and node_iter are equal.
     * @return A bool of 1 if they are equal, 0 if not. 
     * Checks to see if the same graph pointer and also the same index.
     *
     * @post All bool 0 or 1.
     */

    bool operator==(const NodeIterator& node_iter) const {
    	return (this->g == node_iter.g && this->idx == node_iter.idx);
    }

    /** Tests whether this NodeIterator and node_iter are not equal.
     * @return A bool of 1 if not equal, 0 if they are. 
     * Checks to see if the same graph pointer and also the same index.
     *
     * @post All bool 0 or 1.
     */

    bool operator!=(const NodeIterator& node_iter) const {
    	return !(this->g == node_iter.g && this->idx == node_iter.idx);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    unsigned int idx;
    Graph* g;

    NodeIterator(Graph* gpointer, unsigned int index) {
    	idx = index;
    	g = gpointer;
    }

    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

    /** Start of the node iterator.
     * @return A NodeIterator of the first Node to be iterated on.
     *
     * @post For all NodeIterator, the index = 0.
     */

  	node_iterator node_begin() const {
  	  return NodeIterator(const_cast<Graph*>(this), 0);
  	}

    /** End of the node iterator.
     * @return A NodeIterator one past the last Node to be iterated on.
     *
     * @post For all NodeIterator, the index = this->size(), the number of nodes in the graph.
     */

    node_iterator node_end() const {
      return NodeIterator(const_cast<Graph*>(this), this->size());
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

    /** Deference an IncidentIterator.
     * @return The Edge that the iterator is pointing to.
     *
     * @post For all returned Edge, Edge index < size of items in the IncidentIterator.
     */

    Edge operator*() const {
      return Edge(n.idx, n.g->node_idx_key.at(n.idx).at(idx), n.g);
    }

    /** Increment an IncidentIterator.
     * @return The IncidentIterator that has an index one larger than before.
     *
     * @post For all returned Edge, Edge index < size of items in the IncidentIterator.
     */

    IncidentIterator& operator++() {
      idx = idx + 1;
      return (*this);
    }

    /** Tests whether this IncidentIterator and inc_iter are equal.
     * @return A bool of 1 if they are equal, 0 if not. 
     * Checks to see if the same node and also the same index.
     *
     * @post All bool 0 or 1.
     */

    bool operator==(const IncidentIterator& inc_iter) const {
      return (this->n == inc_iter.n && this->idx == inc_iter.idx);
    }

    /** Tests whether this IncidentIterator and inc_iter are not equal.
     * @return A bool of 1 if not equal, 0 if they are. 
     * Checks to see if the same node and also the same index.
     *
     * @post All bool 0 or 1.
     */

    bool operator!=(const IncidentIterator& inc_iter) const {
      return !(this->n == inc_iter.n && this->idx == inc_iter.idx);
    }
  
   private:
    friend class Graph;
    Node n;
    size_type idx;

    // HW1 #3: YOUR CODE HERE
    IncidentIterator(Node node, size_type index) {
      idx = index;
      n = node;
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

    /** Deference an EdgeIterator.
     * @return The Edge that the iterator is pointing to.
     *
     * @post For all returned Edge, Edge index < size of items in the EdgeIterator.
     */

    Edge operator*() const {
      return (Edge(g->edge_key.at(idx)[0], g->edge_key.at(idx)[1], g));
    }

    /** Increment an EdgeIterator.
     * @return The EdgeIterator that has an index one larger than before.
     *
     * @post For all returned Edge, Edge index < size of items in the EdgeIterator.
     */

    EdgeIterator& operator++() {
      idx = idx + 1;
      return (*this);
    }
    
    /** Tests whether this EdgeIterator and edge_iter are equal.
     * @return A bool of 1 if they are equal, 0 if not. 
     * Checks to see if the same node and also the same index.
     *
     * @post All bool 0 or 1.
     */

    bool operator==(const EdgeIterator& edge_iter) const {
      return (this->g == edge_iter.g && this->idx == edge_iter.idx);
    }

    /** Tests whether this EdgeIterator and edge_iter are not equal.
     * @return A bool of 1 if not equal, 0 if they are. 
     * Checks to see if the same node and also the same index.
     *
     * @post All bool 0 or 1.
     */

    bool operator!=(const EdgeIterator& edge_iter) const {
      return !(this->g == edge_iter.g && this->idx == edge_iter.idx);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* g;
    unsigned int idx;

    // HW1 #3: YOUR CODE HERE
    EdgeIterator(Graph* gpointer, unsigned int index) {
      idx = index;
      g = gpointer;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Start of the edge iterator.
     * @return An EdgeIterator of the first Edge to be iterated on.
     *
     * @post For all EdgeIterator, the index = 0.
     */
  edge_iterator edge_begin() const {
    return EdgeIterator(const_cast<Graph*>(this), 0);
  }
  
  /** End of the edge iterator.
     * @return An EdgeIterator one past the last Edge to be iterated on.
     *
     * @post For all EdgeIterator, the index = this->edge_key.size(), the size of all the edges in the graph.
     */
  edge_iterator edge_end() const {
    return EdgeIterator(const_cast<Graph*>(this), this->edge_key.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};
//--functionality_0
//--Great job!
//--END
#endif // CME212_GRAPH_HPP
