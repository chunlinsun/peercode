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

template <typename V, typename E>
class Graph {
 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  /** A struct to hold the attributes of a Node and more information
      Contains uuid, ext_idx, point, value, and g */
  struct InternalNode {
    unsigned uuid;
    unsigned ext_idx;
    Point point;
    V value;
    Graph* g;

    /** Constructor for internal node */
    InternalNode(unsigned ext_id, Point pt, V val, Graph* gpointer) {
      uuid = gpointer->internal_nodes.size();
      ext_idx = ext_id;
      point = pt;
      value = val;
      g = gpointer;
    }

  };

  /** A struct to hold the attributes of an Edge and more information
      Contains uuid, ext_idx, node_1_idx, node_2_idx, value */
  struct InternalEdge {
    unsigned uuid;
    unsigned ext_idx;
    unsigned node_1_idx;
    unsigned node_2_idx;
    E value;

    /** Constructor for internal edge */
    InternalEdge(unsigned ext_id, unsigned a_idx, unsigned b_idx, E val, Graph* gpointer) {
      uuid = gpointer->internal_edges.size();
      ext_idx = ext_id;
      node_1_idx = a_idx;
      node_2_idx = b_idx;
      value = val;
    }
  };

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  using node_value_type = V;
  using edge_value_type = E;
 
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

    size_type ext_idx;  //the index of the node
    Graph* g;  //pointer to the Graph (proxy programming)

    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g->internal_nodes.at(g->node_i_to_u.at(ext_idx)).point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return ext_idx;
    }

    Point& position() {
      // HW0: YOUR CODE HERE
      return g->internal_nodes.at(g->node_i_to_u.at(ext_idx)).point;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the value at a certain node.
     * @return A node_value_type reference with the value at a node.
     *
     * @post For all i, i must be a valid value.
     */

    node_value_type& value() {
    	return g->internal_nodes.at(g->node_i_to_u.at(ext_idx)).value;
    };

    /** Return the const value at a certain node.
     * @return A const node_value_type reference with the value at a node.
     *
     * @post For all i, i must be a valid value and const.
     */

    const node_value_type& value() const {
    	return g->internal_nodes.at(g->node_i_to_u.at(ext_idx)).value;
    };
    
    /** Return the number of incident edges
     * @return A size_type value that is the number of incident edges from a Node.
     *
     * @post For all i, 0 <= i <= g->size().
     */

    size_type degree() const {
      return g->adj_map.at(g->node_i_to_u.at(ext_idx)).size();
    };

    /** Start of the incident iterator.
     * @return An IncidentIterator of the first Edge to be iterated on.
     *
     * @post For all IncidentIterator, the index = 0.
     */

    incident_iterator edge_begin() const {
      return IncidentIterator((*this), g->adj_map.at(g->node_i_to_u.at(ext_idx)).begin());
    };
    
    /** End of the incident iterator.
     * @return An IncidentIterator one past the last Edge to be iterated on.
     *
     * @post For all IncidentIterator, the index = this->degree().
     */

    incident_iterator edge_end() const {
      return IncidentIterator((*this), g->adj_map.at(g->node_i_to_u.at(ext_idx)).end());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */

    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->g == n.g  && this->index() == n.index()) {
        return true; //return true if the indices are equal and in the same graph
      }
      return false; //return false if the indices are not equal or different graphs
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
      if (this->g == n.g) {
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
    Node(size_type n_index, Graph* gpointer){ //private Node constructor
      ext_idx = n_index;
      g = gpointer;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_i_to_u.size();
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
    //*****check for has_edge

    size_type ext_idx = num_nodes();
    size_type uuid = internal_nodes.size();

    Node new_node(ext_idx, this);

    internal_nodes.push_back(InternalNode(ext_idx, position, value, this));
    adj_map.insert(std::pair<size_type, std::map<size_type, size_type>>(ext_idx, {}));

    node_i_to_u.push_back(uuid);

    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.g && this->num_nodes() >= n.ext_idx) {
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
    
    size_type node1_idx_uuid;
    size_type node2_idx_uuid;
    size_type ext_idx;
    Graph* g;
    
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(g->internal_nodes.at(node1_idx_uuid).ext_idx, g);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g->internal_nodes.at(node2_idx_uuid).ext_idx, g);    
    }

    edge_value_type& value() {
      return g->internal_edges.at(g->edge_i_to_u.at(ext_idx)).value;
    };

    const edge_value_type& value() const {
      return g->internal_edges.at(g->edge_i_to_u.at(ext_idx)).value;
    };

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (e.g == this->g && e.node1_idx_uuid == this->node1_idx_uuid && e.node2_idx_uuid == this->node2_idx_uuid && this->ext_idx == e.ext_idx) {
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
        if (std::min(this->node1_idx_uuid, this->node2_idx_uuid) < std::min(e.node1_idx_uuid, e.node2_idx_uuid)) {
          return true; //return true if same graph and the minimum node index of this is less than the minimum node index of e
        }
        else if (std::min(this->node1_idx_uuid, this->node2_idx_uuid) == std::min(e.node1_idx_uuid, e.node2_idx_uuid)) {
            return std::max(this->node1_idx_uuid, this->node2_idx_uuid) < std::max(e.node1_idx_uuid, e.node2_idx_uuid);
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
    Edge(size_type e_n1_idx, size_type e_n2_idx, size_type edg_idx, Graph* gpointer) { //private constructor for Edge
        node1_idx_uuid = e_n1_idx;
        node2_idx_uuid = e_n2_idx;
        ext_idx = edg_idx;
        g = gpointer;
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_i_to_u.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_edges()){
      size_type curr_edge_idx = edge_i_to_u.at(i);
      Edge result_edge(internal_edges.at(curr_edge_idx).node_1_idx, internal_edges.at(curr_edge_idx).node_2_idx, i, const_cast<Graph* >(this));
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
    size_type lower_idx = std::min(node_i_to_u.at(a.ext_idx), node_i_to_u.at(b.ext_idx)); //our edges are added with the lower node index first
    size_type higher_idx = std::max(node_i_to_u.at(a.ext_idx), node_i_to_u.at(b.ext_idx));

    std::map<unsigned, unsigned> map_at_index = adj_map.at(lower_idx);
    if (map_at_index.find(higher_idx) == map_at_index.end()) {
      return false;

    }
    return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE

    if (has_edge(a, b)) { //check to see if current edge already exists
      Edge same_edge(node_i_to_u.at(a.ext_idx), node_i_to_u.at(b.ext_idx), (*(adj_map.at(node_i_to_u.at(a.ext_idx)).find(node_i_to_u.at(b.ext_idx)))).second, this);
      return same_edge;
    }
    else {

      size_type num_edge = num_edges();
      size_type size = internal_edges.size();

      adj_map.at(node_i_to_u.at(a.ext_idx)).insert(std::pair<size_type, size_type>(node_i_to_u.at(b.ext_idx), num_edge));
      adj_map.at(node_i_to_u.at(b.ext_idx)).insert(std::pair<size_type, size_type>(node_i_to_u.at(a.ext_idx), num_edge));
      
      internal_edges.push_back(InternalEdge(num_edge, node_i_to_u.at(a.ext_idx), node_i_to_u.at(b.ext_idx), value, this));
      edge_i_to_u.push_back(size);

      Edge new_edge(node_i_to_u.at(a.ext_idx), node_i_to_u.at(b.ext_idx), num_edge, this);

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
    internal_nodes.clear();
    adj_map.clear();
    internal_edges.clear();
    edge_i_to_u.clear();
    node_i_to_u.clear();
  }

  
  // Node Iterator
  

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
    size_type idx;
    Graph* g;

    NodeIterator(Graph* gpointer, size_type index) {
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
      return NodeIterator(const_cast<Graph*>(this), this->node_i_to_u.size());
    }

  
  // Incident Iterator

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
      return Edge(n.g->node_i_to_u.at(n.ext_idx), map_it->first, n.g->internal_edges.at(map_it->second).ext_idx, n.g);
    }

    /** Increment an IncidentIterator.
     * @return The IncidentIterator that has an index one larger than before.
     *
     * @post For all returned Edge, Edge index < size of items in the IncidentIterator.
    */

    IncidentIterator& operator++() {
      ++map_it;
      return (*this);
    }

    /** Tests whether this IncidentIterator and inc_iter are equal.
     * @return A bool of 1 if they are equal, 0 if not. 
     * Checks to see if the same node and also the same index.
     *
     * @post All bool 0 or 1.
     */

    bool operator==(const IncidentIterator& inc_iter) const {
      return (this->n == inc_iter.n && this->map_it == inc_iter.map_it);
    }

    /** Tests whether this IncidentIterator and inc_iter are not equal.
     * @return A bool of 1 if not equal, 0 if they are. 
     * Checks to see if the same node and also the same index.
     *
     * @post All bool 0 or 1.
     */

    bool operator!=(const IncidentIterator& inc_iter) const {
      return !(this->n == inc_iter.n && this->map_it == inc_iter.map_it);
    }
  
   private:
    friend class Graph;
    Node n;
    std::map<size_type, size_type>::const_iterator map_it;

    // HW1 #3: YOUR CODE HERE
    IncidentIterator(Node node, std::map<size_type, size_type>::const_iterator map_iter) {
      n = node;
      map_it = map_iter;
    }
  };

  
  // Edge Iterator
  

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
      size_type edge_uuid = g->edge_i_to_u.at(idx);
      return(Edge(g->internal_edges.at(edge_uuid).node_1_idx, g->internal_edges.at(edge_uuid).node_2_idx, idx, g));
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
    size_type idx;

    // HW1 #3: YOUR CODE HERE
    EdgeIterator(Graph* gpointer, size_type index) {
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
    return EdgeIterator(const_cast<Graph*>(this), this->edge_i_to_u.size());
  }

  /** Removes the specified edge given two nodes.
     * @return A size_type of 1 if edge was removed, 0 if not removed. 
     * Checks to see if that edge exists in the graph or not.
     *
     * @post A size_type of  0 or 1.
     */
  size_type remove_edge(const Node& n1, const Node& n2) {

    if (has_edge(n1, n2)) {

      size_type curr_edge_idx = adj_map.at(node_i_to_u.at(n1.ext_idx)).at(node_i_to_u.at(n2.ext_idx));

      adj_map.at(node_i_to_u.at(n1.ext_idx)).erase(node_i_to_u.at(n2.ext_idx));
      adj_map.at(node_i_to_u.at(n2.ext_idx)).erase(node_i_to_u.at(n1.ext_idx));
    
      size_type new_ext_idx = internal_edges.at(curr_edge_idx).ext_idx;
      size_type i_to_u_move = edge_i_to_u.size() - 1;

      internal_edges.at(edge_i_to_u.at(i_to_u_move)).ext_idx = new_ext_idx;

      edge_i_to_u.at(new_ext_idx) = edge_i_to_u.at(i_to_u_move);
      edge_i_to_u.pop_back();

      return 1;
    }
    return 0;
  }

  /** Removes the specified edge given an Edge.
     * @return A size_type of 1 if edge was removed, 0 if not removed. 
     * Checks to see if that edge exists in the graph or not.
     *
     * @post A size_type of  0 or 1.
     */
  size_type remove_edge(const Edge& e) {

    size_type node_1_idx = e.node1_idx_uuid;
    size_type node_2_idx = e.node2_idx_uuid;

    if (has_edge(Node(internal_nodes.at(node_1_idx).ext_id, this), Node(internal_nodes.at(node_2_idx).ext_id, this))) {

      adj_map.at(node_1_idx).erase(node_2_idx);
      adj_map.at(node_2_idx).erase(node_1_idx);

      size_type new_ext_idx = internal_edges.at(e.uuid).ext_idx;
      size_type i_to_u_move = edge_i_to_u.size() - 1;

      internal_edges.at(edge_i_to_u.at(i_to_u_move)).ext_idx = new_ext_idx;

      edge_i_to_u.at(e.ext_idx) = edge_i_to_u.at(i_to_u_move);
      edge_i_to_u.pop_back();
      return 1;
    }
    return 0;
  }

   /** Removes the specified edge given an edgeiterator. If no nodes left, return the end Edge Iterator.
     * @return An edgeiterator. 
     * 
     *
     * @post The edge_iterator.
     */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (num_edges() == 0) {
      return edge_end();
    }
    remove_edge(*e_it);
    return e_it;
  }

  /** Removes the specified node given a Node.
     * @return A size_type of 1 if node was removed, 0 if not removed. 
     * Checks to see if that node exists in the graph or not.
     *
     * @post A size_type of 0 or 1.
     */
  size_type remove_node(const Node& n) {
    
    if (has_node(n)) {
      std::map<unsigned, unsigned> curr_map = adj_map.at(node_i_to_u.at(n.ext_idx));

      for (auto it = curr_map.begin(); it != curr_map.end(); ++it) {
        remove_edge(Node(n.ext_idx, this), Node(internal_nodes.at((*it).first).ext_idx, this));
      }

      size_type new_ext_idx = internal_nodes.at(node_i_to_u.at(n.ext_idx)).ext_idx;
      size_type i_to_u_move = node_i_to_u.size() - 1;

      internal_nodes.at(node_i_to_u.at(i_to_u_move)).ext_idx = new_ext_idx;

      node_i_to_u.at(n.ext_idx) = node_i_to_u.at(i_to_u_move);

      node_i_to_u.pop_back();
      return 1;
    }
    
    return 0;
  }

  /** Removes the specified node given a nodeiterator. If no nodes left, return the end Node Iterator.
     * @return A node_iterator. 
     * 
     *
     * @post The node_iterator.
     */
  node_iterator remove_node(node_iterator n_it) {
    if (num_nodes() == 0) {
      return node_end();
    }

    remove_node(*n_it);
    return n_it;
  }

 public:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<InternalNode> internal_nodes; //stores our vector of InternalNode objects
  std::map<unsigned, std::map<unsigned, unsigned>> adj_map; // std::map<unsigned, std::map<unsigned, unsigned>> node_idx_key; //stores a map of map with the node_1 index as the key, node_2 index as key for inner map with edge index as value
  std::vector<InternalEdge> internal_edges; //stores our vector of InternalEdge objects
  std::vector<unsigned> edge_i_to_u;
  std::vector<unsigned> node_i_to_u;

private:

  /** Print helper functions that print out the data structures within the graph class
     */
  void sw_print_nodes() {
    std::cout << "Internal_nodes: " << std::endl;
    for (auto it = internal_nodes.begin(); it != internal_nodes.end(); ++it) {
      std::cout << "Node UUID: " << it->uuid << std::endl;
      std::cout << "Node Ext_ID: " << it->ext_idx << std::endl;
    }

    std::cout << "Node_i_to_u: " << std::endl;
    for (auto it = node_i_to_u.begin(); it != node_i_to_u.end(); ++it) {
      std::cout << (*it) << std::endl;
    }
  }

  void sw_print_edges() {
    std::cout << "Internal_edges: " << std::endl;
    for (auto it = internal_edges.begin(); it != internal_edges.end(); ++it) {
      std::cout << "Edge UUID: " << it->uuid << std::endl;
      std::cout << "Edge Ext_ID: " << it->ext_idx << std::endl;
    }

    std::cout << "Edge_i_to_u: " << std::endl;
    for (auto it = edge_i_to_u.begin(); it != edge_i_to_u.end(); ++it) {
      std::cout << (*it) << std::endl;
    }
  }

  void sw_print_adj_map() {
    std::cout << "Adjacency Map: " << std::endl;
    for (auto it = adj_map.begin(); it != adj_map.end(); ++it) {
      std::cout << "Node 1 UUID: " << (*it).first << std::endl;
      for (auto it_2 = (*it).second.begin(); it_2 != (*it).second.end(); ++it_2) {
        std::cout << "Node 2 UUID: " << (*it_2).first << std::endl;
        std::cout << "Edge UUID: " << (*it_2).second << std::endl;
      }
    }
  }

};

#endif // CME212_GRAPH_HPP
