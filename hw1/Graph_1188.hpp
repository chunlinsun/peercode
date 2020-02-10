#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <utility>
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
public:	
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
  using node_value_type = V;
  
  /** Tupils for storing index info */
  using IndexTupil = typename std::pair<size_type, size_type>;
  
  /** set up a structure to store multiple values required by one node */
  struct NodeValue {
    Point pos;
    node_value_type val;
    NodeValue() {}
    NodeValue(const Point& p): pos(p) {}
    NodeValue( const Point& p, const node_value_type& v ):
      pos(p), val(v) {}
  };

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  std::vector<NodeValue> node_list;
  std::vector< std::vector<IndexTupil> > adjacency_list;
  std::vector<IndexTupil> edge_list;
  // usage: adjacency_list[ start_node_idx ] = (end_node_idx, edge_idx)
  size_type node_num, edge_num;
  

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph():node_num(0), edge_num(0) {
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
  class Node: private totally_ordered<Node> {
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
	 
    Node(): graph_ptr(nullptr), idx(0) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_ptr->node_list[idx].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Query node value */
    node_value_type& value() {
      return graph_ptr->node_list[idx].val;
    }
    const node_value_type& value() const {
      return graph_ptr->node_list[idx].val;
    }
    
    /** Query degree of a node (num of edges connected to it) */
    size_type degree() const {
      return graph_ptr->adjacency_list[idx].size();
    }
    
    /** Get the iterator to the first edge(adjacency) stored in the node's adjacency list
      * Return value is a incident iterator */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_ptr, idx, 0);
    }
    
    /** Get the iterator pointing next to the last edge(adjacency)
      * Return value is a incident iterator */
    incident_iterator edge_end() const {
      return incident_iterator(graph_ptr, idx, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ptr == n.graph_ptr) && (idx == n.idx);
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
      assert( graph_ptr == n.graph_ptr );
      return (idx < n.idx);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    const Graph* graph_ptr;
    size_type idx;
    
    /** a valid constructor */
    Node(const Graph* g, size_type index):
      graph_ptr(g), idx(index) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_num;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node value to specify at the new node (optional)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  //--design_1
  //--This doesn't follow the specification. There should be ONE add_node method
  //--with the signature
  //--  Node add_node(const Point& position, node_value_type& val = node_value_type())
  //--which simply provides a default argument for val if none is given.
  //--START
  Node add_node(const Point& position) {
    node_list.push_back( NodeValue(position) );
    adjacency_list.push_back(std::vector<IndexTupil>());
    return Node(this, node_num++);
  }
  
  Node add_node(const Point& position, node_value_type& val) {
    node_list.push_back( NodeValue(position, val) );
    adjacency_list.push_back(std::vector<IndexTupil>());
    return Node(this, node_num++);
  }
  //--END

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return !(n.graph_ptr != this || n.idx >= node_num);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < node_num);
    return Node(this, i);
  }
  
  //--design_0
  //--This is fine, but i'm not sure why you need it. Instead of calling
  //--  set_node_value(n, v);
  //--you could simply do
  //--  n.value() = v
  //--instead. It should have an identical effect.
  //--START
  /** implement a set_node_value method to update node values, but not node itself */
  void set_node_value(const size_type i, const node_value_type& v) {
    node_list[i].val = v;
  }
  void set_node_value(const node_type& n, const node_value_type& v) {
    set_node_value(n.idx, v);
  }
  //--END
  
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_ptr, node1_idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_ptr, node2_idx);
    }
	
    size_type index() const {
      return idx;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ( (node1_idx == e.node1_idx && node2_idx == e.node2_idx) ||
        (node1_idx == e.node2_idx && node2_idx == e.node1_idx)	);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert( graph_ptr == e.graph_ptr );
      return (idx < e.idx);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const Graph* graph_ptr;
    size_type node1_idx;
    size_type node2_idx;
    size_type idx;
    
    /** a valid constructor */
    Edge(const Graph* g, size_type node1, size_type node2, size_type index):
      graph_ptr(g), node1_idx(node1), node2_idx(node2), idx(index){ }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_num;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < edge_num);
    return Edge(this, edge_list[i].first, edge_list[i].second, i);
  }

  // modification from HW0 regrading feedback
  // added @ 2/1/2020
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    return ( find_edge(a,b) > -1 );
  }

  // modification from HW0 regrading feedback
  // added @ 2/1/2020
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
    const size_type a_idx = a.index(), b_idx = b.index();
    int edge_id = find_edge(a,b);
    // the edge to be added is already found
    if ( edge_id > -1 )
      return Edge(this, a_idx, b_idx, edge_id);
    // add edge info to the adjacency list
    // each adjacency info is stored twice for the convenience of 
    // IncidenceIterrator class and degree() method in Node
    adjacency_list[a_idx]. push_back( IndexTupil(b_idx, edge_num) );
    adjacency_list[b_idx]. push_back( IndexTupil(a_idx, edge_num) );
    // add edge info to the edge index tupil list
    edge_list. push_back( IndexTupil(a_idx, b_idx) );
    return Edge(this, a_idx, b_idx, edge_num++);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    edge_list.clear();
    for (size_type i=0; i < adjacency_list.size(); ++i)
      adjacency_list[i].clear();
    adjacency_list.clear();
    node_list.clear();
    
    // reset node_num and edge_num
    node_num = 0;
    edge_num = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator(): graph_ptr(nullptr), idx(0) {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereferencing operator
     * @return an Node object with the specific graph and index */
    Node operator*() const {
       return Node(graph_ptr, idx);
    }
    
    /** self-incrementing operator */
    NodeIterator& operator++() {
      idx++;
      return *this;
    }
    
    /** equivalence operator for the same type */
    bool operator==(const NodeIterator& it) const {
      return ( graph_ptr == it.graph_ptr && idx == it.idx );
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const graph_type* graph_ptr;
    size_type idx;
    
    /** Valid constructor */
    NodeIterator( const graph_type* g, size_type i ):
      graph_ptr(g), idx(i) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Get the first-indexed node in the graph
   * @return a NodeIterator with n.graph_ptr == this and n.idx == 0 */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }
  
  /** Get the ietrator pointing to the end of node_list
   * @return a NodeIterator with n.graph_ptr == this and n.idx == node_num */
  node_iterator node_end() const {
    return node_iterator(this, node_num);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(): graph_ptr(nullptr), node_id(0), curr_adjacency(0) {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereferencing operator
     * @return an Edge object e with e.node1() == this node */
    Edge operator*() const {
      IndexTupil adjacency_info = graph_ptr->adjacency_list[node_id][curr_adjacency];
      size_type neighbor = adjacency_info.first,
        edge_id = adjacency_info.second;
      return Edge(graph_ptr, node_id, neighbor, edge_id );
    }

    /** self-incrementing operator */
    IncidentIterator& operator++() {
      curr_adjacency++;
      return *this;
    }

    /** equivalence operator for the same type */
    bool operator==(const IncidentIterator& it) const {
      return ( graph_ptr == it.graph_ptr && node_id == it.node_id &&
                 curr_adjacency == it.curr_adjacency );
    }

   private:
    friend class Graph;
    friend class Node;
    // HW1 #3: YOUR CODE HERE
    const graph_type* graph_ptr;
    const size_type node_id;
    size_type curr_adjacency;
    
    /** Valid constructor */
    IncidentIterator( const graph_type* g, const size_type ni, size_type ai ):
      graph_ptr(g), node_id(ni), curr_adjacency(ai) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(): graph_ptr(nullptr), idx(0) {
    }


    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereferencing operator
     * @return an Edge object with the specific graph and index */
    Edge operator*() const {
      IndexTupil edge_info = graph_ptr->edge_list[idx];
      size_type node1_id = edge_info.first, node2_id = edge_info.second;
      return Edge( graph_ptr, node1_id, node2_id, idx );
    }

    /** self-incrementing operator */
    EdgeIterator& operator++() {
      idx++;
      return *this;
    }

    /** equivalence operator for the same type */
    bool operator==(const EdgeIterator& it) const {
      return ( graph_ptr == it.graph_ptr && idx == it.idx );
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const graph_type* graph_ptr;
    size_type idx;
    
    /** Valid constructor */
    EdgeIterator( const graph_type* g, size_type i ):
      graph_ptr(g), idx(i) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Get the first-indexed edge in the graph
   * @return a EdgeIterator with e.graph_ptr == this and e.idx == 0 */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }
  
  /** Get the ietrator pointing to the end of edge_list
   * @return a EdgeIterator with e.graph_ptr == this and e.idx == node_num */
  edge_iterator edge_end() const {
    return edge_iterator(this, edge_num);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
    
  // modification from HW0 regrading feedback
  // added @ 2/1/2020
  /** Find whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return the edge id if there exists a Edge e s.t. e.node1() == @a a and e.node2() == @a b.
   * @return -1 if such an edge do not exists
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  int find_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert( this->has_node(a) && this->has_node(b) );
    const size_type a_idx = a.index(), b_idx = b.index();
    const std::vector<IndexTupil>& a_neighbor = adjacency_list[a_idx];
    for (size_type i=0; i < a_neighbor.size(); ++i)
      if (a_neighbor[i].first == b_idx)
        // edge (a,b) is found
        return a_neighbor[i].second;
        
    // edge (a,b) is not found
    return -1;
  }

};

#endif // CME212_GRAPH_HPP
