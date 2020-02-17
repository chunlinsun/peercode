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
template <typename V, typename E>
class Graph {
public:	
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
  using uid_type  = unsigned;
  using node_value_type = V;
  using edge_value_type = E;
  
  /** Tupils for storing index info */
  using IndexTupil = typename std::pair<size_type, size_type>;
  
  /** set up a structure to store multiple values required by one node */
  struct NodeValue {
    Point pos;
    node_value_type val;
    // add a new attribute as the stored node index
    size_type index;
    // a flag indicates if this node is valid
    bool is_valid;
    
    /** Constructors */
    NodeValue(): pos(0), val(), index(0), is_valid(true) {}
    NodeValue(const NodeValue& nv): 
      pos(nv.pos), val(nv.val), index(nv.index), is_valid(nv.is_valid) {}
    NodeValue( size_type i, const Point& p ): pos(p), index(i), is_valid(true) {}
    NodeValue( size_type i, const Point& p, const node_value_type& v ):
      pos(p), val(v), index(i), is_valid(true) {}
  };
  
  /** set up a structure to store multiple values required by one edge */
  struct EdgeValue {
    size_type node1;
    size_type node2;
    edge_value_type val;
    // add a new attribute as the stored edge index
    size_type index;
    // a flag indicates if this node is valid
    bool is_valid;
    
    /** Constructors */
    EdgeValue(): node1(0), node2(0), val(), index(0), is_valid(true) {}
    EdgeValue(const EdgeValue& ev): node1(ev.node1), node2(ev.node2),
      val(ev.val), index(ev.index), is_valid(ev.is_valid) {}
    EdgeValue( size_type i, size_type n1, size_type n2 ): 
      node1(n1), node2(n2), index(i), is_valid(true) {}
    EdgeValue( size_type i, size_type n1, size_type n2, const edge_value_type& v ): 
      node1(n1), node2(n2), val(v), index(i), is_valid(true) {}
  };

 private:
  std::vector<NodeValue> node_list;
  std::vector< std::vector<IndexTupil> > adjacency_list;
  std::vector<EdgeValue> edge_list;
  // add uid lists for removal operations
  std::vector<uid_type> node_uids;
  std::vector<uid_type> edge_uids;
  

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
  Graph() {}

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
	 
    Node(): graph_ptr(nullptr), uid(0) {
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW2 new features
      // get position of a node via uid
      return graph_ptr-> node_list[uid]. pos;
    }
    /** return a modifiable node position */
    Point& position() {
      return graph_ptr-> node_list[uid]. pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW2 new features
      // get index of a node via uid
      return graph_ptr-> node_list[uid]. index;
    }
    
    /** Query node values */
    node_value_type& value() {
      // HW2 new features
      // get value of a node via uid
      return graph_ptr-> node_list[uid]. val;
    }
    const node_value_type& value() const {
      return graph_ptr-> node_list[uid]. val;
    }
    
    /** Query degree of a node (num of edges connected to it) */
    size_type degree() const {
      return graph_ptr-> adjacency_list[uid]. size();
    }
    
    /** Get the iterator to the first edge(adjacency) stored in the node's adjacency list
      * Return value is a incident iterator */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_ptr, uid, 0);
    }
    
    /** Get the iterator pointing next to the last edge(adjacency)
      * Return value is a incident iterator */
    incident_iterator edge_end() const {
      return incident_iterator(graph_ptr, uid, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ptr == n.graph_ptr) && (uid == n.uid);
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
      return ( std::less<Graph*> {}(graph_ptr, n.graph_ptr) ||
          ( (graph_ptr == n.graph_ptr) && ( index() < n.index() ) ) );
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_ptr;
    size_type uid;
    
    /** a valid constructor */
    Node(const Graph* g, size_type i):
      graph_ptr( const_cast<Graph*>(g) ), uid(i) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_uids.size();
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
  
  Node add_node(const Point& position) {
    // HW2 new features
    // push back a new maping in the node uid list
    // the newly added uid will always point to the back of node_list
    size_type n_idx = node_uids.size(), n_uid = node_list.size();
    node_uids.push_back( n_uid );
    // modified for the new NodeValue structure
    node_list.push_back( NodeValue(n_idx, position) );
    adjacency_list.push_back(std::vector<IndexTupil>());
    return Node(this, n_uid);
  }
  
  Node add_node(const Point& position, const node_value_type& v) {
    Node curr_node = add_node(position);
    curr_node. value() = v;
    return curr_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return node_list[n.uid]. is_valid;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW2 new features
    // get the node index via uid
    if (i < size())
      return Node(this, node_uids[i]);
    else
      return Node();
  }
  
  /** implement a set_node_value method to update node values, but not node itself */
  void set_node_value(const size_type i, const node_value_type& v) {
    node_list[ node_uids[i] ].val = v;
  }
  void set_node_value(const node_type& n, const node_value_type& v) {
    node_list[n.uid].val = v;
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_ptr, node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_ptr, node2_uid);
    }
	
    size_type index() const {
      return graph_ptr-> edge_list[uid]. index;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ( (node1() == e.node1() && node2() == e.node2()) ||
        (node1() == e.node2() && node2() == e.node1())	);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ( std::less<Graph*> {}(graph_ptr, e.graph_ptr) ||
          ( (graph_ptr == e.graph_ptr) && ( index() < e.index() ) ) );
    }
    
    // HW2 new features
    /** Query edge values */
    edge_value_type& value() {
      return graph_ptr-> edge_list[uid]. val;
    }
    const edge_value_type& value() const {
      return graph_ptr-> edge_list[uid]. val;
    }
    
    // HW2 new features
    /** Query the length of an edge according to its end node positions */
    double length() const {
      return norm( node1().position() - node2().position() );
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_ptr;
    size_type node1_uid;
    size_type node2_uid;
    size_type uid;
    
    /** a valid constructor */
    Edge(const Graph* g, size_type n1, size_type n2, size_type i):
      graph_ptr( const_cast<Graph*>(g) ), node1_uid(n1), node2_uid(n2), uid(i) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_uids.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW2 new features
    // Get edge info via edge_uid, then get node info via node_uid
    assert(i < num_edges());
    auto edge_info = edge_list[ edge_uids[i] ];
    return Edge(this, edge_info.node1, edge_info.node2, edge_uids[i]);
  }

  
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return ( get_adjacency_for(a,b) > -1 );
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
    // HW2 new features:
    // The index tupil stored in adjacency list becomes uid!
    int adj_loc = get_adjacency_for(a,b);
    size_type e_idx = edge_uids.size(), e_uid = edge_list.size();
    // the edge to be added is already found
    if ( adj_loc > -1 ) {
      e_uid = adjacency_list[a.uid][adj_loc]. second;
      return Edge(this, a.uid, b.uid, e_uid);
    }
      
    // add a new mapping to the edge_list
    edge_uids. push_back( e_uid );
    // add edge info to the adjacency list
    adjacency_list[a.uid]. push_back( IndexTupil(b.uid, e_uid) );
    adjacency_list[b.uid]. push_back( IndexTupil(a.uid, e_uid) );
    // add edge info to the edge list
    edge_list. push_back( EdgeValue(e_idx, a.uid, b.uid) );
    return Edge(this, a.uid, b.uid, e_uid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    for (size_type i=0; i < adjacency_list.size(); ++i)
      adjacency_list[i].clear();
    adjacency_list.clear();
    node_list.clear();
    edge_list.clear();
    node_uids.clear();
    edge_uids.clear();
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

    /** Dereferencing operator
     * @return an Node object with the specific graph and index */
    Node operator*() const {
       return graph_ptr->node(idx);
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


  /** Get the first-indexed node in the graph
   * @return a NodeIterator with n.graph_ptr == this and n.idx == 0 */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }
  
  /** Get the ietrator pointing to the end of node_list
   * @return a NodeIterator with n.graph_ptr == this and n.idx == node_num */
  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
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
    IncidentIterator(): graph_ptr(nullptr), node_uid(0), curr_adjacency(0) {}

    /** Dereferencing operator
     * @return an Edge object e with e.node1() == this node */
    Edge operator*() const {
      // HW2 new features
      // access incident info via node_uid and edge_uid
      IndexTupil adjacency_info = graph_ptr-> adjacency_list[node_uid] [curr_adjacency];
      return Edge(graph_ptr, node_uid, adjacency_info.first, adjacency_info.second );
    }

    /** self-incrementing operator */
    IncidentIterator& operator++() {
      curr_adjacency++;
      return *this;
    }

    /** equivalence operator for the same type */
    bool operator==(const IncidentIterator& it) const {
      return ( graph_ptr == it.graph_ptr && node_uid == it.node_uid &&
                 curr_adjacency == it.curr_adjacency );
    }

   private:
    friend class Graph;
    friend class Node;
    // HW1 #3: YOUR CODE HERE
    const graph_type* graph_ptr;
    const size_type node_uid;
    size_type curr_adjacency;
    
    /** Valid constructor */
    IncidentIterator( const graph_type* g, const size_type ni, size_type ai ):
      graph_ptr(g), node_uid(ni), curr_adjacency(ai) {}
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
      return graph_ptr->edge(idx);
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
    return edge_iterator(this, num_edges());
  }
  
  
  
  /** remove a node and all its incidents
   * @pre @a n a valid node of this graph
   * @return an index for the target Node object to be removed
   * @post new num_nodes() == old num_nodes() - 1.
   * @post old node_uids[ n.index() ] is invalid
   *
   * Complexity: O(d^2), where d is the maximum degree for all nodes
   */
  size_type remove_node(const Node& n) {
    // if the node to remove is invalid
    if( !has_node(n) )
      return 0;
    // remove all the incidents of the given node
    auto it = n.edge_begin();
    while (it != n.edge_end())
      remove_edge(*it);
    
    // get the target node index
    size_type node_idx = node_list[n.uid].index;
    
    // erase the targeted location in edge_uids
    quick_erase( node_uids, node_idx );
    
    // update index and validity for the item in the same location now
    node_list[ node_uids[node_idx] ]. index = node_idx;
    node_list[n.uid]. is_valid = false;
    
    return 1;
  }
  
  
  /** remove a node and all its incidents
   * @pre @a edge_iterator pointing to a valid edge of this graph
   * @return an index edge_iterator pointing to the same location
   * @post new num_nodes() == old num_nodes() - 1.
   * @post old node_uids[ n.index() ] is inaccessible
   *
   * Complexity: O(d^2), where d is the maximum degree for all nodes
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  /** remove an edge
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an index for the target Edge object to be removed
   * @post has_edge(@a a, @a b) == false
   * @post new num_edges() == old num_edges() - 1.
   *
   * Complexity: O(d), where d is the maximum degree for all nodes
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // find the edge_eid, must be valid
    int adj_loc = get_adjacency_for(a,b);
    // return zero if failure
    if (adj_loc < 0)
      return 0;
      
    size_type e_uid = adjacency_list[a.uid] [adj_loc].second;
    // clear the adjacency info associated with this edge
    quick_erase( adjacency_list[ a.uid ], adj_loc);
    adj_loc = get_adjacency_for(b,a);
    quick_erase( adjacency_list[ b.uid ], adj_loc);
    
    // get the target edge index
    size_type edge_idx = edge_list[e_uid].index;
    
    // erase the targeted location in edge_uids
    quick_erase( edge_uids, edge_idx );
    
    // update index and validity
    edge_list[ edge_uids[edge_idx] ]. index = edge_idx;
    edge_list[e_uid]. is_valid = false;
    
    return 1;
  }


  /** remove an edge
   * @pre @a e a valid edge of this graph
   * @return an index for the target Edge object to be removed
   * @post has_edge(@a e.node1(), @a e.node2()) == false
   * new num_edges() == old num_edges() - 1.
   *
   * Complexity: O(d), where d is the maximum degree for all nodes
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge( e.node1(), e.node2() );
  }

  /** remove an edge
   * @pre @a edge_iterator pointing to a valid edge of this graph
   * @return an index edge_iterator pointing to the same location
   * new num_edges() == old num_edges() - 1.
   *
   * Complexity: O(d), where d is the maximum degree for all nodes
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }


 private:
  /** A helper function to quickly erase the item at a specific location in a vector container
   * @param[in] vector<U>: the target container
   * @param[in] size_type: the location to erase
   *
   * Complexity: O(1)
   */
  template<typename U>
  void quick_erase( std::vector<U>& vec, const size_type loc ) {
    vec[loc] = vec[ vec.size() - 1 ];
    vec.pop_back();
  }
  

  /** Find whether edge (a,b) exists in the adjacency list.
   * @pre @a a and @a b are valid nodes of this graph
   * @return a location in @a a's adjacency list where (a,b) lies
   * @return -1 if such an edge do not exists
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  int get_adjacency_for(const Node& a, const Node& b) const {
    assert( this->has_node(a) && this->has_node(b) );
    const std::vector<IndexTupil>& a_neighbor = adjacency_list[ a.uid ];
    for (size_type i=0; i < a_neighbor.size(); ++i)
      // search in a's adjacency list, compare the neighbor node's uid
      if (a_neighbor[i].first == b.uid)
        // edge (a,b) is found, return the location in the adjacency list
        return i;
        
    // edge (a,b) is not found
    return -1;
  }

};

#endif // CME212_GRAPH_HPP
