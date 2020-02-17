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
template <typename V, typename E>
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

  /** Type of node value. */
  using node_value_type = V;
  /** Type of edge value. */
  typedef E edge_value_type;

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

  // define a pair to store edges
  typedef std::pair<size_type, size_type> npair;

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

    Node() {
      // HW0: YOUR CODE HERE
        graph_ptr = nullptr;
        node_id = 0;
    }
      
    // constructor for valid nodes
    Node(const graph_type* g_ptr, size_type id): graph_ptr(const_cast<graph_type*>(g_ptr)), node_id(id) {};
   
      /** Return this node's position. */
    Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_ptr->nodes[node_id].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_id;
    }

    
      
    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return this node's value */
    node_value_type& value() {
        return graph_ptr->nodes[node_id].second;
    };

    const node_value_type& value() const {
        return graph_ptr->nodes[node_id].second;
    };

    /** Return the number of incident edges of this node, 
      * a number in the range [0, num_edges).
      */
    size_type degree() const {
        return graph_ptr->adj_list[node_id].size();
    }

    /** Return the start of the incident iterator. */
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_ptr, 0, node_id);
    }

    /** Return the end of the incident iterator. 
      * @post incident_edge_id == this->degree(), 
      *       meaning that this iterator does not point to a real edge 
      */
    incident_iterator edge_end() const {
        return IncidentIterator(graph_ptr, this->degree(), node_id);
    };


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
        if (graph_ptr == n.graph_ptr  && node_id == n.index()) {
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
      // HW0: YOUR CODE HERE
        // std function compare two graphs?
        if (node_id < n.index() && graph_ptr == n.graph_ptr ) {
            return true;
        } else if (graph_ptr < n.graph_ptr) {
          return true;
        } 
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_ptr;
    size_type node_id;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return n_nodes;
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
  Node add_node(const Point& position, const node_value_type & node_val= node_value_type ()) {
    // HW0: YOUR CODE HERE
      nodes.push_back(std::make_pair(position, node_val));
      // allocate a new (empty) member in the adjacency list for this node
      adj_list.push_back(std::vector<npair> ());
    return Node(this, n_nodes++);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_ptr == this && n.index() < n_nodes);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i);
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_ptr->node(node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_ptr->node(node2_id);      // Invalid Node
    }
    
    /** Get the private members of the edge */
      
      size_type get_node1_id() const{
          return node1_id;
      }
      
      size_type get_node2_id() const{
          return node2_id;
      }
      
      size_type get_edge_id() const{
          return edge_id;
      }
      
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ptr==e.graph_ptr) {
          if((node1_id==e.get_node1_id()&& node2_id==e.get_node2_id()) ||
             (node2_id==e.get_node1_id() && node1_id==e.get_node2_id()) ) {
              return true;
          }
      }
      return false;
    }
        
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if (graph_ptr==e.graph_ptr ) {
            if (edge_id < e.get_edge_id() ) return true;
        }
        return false;
    }
    
    /* Return the length of the edge, i.e. the distance between the two nodes
     * of the edge.
     */
    double length() const {
        return norm(node1().position() - node2().position());
    }

    /** Return this edge's value */
    edge_value_type & value () {
    	return graph_ptr->edge_value[edge_id];
    }

    const edge_value_type & value () const {
    	return graph_ptr->edge_value[edge_id];
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_ptr;
    size_type node1_id, node2_id;
    size_type edge_id;
      
    // constructor for a valid edge
    Edge(const graph_type* g_ptr, size_type node1, size_type node2, size_type edge): 
        graph_ptr(const_cast<graph_type*>(g_ptr)), node1_id(node1), node2_id(node2),
        edge_id(edge) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, edges[i].first, edges[i].second, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
      // loop over the edge vector:
      for (unsigned i=0; i < n_edges; i++) {
          if ((edges[i].first == a.index() && edges[i].second == b.index() ) ||
              (edges[i].first == b.index() && edges[i].second == a.index() ) )
              return true;
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
  Edge add_edge(const Node& a, const Node& b, 
  	            const edge_value_type& edge_val = edge_value_type()) {
      // check if has the edge
      for (unsigned i=0; i < n_edges; i++) {
          if ((edges[i].first == a.index() && edges[i].second == b.index() ) ||
              (edges[i].first == b.index() && edges[i].second == a.index() ) ) {
              return Edge(this, a.index(), b.index(), i);
          }
      }
      // create new edge
      edge_type new_edge = Edge(this, a.index(), b.index(), n_edges++);
      edges.push_back(std::make_pair(a.index(), b.index()));
      edge_value.push_back(edge_val);

      // add new member for node a and b in the adjacency list
      adj_list[a.index()].push_back(std::make_pair(b.index(), n_edges-1));
      adj_list[b.index()].push_back(std::make_pair(a.index(), n_edges-1));
  
      return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
      nodes.clear();   // all private member
      edges.clear();
      adj_list.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator>  {
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

    NodeIterator(const graph_type* g_ptr, size_type id): 
        graph_ptr(const_cast<graph_type*>(g_ptr)), node_id(id) {};

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference the current node iterator
     * @return The node that the node iterator points to
     */
    Node operator*() const {
        return Node(graph_ptr, node_id);
    }

    /** Increment the current node iterator
     * @pre 0 <= node_id < n_nodes
     * @post The node iterator points to the next node in the graph
     */
    NodeIterator& operator++() {
        node_id ++; 
        return *this;
    }

    /** Determine whether two node iterators are equal
     * @param[in] node_iter The node iterator that the current node iterator compares with
     * @return true if both the graph that the nodes belong to and their id's are the same
     */
    bool operator==(const NodeIterator& node_iter) const {
      return (graph_ptr == node_iter.graph_ptr && node_id == node_iter.node_id);
    }

   private:
    // HW1 #2: YOUR CODE HERE
    friend class Graph;
    graph_type* graph_ptr;
    size_type node_id;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Start of the node iterator 
   * @return An node_iterator that points to the first node of all the nodes
   */
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }

  /** End of the node iterator 
   * @return An node_iterator that points to the last node of all the nodes
   */
  node_iterator node_end() const {
      return NodeIterator(this, n_nodes);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator>{
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
   
    /** Construct a valid IncidentIterator. 
      * @param[in] g_ptr The graph pointer that this node belongs to.
      * @param[in] edge_id_ The index of the current incident edge of this node
      * @param[in] node_id_ The index of this node in the graph   
      * @pre 0<= edge_id < Node(g_ptr, node_id).degree()
      */
    IncidentIterator(const graph_type* g_ptr, size_type edge_id_, 
        size_type node_id_):
        graph_ptr(const_cast<graph_type*>(g_ptr)), incident_edge_id(edge_id_), node_id(node_id_) {
    }

    /** Dereference the edge that the current incident iterator points to. */
    Edge operator*() const {
        return Edge(graph_ptr, node_id, graph_ptr->adj_list[node_id][incident_edge_id].first,
                    graph_ptr->adj_list[node_id][incident_edge_id].second);
    }
    
    /** Increment the incident iterator. */
    IncidentIterator& operator++() {
        incident_edge_id++;
        return *this;
    }

    /** Determine if the current incident iterator is equal to another.
      * @param[in] iter The incident iterator that the current iterator compares which
      * @return true if the graph pointers, the current incident edge ids, and the 
      *        current node id's of the two iterators are all equal
      */
    bool operator==(const IncidentIterator& iter) const {
        return (graph_ptr == iter.graph_ptr && incident_edge_id == iter.incident_edge_id
                && node_id == iter.node_id);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_ptr;
    size_type incident_edge_id;
    size_type node_id; 
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
    EdgeIterator() {
    }

    /** Construct an valid EdgeIterator. 
      * @param[in] g_ptr The graph pointer that the edges belong to.
      *            edge_id The edge index of this graph.
      */
    EdgeIterator(const graph_type* g_ptr, size_type edge_id_):
                 graph_ptr(const_cast<graph_type*>(g_ptr)), edge_id(edge_id_) {};

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference an EdgeIterator. 
      * @pre 0<= edge_id <= num_edge.
      *      If edge_id < num_edge, EdgeIterator points to a valid edge of the graph;
             if edge_id = num_edge, EdgeIterator points to edge_end()
      */          
    Edge operator*() const {
        return Edge(graph_ptr, graph_ptr->edges[edge_id].first, 
                    graph_ptr->edges[edge_id].second, edge_id);
    }

    /** Increment an EdgeIterator. 
      * @pre 0<= edge_id <= num_edge.
      *      If edge_id < num_edge, EdgeIterator points to a valid edge of the graph;
             if edge_id = num_edge, EdgeIterator points to something invalid, i.e it is edge_end()
      * @post 1<= edge_id <= num_edge.
      */  
    EdgeIterator& operator++() {
        if (edge_id < graph_ptr->num_edges()) { // check if reached the end of the edge list
            edge_id ++;
        } 
        return *this;
    }

    /**  Determine if two edge iterators are equal. 
     * @param[in] an edge iterator that the current edge iterator compares with.
     * @return true if both the graph pointers and the edge id's of the two 
     *         iterators are equal.
     */ 
    bool operator==(const EdgeIterator& edge_iter) const {
        return (graph_ptr == edge_iter.graph_ptr && edge_id == edge_iter.edge_id);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_ptr;
    size_type edge_id;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Start of the edge iterator 
   * @return An EdgeIterator that points to the first edge in the edge list
   */
  edge_iterator edge_begin() const {
      return EdgeIterator(this, 0);
  }

  /** End of the edge iterator 
   * @post edge_id == n_edges
   * @return An EdgeIterator that points to something invalid
   */
  edge_iterator edge_end() const {
      return EdgeIterator(this, n_edges);
  }

 /** private members of the Graph class! */
 private:    

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    std::vector<std::pair<Point, node_value_type> > nodes; // vector of struct of point, value, bool(if valid)
    
    // vector of edges, which are pairs of two nodes that are connected
    std::vector<npair> edges;
    
    // adjacency list: a vector (each node) of vectors (its neighbors)
            // outside:id of the node;
            // inside: <node_id, edge_id> of a neighbor
    std::vector<std::vector<npair>> adj_list;

    // vector of the values of edges, accessed by edge id
    std::vector<edge_value_type> edge_value;
    
    // number of nodes and edges
    size_type n_nodes;
    size_type n_edges;
};

#endif // CME212_GRAPH_HPP
