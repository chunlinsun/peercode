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
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;
  struct internal_edge;

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
  using node_value_type = V;

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
  std::vector<internal_node> nodes {}; //default constructor for the std::vector type
  std::vector<internal_edge> edges {};
  std::vector<std::vector<size_type>> adjacency_matrix {};
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
      graph_pointer = nullptr; //does not point to any graph
      n_id = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(graph_pointer != nullptr);
      return graph_pointer->nodes[n_id].position_; //position is stored in internal node in graph
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return n_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value();
    // size_type degree();
    // incident_iterator edge_begin();
    // incident_iterator edge_end();

    /**
    * @brief supplies the value associated with this node
    * 
    * @return the value by reference associated with this node
    * @pre this node belongs to a graph and thus has a value
    */
    node_value_type& value() {
      assert(graph_pointer != nullptr);
      return graph_pointer->nodes[n_id].value_;
    }

    /**
    * @brief supplies the value associated with this node
    * 
    * @return the value associated with this node as a constant reference
    * @pre this node belongs to a graph and thus has a value
    */
    const node_value_type& value() const {
      assert(graph_pointer != nullptr);
      return (const node_value_type) graph_pointer->nodes[n_id].value_;
    }
    
    /**
    * @brief supplies the degree of this node (number of edges connected to it)
    * 
    * @return the degree of this node
    * @pre this node belongs to a graph and thus has a degree
    */
    size_type degree() const {
    assert(graph_pointer != nullptr);
    return graph_pointer->adjacency_matrix[n_id].size();
    }

    /**
    * @brief returns an incident_iterator pointing to the first edge coming out of this node_type
    *
    * @return incident_iterator pointing to fist edge of this node_type
    * @post if this node_type does not have any edges, result may not be dereferenced
    */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_pointer, n_id, graph_pointer->adjacency_matrix[n_id].begin());
    }


    /**
    * @brief returns an incident_iterator pointing to one past the last edge coming out of this node_type
    *
    * @return incident_iterator pointing to one past the last edge of this node_type
    * @post result may not be dereferenced
    */
    incident_iterator edge_end() const {
      return incident_iterator(graph_pointer, n_id, graph_pointer->adjacency_matrix[n_id].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_pointer == n.graph_pointer && n_id == n.n_id);
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
      assert(graph_pointer == n.graph_pointer); //make sure nodes are in same graph before comparison
      return (this->index() < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_pointer;
    size_type n_id;
        
    Node(const graph_type* graph, size_type index){ //constructor used by the graph 
      graph_pointer = const_cast<graph_type*>(graph);
      n_id = index;
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes.size();
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
    adjacency_matrix.push_back(std::vector<size_type>());
    internal_node new_node {position, static_cast<size_type>(nodes.size()), value};
    nodes.push_back(new_node);   
    return node(new_node.n_id_);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.index()+1 <= nodes.size() && n.graph_pointer == this); 
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      graph_pointer = nullptr; //does not belong to any graph 
      node1_id = 0;
      node2_id = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_pointer->node(node1_id);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_pointer->node(node2_id);       
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return ((graph_pointer == e.graph_pointer && node1_id == e.node1_id && node2_id == e.node2_id)
      || (graph_pointer == e.graph_pointer && node1_id == e.node2_id && node2_id == e.node1_id));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      assert(graph_pointer == e.graph_pointer);
      size_type current_id;
      size_type e_id;
      //--functionality_0
      //--I also made changes here
      //--START
      for (size_type i = 0; i < graph_pointer->num_edges(); ++i) {
        if ((graph_pointer->edges[i].node1_id_ == node1_id && graph_pointer->edges[i].node2_id_ == node2_id)
          || (graph_pointer->edges[i].node1_id_ == node2_id && graph_pointer->edges[i].node2_id_ == node1_id)) {
          current_id = i;
        }
        if ((graph_pointer->edges[i].node1_id_ == node2_id && graph_pointer->edges[i].node2_id_ ==node1_id  )
          || (graph_pointer->edges[i].node1_id_ ==node1_id  && graph_pointer->edges[i].node2_id_ ==node2_id  )) {
          e_id = i;
        }
      //--END
      }
      return current_id < e_id;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type node1_id;
    size_type node2_id;
    graph_type* graph_pointer;

    Edge(const graph_type* graph, size_type node1, size_type node2){ //constructor used by the graph
      graph_pointer = const_cast<graph_type*>(graph);
      node1_id = node1;
      node2_id = node2;
      }
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
  edge_type edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, edges[i].node1_id_, edges[i].node2_id_);
  }

  /**
  * @brief helper function I added, to help with incident_iterator 
  * param[in] node1_id, node2_id the id numbers of the nodes connected to the edge_type to be returned
  */
  edge_type edge(size_type node1_id, size_type node2_id) const {
    // HW0: YOUR CODE HERE
    return Edge(this, node1_id, node2_id);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_pointer == b.graph_pointer); //verify pre-conditions
    assert(a.graph_pointer == this);
    return (std::find(adjacency_matrix[a.index()].begin(), adjacency_matrix[a.index()].end(), b.index()) != adjacency_matrix[a.index()].end());
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
  edge_type add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    assert(a.graph_pointer == b.graph_pointer); //verify pre-conditions
    assert(a.graph_pointer == this);
    edge_type return_edge;
    if (has_edge(a,b)){
      return_edge = Edge(this, a.index(), b.index());
    } else {
      adjacency_matrix[a.index()].push_back(b.index());
      adjacency_matrix[b.index()].push_back(a.index());
      internal_edge new_edge {a.index(), b.index(), static_cast<size_type>(edges.size())};
      edges.push_back(new_edge);
      return_edge =  Edge(this, a.index(), b.index());
    }
    return return_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear(); //stl clear method for std::vector class
    edges.clear();
    adjacency_matrix.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator :  private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graph_pointer = nullptr;
      internal_node_iterator = typename std::vector<internal_node>::const_iterator();
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
    * @brief dereferences the NodeIterator to obtain the corresponding node_type
    * 
    * @return node_type behind the NodeIterator
    * @pre this node_iterator points to an actual node_type in a non-empty graph_type
    */
    value_type operator*() const {
      assert(graph_pointer != nullptr && internal_node_iterator!= graph_pointer->nodes.end());
      return graph_pointer->node((*internal_node_iterator).n_id_);
    }

    /**
    * @brief increments node_iterator to point to the next node_type
    * 
    * @return node_iterator reference corresponding to the next node_type
    * @pre this node_iterator points to an actual node_type in a non-empty graph_type
    * @post (*old node_iterator).index() + 1 == (*new node_iterator).index()
    */
    node_iterator& operator++() {
      assert(graph_pointer != nullptr && internal_node_iterator != graph_pointer->nodes.end());
      ++internal_node_iterator;
      return (*this);
    }

    /**
    * @brief checks whether this node_iterator is equal to another node_iterator nodeiter,
    * which is equivalent to checking equality of the underlying node_types
    * 
    * @param[in] nodeiter   other NodeIterator object to compare with this NodeIterator
    * @return true if the NodeIterator's are equal false otherwise
    */
    bool operator==(const NodeIterator& nodeiter) const {
      return (graph_pointer == nodeiter.graph_pointer && internal_node_iterator == nodeiter.internal_node_iterator);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_pointer;
    typename std::vector<internal_node>::const_iterator internal_node_iterator;

    NodeIterator(const graph_type* graphpointer, typename std::vector<internal_node>::const_iterator intni) {
      graph_pointer = const_cast<graph_type*>(graphpointer);
      internal_node_iterator = intni;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // NodeIterator node_begin() const
  // NodeIterator node_end() const

  /**
  * @brief returns a node_iterator corresponding to the first node of this graph
  * 
  * @return node_iterator pointing to the first node of this graph, if graph is empty returns a nullptr
  * @post if graph is non-empty result can be derefenced, otherwise the result cannot be dereferenced
  */
  node_iterator node_begin() const {
    return node_iterator(this, nodes.begin());
  } 

  /**
  * @brief returns a node_iterator corresponding to one past the last node of this graph_type
  * 
  * @return node_iterator pointing to one past the last node of this graph_type
  * @post result cannot be dereferenced
  */
  node_iterator node_end() const {
    return node_iterator(this, nodes.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      graph_pointer = nullptr;
      node1_id = 0;
      node2_id_iterator = typename std::vector<size_type>::const_iterator();
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
    * @brief dereferences the incident_iterator to return the corresponding edge_type
    *
    * @result edge_type that points to the corresponding edge
    * @pre the incident_iterator points to an actual edge_type in a non_empty graph_type
    * @post (*incident_iterator).node1() == graph_pointer->node(node1_id)
    */
    value_type operator*() const {
      return graph_pointer->edge(node1_id, *node2_id_iterator);
    }

    /**
    * @brief increments the incident_iterator to point to the next edge incident from the main node 1
    *
    * @result incident_iterator reference that points to the next edge
    * @pre the incident_iterator points to an actual edge_type in a non_empty graph_type
    * @post (*old incident_iterator).node1() == (*new incident_iterator).node1()
    */
    incident_iterator& operator++() {
      assert(graph_pointer != nullptr && node2_id_iterator != graph_pointer->adjacency_matrix[node1_id].end());
      ++node2_id_iterator;
      return (*this);
    }

     /**
    * @brief checks whether the edges pointed to by the incident_iterator are equal
    *
    * @param[in] incident_iterator iter to which we compare this incident_iterator
    * @result true if equal false otherwise
    */
    bool operator==(const IncidentIterator& iter) const {
      return (graph_pointer == iter.graph_pointer && node1_id == iter.node1_id && node2_id_iterator == iter.node2_id_iterator);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    friend node_type;
    graph_type* graph_pointer;
    size_type node1_id;
    typename std::vector<size_type>::const_iterator node2_id_iterator;

    IncidentIterator(const graph_type* graphpointer, size_type node1, typename std::vector<size_type>::const_iterator node2iter){
      graph_pointer = const_cast<graph_type*>(graphpointer);
      node1_id = node1;
      node2_id_iterator = node2iter;
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      graph_pointer = nullptr;
      internal_edge_iterator = typename std::vector<internal_edge>::const_iterator();
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
    * @brief dereferences the edge_iterator to obtain the corresponding edge_type
    * 
    * @return edge_type behind the edge_iterator
    * @pre this edge_iterator points to an actual edge_type in a non-empty graph_type
    */
    value_type operator*() {
      assert(graph_pointer != nullptr && internal_edge_iterator != graph_pointer->edges.end());
      return graph_pointer->edge((*internal_edge_iterator).e_id_);
    }

    /**
    * @brief increments edge_iterator to point to the next edge_type
    * 
    * @return edge_iterator reference corresponding to the next edge_type
    * @pre this edge_iterator points to an actual edge_type in a non-empty graph_type
    */
    edge_iterator& operator++() {
      ++internal_edge_iterator;
      return (*this);
    }

//--documentation_-1
//--well done on doc
//--END

    /**
    * @brief checks whether this edge_iterator is equal to another edge_iterator nodeiter,
    * which is equivalent to checking equality of the underlying edge_types
    * 
    * @param[in] eiter   other edge_iterator object to compare with this edge_iterator
    * @return true if the edge_iterator's are equal false otherwise
    */
    bool operator==(const edge_iterator& eiter) const {
      //return (**this) == (*eiter); 
      return (graph_pointer == eiter.graph_pointer && internal_edge_iterator == eiter.internal_edge_iterator);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_pointer;
    typename std::vector<internal_edge>::const_iterator internal_edge_iterator;

    EdgeIterator(const graph_type* graphpointer, typename std::vector<internal_edge>::const_iterator iei) {
      graph_pointer = const_cast<graph_type*>(graphpointer);
      internal_edge_iterator = iei;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
  * @brief returns a edge_iterator corresponding to the first edge of this graph
  * 
  * @return edge_iterator pointing to the first edge of this graph, if graph is empty returns a nullptr
  * @post if graph is non-empty result can be derefenced, otherwise the result cannot be dereferenced
  */
  edge_iterator edge_begin() const {
    return edge_iterator(this, edges.begin());
  }

  /**
  * @brief returns a edge_iterator corresponding to one past the last edge of this graph_type
  * 
  * @return edge_iterator pointing to one past the last edge of this graph_type
  * @post result cannot be dereferenced
  */
  edge_iterator edge_end() const {
    return edge_iterator(this, edges.end());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node { //internal node that actually contains all the node data
    Point position_;   
    size_type n_id_;
    node_value_type value_;     
        };

  struct internal_edge { //internal edge that actually contains all the edge data
  //--functionality_1
  //-- you get your variable names mixed up. did you make last minute changes?
  //--START  
    //size_type node_1_id_;
    //size_type node_2_id_;
    size_type node1_id_;
    size_type node2_id_;
  //--END
    size_type e_id_;
  };

  std::vector<internal_node> nodes;
  std::vector<internal_edge> edges;
  std::vector<std::vector<size_type>> adjacency_matrix;
};

#endif // CME212_GRAPH_HPP
