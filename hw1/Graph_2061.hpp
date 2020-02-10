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

  std::vector<internal_node> nodes_;  // node proxy
  std::vector<internal_edge> edges_;  // edge proxy
  std::vector<std::vector<std::pair<unsigned,unsigned>> > neighbors_; 
     // list of nodes (and corresponding edges) adjacent to each node

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type for class template */
  using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph<V>;  //modified

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
  Graph() :  nodes_(),edges_(), neighbors_() {
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
      // HW0: YOUR CODE HERE - unchanged
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[nid_].pos;  
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return nid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for: 

    /**
     * @brief Access the user-specified value stored supported by Node
     *
     * @return Reference to the user-specified value
     */
    node_value_type& value() {
      return graph_->nodes_[nid_].val;
    }

    /**
     * @brief Access the user-specified value stored supported by Node
     *
     * @return Const reference to the user-specified value
     * 
     * Promising not to modify the value
     */
    const node_value_type& value() const {
      return graph_->nodes_[nid_].val;
    }

    /**
     * @brief Get the degree of this Node
     *
     * @return Degree of this Node
     */
    size_type degree() const{
      return graph_->neighbors_[nid_].size();
    }

    /**
     * @brief Obtain start of the iterator over all edges incident to this Node
     
     * @return Start of the incident iterator
     *
     * Complexity: O(1).
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, nid_, 0);
    }

    /**
     * @brief Obtain end of the iterator over all edges incident to this Node
     *
     * @return End of the incident iterator
     *
     * Complexity: O(1). 
     * Iterating over all the edges incident to a Node n has maximum complexity O(n.degree())
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, nid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE - modified
      bool eq_id = nid_==n.nid_;
      bool eq_graph = graph_==n.graph_;
      return (eq_id && eq_graph);
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
      // HW0: YOUR CODE HERE - modified
      return (graph_==n.graph_) && (nid_<n.nid_);
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE - modified
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph 
    Graph* graph_;
    // Node index
    size_type nid_; 
    /** Private Constrctor */
    Node(const graph_type* graph, size_type nid)
     : graph_(const_cast<graph_type*>(graph)), nid_(nid) {  //modified
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE 
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type new_id = num_nodes();

    internal_node new_node(position,new_id,value);
    nodes_.push_back(new_node);  // create an internal node: takes O(1)
    neighbors_.push_back(std::vector<std::pair<size_type,size_type>>()); // insertion to vector is O(1) 

    return Node(this,new_id);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE - modified
    //check if node index within range 
    if (num_nodes()<=n.index()) {  
      return false;
    }
    internal_node node_found = nodes_[n.index()];
    return node_found.pos==n.position();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE - modified
    assert(i<num_nodes());
    return Node(this,i);        
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_,nid1_); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,nid2_); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     *
     * Complexity: O(1).
     */
    bool operator==(const Edge& e) const {
      bool eq_graph = graph_==e.graph_;
      bool eq_nodes = ((nid1_==e.nid1_ && nid2_==e.nid2_) || (nid1_==e.nid2_ && nid2_==e.nid1_));
      return (eq_graph && eq_nodes);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((graph_==e.graph_) && (eid_ < e.eid_));

    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph 
    Graph* graph_;
    // Edge index
    size_type eid_; 
    // Node indices
    size_type nid1_;
    size_type nid2_;
    /** Private Constrctor */
    Edge(const graph_type* graph, size_type eid, size_type nid1, size_type nid2)
     : graph_(const_cast<graph_type*>(graph)), eid_(eid), nid1_(nid1), nid2_(nid2) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i<num_edges()); // the index needs to be within range
    return Edge(this, i, edges_[i].nid1, edges_[i].nid2);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    size_type nid_a = a.index();
    size_type nid_b = b.index();
    // get the list of neighbors (nid, eid) of Node a
    std::vector<std::pair<size_type,size_type>> neighbors_a = neighbors_[nid_a];
    // check if Node b is one of the neighbors
    // complexity at most O(num_edges())
    for(std::vector<std::pair<size_type,size_type>>::const_iterator iter = neighbors_a.begin(); iter != neighbors_a.end();++iter) {
        if (iter->first == nid_b){
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

    size_type nid_a = a.index();
    size_type nid_b = b.index();  
    size_type new_eid = num_edges();  //+1
    if (has_edge(a,b)){ // takes O(num_edges())
      // make no change to the vector of edges
      // find and return the existing edge by locating the specific Node b in the neighbors of Node a
      // complexity at most O(num_edges())
      std::vector<std::pair<size_type,size_type>> neighbors_a = neighbors_[nid_a];
      for(std::vector<std::pair<size_type,size_type>>::const_iterator iter = neighbors_a.begin(); iter != neighbors_a.end();++iter) {
          if (iter->first == nid_b){
            return Edge(this, iter->second, nid_a, nid_b);
          }
      }  
    }else{
      // create an internal edge: takes O(1)
      internal_edge new_edge(nid_a, nid_b, new_eid);
      edges_.push_back(new_edge);  
      // update neighbor list of related nodes
      // searching in vector is O(1) 
      neighbors_[nid_a].push_back({nid_b,new_eid});
      neighbors_[nid_b].push_back({nid_a,new_eid});      
    }
    return Edge(this, new_eid, nid_a, nid_b);        // Edge that accesses the created internal edge


  };

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    neighbors_.clear();
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

    /**
     * @brief Dereference operator of node iterator
     *
     * @return Node object that the current iterator points to
     *
     * Complexity: O(1)
     */
    Node operator*() const 
    {
      return graphPtr_->node(nodeIterIdx_);
    }

    /**
     * @brief Increment operator of node iterator
     *
     * @param[out] nodeIterIdx_ Index of node the iterator currently points to 
     * @return Reference to the current iterator, with pointer pointing to the next
     *
     * Complexity: O(1)
     */
    NodeIterator& operator++()
    {
      ++nodeIterIdx_;
      return *this;
    }

    /**
     * @brief Equality operator of node iterator

     * @param[in] node_iter Another node's iterator
     * @return True if _this_ and _node_iter_ are the same
     *         (belonging to the same graph and pointing to the same node)
     *
     * Complexity: O(1)
     */
    bool operator==(const NodeIterator& node_iter) const
    {
      return (graphPtr_==node_iter.graphPtr_)&&(nodeIterIdx_==node_iter.nodeIterIdx_);
    }

    /**
     * @brief Inequality operator of node iterator

     * @param[in] node_iter Another node's iterator
     * @return True if _this_ and _node_iter_ are not the same
     *         (not belonging to the same graph or not pointing to the same node)
     *
     * Complexity: O(1)
     */
    bool operator!=(const NodeIterator& node_iter) const
    {
      return (graphPtr_!=node_iter.graphPtr_)||(nodeIterIdx_!=node_iter.nodeIterIdx_);
    }

   private:
    friend class Graph;
    // Pointer to graph
    graph_type* graphPtr_;
    // Node index
    size_type nodeIterIdx_;
    /** Private Constrctor */
    NodeIterator(const graph_type* graph, size_type nid) : 
    graphPtr_(const_cast<graph_type*>(graph)), nodeIterIdx_(nid) {}

  };


  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:


  /**
   * @brief Get the start of the iterator over all nodes
   
   * @return Start of the node iterator
   *
   * Complexity: O(1).
   */
  node_iterator node_begin() const {
    return NodeIterator(this,0);
  }

  /**
   * @brief Get the end of the iterator over all nodes
   
   * @return End of the node iterator
   *
   * Complexity: O(1).
   */
  node_iterator node_end() const {
    return NodeIterator(this,num_nodes());
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

    /**
     * @brief Dereference operator of incident iterator
     *
     * @return Edge object that the current iterator points to
     *
     * Complexity: O(1)
     */
    Edge operator*() const {

      std::vector<std::pair<size_type,size_type>> neighbors_1 = graphPtr_->neighbors_[nid1_];
      std::pair<size_type,size_type> cur_edge = neighbors_1[adjIdx_];
      size_type nid2 = cur_edge.first;
      size_type eid = cur_edge.second;

      return Edge(graphPtr_, eid, nid1_, nid2);
    }

    /**
     * @brief Increment operator of incident iterator
     *
     * @param[out] adjIdx_ Index of node in adjacency list that the iterator currently points to 
     * @return Reference to the current incident iterator, with pointer pointing to the next in adjacency list
     *
     * Complexity: O(1)
     */
    IncidentIterator& operator++() {
      ++adjIdx_;
      return *this;
    }

    /**
     * @brief Equality operator of incident iterator

     * @param[in] incidentIter Another incident iterator
     * @return True if _this_ and _incidentIter_ are the same
     *         (belonging to the same graph, incidenting to the same node and pointing to the same neighbor)
     *
     * Complexity: O(1)
     */
    bool operator==(const IncidentIterator& incidentIter) const {
      return ((graphPtr_==incidentIter.graphPtr_) && (nid1_==incidentIter.nid1_) && (adjIdx_==incidentIter.adjIdx_));
    }
    
    /**
     * @brief Inequality operator of incident iterator

     * @param[in] incidentIter Another incident iterator
     * @return True if _this_ and _incidentIter_ are not the same 
     *         (not belonging to the same graph, not incidenting to the same node or not pointing to the same neighbor)
     *
     * Complexity: O(1)
     */
    bool operator!=(const IncidentIterator& incidentIter) const {
      return ((graphPtr_!=incidentIter.graphPtr_) || (nid1_!=incidentIter.nid1_) || (adjIdx_!=incidentIter.adjIdx_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Pointer to graph
    graph_type* graphPtr_;
    // Index of the current node
    size_type nid1_;
    // Index of the adjacent node
    size_type adjIdx_;
    /** Private Constructor */
    IncidentIterator(const graph_type* graph, const size_type nid1, size_type adjIdx): 
      graphPtr_(const_cast<graph_type*>(graph)), nid1_(nid1), adjIdx_(adjIdx) {}

  };

//--design_0
//--good changes to the underlying structure, well done on this homework 
//--END 

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

    /**
     * @brief Dereference operator of edge iterator
     *
     * @return Edge object that the current iterator points to
     *
     * Complexity: O(1)
     */
    Edge operator*() const {
      size_type nid_1 = graphPtr_->edges_[edgeIterIdx_].nid1;
      size_type nid_2 = graphPtr_->edges_[edgeIterIdx_].nid2;
      return Edge(graphPtr_, edgeIterIdx_, nid_1, nid_2);
    }

    /**
     * @brief Increment operator of incident iterator
     *
     * @param[out] edgeIterIdx_ Index of edge the iterator currently points to
     * @return Reference to the current edge iterator, with pointer pointing to the next 
     *
     * Complexity: O(1)
     */
    EdgeIterator& operator++() {
      ++edgeIterIdx_;
      return *this;
    }

//--documentation_0
//--good doxygen comments
//--END

    /**
     * @brief Equality operator of edge iterator

     * @param[in] edgeIter Another edge iterator
     * @return True if _this_ and _edgeIter_ are the same
     *         (belonging to the same graph and pointing to the same edge)
     *
     * Complexity: O(1)
     */
    bool operator==(const EdgeIterator& edgeIter) const {
      return ((graphPtr_==edgeIter.graphPtr_) && (edgeIterIdx_==edgeIter.edgeIterIdx_));
    }

    /**
     * @brief Inequality operator of edge iterator

     * @param[in] edgeIter Another edge iterator
     * @return True if _this_ and _edgeIter_ are not the same
     *         (not belonging to the same graph or not pointing to the same edge)
     *
     * Complexity: O(1)
     */
    bool operator!=(const EdgeIterator& edgeIter) const {
      return ((graphPtr_!=edgeIter.graphPtr_) || (edgeIterIdx_!=edgeIter.edgeIterIdx_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Pointer to graph
    graph_type* graphPtr_;
    // Edge index
    size_type edgeIterIdx_;
    /** Private Constructor */
    EdgeIterator(const graph_type* graph, size_type eid) : 
    graphPtr_(const_cast<graph_type*>(graph)), edgeIterIdx_(eid) {}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /**
   * @brief Get the start of the iterator over all edges
   
   * @return Start of the edge iterator
   *
   * Complexity: O(1).
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this,0);
  }

  /**
   * @brief Get the end of the iterator over all edges
   
   * @return End of the edge iterator
   *
   * Complexity: O(1).
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this,num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    Point pos; // position
    size_type nid; // node id
    node_value_type val; // value to be specified
    /** Constructor */
    internal_node(const Point &position, size_type node_index, node_value_type value):
      pos(position), nid(node_index), val(value) {}
  };

  struct internal_edge {
    size_type nid1;  // first node's id
    size_type nid2;  // second node's id
    size_type eid;   // edge id
    /** Constructor */
    internal_edge(size_type node_index_1, size_type node_index_2, size_type edge_index):
      nid1(node_index_1), nid2(node_index_2), eid(edge_index) {}
  };

};

#endif // CME212_GRAPH_HPP
