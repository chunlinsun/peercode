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

/** Naming rules in this program:
CapitalName: name of a type
calptalName: name of a vector Object
cpital_name: name of a parameter
capital_name_ : name of a class member variable

*/


// =========================== Start: graph class =======================
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
  // To solve the problem of re-index, we give those nodes two kind of
  // indentification, nodeid_(the position in the internalnode) and index.
  // Each time we add or delete a point, we only need to add or del the point in
  // the internalnode and, then, adjust map bewteen index and node_id_.
    
  // Introduce internal private classes:
    class InternalNode;
    class InternalEdge;
    

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
  Graph() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() = default;
    
// -------------------------Start: Node class --------------------------
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
      //Initialize an invalid node with a nullptr
        graph_ = nullptr;
        
    }

    /** Return this node's position. */
    const Point& position() const {
        // HW0: YOUR CODE HERE
        
        //Test if the node is valid
        assert(graph_!=nullptr && node_id_ < graph_->size());
        // If valid return the corresponding point
        return (graph_->graphNodes[node_id_].node_point_);
        
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        // HW0: YOUR CODE HERE
        
        //Test if the node is valid
        assert(graph_!=nullptr && node_id_ < graph_->size());
        
        return (graph_->graphNodes[node_id_].node_idx_);
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
      if (graph_ != n.graph_) {
          return false;
      }
      if (node_id_ != n.node_id_){
          return false;
      }
      return true;
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
      // Use Lexicographic order to decide the "<" relation
        if (graph_ < n.graph_){
            return true;
        }
        if( (graph_ == n.graph_) && (node_id_ < n.node_id_)){
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
    
      // We allow the node class to acess graph
         graph_type* graph_;
      // We have a unique id id for the node
         size_type node_id_;
      // Node constructor initializes a node with a graph and its corresponding internal node id
         Node(const Graph* graph, size_type id) : graph_(const_cast<Graph*>(graph)), node_id_(id) {}
      
      
      
  };
    
// -------------------------  End: Node class --------------------------

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
      
      return this->graphNodes.size();
     
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
    (void) position;      // Quiet compiler warning
      this->graphNodeIdx.push_back(this->size());
      this->graphNodes.push_back(InternalNode(position,this->size()));
      return Node(this,this->size()-1);
  }
    /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
      (void) n;   // Quiet compiler warning
      //Test if the id of the node is contained in the graphNodeId vector
      return ((this == n.graph_) && (n.node_id_ < this->size()));
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
    assert(i < this->size());
    return Node(this, this->graphNodes[i].node_idx_);
  }

// -------------------------Start: Edge class --------------------------

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {     // Begin definition of class Edge
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
        //create an invalid edge with graph pointer as the nullptr
        graph_=nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
        if (graph_ == nullptr){
            return Node();      // Invalid Node
        }
        return Node(graph_, graph_->graphEdges[edge_id_].node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
        if (graph_ == nullptr){
                return Node();      // Invalid Node
            }
        return Node(graph_, graph_->graphEdges[edge_id_].node2_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
      
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0 : YOUR CODE HERE
        bool Condition1 = (node1()==e.node1())&&(node2() ==e.node2());
        bool Condition2 = (node1()==e.node2())&&(node2() ==e.node1());
        return ((graph_ == e.graph_)&&(Condition1 || Condition2));
    }
 
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
        //HW0 : YOUR CODE HERE
          if (graph_ < e.graph_){
              return true;
          }
          if( (graph_ == e.graph_)&&
              (graph_->graphEdges[edge_id_].edge_idx_ <
             e.graph_->graphEdges[edge_id_].edge_idx_))
          {
              return true;
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
    graph_type* graph_;
    size_type edge_id_;
          
    Edge(const Graph* graph, size_type edge_id):
                graph_(const_cast<Graph*>(graph)),
                edge_id_(edge_id){}
    };
    
// ------------------------- End: Edge class --------------------------

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    
    
    
    
    
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
      return this->graphEdges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    assert(i < this->num_edges());
    return Edge(this,graphEdgeIdx[i]);
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
    for (size_type i = 0; i < this->num_edges(); i++){
    if ((a.node_id_ == graphEdges[i].node1_id_ && b.node_id_ == graphEdges[i].node2_id_)||
        (a.node_id_ == graphEdges[i].node2_id_ && b.node_id_ == graphEdges[i].node1_id_)){
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
    // HW0: YOUR CODE HERE'
    (void) a, (void) b;   // Quiet compiler warning
    if(a==b){
        return Edge(); // invalid edge
    }
    if( !has_edge(a,b)){
        graphEdgeIdx.push_back(this->num_edges());
        graphEdges.push_back(InternalEdge(a.node_id_, b.node_id_, this->num_edges()));
        return Edge(this, this->num_edges()-1);
      }
    else{
    for (size_type i = 0; i < this->num_edges(); i++){
    if ((a.node_id_ == graphEdges[i].node1_id_ && b.node_id_ == graphEdges[i].node2_id_)||
        (a.node_id_ == graphEdges[i].node2_id_ && b.node_id_ == graphEdges[i].node1_id_)){
        return Edge(this, this->graphEdgeIdx[i]);
    } // if end
    } // for end

  } //else end
  }//add_edge end

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
            graphNodes.clear();
            graphNodeIdx.clear();
            graphEdges.clear();
            graphEdgeIdx.clear();
 }

  // -------------------------  Start: NodeIterator class --------------------------

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
  // -------------------------  End: NodeIterator class --------------------------

    
  // -------------------------  Start: IncidentIterator class --------------------

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

    // -------------------------  End: IncidentIterator class ----------------------
    
    // -------------------------Start : EdgeIterator class -------------------------
   
    /** @class Graph::EdgeIterator
    @brief Iterator class for edges. A forward iterator. */
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
// -------------------------End: EdgeIterator class ------------------------------

 private:

  // HW0: YOUR CODE HERE
  // Use this space for implementations of your Graph class's internals:
  //   helper functions, data members, and so forth.
    class InternalNode{
        
        // Member variables of class InternalNode
        Point node_point_;
        size_type node_idx_;
        
        
        
        //Constructor of an invalid InternalNode type object
        InternalNode() {}
        //Constructor of a valid InternalNode type object
        InternalNode(Point point, size_type node_idx):
        node_point_(point),//node_point_ (member variabale) is initialized with point (parameter)
        node_idx_(node_idx)// node_idx_ (member variable) is initialized with idx (parameter)
        {}
        
        friend class Graph;
      };

    class InternalEdge{
        // Member variables of class InternalNode
        size_type node1_id_;
        size_type node2_id_;
        size_type edge_idx_;
        
        
        
        //Constructor of an invalid InternalEdge type object
        InternalEdge()  {}
        
        //Constructor of a valid InternalEdge type object
        InternalEdge(size_type id1, size_type id2, size_type edge_idx):
        node1_id_(id1),//node1_id_ (member variabale) is initialized with id1 (parameter)
        node2_id_(id2),//node2_id_ (member variabale) is initialized with id2 (parameter)
        edge_idx_(edge_idx) //edge_idx_ (member variabale) is initialized with edge_idx (parameter)
        {}
        
        friend class Graph;
        
        
      };
    
    // Initialze vectors of InternalNode and InternalEdge type objects
    std::vector<InternalNode> graphNodes;
    std::vector<InternalEdge> graphEdges;
    
    // Initialze vectors of idx's of graphNode and graphEdge
    std::vector<size_type> graphNodeIdx;
    std::vector<size_type> graphEdgeIdx;
    
    
    
    
    
    
    
    

};
 //=====================End: Graph class =======================================

#endif // CME212_GRAPH_HPP



