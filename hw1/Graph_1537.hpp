#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */
 #include <string>
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
  using  node_value_type = V;
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
  Graph() :nodes_(), edges_(), node_size_(0), edge_size_(0), next_node_id_(0), next_edge_id_(0), node_index_vec_(), edge_index_vec_()
   {
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
      //--design_0
      //--Technically no need for anything here. Can be an empty block.
      //--START
      graph_ptr_=nullptr;
      //--END
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      size_type ind = graph_ptr_->node_index_vec_[node_id_];
      return graph_ptr_->nodes_[ind].internal_point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      size_type index = graph_ptr_->node_index_vec_[node_id_];
      return index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    node_value_type& value (){
      size_type ind = graph_ptr_->node_index_vec_[node_id_];
      return graph_ptr_->nodes_[ind].internal_node_value_;
    //--style_0
    //--Blank line between methods please
    //--START
    }
    const  node_value_type& value () const{
    //--END
      size_type ind = graph_ptr_->node_index_vec_[node_id_];

      return graph_ptr_->nodes_[ind].internal_node_value_;
    }

    //--documentation_1
    //--No specifications?
    //--START
    size_type degree() const{
      size_type ind = graph_ptr_->node_index_vec_[node_id_];
      return graph_ptr_->nodes_[ind].incident_edges_.size();

    }

    incident_iterator edge_begin() const{
      size_type ind = graph_ptr_->node_index_vec_[node_id_];
      return IncidentIterator(graph_ptr_, this,graph_ptr_->nodes_[ind].incident_edges_.begin() );
    }

    incident_iterator edge_end() const{
      size_type ind = graph_ptr_->node_index_vec_[node_id_];
      return IncidentIterator(graph_ptr_, this,graph_ptr_->nodes_[ind].incident_edges_.end() );
    }
    //--END




    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph_ptr_==n.graph_ptr_ && this->index() == n.index();
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
      if (graph_ptr_ == n.graph_ptr_){
      return this->node_id_ < n.node_id_;
      }
      return std::less<Graph *>{}(graph_ptr_, n.graph_ptr_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_ptr_;
    size_type node_id_;
    Node(const graph_type* graph_ptr, size_type node_id)
    : graph_ptr_(const_cast<graph_type*>(graph_ptr)), node_id_(node_id){

    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_size_;
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
  Node add_node(const Point& position, const  node_value_type& val = node_value_type ()) {
    // HW0: YOUR CODE HERE
    
    //--functionality_0
    //--internal_node has 4 members, but you only provide 3 arguments to the
    //--constructor here. Please fix. Several other compiler warnings that need
    //--to be fixed.
    //--START
    internal_node new_node{position,next_node_id_,val};
    //--END
    //--style_0
    //--Remove commented out code before submission
    //--START
    //new_node.internal_point_ = position;
    //new_node.internal_node_id_ = next_node_id_;
    //new_node.internal_node_value_ = val;
    //--END
    node_index_vec_.push_back(node_size_);
    nodes_.push_back(new_node);
    ++node_size_;
    ++next_node_id_;
    return Node(this, next_node_id_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //--functionality_1
    //--Missing implementation?
    //--START
    // HW0: YOUR CODE HERE
    //// update!!!!!!!!!!
    
    return this == n.graph_ptr_ ;
    //--END
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i>=0);
    assert(i<size());
    return Node(this,nodes_[i].internal_node_id_);       
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
      graph_ptr_=nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      assert(graph_ptr_);
      //size_type ind = graph_ptr_->edge_index_vec_[edge_id_];
      //return graph_ptr_->edges_[ind].first_node_;
      return Node(graph_ptr_, n1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      assert(graph_ptr_);
      //size_type ind = graph_ptr_->edge_index_vec_[edge_id_];
      //return graph_ptr_->edges_[ind].second_node_;
      return Node(graph_ptr_, n2_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return e.graph_ptr_ == this->graph_ptr_ && e.edge_id_ == this->edge_id_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ptr_ == e.graph_ptr_){
        return this->edge_id_<e.edge_id_;
      }
      return std::less<graph_type*>{}(graph_ptr_, e.graph_ptr_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_ptr_;
    size_type edge_id_;
    size_type n1_id_;
    size_type n2_id_;
    Edge(const graph_type* graph_ptr,  size_type edge_id, size_type n1_id, size_type n2_id)
    : graph_ptr_(const_cast<graph_type*>(graph_ptr)), edge_id_(edge_id), n1_id_(n1_id), n2_id_(n2_id){

    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i>=0);
    assert(i<num_edges());
    //std::cout<<edges_.size()<<"\n";
    return Edge(this, edges_[i].internal_edge_id_, edges_[i].first_node_.node_id_,  edges_[i].second_node_.node_id_);       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(this->has_node(a));
    assert(this->has_node(b));
    //--design_0
    //--This implementation is O(num_edges). Even without changing your data
    //--structures, you can reduce it to O(max(deg(a), deg(b))). Do you see how?
    //--The same thing can be applied to add_edge also.
    //--START
    for (size_type i=0; i<edges_.size(); i++){
      size_type a_id = edges_[i].first_node_.node_id_;
      size_type b_id = edges_[i].second_node_.node_id_;
      if ((a_id == a.node_id_ && b_id == b.node_id_)
      ||(a_id == b.node_id_ && b_id == a.node_id_))
        return true;
    }
    //--END

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
    assert(this->has_node(a));
    assert(this->has_node(b));
    
    // First, loop through the edges to see if the edge already exists between a and b
    // if exists, return that edge
    for (size_type i=0; i<edges_.size(); i++){
      size_type a_id = edges_[i].first_node_.node_id_;
      size_type b_id = edges_[i].second_node_.node_id_;
      if ((a_id == a.node_id_ && b_id == b.node_id_)
      ||(a_id == b.node_id_ && b_id == a.node_id_))
      return Edge(this, edges_[i].internal_edge_id_,a.node_id_, b.node_id_);
    }
    // if the edge does not exist, add the edge to the end of edge_index_vec and edges and return itself
      internal_edge new_edge {a,b,next_edge_id_};
      edge_index_vec_.push_back(edge_size_);
      edges_.push_back(new_edge);
      edge_type edge_for_a = Edge(this, next_edge_id_, a.node_id_, b.node_id_);
      edge_type edge_for_b = Edge(this, next_edge_id_, b.node_id_, a.node_id_);
      this->nodes_[a.node_id_].incident_edges_.push_back(edge_for_a);
      this->nodes_[b.node_id_].incident_edges_.push_back(edge_for_b);
      ++edge_size_;
      ++next_edge_id_;
      return Edge(this, next_edge_id_-1, a.node_id_, b.node_id_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
      nodes_.clear();
      edges_.clear();
      node_size_ = 0;
      edge_size_ = 0;
      next_node_id_ = 0;
      next_edge_id_ = 0;
      node_index_vec_.clear();
      edge_index_vec_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}
    

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const{
      return Node(graph_pt_,node_it_);
    }

    NodeIterator& operator++(){
      ++node_it_;
      return *this;
    }

    bool operator==(const NodeIterator& it) const{
      return (it.graph_pt_ == graph_pt_)&&(it.node_it_==node_it_);
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type node_it_ ;
    const Graph* graph_pt_ ;
    NodeIterator(const Graph* graph_pointer, size_type node_index):
		graph_pt_(graph_pointer), node_it_(node_index) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }

  node_iterator node_end() const{
    return NodeIterator(this,node_size_);

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
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const{
      return *it_;
    }
    IncidentIterator& operator++(){
      it_++;
      return *this;

    }
    bool operator==(const IncidentIterator& inci) const{
      return (graph_pt_ == inci.graph_pt_)&&(node_pt_ == inci.node_pt_)&&(it_==inci.it_);
    }
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph_pt_;
    const Node* node_pt_;
    typename std::vector<edge_type>::iterator it_;
    IncidentIterator(const Graph* graph_pointer,const Node* node_pointer,typename std::vector<edge_type>::iterator it):
		graph_pt_(graph_pointer), node_pt_(node_pointer), it_(it) {}
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
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    //--documentation_1
    //--No documentation/specifications?
    //--START
    Edge operator*() const{
      internal_edge internal_edge_cur = *it_;
      return Edge(graph_pt_,
        internal_edge_cur.internal_edge_id_,
        internal_edge_cur.first_node_.node_id_,
        internal_edge_cur.second_node_.node_id_);
    }
    EdgeIterator& operator++(){
      it_++;
      return *this;

    }
    bool operator==(const EdgeIterator& edge_it) const{
      return (graph_pt_== edge_it.graph_pt_)&&(it_ == edge_it.it_);
    }
    //--END
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_pt_;
    typename std::vector<internal_edge>::const_iterator  it_;
    EdgeIterator(const Graph* graph_pointer, typename std::vector<internal_edge>::const_iterator it):
		graph_pt_(graph_pointer), it_(it) {}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const{
    return EdgeIterator(this, edges_.begin());
}

edge_iterator edge_end() const{
    return EdgeIterator(this, edges_.end());
}
 private:
  struct internal_node{
    Point internal_point_;
    size_type internal_node_id_;
    //added node value and also edges incident to that node
    node_value_type internal_node_value_;
    //--design_1
    //--A graph should not store Edge objects! The Edge object is just an
    //--outward-facing type for the user. You are certainly okay to store
    //--internal_edge structs here, but storing Edge objects breaks the proxy
    //--pattern.
    //--START
    std::vector<edge_type> incident_edges_;
    //--END
  };

  struct internal_edge{
    node_type first_node_;
    node_type second_node_;
    size_type internal_edge_id_;

  };
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  size_type node_size_;
  size_type edge_size_;
  size_type next_node_id_;
  size_type next_edge_id_;
  std::vector<size_type> node_index_vec_;
  std::vector<size_type> edge_index_vec_;
};

//--style_1
//--Very little care paid to spacing, indentation, commented out code, etc.
//--makes the code difficult to read. Please take some time to clean up the
//--formatting and make the code easier to read and understand.
//--END

//--documentation_0
//--It'd be nice to have the README.md file cleaned up; it contains debris from
//--a git merge conflict.
//--END

#endif // CME212_GRAPH_HPP
