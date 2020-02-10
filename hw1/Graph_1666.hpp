#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <tuple>

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
  
  unsigned nnode;
  unsigned nedge;


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using  node_value_type = double;

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
  std::vector<Node> nodelist;
  std::vector<Edge> edgelist;
  std::map<std::tuple<Node,Node>,unsigned> edgemap;
  std::map<Node,std::vector<Edge>> incidentmap;
  
  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    nnode=0;
    nedge=0;
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return realnode->point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return realnode->nid;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
  /** get the value stored
   * @return value stored
   */ 
    node_value_type& value(){
       return realnode->val;
    }
    /** set the value of the node
    *@para[i] value 
    */ 
    void setvalue(node_value_type v){
       realnode->val=v;
    }
  /** return value stored
   *@return the value of the node
   */
    const node_value_type& value() const{
       return realnode->val;
    }
    /** return number of edges
   * @return number of edges
   */
    //--functionality_1
    //--Does not return correct output, perhaps because you are using
    //--Node objects as keys.
    //--START
    size_type degree() const{
      return realnode->graph_->incidentmap[*this].size();
    }
    //--END
  /** return the start of the incident iterator
   * @return  start of the iterator
   */
    incident_iterator edge_begin() const{
      //--functionality_0
      //--Fails test that checks if node1 of incidentiterator == node 
      //--Why not use a constructor to set the members?
      //--START
      incident_iterator it;
      it.index=0;
      it.list=it.list=const_cast<std::vector<Edge>*>(&(realnode->graph_->incidentmap[*this]));
      return it;
      //--END
    }
    /** return the end of the incident iterator
   * @return  end of the iterator
   */
    incident_iterator edge_end() const{
      incident_iterator it;
      it.index=realnode->graph_->incidentmap[*this].size();
      it.list=const_cast<std::vector<Edge>*>(&(realnode->graph_->incidentmap[*this]));
      return it;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // Quiet compiler warning
      if(realnode->graph_==n.realnode->graph_){
        if(realnode->nid==n.realnode->nid){
          return true;
        }
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
      // Quiet compiler warning
      if(realnode->graph_<n.realnode->graph_){
        return true;
      }
      if(realnode->graph_==n.realnode->graph_){
        if(realnode->nid<n.realnode->nid){
          return true;
        }
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    internal_node* realnode;
    Node(Graph* g,const Point& p, unsigned i,node_value_type v){
      realnode=new internal_node(g,p,i,v);
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nnode;
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
  Node add_node(const Point& position, const  node_value_type& v = node_value_type ()) {
    // HW0: YOUR CODE HERE
    nodelist.push_back(Node(this,position,nnode,v));
    nnode+=1;    // Quiet compiler warning
    return nodelist[nnode-1];        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if(n.realnode->nid<nnode){
      if(nodelist[n.realnode->nid]==n){
        return true;
      }
    }           // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE           // Quiet compiler warning
    return nodelist[i];        // Invalid node
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return *n1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return *n2;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        if(graph_==e.graph_&&eid==e.eid){
          return true;
        }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(graph_<e.graph_){
        return true;
      }
      if(graph_==e.graph_){
        if(eid<e.eid){
          return true;
        }
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const Node* n1;
    const Node* n2;
    Graph* graph_;
    unsigned eid;
    Edge(const Node* a1, const Node* a2,Graph* g,unsigned i):n1(a1),n2(a2),graph_(g),eid(i){ }
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
    return nedge;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
           // Quiet compiler warning
    return edgelist[i] ;       // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    std::tuple<Node,Node> ed=std::make_tuple(a,b);
    //std::tuple<Node,Node> ed2=std::make_tuple(b,a);
    if(edgemap.count(ed)){
      return true;
    }
    // for(unsigned i=0;i<nedge;i++){
    //   if((edgelist[i].node1()==a)&&(edgelist[i].node2()==b)){
    //     return true;
    //   }
    //   if((edgelist[i].node2()==a)&&(edgelist[i].node1()==b)){
    //     return true;
    //   }
    // }
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
    std::tuple<Node,Node> ed=std::make_tuple(a,b);
    std::tuple<Node,Node> ed2=std::make_tuple(b,a);
    if(has_edge(a,b)){
       return edgelist[edgemap[ed]];
    }
    if(has_edge(b,a)){
       return edgelist[edgemap[ed2]];
    }
    //--style_1
    //--Please don't turn in code with commented out code blocks.
    //--START

    // for(unsigned i=0;i<nedge;i++){
    //   if((edgelist[i].node1()==a)&&(edgelist[i].node2()==b)){
    //     return edgelist[i];
    //   }
    //   if((edgelist[i].node2()==a)&&(edgelist[i].node1()==b)){
    //     return edgelist[i];
    //   }
    // }
    //--END

    edgelist.push_back(Edge(&a,&b,this,nedge));
    edgemap.insert({ed,nedge});
    nedge+=1;
    if(incidentmap.count(a)){
      incidentmap[a].push_back(edgelist[nedge-1]);
    }
    else{
      std::vector<Edge> newlist;
      newlist.push_back(edgelist[nedge-1]);
      incidentmap.insert({a,newlist});
    }
    return edgelist[nedge-1];        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nnode=0;
    nedge=0;
    nodelist.clear();
    edgelist.clear();
    edgemap.clear();
    incidentmap.clear();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy
    std::vector<Node>* list;
    unsigned int index;
    /** Construct an invalid NodeIterator. */
    
    NodeIterator() {
      index=0;
    }
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const{
      return (*list)[index];
    }
    NodeIterator& operator++(){
      index+=1;
      return *this;
    }
    bool operator==(const NodeIterator& n) const{
      if(index==n.index&&list==n.list){
        return true;
      }
      return false;
    }
    bool operator!=(const NodeIterator& n) const{
      if(index!=n.index || list!=n.list){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    
    // HW1 #2: YOUR CODE HERE
  };
  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** return the start of the node iterator
   * @return  start of the node iterator
   */
  node_iterator node_begin() const{
    node_iterator it;
    it.index=0;
    it.list=const_cast<std::vector<Node>*>(&nodelist);
    return it;
  }
  /** return the end of the node iterator
   * @return  end of the iterator
   */
  node_iterator node_end() const{
    node_iterator it;
    it.index=nodelist.size();
    it.list=const_cast<std::vector<Node>*>(&nodelist);
    return it;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy
    unsigned int index;
    std::vector<Edge>* list;
    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      index=0;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
   
    /** dereference method
   * @return  the edge that the iterator is poiting at 
   */
    Edge operator*() const{
      return (*list)[index];
    }
    /** Increament of iterator
    *@return the updated iterator
    */
    IncidentIterator& operator++(){
      index+=1;
      return *this;
    }
    /** test equality of two iterator
    *para[in] IncidentIterator i
    *return true if the iterator is the same as the current one
    */
    bool operator==(const IncidentIterator& i) const{
      if(list==i.list && index==i.index){
        return true;
      }
      return false;
    }
    /** test inequality of two iterator
    *para[in] IncidentIterator i
    *return true if the iterator is not the same as the current one
    */
    bool operator!=(const IncidentIterator& i) const{
      if(list==i.list && index==i.index){
        return false;
      }
      return true;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy
    std::vector<Edge>* list;
    unsigned int index;
    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      index=0;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** dereference method
   * @return  the edge that the iterator is poiting at 
   */
    Edge operator*() const{
      return (*list)[index];
    }
     /** Increament of iterator
    *@return the updated iterator
    */
    EdgeIterator& operator++(){
      index+=1;
      return *this;
    }
     /** test equality of two iterator
    *para[in] EdgeIterator e
    *return true if the iterator is the same as the current one
    */
    bool operator==(const EdgeIterator& e) const{
      //--functionality_0
      //--Using == on different vector instances will return False
      //--even if they have the same content.
      //--This breaks your == operator.
      //--START
      if(list==e.list&&index==e.index){
        return true;
      }
      //--END
      return false;
    }
     /** test equality of two iterator
    *para[in] EdgeIterator e
    *return true if the iterator is not the same as the current one
    */
    bool operator!=(const EdgeIterator& e) const{
      if(list==e.list&&index==e.index){
        return false;
      }
      return true;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
 /** return the start of the node iterator
   * @return  start of the iterator
   */
  edge_iterator edge_begin() const{
    edge_iterator it;
    it.index=0;
    it.list=const_cast<std::vector<Edge>*>(&edgelist);;
    return it;
  }
   /** return the end of the edge iterator
   * @return  end of the iterator
   */
  edge_iterator edge_end() const{
    edge_iterator it;
    it.index=edgelist.size();
    it.list=const_cast<std::vector<Edge>*>(&edgelist);;
    return it;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node{
     Graph* graph_;
     Point point_;
     unsigned nid;
     node_value_type val;
     internal_node(Graph* g,const Point p, unsigned i,node_value_type v){
      graph_=g;
      point_=p;
      nid=i;
      val=v;
    }
    };

};

#endif // CME212_GRAPH_HPP
