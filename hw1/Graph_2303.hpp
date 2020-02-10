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

using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template  <typename V>
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

  /** Type of this graph. */
  using graph_type = Graph;

  using  node_value_type = V;
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
  class Node: private totally_ordered<Node>  {
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
      return gr->vec_node[nid].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nid;
    }
   /* Update the value of a node*/
   void update_value(const node_value_type& val){
     gr->vec_node[nid].second=val;
   }

/** 
 * @brief get the value of a node.
 *
 * @pre the node index nid must be in vec_node
 * @return  node_value_type of a node.
 */
   
   node_value_type& value (){
     return gr->vec_node[nid].second;
   }	
   const  node_value_type& value () const{
     return gr->vec_node[nid].second;
   } 
    // HW1: YOUR CODE HERE
   
/** 
 * @brief get the degree of a node.
 *
 * @pre the node index nid must be in vec_node
 * @return  size_type associated to the degree of a node.
 */
   

    size_type degree() const{
     return gr->matrix[nid].size();
    }
/** 
 * @brief returns the first incident_iterators of a Node .
 *
 * @pre the node index nid must be in matrix.
 * @return  the incident iterator.
 */
   

    incident_iterator edge_begin() const{
      vector<size_type>::iterator it= gr->matrix[nid].begin();
      return IncidentIterator(gr,nid,it);
    }
/** 
 * @brief returns the last incident_iterators of a Node .
 *
 * @pre the node index nid must be in matrix.
 * @return  the incident iterator.
 */
    /*Ending the incident_iterator of a node*/
    incident_iterator edge_end() const{
      vector<size_type>::iterator it= gr->matrix[nid].end();
      return IncidentIterator(gr,nid,it);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
               // Quiet compiler warning
      return (n.nid==nid)&&(n.gr==gr);
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
      return (n.nid>nid)&&(n.gr==gr);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* gr;
    size_type nid;

    Node(const Graph* graph, size_type id)
    : gr(const_cast<Graph*>(graph)), nid(id) {}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return vec_node.size();
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
  //--functionality_0
  //--add_node should take in a value as specified in the starter code.
  //--START
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
   vec_node.push_back(make_pair(position,value));
   vector<size_type> ve;
   matrix.push_back(ve);
    return Node(this,vec_node.size()-1);        // Invalid node
  }
  //--END

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this==n.gr;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
   
    return Node(this,i);        // Invalid node
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
  class Edge : private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return n1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return n2;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.gr==gr)&&((e.n1==n1 && e.n2==n2)||(e.n1==n2 && e.n2==n1));
              // Quiet compiler warning
      //HW0: YOUR CODE HERE
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
     
    
      if (e.nedge<nedge){ return true;}
      else if (e.nedge==nedge){
           if (max(n1.nid,n2.nid)>max(e.n1.nid,e.n2.nid)){ return true;}
           else if (max(n1.nid,n2.nid)==max(e.n1.nid,e.n2.nid)){
                 return min(n1.nid,n2.nid)>min(e.n1.nid,e.n2.nid);}
      }
      return false;
      
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
   //--design_1
   //--Edge class not meet proxy size requirements. 
   //--It will be penalized again if not fixed because it does not meet the proxy desgin requirements.
   //--START
   Node n1;
   Node n2;
   size_type nedge;
   Graph* gr;
   //--END 
   Edge(const Graph* graph, Node node1,Node node2,size_type ned)
    : gr(const_cast<Graph*>(graph)), n1(node1), n2(node2),nedge(ned) {}


    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return vec_edge.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
                 // Quiet compiler warning
    return vec_edge[i];        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
   bool val=false;
   if (a.gr==this && this==b.gr){
     for (int i=0;i<vec_edge.size();i=i+1){
        if ((vec_edge[i].n1==a && vec_edge[i].n2==b)||(vec_edge[i].n1==a && vec_edge[i].n2==b)){val=true;}
      }
   }
   return val;
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
  //--functionality_0
  //--add_edge has to check if an edge already exists in the graph.
  //--START
  Edge add_edge(const Node& a, const Node& b) {
   Edge e=Edge(this,a,b,vec_edge.size());
   matrix[a.nid].push_back(vec_edge.size());
   matrix[b.nid].push_back(vec_edge.size());
   vec_edge.push_back(e);

   // Quiet compiler warning
   return e;        // Invalid Edge
  }
  //--END

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    vec_node.clear();
    vec_edge.clear();
    matrix.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator(const Graph* g,const size_type i) 
    : gr(const_cast<Graph*>(g)), it(i){};
    
    NodeIterator() {
    } 
    

   /** 
 * @brief Accessing the node from the iterator.
 * @return  the Node Iterator.
 */ 
    // Accessing the Node from the iterator
    Node operator*() const{
	return Node(gr,it);
    }
   /** 
 * @brief moves to the next iterator.
 * @return  the next iterator.
 */ 

    NodeIterator& operator++(){
       it++;
       return *this;
    }
  
   /** 
 * @brief Testing the equality with another node.
 * @param ni the other node iterator.
 * @return  the Node Iterator.
 */ 

    bool operator==(const NodeIterator& ni) const{
       return (gr==ni.gr && it==ni.it);
   }

   private:
    friend class Graph;
    Graph* gr;
    size_type it;
   
  };

  
    /** 
 * @brief get the first node_iterator.
 * @return  the node_iterator.
 */ 

  node_iterator node_begin() const{
     return NodeIterator(this,0);
}
    /** 
 * @brief get the first node_iterator.
 * @return  the node_iterator.
 */ 
  node_iterator node_end() const{
     return NodeIterator(this,this->vec_node.size());
}


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(const Graph* g,size_type i,vector<size_type>::iterator itt)
    : gr(const_cast<Graph*>(g)),in(i),it(itt){};

    IncidentIterator() {
    }
   
 // HW1 #3: YOUR CODE HERE
      /** 
 * @brief get the edge from the incident_iterator associated.
 * @pre *it must be a cell in vec_edge.
 * @return  the edge.
 */ 
   Edge operator*() const{
      //--functionality_1
      //--node1 of the returned edge should be the node corresponding to the incident iterator.
      //--Returning an edge like this makes the ordering arbitrary.
      //--START
      return gr->vec_edge[*it];
      //--END
    }
     /** 
 * @brief moves to the next incident_iterator.
 * @return  the incident_iterator.
 */ 
    IncidentIterator& operator++(){
      it++;
      return *this;
   }
   /** 
 * @brief test equality with an other iterator
 * @param other the incident iterator to compare
 * @return  a boolean.
 */ 
   bool operator==(const IncidentIterator& other) const{
      return (gr==other.gr && in==other.in && it==other.it);
   }

   private:
    friend class Graph;
     Graph* gr;
     size_type in;
     vector<size_type>::iterator it;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(const Graph*  g, size_type ind):
   gr(const_cast<Graph*>(g)), index(ind){};

   EdgeIterator() {
    };
   

      /** 
 * @brief get the edge associated to an edge iterator.
 * @pre the index must be in vec_edge.
 * @return  the edge.
 */ 
    Edge operator*() const{
     return gr->vec_edge[index];
   }
   /** 
 * @brief moves to the next  edge_iterator.
 * @return  the next edge_iterator.
 */ 
    EdgeIterator& operator++(){
     index++;
     return *this;
   }
   /** 
 * @brief test the equality between two EdgeIterators.
 * @return  the EdgeIterator.
 */ 
    bool operator==(const EdgeIterator &ei ) const{
      return (gr==ei.gr && index==ei.index);
   }


   private:
    friend class Graph;
     Graph* gr;
     size_type index;
  };
  

  
     /** 
 * @brief get the first edge_iterator.
 * @return  the edge_iterator.
 */ 
  edge_iterator edge_begin() const{
     return EdgeIterator(this,0);
}
     /** 
 * @brief get the last node_iterator.
 * @return  the last_iterator.
 */ 
   edge_iterator edge_end() const{
     return EdgeIterator(this,this->vec_edge.size());
}

 private:

  
  // vec_node stores the value and the coordinate of each Node according to its index
  vector<pair<Point, node_value_type>> vec_node;
 // vec_edge stores the value and the coordinate of each Node according to its index
  vector<Edge> vec_edge;
 //This vector stores a vector associated to the edge index of each node/
  vector<vector<size_type>> matrix;
};

#endif // CME212_GRAPH_HPP
