#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */



#include <algorithm>
#include <vector>
#include <cassert>
#include <queue>
#include <set>
#include <unordered_map>

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
     
   struct nodeindex;
   struct internal_edges;
   struct internal_nodes;  
  
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
  using node_value_type = V;
  using edge_value_type = E;

  

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
  using ne_data=std::unordered_map<size_type, std::pair<size_type, bool>>;


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph.  
  Here next_nid, next_eid are used in the construction of nodes and edges.  node_size, edge_size are used to record the size of the whole graph
  in terms of nodes and edges, rather than the size of a single node or edge.
   */
 Graph() :  node_idx(), edge_idx(), graph_nodes(), graph_edges(), next_nid(0), next_eid(0), node_size(0), edge_size(0)  { 
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
    /** Construct an i
    id node.
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

    //Construct an invalid node with nullptr as the graph pointer.
    // HW0: YOUR CODE HERE
       graph_  = nullptr;  
    }

    /** Return this node's position. */
    
    const Point& position() const {
        //Check whether the node is a `valid` node and return the corresponding point.
   if (graph_->node_idx[node_id]==-1) {
            std::cout<<node_id<<"\n";
        }
        assert(graph_->node_idx[node_id]!=-1);
        return graph_->graph_nodes[(unsigned)graph_->node_idx[node_id]].point;
    }

//Modifiable position.
Point& position() {
        //Test if the node is valid and return the corresponding point
        if (graph_->node_idx[node_id]==-1) {
            std::cout<<node_id<<"\n";
        }
        assert(graph_->node_idx[node_id]!=-1);
        return graph_->graph_nodes[(unsigned)graph_->node_idx[node_id]].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //Check 
      assert(graph_!=nullptr);
      assert(graph_->node_idx[node_id]!=-1);
      return (unsigned)graph_->node_idx[node_id];
      //return size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

         /** Find the node's data/value.
         * @return A the reference to the value of this node.
         * @pre The node is valid.
         */
        node_value_type& value() {
       //Check.
        assert(graph_->node_idx[node_id]!=-1);
        return graph_->graph_nodes[index()].data;
        }
        /** Find the node's value.
         * @return A constant reference to the value of this node.
         *@post The node address is unchanged.
         */
        const node_value_type& value() const {
            //Check.
            assert(graph_->node_idx[node_id]!=-1);
            return graph_->graph_nodes[index()].data;
        }
        /** Find the node's degree.
         * @return The number of its neighbors.
         */
        size_type degree() const {
            return graph_->graph_nodes[index()].neighbors_edges.size();
        }
     

     /** Find the beginning iterator of the collection of incident edges.
         * @return The iterator pointing to the first incident edge to this node.
         */
        incident_iterator edge_begin() {
            typename ne_data::iterator it_temp =graph_->graph_nodes[this->index()].neighbors_edges.begin();
            return IncidentIterator(graph_,this, it_temp);
        }
      /** Find the las iterator of the collection of incident edges.
         * @return The iterator pointing to the last incident edge to this node.
         */
        incident_iterator edge_end() {
            typename ne_data::iterator it_temp=graph_->graph_nodes[this->index()].neighbors_edges.end();
            return IncidentIterator(graph_,this, it_temp);
        }
        




    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
        //check equality through checking pointers.
        if ( (this->graph_ == n.graph_) && (this->node_id == n.node_id) ) {
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
      //Determine the order of nodes according to their indices.
      if((this->graph_ < n.graph_) || (this->node_id < n.node_id))
    { return true;
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
    
    //accessor to the graph from the node.
    graph_type* graph_ ;
    //The id of this node
    size_type node_id;
    //Construct a node and assign its id.
    Node(const Graph* graph, size_type id) : graph_(const_cast<Graph*>(graph)), node_id(id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    //The way to count node_size will be specified below.
    return node_size;
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
    
     Node add_node(const Point& position,const node_value_type& n_val = node_value_type())  {
        node_idx.push_back(node_size);
 
        next_nid++;
        
        internal_nodes n1;

        n1.nid = next_nid-1;
        n1.point = position;
        n1.data = n_val;
        
        if (graph_nodes.size() > node_size) {
            graph_nodes[node_size]=n1;
        }
        else {
            graph_nodes.push_back(n1);
        }
        //Increase the number of nodes by 1.
        node_size++;
        //Return this node.
        return Node(this, next_nid-1);
    }



  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
   // HW0: YOUR CODE HERE
    //Check validity of this node, and return the result of checking.
   return (node_idx[n.node_id]!=-1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // Check whether the index is within the range, and then 
    // return the corresponding node.
    assert(i<node_size);
    return Node(this, graph_nodes[i].nid);   
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
      graph_=nullptr;
    }

    /** Return a node of this Edge. The `false` orientation is equivalent to exchanging these two nodes. */
    Node node1() const {
      // HW0: YOUR CODE HERE
      //Determine the order of endpoints according to the bool orientation.
        if (orientation) {
            return graph_->graph_edges[index()].end_1st;
        }
        else {
            return graph_->graph_edges[index()].end_2nd;
        }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
        if (orientation) {
            return graph_->graph_edges[index()].end_2nd;
        }
        else {
            return graph_->graph_edges[index()].end_1st;
        }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((this->edge_id==e.edge_id)&&(this->graph_ == e.graph_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
    //Compare two edges according to their ids.
      return ((this->graph_ < e.graph_) || (this->edge_id<e.edge_id));
    }

/** Find the edge's value.
     * @return A reference to the value of this edge.
     */
    edge_value_type& value() {
        return graph_->graph_edges[index()].data;
    }
    /** Find the node's value.
     * @return A const reference to the value of this node.
     */
    const edge_value_type& value() const {
        return graph_->graph_edges[index()].data;
    }




   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // Accessor to the graph from the edges.
    graph_type* graph_;
    //The index of this edge.
    size_type  edge_id;
    //Use a bool to determine the order of nodes in an edge.
    bool orientation;
    //Return the index of this edge.
    size_type index() const {
      //Check validity.
    assert(graph_!=nullptr);
    assert(graph_->edge_idx[edge_id]!=-1);
    return ((unsigned)graph_->edge_idx[edge_id]);
  }

  //Initialize an edge with given data. The initial orientation is set to be true.
  Edge(const Graph* graph, size_type eid, bool o1=true): graph_(const_cast<Graph*>(graph)), edge_id(eid), orientation(o1) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_size;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //Check whether this index is within the desired range.
    assert(i<edge_size);
    return Edge(this, graph_edges[i].eid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //Check whether b is in node a's neighbor.
    return ( graph_nodes[a.index()].neighbors_edges.count(b.node_id) > 0);
  }





  /** Add an edge to the graph, or return the current edge if it already exists. The default initialization of value enables
  this method to take only two nodes as variables to construct an edge.
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

  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val=edge_value_type()) {
    // HW0: YOUR CODE HERE
    
   //Check whether this edge exists already.
    if (graph_nodes[a.index()].neighbors_edges.count(b.node_id) > 0) {
          return Edge(this, graph_nodes[a.index()].neighbors_edges[b.node_id].first);
      }
      
    //Set up the data of this edge.  
    internal_edges e1;
    e1.eid = next_eid;
    e1.end_1st = a;
    e1.end_2nd = b;
    e1.data = val;
    
    //Record this edge in the data of neighbors.
    graph_nodes[a.index()].neighbors_edges[b.node_id]=(std::pair<size_type,bool>(e1.eid, true));
    graph_nodes[b.index()].neighbors_edges[a.node_id]=(std::pair<size_type,bool>(e1.eid, false));

    edge_idx.push_back(edge_size);
    if (graph_edges.size() > edge_size) {
        graph_edges[edge_size]=e1;
    }
    else {
      //Add this edge to the end of the collection of edges.
        graph_edges.push_back(e1);}
    next_eid++;
    edge_size++;
    return Edge(this, next_eid-1);
}


//HW2 Part 7: Removing edges and nodes.

  /** A method that is used to remove an edge in a given graph. Return an integer indicating whether the deletion is successful.
     * @pre @a e is an edge of this graph.
     * @return 1 if the edge is sucessfully deleted.
     *         0 if the edge is invalid.
     * @post has_edge(@a e.node1(), @a e.node2()) == false
     * @post if return 1 , new num_edges() == old num_edges() - 1.
     *       Else, new num_edges() == old num_edges().
     *
     *
     * Time Complexity: O(1) on average, O(max degree) in the worst case.
     */
    size_type remove_edge(const Edge& e) {
        if (!has_edge(e.node1(),e.node2())) {
            return 0;
        }
        Node a=e.node1();   
        Node b=e.node2();
        if (edge_size > 1)   {
            graph_edges[e.index()]=graph_edges[edge_size-1];
            edge_idx[graph_edges[e.index()].eid]=e.index();
        }
        graph_nodes[a.index()].neighbors_edges.erase(b.node_id);
        graph_nodes[b.index()].neighbors_edges.erase(a.node_id);
        edge_idx[e.edge_id]=-1;
        edge_size--;
        return 1;
    }

  /** A method that is used to remove an edge (given by its end points) in a given graph. 
  Return an integer indicating whether the deletion is successful. 
     * @pre @a e is an edge of this graph.
     * @return 1 if the edge is sucessfully deleted.
     *         0 if the edge is invalid.
     * @post has_edge(@a e.node1(), @a e.node2()) == false
     * @post if return 1 , new num_edges() == old num_edges() - 1.
     *       Else, new num_edges() == old num_edges().
     *
     *
     * Time Complexity: O(1) on average, O(d_max) in the worst case. d_max is the maximum degree.
     */
    
    size_type remove_edge(const Node& a, const Node&b) {
        if (!has_edge(a,b)) {
            return 0;
        }
        return remove_edge(add_edge(a,b));
    }

// Iterator version of above methods. We do not return the integer for indication now.
  edge_iterator remove_edge(edge_iterator e_iter) {
        remove_edge(*e_iter);
        return e_iter;
    }


//Methods to remove nodes.

 /** A method that is used to remove a node in a given graph. Return an integer indicating whether the deletion is successful.
     * @pre @a n is a valid or invalid node of this graph.
     * @return 1 if the node is sucessfully deleted.
     *         0 if the node is invalid.
     * @post has_node(@a n) == false
     * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
     *       Else,                  new num_nodes() == old num_nodes().
     *
     *
     * Complexity: O(d_max * T), where T is the time complexity for remove_edge and d_max is the maximum degree.
     */
    
    size_type remove_node(const Node& n) {
        //When the node is invalid.
        if (!has_node(n)) {
            return 0;
        }
        //remove the edges incident to this node.
        while (graph_nodes[n.index()].neighbors_edges.size()>0) {
            ne_data::iterator it=graph_nodes[n.index()].neighbors_edges.begin();
            remove_edge(Edge(this,(*it).second.first));
        }
        if (node_size>1)    {
            graph_nodes[n.index()]=graph_nodes[node_size-1];
            node_idx[graph_nodes[n.index()].nid]=n.index();
        }
        node_idx[n.node_id]=-1;
        node_size--;
        return 1;
    }
    
    /** Delete a node from the graph. We do not return the integer for indication now.
     * @pre @a n_it is a vlid node iterator of this graph.
     * @return the iterator pointing to the next node.
     * @post has_node(@a (*n_it)) == false
     * @post If old has_node(@a (*n_it)), new num_nodes() == old num_nodes() - 1.
     *       Else,                        new num_nodes() == old num_nodes().
     *
     *
     * Complexity: O(d_max *　T), where T is the time complexity for remove_edge and d_max is the maximum degree.
     */
    
    node_iterator remove_node(node_iterator n_iter)  {
        remove_node(*n_iter);
        return n_iter;
    }
    









  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    //Delete all data and set sizes to 0.
    node_idx.clear();
    edge_idx.clear();
    graph_nodes.clear();
    graph_edges.clear();
    next_nid=0;
    next_eid=0;
    node_size=0;
    edge_size=0;
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
    NodeIterator() {
    	graph_  = nullptr;
      iter_node = nullptr;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return the node which the iterator is pointing to 
     *@pre The node is a valid node.
    */
      Node operator*() const {
          //Check if the iterator is valid.
          assert(graph_ != nullptr);
          return Node(graph_, iter_node->nid);
      }
      
    /** Return the incremented iterator, node version.*/
      NodeIterator& operator++() {
          iter_node ++;
          return *this;
      }
    // bool operator==(const IncidentIterator&) const
    /** Check if two iterators are equal
     * @return True if the iterators are the same.
     * 'iter_temp' is the other object being compared with.
     */
      bool operator==(const NodeIterator& it_temp) const{
          return (iter_node == it_temp.iter_node);
      }



   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
     Graph* graph_;
     typename std::vector<internal_nodes>::const_iterator iter_node;
/** Constructor */
NodeIterator(const Graph* graph, typename std::vector<internal_nodes>::iterator it_initial): graph_(const_cast<Graph*>(graph)), iter_node(it_initial) {}
/** Constructor, const_iterator version.*/
NodeIterator(const Graph* graph, typename std::vector<internal_nodes>::const_iterator it_initial): graph_(const_cast<Graph*>(graph)), iter_node(it_initial) {}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

    
    /** Return the beginning iterator of the collection of nodes.*/
    node_iterator node_begin() const {
        return NodeIterator(this, graph_nodes.begin());
    }
    /** Return the last iterator, by adding the increment of the address.*/
    node_iterator node_end() const {
        typename std::vector<internal_nodes>::const_iterator it_temp=graph_nodes.begin()+node_size;
        return NodeIterator(this, it_temp);
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
    	graph_  = nullptr;
      t_node  = nullptr;
      iter_incident = nullptr;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const



    /** Return the edge which the iterator is pointing to 
     *@pre The node is a valid node.
    */
      Edge operator*() const {
          //Check if the iterator is valid.
          assert(graph_ != nullptr);
          std::pair<size_type, bool> edge_temp=iter_incident->second;
          return graph_->edge(edge_temp.first,edge_temp.second);
      }
    // IncidentIterator& operator++()
    /** Return the incremented iterator*/
      IncidentIterator& operator++() {
          iter_incident ++;
          return *this;
      }
    // bool operator==(const IncidentIterator&) const
    /** Check if two iterators are equal
     * @return True if the iterators have the same graph, nodes and corresponding incident edges.
     *           False otherwise.
     */
      bool operator==(const IncidentIterator& it_temp) const {
          return ((graph_ == it_temp.graph_) && (t_node == it_temp.t_node) && (iter_incident == it_temp.iter_incident));
      }





   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

     //The parent graph.
      Graph* graph_;
      //The node to which this edge is incident.
      Node* t_node;
      //Iiterator of this edge.
      typename ne_data::iterator iter_incident;
      /** Constructor */
        IncidentIterator(const Graph* graph, const Node* node, typename ne_data::iterator it_initial): graph_(const_cast<Graph*>(graph)), 
        t_node(const_cast<Node*>(node)), iter_incident(it_initial) {}     
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
      iter_edge=nullptr;
      graph_=nullptr;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Evaluation operator of the iterator.
     * @return The edge pointed to by the iterator.
     */

      Edge operator*() const {
          return Edge(graph_, iter_edge->eid);
      }
      
    /** `Next` operaotr of the iterator.
    @return the next iterator.*/
      EdgeIterator& operator++() {
          iter_edge++;
          return *this;
      }
    /** Check whether two iterators are the same.
     * @return True: if it_edge of these two edges are equal.
     *         False: if they are not equal.
     */
      
      bool operator==(const EdgeIterator& other_edge) const {
          return (iter_edge == other_edge.iter_edge);
      }


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // The parent graph.
    Graph* graph_;
    // The iterator that goes through graph's vector containing data of the edges
      typename std::vector<internal_edges>::const_iterator iter_edge;
      /** Constructor */
      EdgeIterator(const Graph* graph, typename std::vector<internal_edges>::iterator it_initial): graph_(const_cast<Graph*>(graph)), iter_edge(it_initial) {}
      /** Constructor, const version.*/
      EdgeIterator(const Graph* graph, typename std::vector<internal_edges>::const_iterator it_initial): graph_(const_cast<Graph*>(graph)), iter_edge(it_initial) {}
  };

 
  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

    /**
    @return The beiginning iterator of the collection of edges.
   */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, graph_edges.begin());
    }
    /**
    @return The last iterator of the collection of edges.
   */
    edge_iterator edge_end() const {
        typename std::vector<internal_edges>::const_iterator iter_temp=graph_edges.begin()+edge_size;
        return EdgeIterator(this, iter_temp);
    }




private: 

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

//The struct of nodes, which include its index and position.
   struct internal_nodes {
    //Use nid to identify all nodes that we have ever addded.
    size_type nid;
    Point point;
    ne_data neighbors_edges;
    node_value_type data;
   };
  
//The struct of edges, which include its end points.
   struct internal_edges {
     size_type eid; 
     node_type end_1st;
     node_type end_2nd;
     edge_value_type data;
   };
    


    
    /** Return the edge with the desired orientation and designated edge internal id */
    Edge edge(size_type id, bool ott) {
        return Edge(this, id, ott);
    }
    


   // node_idx records the correspondence between node internal id and node's index
   // the same applies to edges
   // the node's and edge's index are then used to access points and nodes in the corresponding vector 
   std::vector<int> node_idx;
   std::vector<int> edge_idx;

     //The collection of nodes.
   std::vector<internal_nodes> graph_nodes;
   //The collection of edges.
   std::vector<internal_edges> graph_edges;
     // next_nid and next_eid determine the next node id and the next edge id respecrtively
   size_type next_nid;
   size_type next_eid;

     size_type node_size;
     size_type edge_size;

};

#endif // CME212_GRAPH_HPP
