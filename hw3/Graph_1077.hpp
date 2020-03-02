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

/** @pre/@post
 * Jialu's naming rules in this program:
 * CapitalName: name of a type
 * calptalName: name of a vector Object
 * cpital_name: name of a parameter
 * capital_name_ : name of a class member variable
*/


// =========================== Start: graph class =======================
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */


  // Using templates here
template <typename V,  typename E>


/** @class Graph
 * @brief A class for 3D undirected graphs.
 * Users can add and check information for nodes and edges. Edges are undirected.
 */
class Graph {
 private:

  // Introduce internal private classes:
    class InternalNode;
    class InternalEdge;
    

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
   

  /** Type of this graph. */
  using graph_type = Graph;
   
  /** Type of node_value and edge_value. */
  using node_value_type = V;
  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Type of node_value and edge_value. */
  using edge_value_type = E;
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
    
// -------------------------Start: Node class --------------------------
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
      //Initialize an invalid node with a nullptr
        graph_ = nullptr;
        
    }

    /**  Return the Node's position by reference to make it modifiable.
           * @return The reference of the Node position.
           */
      Point& position(){
          //Test if the node is valid
          assert(graph_!=nullptr && node_id_ < graph_->size());
          // If valid return the corresponding point
          return (graph_->graphNodes[node_id_].node_point_);

      }
    /** Return this node's position. */
    const Point& position() const {
        //Test if the node is valid
        assert(graph_!=nullptr && node_id_ < graph_->size());
        // If valid return the corresponding point
        return (graph_->graphNodes[node_id_].node_point_);
        
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
       //Test if the node is valid
        assert(graph_!=nullptr );
        assert( node_id_ < graph_->size());
        
        return (graph_->graphNodes[node_id_].node_idx_);
    }

      /** Return the node's value by reference to make it modifiable.
      * @return A reference to the value of this node.
      */
     node_value_type& value(){
        //Test if the node is valid
        assert(graph_!=nullptr && node_id_ < graph_->size());
        // If valid return the corresponding value
        return graph_->graphNodes[node_id_].node_value_;
     }

      /** Return the node's value.
          * @return A constant reference to the value of this node.
          */
     const node_value_type& value() const{
        //Test if the node is valid
        assert(graph_!=nullptr && node_id_ < graph_->size());
          // If valid return the corresponding value
        return graph_->graphNodes[node_id_].node_value_;
     }

      /** Find the node's degree.
          * @return The number of edges incident to this node.
          */
     size_type degree() const{
       return graph_->graphEdgesIncidentToNode[node_id_].size();
     }


      /** Find the beginning iterator of all the incident edges
          * @return The iterator pointing to the first incident edge to this node.
          *           return this->edge_end() if this->degree()==0.
          */
     incident_iterator edge_begin() const{
        return IncidentIterator(graph_,
                graph_->graphEdgesIncidentToNode[node_id_].begin(), node_id_);
     }

      /** Find the end iterator of all the incident edges
          * @return The iterator of the last incident edge to the point
          */
     incident_iterator edge_end() const{
        return IncidentIterator(graph_,
         graph_->graphEdgesIncidentToNode[node_id_].end(), node_id_);
     }

     
    /** Test whether this node and @a n are equal.
     *
     * @return True if the two nodes are equal.
     *           False otherwise.
     */
      
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (node_id_ == n.node_id_);
    }

      /** Test whether this node is less than @a n in a global order.
       *  @return True if the index of the node is less than the index of new node.
       *           False otherwise.
       */
    bool operator<(const Node& n) const {
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
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
      // We allow the node class to acess graph
         graph_type* graph_;
      // We have a unique id id for the node
         size_type node_id_;
      // Node constructor initializes a node with a graph
      // and its corresponding internal node id
      /** Constructor*/
         Node(const Graph* graph, size_type id) :
         graph_(const_cast<Graph*>(graph)),
         node_id_(id) {}
      
      
      
  };
    
// -------------------------  End: Node class --------------------------

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {

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
  Node add_node(const Point& position,
                const node_value_type& nval = node_value_type()) {
     
      graphNodeIdx.push_back(this->size());
      InternalNode n;
      n.node_idx_ = this->size();
      n.node_point_ = position;
      n.node_value_ = nval;
      
     
      graphNodes.push_back(n);
     // Update Edges incident to the new node with index graphNodes.size()-1
      std::vector<size_type> incident_edges_to_new_node;
     // Save Edges incident to the new node with in the (graphNodes.size()-1)th
     // member vector in the graphEdgesIncidentToNode vector
      graphEdgesIncidentToNode.push_back(incident_edges_to_new_node);
       
      //***********
      return Node(this,this->size()-1);
  }
    /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
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
  class Edge: private totally_ordered<Edge> {
     // Begin definition of class Edge
   public:
    /** Construct an invalid Edge. */
    Edge() {
        //create an invalid edge with graph pointer as the nullptr
        graph_=nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
        if (graph_ == nullptr){
            return Node();      // Invalid Node
        }
        return Node(graph_, node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
        if (graph_ == nullptr){
                return Node();      // Invalid Node
            }
        return Node(graph_, node2_id_);
    }

    /** Return the length of an Edge under the Euclidean norm */
    double length() const{
          return norm(node1().position() - node2().position());
    }

    /** Return the value of the edge by reference to make it modifiable. */
    edge_value_type& value(){
        //Test if the edge is valid
        assert(graph_!=nullptr && edge_id_ < graph_->num_edges());
        return graph_->graphEdges[edge_id_].edge_value_;
    }

    /** Return the constant value of the edge by reference. */
    const edge_value_type& value() const{
        //Test if the edge is valid
        assert(graph_!=nullptr && edge_id_ < graph_->num_edges());
        // If valid return the corresponding value
        return graph_->graphEdges[edge_id_].edge_value_;
    }

    /** Test whether this edge and @a e are equal.
     * Equal edges represent the same undirected edge between two nodes.
     * @return True if the two iterators are equal.
     *           False otherwise.
     */
    bool operator==(const Edge& e) const {

        bool Condition1 = (node1()==e.node1())&&(node2() ==e.node2());
        bool Condition2 = (node1()==e.node2())&&(node2() ==e.node1());
        return ((graph_ == e.graph_)&&(Condition1 || Condition2));
    }
 
    /** Test whether this edge is less than @a e in a global order.
     *  @return True if the index of the edge is less than the index of new edge.
     *           False otherwise.
     */
    bool operator<(const Edge& e) const {

          if (graph_ < e.graph_){
              return true;
          }
          if( (graph_ == e.graph_)&&
              (graph_->graphEdges[edge_id_].edge_idx_ <
             e.graph_->graphEdges[e.edge_id_].edge_idx_))
          {
              return true;
                   }
        return false;
    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type edge_id_;
    size_type node1_id_;
    size_type node2_id_;
      /** Constructor*/
    Edge(const Graph* graph, size_type edge_id, size_type node1_id,
         size_type node2_id):
                graph_(const_cast<Graph*>(graph)),
                edge_id_(edge_id),
                node1_id_(node1_id),
                node2_id_(node2_id)
                
     {}
                       
   };
    
// ------------------------- End: Edge class --------------------------

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
   
  size_type num_edges() const {

      return this->graphEdges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < this->num_edges());
    return Edge(this, graphEdgeIdx[i], graphEdges[i].node1_id_, graphEdges[i].node2_id_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

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
      if(a==b){
        return Edge(); // invalid edge
      }
      for (size_type i = 0; i < this->num_edges(); i++){
          if ((a.node_id_ == graphEdges[i].node1_id_ && b.node_id_ == graphEdges[i].node2_id_)||
           (a.node_id_ == graphEdges[i].node2_id_ && b.node_id_ == graphEdges[i].node1_id_)){
          return Edge(this, this->graphEdgeIdx[i], a.node_id_, b.node_id_);
          } // if end
       } // for end
       graphEdgeIdx.push_back(num_edges());
      // num_edges count the amount of members in graphEdges but not graphEdgeIdx
       graphEdges.push_back(InternalEdge(a.node_id_, b.node_id_, this->num_edges()));

       graphEdgesIncidentToNode[a.index()].push_back(num_edges()-1);
       graphEdgesIncidentToNode[b.index()].push_back(num_edges()-1);
       return Edge(this,num_edges()-1, a.node_id_, b.node_id_);

  }//add_edge end



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


   Edge add_edge(const Node& a, const Node& b, edge_value_type edge_data) {
        if(a==b){
            return Edge(); // invalid edge
        }
        for (size_type i = 0; i < this->num_edges(); i++){
            if ((a.node_id_ == graphEdges[i].node1_id_ && b.node_id_ == graphEdges[i].node2_id_)||
                (a.node_id_ == graphEdges[i].node2_id_ && b.node_id_ == graphEdges[i].node1_id_)){
                return Edge(this, this->graphEdgeIdx[i], a.node_id_, b.node_id_);
            } // if end
        } // for end
        graphEdgeIdx.push_back(num_edges());
        // num_edges count the amount of members in graphEdges but not graphEdgeIdx
        graphEdges.push_back(InternalEdge(a.node_id_, b.node_id_, num_edges(), edge_data));


        graphEdgesIncidentToNode[a.index()].push_back(num_edges()-1);
        graphEdgesIncidentToNode[b.index()].push_back(num_edges()-1);
        return Edge(this,num_edges()-1, a.node_id_, b.node_id_);

   }//add_edge end

   /** Remove an edge from the graph.
   * @param[in] A Node object passed by reference which is the first node
   * of the edge.
   * @param[in] A Node object passed by reference which is the second node
   * ofthe edge.
   * @return This function returns 0 if the edge is not in the graph
   * and 1 if it is removed from the graph.
   * @post has_node(@a node1, node2) is false.
   *
   * Complexity: O(num_nodes() + num_edges())
   */
    size_type remove_edge(const Node& node1, const Node& node2){
       bool has_edge = false;
       size_type edge_id = 0;

       // See if the edge is in this graph
       // If yes, give has_edge = true an its edge_id;
       // If not, has_edge = false and return 0.
       for (size_type i = 0; i < this->num_edges(); i++){
           if ((node1.node_id_ == graphEdges[i].node1_id_ && node2.node_id_ == graphEdges[i].node2_id_)||
               (node1.node_id_ == graphEdges[i].node2_id_ && node2.node_id_ == graphEdges[i].node1_id_)){
               has_edge = true;
               edge_id = i;
           }
       }

       if(!has_edge){
           return 0;
       }
       //Go through all edges incident to every nodes and delete the edge we want to remove
       size_type num_nodes = size();
       for(size_type i=0; i<num_nodes; ++i){
           std::vector<size_type> edges_adj_node_i = graphEdgesIncidentToNode[i];
           size_type node_i_degree = edges_adj_node_i.size();
           for(size_type j = 0; j < node_i_degree;){
               if (graphEdgesIncidentToNode[i][j] == edge_id){
                   graphEdgesIncidentToNode[i].erase(graphEdgesIncidentToNode[i].begin()+j);
                   --node_i_degree;
               }
               else if(graphEdgesIncidentToNode[i][j] > edge_id){
                   --graphEdgesIncidentToNode[i][j];
                   ++j;
               }
               else{
                   ++j;
               }
           }
       }

       // Delete the edge in the graph
       graphEdges.erase(graphEdges.begin()+edge_id);
       // Size of graphEdges-1 and update the index
       graphEdgeIdx.pop_back();
       return 1;

    }

   /** Remove an edge from the graph.
   * @param[in] A Node object passed by reference which is the first node
   * of the edge.
   * @param[in] A Node object passed by reference which is the second node
   * ofthe edge.
   * @return This function returns 0 if the edge doesn't exist in the graph
   * and 1 if it does and got removed from the graph.
   * @post num_edges() decreases by 1.
   *
   * Complexity: O(num_nodes() + num_edges())
   */
    size_type remove_edge(const Edge& edge){
       return remove_edge(edge.node1(), edge.node2());


    }


    /** Remove an edge from the graph.
    * @param[in] edge_iter An edge_iterator=of the edge we want to remove
     * from the graph.
    * @return This function returns an edge_iterator pointing to the first Edge.
    * @post num_edges() decreases by 1.
    *
    * Complexity: O(num_nodes() + num_edges())
    */
    edge_iterator remove_edge(edge_iterator edge_iter){
        remove_edge((*edge_iter).node1(), (*edge_iter).node2());
        return edge_begin();

    }


   /** Remove a node from the graph.
     * @param[in] node A Node object passed by reference.
     * @return Return 0 if the node is not in the graph
     * and 1 if the node removed from the graph.
     * @post has_node(@a node) is false.
     *
     * Complexity: O(num_nodes() )
     */
   size_type remove_node(const Node& node){
       //Return 0 if the node is not in the graph
       if (!has_node(node))
       {
           return 0;
       }
       //First remove the edges that are incident to this node
       std::vector<size_type> edges_adj_node = graphEdgesIncidentToNode[node.index()];
       size_type node_degree = edges_adj_node.size();

       for(size_type i= 0; i <node_degree ; ){
           Edge e = edge(graphEdgesIncidentToNode[node.index()][i]);
           remove_edge(e);
           node_degree--;
       }

       // For nodes with index greater than the to-be-removed node:
       // reduce their index by one!
       for(size_type i=node.node_id_ +1; i < num_nodes();i++)
       {
           graphNodes[i].node_idx_ -= 1;
       }

       // Update the index of node1 and node2 for all the edges that
       // have node index greater than the index of the to-be-removed node
       for(size_type i=0; i < num_edges(); ++i)
       {
           if(graphEdges[i].node1_id_ > node.node_id_){
               graphEdges[i].node1_id_ -= 1;}
           if(graphEdges[i].node2_id_> node.node_id_){
               graphEdges[i].node2_id_ -= 1;}
       }

       graphEdgesIncidentToNode.erase(graphEdgesIncidentToNode.begin()+node.node_id_);
       graphNodes.erase(graphNodes.begin()+ node.node_id_);
       graphNodeIdx.pop_back();
       return 1;
  }


  /** Remove a node from the graph.
    * @param[in] node_iter A node_iterator of the node to remove from the graph.
    * @return This function returns a node_iterator pointing to the first Node.
    * @post num_nodes() decreases by 1.
    *
    * Complexity: O(num_nodes())
    */
  node_iterator remove_node(node_iterator node_iter){

      remove_node(*node_iter);

      return node_begin();

  }







  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

            graphNodes.clear();
            graphNodeIdx.clear();
            graphEdges.clear();
            graphEdgeIdx.clear();
            graphEdgesIncidentToNode.clear();
            
 }

  // -------------------------  Start: NodeIterator class --------------------------

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }


    // I use the graphNodeIdx to construct NodeIterator for a graph object
      /** Dereference the node iterator.
       * @return The node pointed to by the iterator.
       */
     Node operator*() const {

        return graph_->node(graph_->graphNodeIdx[*node_iterator_]);
     }
    /** Return the node iterator incremented by one */
     NodeIterator& operator++(){
        ++node_iterator_;
        return *this;

     }

      /** Test if two iterators are equal.
        * @return True if the two iterators are equal.
        *           False otherwise.
        */
     bool operator==(const NodeIterator& iter) const{
        return (node_iterator_ == iter.node_iterator_);
     }



   private:
    friend class Graph;
    const Graph* graph_;
    typename std::vector<size_type>::const_iterator node_iterator_;
      /** Constructor*/
    NodeIterator(const Graph* graph,
                 typename std::vector<size_type>::const_iterator iter):
     graph_(const_cast<Graph*>(graph)),
     node_iterator_(iter) {
     }


  };

    /** Return the beginning iterator of all the nodes of the graph */
   node_iterator node_begin() const{
      
      return NodeIterator(this, this->graphNodeIdx.begin());
   }


   /** Return the end iterator of all the nodes of the graph */
   node_iterator node_end() const{
      return NodeIterator(this, this->graphNodeIdx.end());
   }
   
   
   
  // -------------------------  End: NodeIterator class --------------------------

    
  // -------------------------  Start: IncidentIterator class --------------------

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
     
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }
    /** Dereference the incident iterator.
     * @return The edge incident to a node pointed to by the iterator.
     */
     Edge operator*() const{
        size_type curr_edge_id = *incident_iterator_;
        value_type curr_edge = graph_->edge(curr_edge_id);
        size_type n1_id = curr_edge.node1().index();
        size_type n2_id = curr_edge.node2().index();
        
        if(node_index_ == n1_id)
        {
         return Graph::Edge(graph_, graph_->graphEdgeIdx[curr_edge_id],
                               node_index_, n2_id);
        }
        else
        {
         return Graph::Edge(graph_, graph_->graphEdgeIdx[curr_edge_id],
                            node_index_, n1_id);
        }
        
     }
      /** Return the incident iterator incremented by one */
     IncidentIterator& operator++(){
        ++incident_iterator_ ;
        return *this;
     }

      /** Test if two iterators are equal.
        * @return True if the two iterators are equal.
        *           False otherwise.
        */
     bool operator==(const IncidentIterator& iter) const{
     bool flag =((graph_ == iter.graph_)&&
                 (incident_iterator_ == iter.incident_iterator_)&&
                 (node_index_ == iter.node_index_));
        return flag;
     }

   private:
    friend class Graph;

     Graph* graph_;
     typename std::vector<size_type>::const_iterator incident_iterator_;
     size_type node_index_;

      /** Constructor*/
     IncidentIterator(const Graph* graph,
     typename std::vector<size_type>::const_iterator incident_iterator,
                  size_type node_index):
     graph_(const_cast<Graph*>(graph)),
     incident_iterator_(incident_iterator),
     node_index_(node_index){
     }
     

  };

    // -------------------------  End: IncidentIterator class ------------------
    
    // -------------------------Start : EdgeIterator class ---------------------
   
    /** @class Graph::EdgeIterator
    @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // I use the graphEdgeIdx to construct NodeIterator for a graph object
    /** Dereference the edge iterator.
     * @return The edge pointed to by the iterator.
     */
     Edge operator*() const{
        size_type curr_edge_id = *edge_iterator_;
      
        return graph_->edge(graph_->graphEdgeIdx[ curr_edge_id]);
      }
      /** Return the edge iterator incremented by one */
     EdgeIterator& operator++(){
        ++ edge_iterator_;
        return *this;
      }

      /** Test if two iterators are equal.
       * @return True if the two iterators are equal.
       *           False otherwise.
       */
     bool operator==(const EdgeIterator& iter) const{
        return (edge_iterator_ == iter.edge_iterator_);
      }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    typename std::vector<size_type>::const_iterator edge_iterator_;
    /** Constructor*/
    EdgeIterator(const Graph* graph,
                 typename std::vector<size_type>::const_iterator iter):
                  graph_(const_cast<Graph*>(graph)),
                  edge_iterator_(iter) {
    }
                   
                   
                   
  };
    // HW1 #5: YOUR CODE HERE
      /** Return the beginning iterator of all the edges of the graph */
     edge_iterator edge_begin() const {
          return EdgeIterator(this, this->graphEdgeIdx.begin());
      }
      /** Return the end iterator of all the edges of the graph */
     edge_iterator edge_end() const{
        return EdgeIterator(this, this->graphEdgeIdx.end());
                }
   
   
   
   
   
   
// -------------------------End: EdgeIterator class ----------------------------



private:

  // Use this space for implementations of your Graph class's internals:
  //   helper functions, data members, and so forth.
    class InternalNode{

        // Member variables of class InternalNode
        Point node_point_;
        size_type node_idx_;
        node_value_type node_value_;
        
        
        /**Constructor of an invalid InternalNode type object*/
        InternalNode() {}
        /** Constructor*/
        InternalNode(Point point, size_type node_idx):
        //node_point_ (member variabale) is initialized with point (parameter)
        node_point_(point),
        // node_idx_ (member variable) is initialized with idx (parameter)
        node_idx_(node_idx)
        {}
        
        friend class Graph;
      };

    class InternalEdge{

        // Member variables of class InternalEdge
        size_type node1_id_;
        size_type node2_id_;
        size_type edge_idx_;
        edge_value_type edge_value_;



        /**Constructor of an invalid InternalEdge type object*/
        InternalEdge()  {}

        /** Constructor*/
        //Constructor of a valid InternalEdge type object
        InternalEdge(size_type id1, size_type id2, size_type edge_idx):
        node1_id_(id1),//node1_id_ (member variabale) is initialized with id1 (parameter)
        node2_id_(id2),//node2_id_ (member variabale) is initialized with id2 (parameter)
        edge_idx_(edge_idx) //edge_idx_ (member variabale) is initialized with edge_idx (parameter)
        {}

        /** Constructor*/
        // !!!
        //Constructor of a valid InternalEdge type object with initial edge value
        InternalEdge(size_type id1, size_type id2, size_type edge_idx, edge_value_type edge_value):
                node1_id_(id1),//node1_id_ (member variabale) is initialized with id1 (parameter)
                node2_id_(id2),//node2_id_ (member variabale) is initialized with id2 (parameter)
                edge_idx_(edge_idx), //edge_idx_ (member variabale) is initialized with edge_idx (parameter)
                edge_value_(edge_value) //edge_value_ (member variabale) is initialized with edge_value (parameter)
        {}

        friend class Graph;


      };
    
    // Initialze vectors of InternalNode and InternalEdge type objects
    std::vector<InternalNode> graphNodes;
    std::vector<InternalEdge> graphEdges;
    
    // Initialze vectors of idx's of graphNode and graphEdge
    std::vector<size_type> graphNodeIdx;
    std::vector<size_type> graphEdgeIdx;
    
    // Matches the indices of nodes and their incident edges' index
    std::vector<std::vector<size_type> >  graphEdgesIncidentToNode;
    
    

};
 //=====================End: Graph class =======================================

#endif // CME212_GRAPH_HPP



