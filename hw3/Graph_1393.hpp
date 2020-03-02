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
 
 // Intern_Node is defined at the end of this code. 
 struct Intern_Node ;  

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
 
  /** Node value type - Graph is template */
  using node_value_type = V ; 

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
    
  /** Edge value type */
  using edge_value_type = E;

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
    Graph(): nodes_(), edge_adj_(), num_edges_(0), index2nodeid_()
            {
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
    }

    /**
     * @brief Return this node's position.
     */
    const Point& position() const {
      return node_graph_->nodes_[node_idx_].position_ ;
    }
      
    /** @brief Return this node's position.
     *  The position is returned as a reference that ca be updated.
     *
     */
    Point& position() {
        return node_graph_->nodes_[node_idx_].position_;
      }

    /**
     * @brief Return this node's index, a number in the range [0, graph_size).
     */
    size_type index() const {
      // We return the attribute `node_idx_` from the `Node` class.
      return node_graph_->nodes_[node_idx_].node_idx_ ; 
    }
      
     /**
      *  @brief Returns the node value as a reference that can be modified.
      */
       node_value_type& value(){
           return node_graph_->nodes_[node_idx_].val_ ; 
       }
      
     /**
      * @brief Returns the node value as a reference that cannot be modified.
      */
       const node_value_type& value() const{
           return node_graph_->nodes_[node_idx_].val_; 
       }
    /**
     * @brief Returns the degree of a node.
     */
      size_type degree() const{
          return node_graph_->edge_adj_[node_idx_].size() ;
      }
    /**
     * @Brief Methods to initialize the start of EdgeIncident iterator.
     */
      incident_iterator edge_begin() const{
          return IncidentIterator(node_graph_, node_idx_, 0);
      }
      
     /**
      * @Brief Methods to initialize the end EdgeIncident iterator.
      */
      incident_iterator edge_end() const{
          return IncidentIterator(node_graph_, node_idx_, degree());
      }

     /**
      * @brief Test whether this node and _n_ are equal.
      *
      * Equal nodes have the same graph and the same index.
      */
      bool operator==(const Node& n) const {
          (void) n;          // Quiet compiler warning
          bool res = false ;
          // We compare both the pointers value and the index of the nodes.
          if ((n.node_graph_ == node_graph_) && (n.node_idx_ == node_idx_)){res = true ;}
            return res;
      }

    /**
     * @brief Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      (void) n;           // Quiet compiler warning
      // Two cases are possible : if the nodes belong to the same graph,
      // we compare their IDs. 
      if (n.node_graph_ == node_graph_){
          return node_idx_ < n.node_idx_ ; }
      // Otherwise we compare the address of the graphs.
      return node_graph_ < n.node_graph_ ;
    }
      
    /** @brief Boolean function to determine if a node is valid.
     *
     * @returns true      if 0 <= @a node_id_ < @a nodes_.size() (less than the
     *                    total number of nodes ever) and
     *                    0 <= @a node_idx_ < @a index2nodeid_.size() (node
     *                    index is in range)
     *          false     else
     */
    bool valid() const {
      return node_idx_ >= 0 && node_idx_ < node_graph_->nodes_.size() &&
             node_graph_->nodes_[node_idx_].node_idx_ < node_graph_->index2nodeid_.size();
    }
      
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
   
    // We declare a pointer to the graph associated with the node.
    Graph* node_graph_ ; 
   
    // Node ID for the node.  
    size_type node_idx_ ; 
   
    // A private constructor for Node class. 
    Node(const Graph* node_graph, size_type node_idx)
     : node_graph_(const_cast<Graph*>(node_graph)), node_idx_(node_idx){
     }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // We return the size of the container `nodes_`. 
    return index2nodeid_.size();
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
   * Modified for the purpose of HW2.
   */
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
      // New node index equals the current total number of nodes the graph has
      // ever had, which equals the number of elements of edge_adj_
      size_type new_node_id = edge_adj_.size();
      // New node index equals the current number of nodes remaining in the graph,
      // which is the size of index2nodeid
      size_type new_node_idx = index2nodeid_.size();
      index2nodeid_.push_back(new_node_id);
      // We use the Intern_Node struct to instantiate the new node.
      Intern_Node newNode = Intern_Node(new_node_idx, position, value);
      nodes_.push_back(newNode);
      // We add an empty vector because the new node has no neighbors for the moment
      edge_adj_.push_back(std::vector<std::tuple<size_type,edge_value_type>>());
      return Node(this, new_node_id);
    }
    
    /** @brief Given a node @a a, remove it from the graph (if it exists in the
     * graph), and return @a node_id_ OR 0 otherwise.
     *
     * @param[in]  a            Node to attempt to remove
     *
     * Complexity: O(num_nodes() + num_edges_), hopefully less.
     */
     size_type remove_node(const Node& a){

       if(!a.valid()){
           std::cout << "Node is invalid" << std::endl ;
           return size_type(0);}

       unsigned int num_edges_del = 0;

       // Make a **copy** of a's neighbors to iterate over as we delete them
       std::vector<std::tuple<size_type,edge_value_type>> neighbors =
                                                          edge_adj_[a.node_idx_];
       // Remove all edges incident to this node via remove_edge() method
       for(unsigned int k = 0; k < neighbors.size(); k++){
         
         size_type neighbor_id = std::get<0>(neighbors[k]);
         // Attempt to delete an edge
         remove_edge(a, Node(this, neighbor_id));
         // A new edges is deleted.
         num_edges_del++;
       }
        
       // Clear all pointers to this node in edge_ajd.
       edge_adj_[a.node_idx_].clear();
         
       // Erase the node ID, and return a pointer to the next position in
       // index2nodeid_ if the last node ID was
       // deleted)
       auto iter_to_erase = index2nodeid_.erase(index2nodeid_.begin() + a.index());
         
       // Cost of while loop: O(num_nodes) since in the worst case, the first node was
       // deleted and all other node indices need to be updated
       while(iter_to_erase != index2nodeid_.end()){
         nodes_[*iter_to_erase].node_idx_ = iter_to_erase - index2nodeid_.begin();
         ++iter_to_erase;
       }
       return a.node_idx_;
     }
    
    /** @brief Given a NodeIterator @a n_it, remove the associated node from the graph (if it exists in the
     *  graph).
     *
     * @param[in]  n_it            NodeIterator
     *
     * It uses the -- operator from NodeIterator class (written for the purpose of HW2).
     * Complexity: O(num_nodes() + num_edges_), hopefully less.
     */
    NodeIterator remove_node(NodeIterator n_it){
      node_type node_to_remove = *n_it;
      remove_node(node_to_remove);
      return n_it.node_iter_idx_ != node_end().node_iter_idx_ - 1
             ? n_it : --n_it;}
    
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;            // Quiet compiler warning
    // We test for the equality of the pointers value. 
    return n.node_graph_ == this ;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    (void) i;             // Quiet compiler warning
    return Node(this, index2nodeid_[i]);   // Invalid node
  }
 
  //
  // EDGES
  //

  /** @class Graph::Edge
   *  @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(edge_graph_, node1_idx_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(edge_graph_, node2_idx_);      // Invalid Node
    }

    /**
     * @brief Test whether this edge and _e_ are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      (void) e;           // Quiet compiler warning
      // The edge is assumed to be undirected. 
      // We test for the nodes index equality. 
      bool nodes_match = (e.node1_idx_ == node1_idx_ && e.node2_idx_ == node2_idx_)
                      || (e.node1_idx_ == node2_idx_ && e.node2_idx_ == node1_idx_);
      // We test for the pointers value equality. 
      bool graphs_match = e.edge_graph_ == edge_graph_;

      return (nodes_match && graphs_match);

    }

    /**
     * @brief Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {

      (void) e ; // Quiet compiler warning.

      // Do the edges belong to the same graph?
      if (edge_graph_ ==e.edge_graph_){
      // Test if the nodes are equal. 
      if (e.node1_idx_ == node1_idx_){
      // then we compare the other pair of nodes. 
          return node2_idx_ < e.node2_idx_;
        }
        // Then compare first nodes
        return node1_idx_ < e.node1_idx_;
      }
      // Otherwise, let us compare the pointer value. 
      return edge_graph_ < e.edge_graph_;
    }
      
    /* Return the length of this edge, which is the Euclidean distance between
     * its nodes. */
      double length() const {
        return norm(node1().position() - node2().position());
      }
      
    /** @brief Return this edge's value as a reference that CAN be updated or
     *  modified.
     */
    edge_value_type& value() {
        size_type min_idx = std::min(node1_idx_, node2_idx_);
        size_type max_idx = std::max(node1_idx_, node2_idx_);
        // Value is stored in the entry of edge_adj_ corresponding to edge
        // representation (min_id, max_id)
        std::vector<std::tuple<size_type,edge_value_type>> neighbordata = edge_graph_->edge_adj_[min_idx];
        size_type edge_idx = neighbordata.size();
        // Locate edge_idx such that
        // edge_adj_[min_id][edge_idx] == tuple(max_id, edge_value)
        for(unsigned int k = 0; k < neighbordata.size(); k++){
            size_type neighbor_id = std::get<0>(neighbordata[k]);
            if(neighbor_id == max_idx){
                edge_idx = k;
            }
        }
        // Return the edge's value (second entry in the tuple at
        // edge_adj_[min_id][edge_idx])
        return(std::get<1>(edge_graph_->edge_adj_[min_idx][edge_idx]));
    }
      
    /** @brief Return this edge's value as a reference that CANNOT be updated or
     *  modified.
     */
    const edge_value_type& value() const {
        size_type min_idx = std::min(node1_idx_, node2_idx_);
        size_type max_idx = std::max(node1_idx_, node2_idx_);
        // Value is stored in the entry of edge_adj_ corresponding to edge
        // representation (min_id, max_id)
          std::vector<std::tuple<size_type,edge_value_type>> neighbordata = edge_graph_->edge_adj_[min_idx];
        size_type edge_idx = neighbordata.size();
        // Locate edge_idx such that
        // edge_adj_[min_id][edge_idx] == tuple(max_id, edge_value)
        for(unsigned int k = 0; k < neighbordata.size(); k++){
            size_type neighbor_id = std::get<0>(neighbordata[k]);
            if(neighbor_id == max_idx){
                edge_idx = k;
            }
        }
        // Return the edge's value (second entry in the tuple at
        // edge_adj_[min_id][edge_idx])
          return(std::get<1>(edge_graph_->edge_adj_[min_idx][edge_idx]));
      }
      
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
   
    // We define a pointer to the Graph to which the edge belongs.
    Graph* edge_graph_ ; 
   
   // Two nodes are asociated t an edge, we identify them with their ID.
    size_type node1_idx_ ; 
    size_type node2_idx_ ; 
   
    // Constructor for `Edge` class.
    Edge(const Graph* edge_graph, size_type node1_idx, size_type node2_idx)
     : edge_graph_(const_cast<Graph*>(edge_graph)), node1_idx_(node1_idx),
       node2_idx_(node2_idx){
       }
  };

  /**
   * @brief Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // Returns the attribute of Graph `num_edges_`.
    return num_edges_;
  }

  /**
   * @brief Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Locate ith edge by index with EdgeIterator
    auto e_itr = edge_begin() + i;
    return *e_itr;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    std::vector<std::tuple<size_type,edge_value_type>> neighbors =
                                                       edge_adj_[a.node_idx_];
    // Iterate through neighboring nodes of a
    for(auto k = neighbors.begin(); k != neighbors.end(); ++k){
      // Check if b is adjacent to a. This is sufficient b/c of the
      // pre-condition (i.e. we assume edges are unordered)
      if (std::get<0>(*k) == b.node_idx_){return true;}
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
   * Use the edge iterator defined in HW1.
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    Edge add_edge(const Node& a, const Node& b,
                  const edge_value_type& value = edge_value_type()) {

      // Iterate through the neighbors of a
      for(auto e_it = a.edge_begin(); e_it != a.edge_end(); ++e_it){
        Edge e = *e_it;
        if(e.node2() == b){
          // This means that the edge already exists, we can return it.
           return Edge(this, a.node_idx_, b.node_idx_);
        }
      }

      // Otherwise (a,b) is not in the graph yet

      // Get ordered edge pair.
      size_type min_id = std::min(a.node_idx_, b.node_idx_);
      size_type max_id = std::max(a.node_idx_, b.node_idx_);

      // We must update _edge_adj_ because two news nodes are neighbors.
      edge_adj_[min_id].push_back(std::make_tuple(max_id,value));
      edge_adj_[max_id].push_back(std::make_tuple(min_id,edge_value_type()));

      num_edges_ ++;

      return Edge(this, a.node_idx_, b.node_idx_);
    }
    
    /** @brief Given two nodes, remove the edge between them from a graph (if the
     * edge exists), returning the number of edges successfully deleted (0 or 1).
     *
     * @param[in,out] a   One node of the edge
     * @param[in,out] b   Other node of the edge
     *
     * @returns         Number of edges deleted (0 or 1)
     *
     * @post            IF has_edge(@a a, @a b) == true,
     *                       old @a num_edges_ == new @a num_edges_ + 1
     *                  ELSE
     *                       old @a num_edges_ == new @a num_edges_
     *
     * Runtime: The for loops each have at most O(d) iterations. Within the
     * for loop, calling @a edge_adj_[@a a.node_id_].erase() happens at most once (from
     * the break statement). Since calling @a edge_adj_[@a a.node_id_].erase() is O(d)
     * the runtime is O(d) (=O(d+d)).
     */
     size_type remove_edge(const Node& a, const Node& b) {
        // First, we need to check if the graph has the edges.
        if(!has_edge(a,b)){
            // Otherwise we return 0.
            return 0;
        }
        // For sanity checks purposes
        size_type adegree_before = a.degree();
        size_type bdegree_before = b.degree();

        // Iterate through the neighbors of node a
        for(auto ni = edge_adj_[a.node_idx_].begin(); ni != edge_adj_[a.node_idx_].end(); ++ni){
          size_type neighbor_idx = std::get<0>(*ni);
          if(neighbor_idx == b.node_idx_){
            // Delete edge data via pointer
            edge_adj_[a.node_idx_].erase(ni);
            // No need to go further.
            break;
          }
        }
        // Do the same for node b.
        for(auto ni = edge_adj_[b.node_idx_].begin(); ni != edge_adj_[b.node_idx_].end(); ++ni){
          size_type neighbor_idx = std::get<0>(*ni);

          if(neighbor_idx == a.node_idx_){
            // Delete edge
            edge_adj_[b.node_idx_].erase(ni);
            break;
          }
        }

        // SANITY CHECKS : Ensure edges were properly deleted from adj_
        assert(a.degree() == adegree_before - 1);
        assert(b.degree() == bdegree_before - 1);
        // We update the number of edges in the graph.
        num_edges_ --;
        // Return that 1 edge was deleted
        return size_type(1);
      }
      
    /** @brief A variant of remove_edge function  where the argument is now an edge.*/
    size_type remove_edge(const Edge& edge){
        Node node_a = edge.node1() ;
        Node node_b = edge.node2() ;
        return remove_edge(node_a, node_b) ;
    }
      
    /** @brief Given an EdgeIterator @ ei, remove the edge pointed by the iterator (if the
     * edge exists), returning the number of edges successfully deleted (0 or 1).
     *
     * @param[in,out] ei   EdgeIterator
     *
     * Runtime : Same as remove_edge function.
     */
    EdgeIterator remove_edge(EdgeIterator e_itr) {
        Edge edge = *e_itr;
        // If this edge is the last edge in the graph
        // We decrement the iterator so it can point to the new last iterable edge.
        if(e_itr.root_node_idx_  == num_nodes() - 1 &&
           e_itr.num_edges_traversed_ == node(e_itr.root_node_idx_).degree() - 1){
          --e_itr;
        }
        // Removing an edge will cause @a ei to point to the next edge in the graph
        remove_edge(edge);

        // Return the iterator
        return e_itr;
      }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
    // We clear the atributes of the clas `Graph`. 
    void clear() {
      edge_adj_.clear();
      index2nodeid_.clear();
      nodes_.clear();
      num_edges_ = 0;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }
    /**
     * @brief Returns the dereferenced value of the iterator using _node_iter_idx_.
     */
    Node operator*() const {
        return graph_->node(node_iter_idx_);
    }
    /**
     * @brief Returns the NodeIterator after one node incrementation.
     */
    NodeIterator& operator++(){
        this->node_iter_idx_++ ;
        return *this ;
    }
      
    /** @brief Return the iterator pointing to the previous node OR the
     * NodeIterator object with @a nodeiterator_idx == 0 i.e. @a node_begin().
     */
    NodeIterator& operator--() {
      // Increment current iteration index
      this->node_iter_idx_--;
      // Dereference & return NodeIterator object
      return *this;
    }
      
    /**
     * @brief Tests the equality of two NodeIterator objects.
     */
    bool operator==(const NodeIterator& nodeiterator) const {
        bool graph = graph_ == nodeiterator.graph_ ;
        bool idx = node_iter_idx_ == nodeiterator.node_iter_idx_ ;
        return graph && idx ;
    }

   private:
    friend class Graph;
   
    // We declare a pointer to a graph object. 
    Graph* graph_ ;
   
    // Node ID.
    size_type node_iter_idx_ ; 
    
    /**
     * @brief Constructor for the _NodeIterator_ class.
     *
     * @param[in] graph                                A pointer to _Graph._
     * @param[in] node_iter_ix                The index of  node.
     *
     */
    NodeIterator(const Graph* graph, const size_type node_iter_idx = 0 )
     : graph_(const_cast<Graph*>(graph)), node_iter_idx_(node_iter_idx) {
    }
  };
    
    /**
     * @brief Method to initialize the first NodeIterator.
     */
    node_iterator node_begin() const {
        return NodeIterator(this,0);
    }
    /**
     * @brief Method to initialize the last NodeIterator.
     */
    node_iterator node_end() const {
        return NodeIterator(this, num_nodes()); // Error correction from HW1. 
    }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
      /**
       * @brief Returns the dereferenced value of the iterator..
       */
      Edge operator*() const{
          size_type a_edge_idx = std::get<0>(graph_->edge_adj_[root_node_idx_][traversed_edges_num_]) ;
          return Edge(graph_, root_node_idx_, a_edge_idx) ;
      }
      
      /**
       * @brief Returns the iterator after one incrementation..
       */
      IncidentIterator& operator++(){
          ++traversed_edges_num_;
          return *this ;
      }
      
      /** @brief operator -- for edge and node removal. */
      IncidentIterator& operator--() {
        --traversed_edges_num_;
        // Return the (dereferenced) IncidentIterator object
        return *this;
      }
      
      /**
       * @brief Tests the equality between _incident_ite_ and this iterator.
       */
      bool operator==(const IncidentIterator& incident_ite) const {
          return (graph_ == incident_ite.graph_) &&
                 (root_node_idx_ == incident_ite.root_node_idx_) &&
                 (traversed_edges_num_ == incident_ite.traversed_edges_num_) ;
      }
      
      /** @brief Increment the IncidentIterator forward by some @a diff.
       * Used for the purpose of HW2.
       * Runtime: O(@a diff)
       */
       IncidentIterator& operator+(const int diff) {
         auto ei_out = this;
         for (int i = 0; i < diff; ++i)
           ++(*ei_out);
         return *ei_out;
       }

   private:
    friend class Graph;
      
      // We define a pointer to _Graph_.
      Graph* graph_ ;
      
      // Index of the root node.
      size_type root_node_idx_ ;
      
      // Number of edges traversed since the beginning of the iterations.
      size_type traversed_edges_num_ ;
      
      /**
       * @brief Constructor for the _IncidentIterator_ class.
       *
       * @param[in] graph                                A pointer to _Graph._
       * @param[in] node_ix                           The index of the root node.
       * @param[in] traversed_edges_num Number of edges traversed (default 0).
       *
       */
      IncidentIterator(const Graph* graph, const size_type node_idx,
                       const size_type traversed_edges_num=0)
       : graph_(const_cast<Graph*>(graph)), root_node_idx_(node_idx),
         traversed_edges_num_(traversed_edges_num) {}
      
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }
      /**
       * @brief Returns the dereferenciated value of the iterator.
       */
      Edge operator*() const {
        size_type node_id = graph_->index2nodeid_[root_node_idx_];
        // Node ID of the current neighbor of the root node
        size_type neigh_id = std::get<0>(graph_->edge_adj_[node_id]
                                                     [num_edges_traversed_]);
        // Return edge defined by (node_id, neigh_id)
        return Edge(graph_, node_id, neigh_id);
      }
      
      /**
       * @brief Returns the  value of the iterator after one increment.
       */
      EdgeIterator& operator++(){
          num_edges_traversed_ ++ ;
          for( ; root_node_idx_ != graph_->num_nodes(); ++root_node_idx_) {
            // For each node, we look at its degree.
            size_type node_degree = graph_->node(root_node_idx_).degree();
            size_type current_node_idx = graph_->index2nodeid_[root_node_idx_];
            // Knowing its degree, we can iterate throuh the edges.
            for ( ; num_edges_traversed_ != node_degree; ++num_edges_traversed_) {
              size_type edge_id = std::get<0>(graph_->edge_adj_[current_node_idx][num_edges_traversed_]);

              
              if (current_node_idx < edge_id)
                // If this condition is not satisfied, it means that
                // the edge has already been visited.
                // Return the (dereferenced) EdgeIterator object.
                return *this;
            }
            // All edges have been visited so we can set the number of edges
            // visited to 0 and move to the next node.
            num_edges_traversed_  = 0;
               }
               num_edges_traversed_  = 0;
               // Return the (dereferenced) EdgeIterator object
               return *this;
             }
          
      bool operator==(const EdgeIterator& edge_ite) const{
          
          return (graph_ == edge_ite.graph_) &&
                 (root_node_idx_ == edge_ite.root_node_idx_) &&
                 (num_edges_traversed_ == edge_ite.num_edges_traversed_) ;
      }
        
    /** @brief Increment the EdgeIterator forward by some @a diff.
     * Used for the purpose of HW2.
     * Runtime: O(@a diff)
     */
    EdgeIterator& operator+(const int diff) {
      auto ei_out = this;
      for (int i = 0; i != diff; ++i)
        ++(*ei_out);
      return *ei_out;
    }

   private:
    friend class Graph;

      // We declare a pointer to _Graph_.
      Graph* graph_ ;
      
      // Root node index.
      size_type root_node_idx_;
      
      // Num of edges traversed.
      size_type num_edges_traversed_ ;
        
      /**
       * @brief Constructor for the _EdgeIterator_ class.
       *
       * @param[in] graph                                A pointer to _Graph._
       * @param[in] root_node_ix                The index of the root node.
       * @param[in] num_edges_traversed Number of edges traversed (default 0).
       *
       */
      EdgeIterator(const Graph* graph, const size_type root_node_idx=0,
                   const size_type num_edges_traversed=0)
       : graph_(const_cast<Graph*>(graph)), root_node_idx_(root_node_idx),
         num_edges_traversed_(num_edges_traversed) {}
      
  };
     /**
      * @brief Method to initialize the first EdgeIterator.
      */
      edge_iterator edge_begin() const{
          return EdgeIterator(this,0,0);
      }
      /**
       * @brief Method to initialize the last EdgeIterator.
       */
      edge_iterator edge_end() const{
          return EdgeIterator(this, num_nodes(), 0);
      }

 private:

 // A struct that enables the storing of the data related to each node.
 struct Intern_Node {
  
     // Node Index
     size_type node_idx_ ; 
  
     // Space position 
     Point position_ ; 
  
     // Node Value 
     node_value_type val_ ; 
 
     // Private Constructor for Intern_Node
     Intern_Node(const size_type node_idx, const Point& position,
                 node_value_type val)
       : node_idx_(node_idx), position_(position), val_(val) {
       }
  } ; 
  
 // A container from STL that stocks every node in the graph.
 std::vector<Intern_Node> nodes_ ;  
 
 // A container for edges adjacency data that now takes edge value into account.
 // Modified for the purpose of HW2 (node value is added to each neighbor).
 std::vector<std::vector<std::tuple<size_type,edge_value_type>>> edge_adj_;

 // The total number of edges in the graph.
 size_type num_edges_ ;
    
// Vector indexed by node_idx_, where index2nodeid_[i] = node_idx_ for Node i
// Used to track Nodes removal.
std::vector<size_type> index2nodeid_;
};

#endif // CME212_GRAPH_HPP
