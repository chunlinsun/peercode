#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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
    // The Graph class contains all the information which the proxy classes
    // Node() and Edge() can map to, and therefore be lightweight.
	
	// Predeclare the internal node struct
	struct internal_node;
	
	// Predeclare the internal edge struct
	struct internal_edge;
 
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  /** Type of nodes. */
  using node_value_type = V;
  /** Type of edges. */
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
	
  /** Type of NodeIteratior*/
  using node_iterator_type = typename std::vector<size_type>::const_iterator;
	
  /** Type of EdgeIteratior*/
  using edge_iterator_type = typename std::vector<size_type>::const_iterator;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
      // Nothing required here as the private variables are either
      // int, map, or vector and are correctly initialized by their
      // default constructors.

      // Make sure Node and Edge classes meet size restrictions
      assert(sizeof Graph::Node() <= 16);
      assert(sizeof Graph::Edge() <= 32);
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
        // Nothing needs to go in here because we don't need to
        // do anything with a Node defined in this way. Calling
        // Graph::Node::Node() will just create an empty node with
        // no information, an invalid Node.
        
        // See the private variables and constructor.
    }

      
    /** Return this node's position. */
    const Point& position() const {
        // Access the parent graph by deferencing the pointer,
        // use the vector of points and the index to get the
        // position
		assert(valid());
		const Point& p = gp_->points_[index_].p_;
		return p;
    }
	  
	/** Return this node's position. Pass by reference so it can be modidfied
	 *  @pre this is a valid node of this Graph.*/
	  Point& position () {
		  assert(valid());
		  Point& p = gp_->points_[index_].p_;
		  return p;
	  }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
		size_type result = gp_->points_[index_].index_;
		assert(result < gp_->num_active_points_);
      return result;
    }
      /** Return this node's graph pointer. */
      const Graph* gp() const {
        return gp_;
      }
 
	  /** Return this node's value of type node_value_type. Pass by reference
	   *  so it can be set by node.value() = val
	   */
	  node_value_type& value() {
		  return gp_->points_.at(index_).val_;
	  }
	  
	  /** Return this node's value of type node_value_type. Read only. */
	  const node_value_type& value() const {
		  return gp_->points_.at(index_).val_;
	  }
	  
	  /** Return the number of nodes this node is connected to via
	   *  a valid edge
	   */
	  size_type degree() const {
		  size_type result = (gp_->adj_map_.at(index_)).size();
		  return result;
	  };
	  
	  /** Returns an IncidentIterator object for the beginning of the STL
	   *  container which holds the adjacency information for this node
	   */
	  incident_iterator edge_begin() const {
		  assert(valid());
		  typename std::map<int, int>::const_iterator it =
			(gp_->adj_map_.at(index_)).begin();
		  return IncidentIterator(it,*this,gp_);
	  }
	  
	  /** Returns an IncidentIterator object for the end of the STL
	   *  container which holds the adjacency information for this node
	   */
	  incident_iterator edge_end() const {
		  assert(valid());
		  typename std::map<int, int>::const_iterator it =
			(gp_->adj_map_.at(index_)).end();
		  return IncidentIterator(it,*this,gp_);
	  }
	

	/** Test whether this node and @a n are equal.
	 * @pre n is a valid node, and so is this one.
	 *
	 * Equal nodes have the same graph and the same index.
	 */
    bool operator==(const Node& n) const {
        // check if index of this node is same as index
        // of node n and they are part of the same
        // graph
		assert(valid());
		return (gp_->points_[index_].index_ == n.index() && (gp_ == n.gp()));
    }

	/** Test whether this node is less than @a n in a global order.
	 * @pre n is a valid node, and so is this one.
	 *
	 * This ordering function is useful for STL containers such as
	 * std::map<>. It need not have any geometric meaning.
	 *
	 * The node ordering relation must obey trichotomy: For any two nodes x
	 * and y, exactly one of x == y, x < y, and y < x is true.
	 */
    bool operator<(const Node& n) const {
        // check if index of this node is less than index
        // of node n. Needs to work even if they belong
        // to different graphs.
		 assert(valid());
        if (gp_->points_[index_].index_< n.index())
            return true;
		else
			return (std::less<const Graph*>{}(gp_, n.gp_));
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Define constructor which takes in the index of a Node. Must exist
    // here as we need the graph it has to point to in order to add the point,
    // and to fulfill requirement that a valid node can only be created within
    // the Graph class.
      
      // Pointer back to the parent Graph
      Graph* gp_;
      // int to store the index of the node in the graph's vector of Points
      size_type index_; // UID
      /** Private Constructor. Uses UID as index. */
      Node(const Graph* graph, int index)
      : gp_(const_cast<Graph*>(graph)), index_(index) {}
	  
	  
	  /** Function to test for invariants for Nodes
	   * @brief Tests if this node has active index in the correct range, UID in the correct range,
	   * and that the node is in sync for the containers @a points_ and @a i2u_nodes_
	   * @return true if all invariants are passed, false if not, along with a printout of useful
	   * diagnostic information for edbugging.
	   */
	  bool valid() const {
		  if (!(index_ >=0 && index_ < gp_->points_.size())) {
			  std::cout<<"UID out of range"<<std::endl;
		  }
		  else if (!(gp_->points_[index_].index_ < gp_->i2u_nodes_.size())) {
			  std::cout<<"Active out of range"<<std::endl;
		  }
					else if (!(gp_->i2u_nodes_[(gp_->points_)[index_].index_] == index_)) {
			  std::cout<<"Node UID "<<index_<<" out of sync; UID via points_ = "<<(gp_->i2u_nodes_)[(gp_->points_)[index_].index_]<<std::endl;
						std::cout<<"Node active ID "<<(gp_->points_)[index_].index_<<std::endl;
		  }
		  return (index_ >=0 && index_ < gp_->points_.size())
		  && (gp_->points_[index_].index_ < gp_->i2u_nodes_.size())
		  && (gp_->i2u_nodes_[(gp_->points_)[index_].index_] == index_);
	  }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_active_points_;
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
	// Add the point to the vector of points which
	// exists in the graph, and create an instance
	// of the Node() class to return.
	  
	  // Declare internal_node containing the information for
	  // the Graph class.
	  internal_node in = {
		.p_ = position,
		.index_ = num_active_points_, // next active ID
		.gp_ = this
	  };
      points_.push_back(in);
	  
	  // Add to the set of currently active nodes
	  i2u_nodes_.push_back(num_points_); // i2u_nodes_[num_active_points_] == UID
	  
	  // Need to add an empty map to the adjacency map for this
	  // node in case we try to invoke incident iterator on a
	  // Node with no edges
	  std::map<int, int> empty_map {};
	  adj_map_[num_points_] = empty_map;
	  assert(adj_map_[num_points_].empty());
	  
      Node n = Node(this, num_points_);
      ++num_points_;
	  ++num_active_points_;
    return n;
  }

	Node add_node(const Point& position, const node_value_type& val) {
	// Add the point to the vector of points which
	// exists in the graph, and create an instance
	// of the Node() class to return.
	  
	  // Declare internal_node containing the information for
	  // the Graph class.
	  internal_node in = {
		.p_ = position,
		.index_ = num_active_points_, // next active ID
		.gp_ = this,
		.val_ = val
	  };
      points_.push_back(in);
	  
	  // Add to the set of currently active nodes
	  i2u_nodes_.push_back(num_points_); // i2u_nodes_[num_active_points_] == UID
	  
	  // Need to add an empty map to the adjacency map for this
	  // node in case we try to invoke incident iterator on a
	  // Node with no edges
	  std::map<int, int> empty_map {};
	  adj_map_[num_points_] = empty_map;
	  assert(adj_map_[num_points_].empty());
		
	  Node n = Node(this, num_points_);
      ++num_points_;
	  ++num_active_points_;
    return n;
	}
	

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      // Is this node's active index in range,
      // and are they pointing to the same graph?
    return (n.index() >= 0) && (n.index() < num_active_points_) && (this==n.gp());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
      // Construct a new Node pointing to this graph
	  assert(points_[i2u_nodes_[i]].index_==i);
    return Node(this, i2u_nodes_[i]);
  }
	
	/** Remove this node from this graph
	 * @param[in] @a n node to remove
	 * @brief Remove this node from the graph, along with all of its incident
	 * edges. Complexity is O(num_edges)
	 * @return false if this is not a valid node of this graph, otherwise the
	 * active index of this node
	 * @post this node, and all its incident edges, are no longer part of this graph
	 * @post graph.num_nodes() has decreased by 1
	 * @post graph.num_edges() has decreased by the number of edges
	 * incident to this node
	 * @post Active index for the node with the largest active index becomes
	 * this node's acive index
	 * @post NodeIterator objects are invalidated
	 */
  size_type remove_node(const Node& n){
	  // check the node we are trying to remove is in the graph
	  if (!has_node(n)) {
		  return 0;
	  }
	  
	  // active index of the node being removed
	  size_type result = n.index();
	  
	  // Remove all edges incident to this node via incident iterator
	  while (n.edge_begin()!=n.edge_end()) {
		  remove_edge(*n.edge_begin());
	  }
	  
	  // Also update active index for edges incident to the node
	  // you are swapping in
	  const size_type swap_uid = i2u_nodes_.back();
	  Node n_last = Node(this, swap_uid);
	
	  for (auto it = n_last.edge_begin(); it!=n_last.edge_end(); ++it){
		  Edge e = *it;
		  size_type e_UID = i2u_edges_[e.index()];
		  // Change active index for the node which is not n
		  if (index_edge_map_[e_UID].node_idx_1_==num_active_points_-1) {
			  index_edge_map_[e_UID].node_idx_1_=result;
		  }
		  else {
			  index_edge_map_[e_UID].node_idx_2_=result;
		  }
	  }
	  
	  // Remove the now empty map for this edgeless node
	  size_type n_UID = i2u_nodes_[n.index()];
	  assert(adj_map_[n_UID].empty());
	  adj_map_.erase(n_UID);
	  
	  // Change active index for the node being swapped in
	  points_[swap_uid].index_ = result;
	  
	  // Remove node from containter of active nodes. The node with the
	  // highest index (==num_active_points_-1) will now have this node's
	  // index
	  std::swap(i2u_nodes_[result], i2u_nodes_.back());
	  i2u_nodes_.pop_back();
	  
	  // Update the node moved into position @result has correct
	  // active / external index.
	  --num_active_points_;
	  assert(num_active_points_==i2u_nodes_.size());
	  return 1;
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
        // No code required to create an invalid Edge,
        // Graph::Edge::Edge() will create an empty Edge
        
        // See private variables and functions
    }

    /** Return a node of this Edge */
    Node node1() const {
		assert(valid());
      return Node(gp_, node_a_index_);;      // node_a_index_ is the UID
    }

    /** Return the other node of this Edge */
    Node node2() const {
		assert(valid());
      return Node(gp_, node_b_index_);      // node_b_index_ is the UID
    }
      
      /** Return this edge's index, a number in the range [0, num_nodes). */
      size_type index() const {
        return gp_->index_edge_map_[index_].index_;
      }
      
      /** Return this edge's graph pointer. */
      const Graph* gp() const {
		  assert(valid());
        return gp_;
      }
	  
	  /** Return this edge's value of type edge_value_type. Pass by reference
	   *  so it can be set by edge.value() = val
	   */
	  edge_value_type& value() {
		  assert(valid());
		  return gp_->index_edge_map_.at(index_).e_val_;
	  }
	  
	  /** Return this edge's value of type edge_value_type. Read only. */
	  const edge_value_type& value() const {
		  assert(valid());
		  return gp_->index_edge_map_.at(index_).e_val_;
	  }

	/** Test whether this edge and @a e are equal.
	 * @pre @a e is a valid edge, and so is this one.
	 *
	 * Equal edges represent the same undirected edge between two nodes.
	 */
    bool operator==(const Edge& e) const {
        // Check if index of this edge is equal to index
        // of edge e. Needs to work even if they belong
        // to different graphs.
		assert(valid());
		return ((gp_->index_edge_map_)[index_].index_ == e.index() && gp_==e.gp());
    }

	/** Test whether this edge is less than @a e in a global order.
	 * @pre @a e is a valid edge, and so is this one.
	 *
	 * This ordering function is useful for STL containers such as
	 * std::map<>. It need not have any interpretive meaning.
	 */
    bool operator<(const Edge& e) const {
        // Check if index of this edge is less than index
        // of edge e. Needs to work even if they belong
        // to different graphs.
		assert(valid());
		if ((gp_->index_edge_map_)[index_].index_ < e.index()) {
			return true;
		}
		return (std::less<const Graph*>{}(gp_, e.gp_));
    }
	  
	  /** Return the L2 distance between this Edge's two nodes */
	  double length() const {
		  assert(valid());
		  Point diff = node1().position() - node2().position();
		  return norm_2(diff);
	  }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
      
      // Pointer back to the Graph class
      Graph* gp_;
 
      // ints to store the index of the two nodes
      // nodes are strictly orederd in the container in the
      // Graph class, but can be either order or equal here.
      size_type node_a_index_; //UID
      size_type node_b_index_; //UID
      
      // int to store the index of the edge
	  size_type index_; //UID
      
      // Valid Edge constructor. Private to fulfill requirement that valid nodes
      // can only be construcrted within the Graph class.
      Edge(const Graph* graph, const Node& a, const Node& b, size_type index)
      : gp_(const_cast<Graph*>(graph)), index_(index) {
		  // No self edges
		  assert(a.index()!=b.index());
		  
		  size_type a_UID = gp_->i2u_nodes_[a.index()];
		  size_type b_UID = gp_->i2u_nodes_[b.index()];
		  //declare lower and upper index, l < u
		  size_type lower = (a_UID < b_UID) ? a_UID : b_UID;
		  size_type upper = (a_UID> b_UID) ? a_UID : b_UID;
		  
		   // Ensure node_a_index_ < node_b_index_
		  node_a_index_ = lower;
		  node_b_index_ = upper;
	  }
	  /** Function to test for invariants for Edges
	   * @brief Tests if this Edge has active index in the correct range, UID in the correct range,
	   * the edge is in sync for the containers @a index_edge_map_ and @a i2u_edges_, and
	   * the two nodes this edge connects are in sync in the containers @a index_edge_map_
	   * and @a i2u_nodes_
	   * @return true if all invariants are passed, false if not, along with a printout of useful
	   * diagnostic information for edbugging.
	   */
	  bool valid() const {
		  bool result =
		  // Edge UID is in range
		  index_ >=0 && index_ < gp_->index_edge_map_.size()
		  // Edge Active ID is in range
		  && gp_->index_edge_map_[index_].index_ < gp_->i2u_edges_.size()
		  // node UIDs are in correct order
		  && node_a_index_ < node_b_index_
		  // Edge UID in sync
		  && gp_->i2u_edges_[(gp_->index_edge_map_)[index_].index_] == index_
		   // UID for node a of this edge in sync
		  && (gp_->i2u_nodes_[(gp_->index_edge_map_)[index_].node_idx_1_]==
				node_a_index_ || gp_->i2u_nodes_[(gp_->index_edge_map_)[index_].node_idx_2_]==
				node_a_index_)
		  // UID for node b of this edge in sync
		  && (gp_->i2u_nodes_[(gp_->index_edge_map_)[index_].node_idx_1_]==
				node_b_index_ || gp_->i2u_nodes_[(gp_->index_edge_map_)[index_].node_idx_2_]==
				node_b_index_);
		  if (!result) {
			  auto in = (gp_->index_edge_map_)[index_];
			  std::cout<<"Validation failed at Edge UID: "<<index_<<std::endl;
			  std::cout<<"Validation failed at Edge active ID: "<<in.index_<<std::endl;
			  			  std::cout<<"Validation failed at Node1 active ID: "<<in.node_idx_1_<<std::endl;
			  			  std::cout<<"Validation failed at Node2 active ID: "<<in.node_idx_2_<<std::endl;
		  }
		  return result;
	  }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
	  return num_active_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Use index_edge_map_ to get the indices of the nodes,
    // then return a new instance of that edge.
      
      // check i is a valid edge index
      assert(i<i2u_edges_.size());
	  
	  size_type e_UID = i2u_edges_[i];
	  
	  // get both node UIDs.
	  size_type ia = i2u_nodes_[index_edge_map_[e_UID].node_idx_1_];
	  size_type ib = i2u_nodes_[index_edge_map_[e_UID].node_idx_2_];
      
      // Construct valid Edge. ia and ib do not need to be ordered
      Edge e = Edge(this, Node(this, ia), Node(this, ib), i);
          
    return e;
  }
    
    /** Helper function to check if two nodes belong to this graph
     */
    void valid_nodes(const Node& a, const Node& b) const {
        // Check nodes are valid nodes of this graph
        assert(has_node(a));
        assert(has_node(b));
        assert(a.gp()==b.gp());
    }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {
      
      // Check nodes are valid distinct nodes of this graph
      valid_nodes(a, b);
	  
	  // Self edges are disallowed
	  if (a.index()==b.index()){
		  return false;
	  }

	  // Get UID for each node
	  size_type a_UID = i2u_nodes_[a.index()];
	  size_type b_UID = i2u_nodes_[b.index()];

      // map with node a UID as its key
      std::map<int, int> ma = adj_map_[a_UID];
	  // map with node b UID as its key
	  std::map<int, int> mb = adj_map_[b_UID];
	  assert((ma.count(b_UID)>0)==(mb.count(a_UID)>0));
      // if and only if map[UID] is not empty, edge exists
      if (ma.count(b_UID)>0||mb.count(b_UID)>0)
          return true;
	  else
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
      
      // Check nodes are valid distinct nodes of this graph
	  assert(a.index()!=b.index());
      valid_nodes(a, b);
      
      // If the Edge does not exist its information needs to be added
      // to the Graph class containers adj_map_ and index_edge_map_
      
      // If the Edge already exists create a new Edge object to return
      // but don't do anything in the Graph class' containers. Need to
      // use its correct index.
      
	  size_type i;
	  
	  // Get UID of the 2 nodes
	  size_type a_uid = i2u_nodes_[a.index()];
	  size_type b_uid = i2u_nodes_[b.index()];
	  
	  if (has_edge(a, b)) {
		  // if it already exists use its index i which can be found
		  // from the adjancy nested map. Also i==adj_map_[b_uid][a_uid] valid.
		  assert(adj_map_[a_uid][b_uid]==adj_map_[b_uid][a_uid]);
		  i = adj_map_[a_uid][b_uid];
	  }

      else {
          // update information in Graph class for new edge
          adj_map_[a_uid][b_uid] = num_edges_;
		  adj_map_[b_uid][a_uid] = num_edges_;

		  // Declare internal_edge containing the information for
		  // the Graph class.
		  internal_edge ie = {
			.node_idx_1_=a.index(), // active index
			.node_idx_2_=b.index(), // active index
			.index_=num_active_edges_, // active index
			.gp_=this
		  };
		  
		  // add to container of all edges which have ever been
		  // part of this Graph
          index_edge_map_.push_back(ie);
		  
		  // Add to the set of currently active edges
		  // i2u_nodes_[num_active_points_] == UID
		  i2u_edges_.push_back(num_edges_);
		  
          // set UID of Edge we are returning
          i = num_edges_;
		  
		  num_edges_++;
		  num_active_edges_++;
		  
		  assert(num_edges_==index_edge_map_.size());
		  assert(num_active_edges_==i2u_edges_.size());
      }
  
    return Edge(this, a, b, i);
  }
	
	/** Remove edge connecting these two nodes from this graph
	 * @param[in] @a n1 node 1 of edge to remove
	 * @param[in] @a n2 node 2 of edge to remove
	 * @brief Remove this edge from the graph. Complexity is
	 * O(num_edges)
	 * @return false if this is not a valid edge of this graph, otherwise the
	 * active index of this edge
	 * @post this edge is no longer part of this graph
	 * @post graph.num_edges() has decreased by 1
	 * @post Active index for the edge with the largest active index becomes
	 * this edge's acive index
	 * @post EdgeIterator objects are invalidated
	 */
size_type remove_edge(const Node& n1, const Node& n2){
	
	if (!has_edge(n1, n2)) {
		return 0;
	}
	
	// Get UID for the edge and each node
	size_type n1_UID = i2u_nodes_[n1.index()];
	size_type n2_UID = i2u_nodes_[n2.index()];
	size_type e_UID = adj_map_[n1_UID][n2_UID];
	
	// active index of the edge being removed
	size_type result = index_edge_map_[e_UID].index_;

	// Remove from the adjacency map
	adj_map_[n1_UID].erase(n2_UID);
	adj_map_[n2_UID].erase(n1_UID);
	
	// Remove node from containter of active edges. The edge with the
	// highest index (==num_active_edges_-1) will now have this node's
	// index
	
	//size_type swap_uid = i2u_edges_[(num_active_edges_ -1)];
	size_type swap_uid = i2u_edges_.back();
	std::swap(i2u_edges_[result], i2u_edges_.back());
	i2u_edges_.pop_back();
	
	// Update the edge moved into position @result has correct
	// active / external index. No need to change the active
	index_edge_map_[swap_uid].index_ = result;
	
	--num_active_edges_;
	assert(num_active_edges_==i2u_edges_.size());
	
	return 1;
	}
	
	/** Remove this edge from this graph
	 * @param[in] @a e edge to remove
	 * @brief Remove this edge from the graph. Complexity is
	 * O(num_edges)
	 * @return false if this is not a valid edge of this graph, otherwise the
	 * active index of this edge
	 * @post this edge is no longer part of this graph
	 * @post graph.num_edges() has decreased by 1
	 * @post Active index for the edge with the largest active index becomes
	 * this edge's acive index
	 * @post EdgeIterator and IncidentIterator objects are invalidated
	 */
size_type remove_edge(const Edge& e) {

	size_type e_UID = i2u_edges_[e.index()];
	size_type n1_UID = i2u_nodes_[index_edge_map_[e_UID].node_idx_1_];
	size_type n2_UID = i2u_nodes_[index_edge_map_[e_UID].node_idx_2_];
	
	return remove_edge(Node(this, n1_UID), Node(this, n2_UID));
	}


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      
      // Set containers to be empty
      points_.clear();
	  i2u_nodes_.clear();
	  
      adj_map_.clear();
      index_edge_map_.clear();
	  i2u_edges_.clear();
	  
      num_points_=0;
      num_edges_=0;
	  
	  num_active_points_=0;
	  num_active_edges_=0;
      
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

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }
	  /** Dereference the NodeIterator
	   *
	   * @return the Node object we are currently iterating over
	   */
      Node operator*() const {
		  return Node(gp_, *it_);
      }
	  
	  /** Increment NodeIterator */
	  NodeIterator& operator++() {
		  ++it_;
		  return *this;
	  }
	  
	  /** @return true if NodeIterator is the same node as this one */
	  bool operator==(const NodeIterator& node_iter) const {
		  return (it_==node_iter.it_);
	  }
	  
	  /** @return false if a NodeIterator is the same node as this one */
	  bool operator!=(const NodeIterator& node_iter) const {
		  return (it_!=node_iter.it_);
	  }

   private:
    friend class Graph;
	  
	  // Utilize the iterator induced on the STL container for nodes
	  // typename std::vector<internal_node>::const_iterator it_;
		node_iterator_type it_;
	  
	  // Keep track of which graph it belongs to TODO CONST?
		const Graph* gp_;
	  
	  // Private constructor that can be accesed by the Graph class (not Node)
	  // Iterate over vector of index->UID @i2u_nodes_
	  NodeIterator(node_iterator_type it,
				   const Graph* gp)
	  : it_(it), gp_(const_cast<Graph*>(gp)) {}
	  
  };
	
	/** Returns a NodeIterator object for the beginning of the STL
	 *  container which holds the information of the node, edge mapping
	 */
	node_iterator node_begin() const {
		return NodeIterator(i2u_nodes_.begin(), this);
	}
	
	/** Returns a NodeIterator object for the end of the STL
	 *  container which holds the information of the node, edge mapping
	 */
	node_iterator node_end() const {
		return NodeIterator(i2u_nodes_.end(), this);
	}
	
	/** Remove node and return node iterator
	 * @return A node iterator pointing to the same position in
	 * @a i2u_nodes_ as was input.
	 * @post The node whcih @a n_it was pointing to has been removed
	 */
	node_iterator remove_node(node_iterator n_it) {
		remove_node(*n_it);
		return n_it;
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

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }
	  
	  /** Dereference the IncidentIterator
	   *
	   * @return the Edge object we are currently iterating over
	   *  out of  the Edges incident to this Node.
	   * @post e.node1() is the Node this IncidentIterator belongs to
	   */
	  Edge operator*() const {
		  size_type e_UID = it_->second;
		  size_type node_b_UID = it_->first;
		  Node node_b = Node(gp_, node_b_UID);
		  return Edge(gp_, node_a_, node_b, e_UID);
	  }
	  
	  /** Increment IncidentIterator */
	  IncidentIterator& operator++() {
		  ++it_;
		  return *this;
	  }
	  
	  /** Return true if a IncidentIterator is the same node as this one */
	  bool operator==(const IncidentIterator& incident_iter) const {
		  return (it_==incident_iter.it_);
	  }
	  
	  /** Return false if a IncidentIterator is the same node as this one */
	  bool operator!=(const IncidentIterator& incident_iter) const {
		  return (it_!=incident_iter.it_);
	  }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
	  // Utilize the iterator induced on the STL container for node,
	  // edge mapping
	  typename std::map<int, int>::const_iterator it_;
	  
	  // Also need to keep track of which node this is and which
	  // Graph it belongs to
	  const Node node_a_;
	  const Graph* gp_;
	  
	  /** Private Constructor
	   *  @param[in] it an interator induced from this node's map of node index, edge index
	   *  @param[in] node_a the Node over which's Edge s we are iterating
	   *  @param[in] gp the Graph these Nodes and Edge belong to */
	  IncidentIterator(typename std::map<int, int>::const_iterator it,
					   const Node node_a, const Graph* gp)
	  : it_(it), node_a_(node_a), gp_(const_cast<Graph*>(gp)) {}
	  
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

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }
	  /** Dereference the EdgeIterator
	   *
	   * @return the Edge object we are currently iterating over
	   */
	  Edge operator*() const {
		  size_type	e_UID = *it_;
		  
		  Node a = Node(gp_, gp_->i2u_nodes_[gp_->index_edge_map_[e_UID].node_idx_1_]);
		  Node b = Node(gp_, gp_->i2u_nodes_[gp_->index_edge_map_[e_UID].node_idx_2_]);
		  
		  return Edge(gp_, a, b, *it_);
	  }
	  
	  /** IncrementEdgeIterator */
	  EdgeIterator& operator++() {
		  ++it_;
		  return *this;
	  }
	  
	  /** @return true if NodeIterator is the same node as this one */
	  bool operator==(const EdgeIterator& edge_iter) const {
		  return (it_==edge_iter.it_);
	  }
	  
	  /** @return false if NodeIterator is the same node as this one */
	  bool operator!=(const EdgeIterator& edge_iter) const {
		  return (it_!=edge_iter.it_);
	  }
	

   private:
    friend class Graph;
	  // Utilize the iterator induced on the STL container for internal_edges
	  // This contains each edge exactly once; we want to iterate over
	  // each edge exactly once
	  edge_iterator_type it_;
	  
	  // Keep track of which graph it belongs to TODO CONST?
	  const Graph* gp_;
	  
	  /** Private Constructor
	   *  @param[in] it an interator induced from this node's container of internal_node s */
	  EdgeIterator(edge_iterator_type it, const Graph* gp)
	  : it_(it), gp_(const_cast<Graph*>(gp)) {}
  };
	
	/** Returns a EdgeIterator object for the beginning of the STL
	 *  container which holds the information of the node, edge mapping
	 */
	edge_iterator edge_begin() const {
		return EdgeIterator(i2u_edges_.begin(), this);
	}
	
	/** Returns a EdgeIterator object for the end of the STL
	 *  container which holds the information of the node, edge mapping
	 */
	edge_iterator edge_end() const {
		return EdgeIterator(i2u_edges_.end(), this);
	}
	
	/** Remove edge and return edge iterator
	 * @return A node iterator pointing to the same position in
	 * @a i2u_edges_ as was input.
	 * @post The edge whcih @a e_it was pointing to has been removed
	 */
	edge_iterator remove_edge(edge_iterator e_it) {
		remove_edge(*e_it);
		return e_it;
	}

 private:
	
	// struct containing any information we need to keep about the nodes
	struct internal_node {
		Point p_;
		size_type index_; // active index
		Graph* gp_;
		node_value_type val_;
	};
	
	// struct containing any information we need to keep about the edges
	// @node_idx_1_ < @node_idx_2_ when added
	struct internal_edge {
		size_type node_idx_1_; // active index
		size_type node_idx_2_; // active index
		size_type index_; // active index
		Graph* gp_;
		edge_value_type e_val_;
	};
	
	 // Container which holds the information of the Nodes
	 typename std::vector<internal_node> points_;
	
	 // Store the currently "active" set of nodes.
	 std::vector<size_type> i2u_nodes_;   // Indexed by node idx
	 
	 // Counter to keep track of the size of the vector of Nodes which
	 // have been part of the graph at any point.
	 // num_points_ = max{UID_nodes} + 1 == points_.size()
	 size_type num_points_ = 0;
	
	 // Counter to keep track of number of active nodes
	 // num_active_points_ == i2u_nodes_.size()
	 size_type num_active_points_ = 0;
	 
	 // For Edge, it contains an adjacency map @adj_map_ which holds the
	 // information for which pairs of Nodes have an edge. Both orderings
	 // of nodes exist. This structure can
	 // thus both check with pairs of nodes have edges, and also return the
	 // edge index for a given pair of nodes.
	
	 // A container for easy mapping between edge index and edge itself is
	 // @index_edge_map_. This is a vector of containing internal_edge
	 // struct s.
	
	 // We also keep a counter of number of edges @num_edges_.
	 
	 
	 // Adjacency map
	 // map<node1_idx, map<node2_idx, edge_inx>>
	 // Each index in this container is a UID, is added when
	 // an edge is added, and never removed even if an edge is
	 // removed from the graph
	 std::map<int, std::map<int, int>> adj_map_;
	 
	 // Vector holding <Index, Edge> mapping via the internal_edge Struct
	// Each index in this container is a UID
	 std::vector<internal_edge> index_edge_map_;
	
	// Store the currently "active" set of edges.
	std::vector<size_type> i2u_edges_;   // Indexed by active edge idx

	 // Counter to keep track of the size of the vector of Edges which
	 // have been part of the graph at any point.
	 // num_edges_ = max{UID_edges} + 1 == index_edge_map_.size()
	 size_type num_edges_ = 0;
	
	// Counter to keep track of number of active edges
	// num_active_edges_ == i2u_edges_.size()
	size_type num_active_edges_ = 0;
	
};

#endif // CME212_GRAPH_HPP
