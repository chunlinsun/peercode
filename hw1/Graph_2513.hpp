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

//--documentation_0
//--docs could benefit from a few pre/post conditions
//--END

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
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
    // ADDED THIS
  using node_value_type = V;

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

        const Point& p = gp_->points_[index_].p_;
      return p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
		assert(index_==gp_->points_[index_].index_);
      return index_;
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
//--functionality_1
//--error when node has degree 0 (not in adj_map_)
//--START
	  size_type degree() const {
		  size_type result = (gp_->adj_map_.at(index_)).size();
		  return result;
	  };
//--END
	  
	  /** Returns an IncidentIterator object for the beginning of the STL
	   *  container which holds the adjacency information for this node
	   */
	  incident_iterator edge_begin() const {
		  typename std::map<int, int>::const_iterator it =
			(gp_->adj_map_.at(index_)).begin();
		  return IncidentIterator(it,*this,gp_);
	  }
	  
	  /** Returns an IncidentIterator object for the end of the STL
	   *  container which holds the adjacency information for this node
	   */
	  incident_iterator edge_end() const {
		  typename std::map<int, int>::const_iterator it =
			(gp_->adj_map_.at(index_)).end();
		  return IncidentIterator(it,*this,gp_);
	  }
	

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        // check if index of this node is same as index
        // of node n and they are part of the same
        // graph

        
      return ((index_ == n.index()) && (gp_ == n.gp()));
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
        // check if index of this node is less than index
        // of node n. Needs to work even if they belong
        // to different graphs.
        if (index_ < n.index())
            return true;
            
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
      size_type index_;
      /** Private Constructor*/
      Node(const Graph* graph, int index)
      : gp_(const_cast<Graph*>(graph)), index_(index) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_points_;
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
//--this method needs to be able to add a node with value equal to the second input.
//--i.e. give the second input a name and assign in.val_ to it.
//--START
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
	// Add the point to the vector of points which
	// exists in the graph, and create an instance
	// of the Node() class to return.
	  
	  // Declare internal_node containing the information for
	  // the Graph class.
	  internal_node in = {
		.p_ = position,
		.index_ = num_points_,
		.gp_ = this
	  };
      
      points_.push_back(in);
      Node n = Node(this, num_points_);
      ++num_points_;
    return n;        // Invalid node
  }
//--END
    
  //Node add_node(const Point&, );

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      // Does the index of this node exist in the Vector,
      // and are they pointing to the same graph?
    return ((n.index() <= num_points_ - 1) && (this==n.gp_));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
      // Construct a new Node pointing to this graph
    return Node(this, i);        // Invalid node
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
        // No code required to create an invalid Edge,
        // Graph::Edge::Edge() will create an empty Edge
        
        // See private variables and functions
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(gp_, node_a_index_);;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(gp_, node_b_index_);      // Invalid Node
    }
      
      /** Return this edge's index, a number in the range [0, num_nodes). */
      size_type index() const {
        return index_;
      }
      
      /** Return this edge's graph pointer. */
      const Graph* gp() const {
        return gp_;
      }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        // Check if index of this edge is equal to index
        // of edge e. Needs to work even if they belong
        // to different graphs.

      return (index_ == e.index() && gp_==e.gp());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        // Check if index of this edge is less than index
        // of edge e. Needs to work even if they belong
        // to different graphs.
		
		if (index_ < e.index()) {
			return true;
		}
		return (std::less<const Graph*>{}(gp_, e.gp_));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
      
      // Pointer back to the Graph class
      Graph* gp_;
 
      // ints to store the index of the two nodes
      // nodes are strictly orederd in the container in the
      // Graph class, but can be either order or equal here.
      size_type node_a_index_;
      size_type node_b_index_;
      
      // int to store the index of the edge
      size_type index_;
      
      // Valid Edge constructor. Private to fulfill requirement that valid nodes
      // can only be construcrted within the Graph class.
      Edge(const Graph* graph, const Node& node_a, const Node& node_b, size_type index)
      : gp_(const_cast<Graph*>(graph)), node_a_index_(node_a.index()), node_b_index_(node_b.index()), index_(index) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
	  return num_edges_;
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
      assert(i<index_edge_map_.size());
	  
	  // get both node indices.
	  size_type ia = index_edge_map_[i].node_idx_1_;
	  size_type ib = index_edge_map_[i].node_idx_2_;
      
      // Construct valid Edge
      Edge e = Edge(this, Node(this, ia), Node(this, ib), i);
          
    return e;        // Invalid Edge
  }
    
    /** Helper function to check if two nodes are distinct and belong to this graph
     */
    void distinct_valid_nodes(const Node& a, const Node& b) const {
        // Check nodes are valid nodes of this graph
        assert(has_node(a));
        assert(has_node(b));
        assert(a.gp()==b.gp());
        
        // Check nodes are distinct
        assert(a.index()!=b.index());
    }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {
      
      // Check nodes are valid distinct nodes of this graph
      distinct_valid_nodes(a, b);
      
      // declare lower and upper node index, l < u
      size_type l = (a.index() < b.index()) ? a.index() : b.index();
      size_type u = (a.index() > b.index()) ? a.index() : b.index();
      
      // map with upper node index as its key
      std::map<int, int> m1 = adj_map_[l];
	  // and lower too
	  std::map<int, int> m2= adj_map_[u];
      
      // if and only if u exists as a key, edge exists
      if (m1.count(u)==1||m2.count(u)==1)
          return true;

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
      distinct_valid_nodes(a, b);
      
      // If the Edge does not exist its information needs to be added
      // to the Graph class containers adj_map_ and index_edge_map_
      
      // If the Edge already exists create a new Edge object to return
      // but don't do anything in the Graph class' containers. Need to
      // use its correct index.
      
      size_type i;
      
      // declare lower and upper index, l < u
      size_type l = (a.index() < b.index()) ? a.index() : b.index();
      size_type r = (a.index() > b.index()) ? a.index() : b.index();

      
      if (!has_edge(a, b)) {
          // update information in Graph class for new edge
          adj_map_[l][r] = num_edges_;
		  adj_map_[r][l] = num_edges_;

		  // Declare internal_edge containing the information for
		  // the Graph class.
		  internal_edge ie = {
			.node_idx_1_=l,
			.node_idx_2_=r,
			.index_=(num_edges_),
			.gp_=this
		  };
		  
		  // add to container in Graph class
          index_edge_map_.push_back(ie);
		  
          // set index of Edge we are returning
          i = num_edges_;
		  
		  num_edges_++;
      }
      
      // if it already exists use its index i which can be found
      // from the adjancy nested map. Also i==adj_map_[r][l] valid.
      i = adj_map_[l][r];
      
      // Construct edge e
      Edge e = Edge(this, a, b, i);
  
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      
      // Set containers to be empty - are objects inside released from
      // memory?
      points_.clear();
      adj_map_.clear();
      index_edge_map_.clear();
      num_points_=0;
      num_edges_=0;
      
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
	  // Utilize the iterator induced on the STL container for nodes
	  typename std::vector<internal_node>::const_iterator it_;
      
	  /** Dereference the NodeIterator
	   *
	   * @return the Node object we are currently iterating over
	   */
      Node operator*() const {
		  internal_node in = *it_;
		  return Node(in.gp_,in.index_);
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
	  // Private constructor that can be accesed by the Graph class (not Node)
	  NodeIterator(typename std::vector<internal_node>::const_iterator it)
	  : it_(it) {}
  };
	
	/** Returns a NodeIterator object for the beginning of the STL
	 *  container which holds the information of the node, edge mapping
	 */
	node_iterator node_begin() const {
		typename std::vector<internal_node>::const_iterator it =
			points_.begin();
		NodeIterator node_iter = NodeIterator(it);
		return node_iter;
	}
	
	/** Returns a NodeIterator object for the end of the STL
	 *  container which holds the information of the node, edge mapping
	 */
	node_iterator node_end() const {
		typename std::vector<internal_node>::const_iterator it =
			points_.end();
		NodeIterator node_iter = NodeIterator(it);
		return node_iter;
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
		  size_type edge_index = it_->second;
		  size_type node_b_index = it_->first;
		  Node node_b = Node(gp_, node_b_index);
		  Edge e = Edge(gp_, node_a_, node_b, edge_index);
		  return e;
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
		  internal_edge ie = *it_;
		  Edge e = (*(ie.gp_)).edge(ie.index_);
		  return e;
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
    // HW1 #5: YOUR CODE HERE
	  
	  // Utilize the iterator induced on the STL container for internal_edges
	  // This contains each edge exactly once; we want to iterate over
	  // each edge exactly once
	  typename std::vector<internal_edge>::const_iterator it_;
	  
	  /** Private Constructor
	   *  @param[in] it an interator induced from this node's container of internal_node s */
	  EdgeIterator(typename std::vector<internal_edge>::const_iterator it)
	  : it_(it) {}
	
  };
	
	/** Returns a EdgeIterator object for the beginning of the STL
	 *  container which holds the information of the node, edge mapping
	 */
	edge_iterator edge_begin() const {
		return EdgeIterator(index_edge_map_.begin());
	}
	
	/** Returns a EdgeIterator object for the end of the STL
	 *  container which holds the information of the node, edge mapping
	 */
	edge_iterator edge_end() const {
		return EdgeIterator(index_edge_map_.end());
	}

 private:
	
	// struct containing any information we need to keep about the nodes
	struct internal_node {
		Point p_;
		size_type index_;
		Graph* gp_;
		V val_;
	};
	
	// struct containing any information we need to keep about the edges
	// @node_idx_1_ < @node_idx_2_
	struct internal_edge {
		size_type node_idx_1_;
		size_type node_idx_2_;
		size_type index_;
		Graph* gp_;
	};

	 // Container which holds the information of the Nodes
	 typename std::vector<internal_node> points_;
	 
	 // Counter to keep track of the size of the vector of Nodes - optional
	 // size_type not yet declared so used unsigned.
	 size_type num_points_ = 0;
	 
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
	 std::map<int, std::map<int, int>> adj_map_;
	 
	 // Vector holding <Index, Edge> mapping via the internal_edge Struct
	 std::vector<internal_edge> index_edge_map_;

	 // Counter to keep track of the size of the vector of Edges - optional
	 size_type num_edges_ = 0;

};

#endif // CME212_GRAPH_HPP
