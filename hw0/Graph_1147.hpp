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
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
    
    // The Graph class contains all the information which the proxy classes
    // Node() and Edge() can map to, and therefore be lightweight.
    
    // For Node, it contains a vector of Point() @points_ with index
    // of the vector equal to index of the Node() and a counter num_points_
    
    
    // Container which holds the information of the Nodes
    std::vector<Point> points_;
    
    // Counter to keep track of the size of the vector of Nodes - optional
    // size_type not yet declared so used unsigned.
    unsigned num_points_;
    
    // For Edge, it contains an adjacency map @ adj_map_ which holds the
    // information for which pairs of Nodes have an edge. The Nodes are
    // ordered such that the the index of the map is the node with lower
    // index, and the value is itself a map with key the node with higher
    // value, and value being the index of that edge. This structure can
    // thus both check with pairs of nodes have edges, and also return the
    // edge index for a given pair of nodes.
    // A container for easy mapping between edge index and edge itself is
    // @index_edge_map_. This is a vector, with index equal to edge index,
    // containing tuples of ordered pairs of int's which are the indices
    // of the respective nodes.
    // We also keep a counter of number of edges @num_edges_.
    
    
    // Adjacency map
    // map<node1_idx, map<node2_idx, edge_inx>>
    std::map<int, std::map<int, int>> adj_map_;
    
    // Vector holding <Index, Edge> mapping
    std::vector<std::tuple<int, int>> index_edge_map_;

    // Counter to keep track of the size of the vector of Edges - optional
    unsigned num_edges_;
    


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
        
        // Nothing needs to go in here because we don't need to
        // do anything with a Node defined in this way. Calling
        // Graph::Node::Node() will just create an empty node with
        // no information, an invalid Node.
        
        // See the private variables and constructor.
        
        
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
        
        // Access the parent graph by deferencing the pointer,
        // use the vector of points and the index to get the
        // position

        const Point& p = gp_->points_[index_];
      return p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }
      /** Return this node's graph pointer. */
      Graph* gp() const {
        // HW0: YOUR CODE HERE
        return gp_;
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
        
        // check if index of this node is same as index
        // of node n. Needs to work even if they belong
        // to different graphs.

        
      (void) n;          // Quiet compiler warning
      return (index_ == n.index());
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
        
        // check if index of this node is less than index
        // of node n. Needs to work even if they belong
        // to different graphs.
        
      (void) n;           // Quiet compiler warning
      return (index_ < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
      
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
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
      
    // Add the point to the vector of points which
    // exists in the graph, and create an instance
    // of the Node() class to return.
      
      points_.push_back(position);
      Node n = Node(this, num_points_);
      ++num_points_;
      
    (void) position;      // Quiet compiler warning
    return n;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
      // Does the index of this node exist in the Vector,
      // and are they pointing to the same graph?
    
      
    (void) n;            // Quiet compiler warning
    return ((n.index() <= num_points_ - 1) && (this==n.gp_));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
      
      // Construct a new Node pointing to this
      // graph.
      
    (void) i;             // Quiet compiler warning
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
        
        // No code required to create an invalid Edge,
        // Graph::Edge::Edge() will create an empty Edge
        
        // See private variables and functions
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
        
      return Node(gp_, node_a_index_);;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
        
      return Node(gp_, node_b_index_);      // Invalid Node
    }
      
      /** Return this edge's index, a number in the range [0, num_nodes). */
      size_type index() const {
        // HW0: YOUR CODE HERE
        return index_;
      }
      
      /** Return this edge's graph pointer. */
      Graph* gp() const {
        // HW0: YOUR CODE HERE
        return gp_;
      }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        // ADDED HW0: YOUR CODE HERE
        
        // Check if index of this edge is equal to index
        // of edge e. Needs to work even if they belong
        // to different graphs.
        
      (void) e;           // Quiet compiler warning
      return (index_ == e.index());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        // ADDED HW0: YOUR CODE HERE
        
        // Check if index of this edge is less than index
        // of edge e. Needs to work even if they belong
        // to different graphs.
        
      (void) e;           // Quiet compiler warning
      return (index_ < e.index());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
      
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
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
          
    // Use index_edge_map_ to get the indices of the nodes,
    // then return a new instance of that edge.
      
      // check i is a valid edge index
      assert(i<index_edge_map_.size());
      
      // tuple of ordered node indices
      std::tuple<int, int> t = index_edge_map_[i];
      
      size_type ia = std::get<0>(t);
      size_type ib = std::get<1>(t);
      
      // Construct valid Edge
      Edge e = Edge(this, Node(this, ia), Node(this, ib), i);
          
    (void) i;             // Quiet compiler warning
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
    // HW0: YOUR CODE HERE
      
      // Check nodes are valid distinct nodes of this graph
      distinct_valid_nodes(a, b);
      
      // declare lower and upper node index, l < u
      size_type l = (a.index() < b.index()) ? a.index() : b.index();
      size_type u = (a.index() > b.index()) ? a.index() : b.index();
      
      // map with upper node index as its key
      std::map<int, int> m = adj_map_[l];
      
      // if and only if u exists as a key, edge exists
      if (m.count(u)==1)
          return true;

    (void) a; (void) b;   // Quiet compiler warning
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
          num_edges_++;
          adj_map_[l][r] = num_edges_;
          index_edge_map_.push_back ({l,r});
          // set index of Edge we are returning
          i = num_edges_;
      }
      
      // if it already exists use its index i which can be found
      // from the adjancy nested map
      i = adj_map_[l][r];
      
      // Construct edge e
      Edge e = Edge(this, a, b, i);
      
    (void) a, (void) b;   // Quiet compiler warning
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
      
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

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

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
