#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * 
 * The code from HW0 has been copied and modified from Graph_1152.hpp from
 * the peer code
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <functional>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

#define print std::cout << "File: " << __FILE__ << " Function: " << __func__ << " Line: " << __LINE__ << std::endl


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {

  /** Private predeclaration of the struct internal_node */
  // NEEDS TO BE PRIVATE SINCE IT IS DECLARED AS PRIVATE IN GRAPH CLASS
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Templatizing the Graph class for nodes (HW1) */
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

  /** Alias for stored edge element */
  using index_pair = std::pair<size_type, size_type>;
  
  /** Invalid node index value */
  static size_type constexpr invalid_index = size_type(-1);

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //


  Graph() {
    // construct an empty graph
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
    Node()
    : node_ptr(nullptr),
      node_idx(invalid_index) {
        // no further initialization
    }

    /** Return this node's position. */
    const Point & position() const {
      return node_ptr->nodes_vec.at(node_idx).coords;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Returns a node_value
    * 
    * This function returns the node_value variable for each node
    * This may be an int like the distance, as measured in 
    * number of edges, from a given node 
    */
    node_value_type& value() {
      //--style_1
      //--Don't use const_cast unless you absolutely have to! There is a way to
      //--get around it here.
      //--START
      return const_cast<internal_node&>(node_ptr->nodes_vec.at(node_idx)).property;
      //--END
    }

    /** Returns a const node_value
    * 
    * This function returns the node_value variable for each node
    * This may be an int like the distance, as measured in 
    * number of edges, from a given node 
    */
    const node_value_type& value() const {
      return node_ptr->nodes_vec.at(node_idx).property;
    }

    /** Returns the degree of the node
    * 
    * This function returns the degree for each node
    * The degree of a node is the number of edges connected to the node
    */
    size_type degree() const {
      return node_ptr->nodes_vec.at(node_idx).degree;
    }

    /** Returns a begin iterator to an edge that connects to the given node
    * 
    * This function returns the begin incident iterator for the node
    * This can be used to iterate over all the edges that connect to the node
    * Such as when finding the distance of the node from a root node
    * (see shortest_path.cpp)
    */
    incident_iterator edge_begin() const {
      //--style_0
      //--Remove commented-out code before submission
      //--START
      // print;
      // for(auto it = mp.begin(); it != mp.end(); ++it) {
      //   std::cout << "Edge " << it->second << " - (" 
      //             << node_idx << ", " << it->first << ")" << std::endl;
      // }
      //--END

      return incident_iterator(node_idx, node_ptr, node_ptr->nested_map.at(node_idx).begin());
    }


    /** Returns an end iterator to an edge that connects to the given node
    * 
    * This function returns the end incident iterator for the node
    * This can be used to iterate over all the edges that connect to the node
    * Such as when finding the distance of the node from a root node
    * (see shortest_path.cpp)
    */
    incident_iterator edge_end() const {
      return incident_iterator(node_idx, node_ptr, node_ptr->nested_map.at(node_idx).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node & n) const {
      return node_ptr == n.node_ptr && node_idx == n.node_idx;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node & n) const {
      return std::less<const Graph *>{}(node_ptr, n.node_ptr) || 
             (node_ptr == n.node_ptr && node_idx < n.node_idx);
    }

    /** constructor used by Graph and NodeIterator to create Node objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] index of the node in owning Graph
     */
    Node(const graph_type* graph, const size_type index)
      : node_ptr(graph),
        node_idx(index) {
      if(!is_valid()) {
        std::cout << "Node - sanity check - this should never occur." << std::endl;
      }
    }

    /** Check node for validity.
     * @return @p true if node is valid (i.e. represents an actual graph node)
     * or @p false otherwise.
     */ 
    bool is_valid() const {
      return node_ptr != nullptr &&
             node_idx != invalid_index;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /** Pointer to the containing Graph */
    const graph_type* node_ptr;
    
    /** Index of the node in the containing Graph */
    size_type node_idx;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    //--style_1
    //--Why the static_cast? Not necessary.
    //--START
    return static_cast<size_type>(nodes_vec.size());
    //--END
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val The new node's node_value_type
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point & position, const node_value_type& val = node_value_type()) {
    nodes_vec.emplace_back(position, val, 0);
    // for now the node has 0 degree since an edge has not been added
    return Node(this, num_nodes() - 1);
  }

//   Node add_node(const Point & position, const node_value_type& v = node_value_type())
//   {
// ////////////////////////////////////////////////////////////////////////////////////////////
//     nodes_vec.push_back(position);
//     return Node(this, num_nodes() - 1);
//   }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node & n) const {
    return n.node_ptr == this && n.index() < num_nodes();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
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
    Edge()
    : edge_ptr(nullptr),
      edge_idx(invalid_index),
      first_node(invalid_index),
      second_node(invalid_index) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(edge_ptr, first_node);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(edge_ptr, second_node);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge & e) const
    {
      // Graph class maintains the invariant first_node <= second_node.
      // Therefore we don't have to check the other combination.
      return edge_ptr == e.edge_ptr &&
             first_node == e.first_node &&
             second_node == e.second_node;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge & e) const
    {
      return std::less<const Graph *>{}(edge_ptr, e.edge_ptr) ||
             (edge_ptr == e.edge_ptr && (first_node < e.first_node || 
                                      (first_node == e.first_node && second_node < e.second_node)));
    }

    size_type index() {
      return edge_idx;
    }

    /** constructor used by Graph, IncidentIterator, EdgeIterator to create Edge objects.
     * @param[in] graph pointer to the owning Graph
     * @param[in] node1 index of first node in owning Graph
     * @param[in] node2 index of second node in owning Graph
     */
    Edge(const graph_type* const graph, 
         const size_type index,
         const size_type node1, 
         const size_type node2)
    : edge_ptr(graph),
      edge_idx(index),
      first_node(node1),
      second_node(node2) {
        // do nothing
    }

   private:

    /** Check if the edge is valid.
     * @return @p true if edge is valid (i.e. represents an actual graph edge)
     * or @p false otherwise.
     */ 
    bool is_valid() const
    {
      return edge_ptr != nullptr && 
             edge_idx != invalid_index &&
             first_node != invalid_index && 
             second_node != invalid_index &&
             first_node <= second_node;
    }
    
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    /** Pointer to the containing Graph */
    const graph_type* edge_ptr;
    
    /** Index of the edge */
    size_type edge_idx;
    
    /** Index of first node */
    size_type first_node;
    
    /** Index of second node */
    size_type second_node;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    //--style_0
    //--What's going on with all the casting? edges_vec.size() is already
    //--a size_type or convertible to size_type.
    //--START
    return (size_type) edges_vec.size();
    //--END
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    index_pair p = edges_vec.at(i);
    return Edge(this, i, p.first, p.second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return edges_map.count(get_node_indices(a, b)) > 0;
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
  Edge add_edge(const Node & a, const Node & b) {
    const index_pair nodes = get_node_indices(a, b);
    const size_type index = num_edges();
    
    auto it = edges_map.emplace(nodes, index); 
    // pair of iterator to position and bool
      // bool is true if inserted
      // bool is false if element already exists

    if(it.second) {
      // enters this block only if a new key-val pair have been added to edges_map
      edges_vec.push_back(nodes);

      // incrementing the degree of each node if the edge was added
      nodes_vec.at(a.node_idx).degree++;
      nodes_vec.at(b.node_idx).degree++;


      nested_map.insert(std::make_pair(nodes.first, std::map<size_type, size_type>()));
      nested_map.at(nodes.first).insert(std::make_pair(nodes.second, index));


      nested_map.insert(std::make_pair(nodes.second, std::map<size_type, size_type>()));
      nested_map.at(nodes.second).insert(std::make_pair(nodes.first, index));
    }

    else {
      nested_map.at(nodes.first).insert(std::make_pair(nodes.second, index));
      nested_map.at(nodes.second).insert(std::make_pair(nodes.first, index));
    }
    
    return Edge(this, index, nodes.first, nodes.second);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_map.clear();
    edges_vec.clear();
    nodes_vec.clear();
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
      // no code added here
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Returns the node that the node iterator points to
    * 
    * This function returns the Node that the node iterator points to
    * This can be used to obtain information about the node of interest
    * If the iterator is at the end, an invalid Node is returned
    */
    node_type operator*() const {
      if(!(node_iter_idx < node_iter->size())) {
        return Node();
      }
      
      else {
        return Node(nodeiter_ptr, node_iter_idx);
      }
    }

    /** Increments the node iterator
    * 
    * This function returns the Node iterator that points to the next Node
    * in the graph. This can be useful when one needs to iterate over all 
    * the nodes in a given graph.
    */
    node_iterator& operator++() {
      if(node_iter_idx < node_iter->size()) {
        node_iter_idx++;
        return *this;
      }

      else {
        node_iter_idx = node_iter->size();
        return *this;
      }
    }

    /** Tests the equality of the current instance of the node iterator to
    * another instance
    * 
    * Two node iterators are the same when the node_iter's are equal
    * the index of the node, being pointed to, are also equal
    */
    bool operator==(const NodeIterator& NodeIt) const {
      bool pred1 = node_iter == NodeIt.node_iter;
      bool pred2 = node_iter_idx == NodeIt.node_iter_idx;

      return pred1 && pred2;
    }

    /** Constructor of a node iterator
    * 
    * allows for the construction of an node iterator
    */
    NodeIterator(std::vector<internal_node>* p, size_type i, const graph_type* g)
    : node_iter(const_cast<std::vector<internal_node>*>(p)),
      node_iter_idx(i),
      nodeiter_ptr(const_cast<graph_type*>(g)) {
        // do nothing
      }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    std::vector<internal_node>* node_iter;
    size_type                   node_iter_idx;
    const graph_type*           nodeiter_ptr;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns a node iterator that points to the beginning node within a graph
  * 
  * This can be used to iterate over all the nodes in a graph.
  */
  node_iterator node_begin() const {
    return NodeIterator(const_cast<std::vector<internal_node>*>(&nodes_vec),
                        0, this);
  }

  /** Returns a node iterator that points to the end node within a graph
  * 
  * This can be used to iterate over all the nodes in a graph.
  */
  node_iterator node_end() const {
    return NodeIterator(const_cast<std::vector<internal_node>*>(&nodes_vec),
                        nodes_vec.size(), this);
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
      // do nothing
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Returns an edge that connects to the node of interest
    * 
    * When dereferencing the iterator, the user obtains an Edge object
    * @post The first node of the edge is always the node of interest
    */
    edge_type operator*() {
      if(map_iter == incid_ptr->nested_map.at(first_idx).end()) {
        return Edge();
      }

      return Edge(incid_ptr, map_iter->second, first_idx, map_iter->first);
    }

    /** Increments the incident iterator
    * 
    * This results in the incident iterator pointing to the next edge that
    * connects to the node of interest
    */
    incident_iterator& operator++() {
      //--style_0
      //--What does this do?
      //--START
      if(map_iter == incid_ptr->nested_map.at(first_idx).end()) {
        map_iter = incid_ptr->nested_map.at(first_idx).end();
      }
      //--END

      ++map_iter;
      return *this;
    }

    /** Tests equality of two incident iterators
    * 
    * Two incident iterators are equal when the two iterators point
    * to the same edge 
    */
    bool operator==(const incident_iterator& it) const {
      return map_iter == it.map_iter;
    }

    /** Constructor for the incident iterator
    * 
    */
    IncidentIterator(size_type index,
                     const graph_type* g,
                     std::map<size_type, size_type>::const_iterator p)
    : first_idx(index),
      incid_ptr(const_cast<graph_type*>(g)),
      map_iter(p) {
        // do nothing
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type                                      first_idx;
    graph_type*                                    incid_ptr;
    std::map<size_type, size_type>::const_iterator map_iter;
    // 1st node idx = first_idx
    // 2nd node idx = map_iter->first
    // edge idx     = map_iter->second
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

    /** Returns an edge that the edge iterator points to
    * 
    * @post Dereferencing the end iterator returns an invalid Edge
    */
    Edge operator*() const {
      if(!(edgeiter_idx < edgeiter->size())) {
        return Edge();
      }

      else {
        index_pair p = edgeiter->at(edgeiter_idx);

        return Edge(edgeiter_ptr, edgeiter_idx, p.first, p.second);
      }
    }

    /** Increments the edge iterator
    * 
    * The iterator now points to the next edge in the graph.
    */
    edge_iterator& operator++() {
      if(edgeiter_idx < edgeiter->size()) {
        edgeiter_idx++;
        return *this;
      }

      else {
        edgeiter_idx = edgeiter->size();
        return *this;
      }
    }

    /** Tests the equality of two edge iterators
    * 
    * Two edge iterators are equal when they point to the same graph,
    * the same edge 
    */
    bool operator==(const EdgeIterator& it) const {
      // print;
      bool pred1 = edgeiter_idx == it.edgeiter_idx;
      bool pred2 = edgeiter_ptr == it.edgeiter_ptr;
      bool pred3 = edgeiter     == it.edgeiter;

      return pred1 && pred2 && pred3;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const graph_type*              edgeiter_ptr;
    //--design_0
    //--Why is this a pointer to const vector rather than simply a
    //--std::vector<index_pair>::const_iterator?
    //--START
    const std::vector<index_pair>* edgeiter;
    //--END
    size_type                      edgeiter_idx;

    /** Private constructor 
    * 
    */
    EdgeIterator(std::vector<index_pair>* ptr, size_type i)
    : edgeiter(ptr),
      edgeiter_idx(i) {
        // doing nothing
      }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns an edge iterator that points to the first edge of the graph
  * 
  * This can be used to iterate over all edges within a graph
  */
  edge_iterator edge_begin() const {
    EdgeIterator it(const_cast<std::vector<index_pair>*>(&edges_vec), 0);
    it.edgeiter_ptr = this;

    return it;
  }

  /** Returns an edge iterator that points to one past the last edge of the graph
  * 
  * This can be used to iterate over all edges within a graph
  */
  edge_iterator edge_end() const {
    EdgeIterator it(const_cast<std::vector<index_pair>*>(&edges_vec), edges_vec.size());
    it.edgeiter_ptr = this;

    return it;
  }

 private:

  // Enforce the size constraint
  static_assert(sizeof(Node) <= 16, "Node class exceeds size limit");
  static_assert(sizeof(Edge) <= 32, "Edge class exceeds size limit");
              
  /** Get node indices in correct order for an edge
   * @pre @a a and @a b are valid nodes of this graph
   * @param[in] a first node
   * @param[in] b second node
   * @return pair of node indices in correct order
   */
  index_pair get_node_indices(const Node & a, const Node & b) const {
    //--style_0
    //--Use std::minmax here instead.
    //--START
    index_pair p(std::min(a.index(), b.index()),
                 std::max(a.index(), b.index()));

    return p;
    //--END
  }

  struct internal_node {
    Point coords;
    node_value_type property;
    //--design_0
    //--If you think about it, there's really no reason to store the node degree
    //--here. You can always derive it by looking up its neighbors in nested_map
    //--and getting the size of the set of neighbors.
    //--START
    size_type degree;
    //--END

    internal_node(const Point p, const node_value_type val, const size_type d)
    : coords(p),
      property(val),
      degree(d) {
        // doing nothing here
      }
  };
  
  /** Storage for node points */
  std::vector<internal_node> nodes_vec;
  
  /** Storage for edges, smaller node number first */
  std::vector<index_pair> edges_vec;
  
  /** Fast edge lookup by index nodes */
  //--design_1
  //--There is no need for edges_map. Edge lookup can already be accomplished in
  //--O(log(max node degree)) by using the adjacencies in nested_map. Find a way
  //--to remove this data structure. I think it will also greatly simplify the
  //--code elsewhere.
  //--START
  std::map<index_pair, size_type> edges_map;
  // <<1st node index, 2nd node index>, edge index>
  //--END

  /** Fast iteration over all edges */
  std::map<size_type, std::map<size_type, size_type>> nested_map;
  // <1st node index, <2nd node index, edge index>>
};

//--style_0
//--There are a large number of unnecessary const_casts throughout the code that
//--indicate a lack of clear understanding. I would suggest removing all
//--const_casts and then seeing what compilation failures you get, then only
//--slowly reintroducing the necessary ones as you understand why they are
//--needed.
//--END

#endif
