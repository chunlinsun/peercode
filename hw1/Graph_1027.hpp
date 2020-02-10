#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <tuple>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// print edge hashmap
template<typename K, typename E>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<K,std::vector<E>> &m) {
  os << "{";
  for (const std::pair<K,std::vector<E>>& p:m) {
    os << "(" << p.first << ", " ;
    os << "[";
    for (const E& e:p.second) {
      os << e << ",";
    }
    os << "]";
    os << ")" << std::endl;
  }
  os << "}";
  return os;
}

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph<V>;

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
    // Rely on default initialization
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
      graph = nullptr;
      uid = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph->get_point_by_uid(uid);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph->get_index_by_uid(uid);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return the reference to value of this node. 
     * @pre: This node is a valid node.
    */
    node_value_type& value() {
      return graph->get_value_by_uid(uid);
    }
    /** Return constant reference to value of this node. 
     * @pre: This node is a valid node.
     //--documentation_0
     //--This postcondition is already implied because the method is marked const;
     //--no need to explicitly call it out.
     //--START
     * @post: value_old == value_old 
     //--END
    */
    const node_value_type& value() const {
      return graph->get_value_by_uid(uid);
    }
    
    /**
     * Return the number of incident edges.
     * @post: Node is unchanged.
     */
    size_type degree() const {
      return graph->edges[uid].size();
    }
    /** Node constructs an incident iterator to firs edge incident
     * to itself.
     * @pre: This is a valid node.
     * @post: The incident_iterator returned is valid and points to
     * first incident edge in some ordering IF there are any incident edges
     * OTHERWISE returns an INVALID incident_iterator
     */
    incident_iterator edge_begin() const {
      if (degree() > 0) {
        std::cout << "deg > 0" << std::endl;
        return IncidentIterator(this, 0); // returns pointer to this node and index
      } else {//degree == 0 
        std::cout << "deg == 0" << std::endl;
        return IncidentIterator();
      }
    }
    /** Node constructs an incident iterator to last edge incident
     * to itself.
     * @pre: This is a valid node.
     * @post: The incident_iterator returned is valid and points to
     * last incident edge in some ordering IF there are any incident edges
     * OTHERWISE returns an INVALID incident_iterator
     */
    incident_iterator edge_end() const{
      if (degree() > 0) {
        return IncidentIterator(this, degree());
      } else {//degree == 0 
        return IncidentIterator();
      }
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (uid == n.uid) && (graph == n.graph);
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
      return (uid < n.uid);
    }

    size_type get_uid() const {
      return uid;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // data members
    // Pointer to the Graph
    Graph* graph;
    // uid of the Node (corresponds to internal_node.uid in graph)
    size_type uid;

    // constructor of valid nodes
    Node(Graph* g, size_type id) {
      graph = g;
      uid = id;
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes.size();
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
  Node add_node(const Point& position, const node_value_type& value= node_value_type()) {
    // HW0: YOUR CODE HERE
    // add new node to internal graph
    size_type uid = next_node_uid;
    internal_node inode {position, size(), value};
    nodes.push_back(uid);
    node_data[uid] = inode;
    next_node_uid = next_node_uid + 1;
    // add node to edge adjacency list
    edges[uid]; // = std::vector<size_type>();
    // construct accessor to this node as Node object and return it
    Node n = Node(this, uid);
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph == this) && has_uid(n.uid);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // retrive data for index i node
    // construct Node
    return Node(const_cast<Graph*>(this), nodes[i]);
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
      graph = nullptr;
      std::tuple<size_type, size_type> d {0,0};
      uids = d;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph, std::get<0>(uids));     
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph, std::get<1>(uids)); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE 
      // test: same graph, same uids
      bool same_graph = (graph == e.graph);
      bool uid_match = (std::get<0>(uids) == std::get<0>(e.uids)) && (std::get<1>(uids) == std::get<1>(e.uids));
      bool uid_reverse_match = (std::get<0>(uids) == std::get<1>(e.uids)) && (std::get<1>(uids) == std::get<0>(e.uids));
      bool same_uids = uid_match || uid_reverse_match;
      return same_graph && same_uids;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // compare sum of edge uids
      size_type this_e = std::get<0>(uids) + std::get<1>(uids);
      size_type other_e = std::get<0>(e.uids) + std::get<1>(e.uids);
      return this_e < other_e; 

    }

    std::string print() {
      return "("+ std::to_string(std::get<0>(uids)) + " , " + std::to_string(std::get<1>(uids)) + ")";
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    
    // data members of Edge
    Graph* graph;
    std::tuple<size_type, size_type> uids;
  
    // constructor of valid Edges
    Edge(Graph* g, std::tuple<size_type, size_type> ids){
      graph = g;
      uids = ids;
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return number_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
      // use edge vector for fast access
      return Edge(const_cast<Graph*>(this), edges_i[i]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // walk through incident edges of a node  in adjacency list
    size_type num_incident_a = edges.at(a.uid).size();
    for (size_type i = 0; i < num_incident_a; ++i) {
      if (edges.at(a.uid)[i] == b.uid){
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
    // HW0: YOUR CODE HERE
    if (!has_edge(a, b)) {
      // add new edge to adjacency list in 2 places
      edges[a.uid].push_back(b.uid);
      edges[b.uid].push_back(a.uid);
      // also add to edge vector
      edges_i.push_back(std::make_tuple(a.uid, b.uid));
      number_edges += 1;
    } 
    // return edge object
    return Edge(this, std::make_tuple(a.uid, b.uid)); 
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // clear nodes
    nodes.clear();
    node_data.clear();
    //next_node_uid = 0; < not strictly necessary
    // clear edges
    edges.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graph=nullptr;
      index=-1;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference this NodeIterator to a Node.
     * Does not alter any arguments during dereference.
     * @pre: This is a valid interator to valid graph.
    */
    Node operator*() const {
      return Node(const_cast<Graph*>(graph), graph->nodes[index]);
    }
    /** Increment this NodeIterator
     * @pre: This is a valid interator to valid graph.
     * @pre: Incrementing this iterator will not put it
     * out of bounds for the corresponding graph. e.g. 0<= index <= num_nodes
    */
    NodeIterator& operator++() {
      assert(index < graph->num_nodes());
      index += 1;
      return *this;
    }
    /** Compare another NodeIterator to this one.
     * @post: Returns equal only if graph and index are the same.
    */
    bool operator==(const NodeIterator& ni) const {
      return ((graph == ni.graph) && (index == ni.index));
    }

   private:
    //--design_0
    //--There really should be a private constructor here so that node_begin()
    //--and node_end() can directly invoke NodeIterator(const Graph<V>*, size_type).
    //--START
    //--END
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph<V>* graph; // pointer to a constant graph
    size_type index; // index into vector of nodes inside graph
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return a NodeIterator pointing to "beginning" of the nodes.
   * @post: Nothing about graph has changed.
   */
  node_iterator node_begin() const {
    NodeIterator ni;
    ni.graph = this;
    ni.index = 0;
    return ni;
  }
  /** Return a NodeIterator pointing to "end" of the nodes.
   * @post: Nothing about graph has changed.
   */
  node_iterator node_end() const {
    NodeIterator ni;
    ni.graph = this;
    ni.index = num_nodes();
    return ni;
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      node=nullptr;
      index=-1;
    }
    //--design_1
    //--This constructor should be private! If you make this constructor public,
    //--you lose the tight control over constructing valid incident iterators.
    //--The only way i should be able to obtain a valid IncidentIterator from
    //--the graph is by calling n.edge_begin() or n.edge_end() on a valid Node n.
    //--START
    /** Construct a valid incident iterator. */
    IncidentIterator(const Node* n, int i) {
      node=n;
      index=i;
    }
    //--END

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference an IncidentIterator to get an edge incident to this Node.
     * @pre: This is a valid IncidentIterator.
     * @post: IncidentIterator is unchanged.
     */
    Edge operator*() const {
      size_type node2uid = node->graph->edges[node->uid][index];
      return Edge(node->graph, std::make_tuple(node->uid, node2uid));
    }
    /** Increment this incident iterator.
     * @pre: Incrementing will not cause the incident iterator to go out of bounds.
     */
    IncidentIterator& operator++() {
      assert(index < node->graph->edges[node->uid].size());
      index += 1;
      return *this;
    }
    /** Compare another inicident iterator to this one. Return true only if the node
     * and index are the same.
     * @post: Nothing has changed about this or the other IncidentIterator.
     */
    bool operator==(const IncidentIterator& ii) const {
      return ((node->uid == ii.node->uid) && (index == ii.index)); 
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // things for accessinig graph.edges adjacency list
    const Node* node; // has uid and graph pointer
    size_type index; // index into vector of incident nodes/edges
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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

    //--design_0
    //--As with IncidentIterator, this constructor should be private.
    //--START
    /** Construct a VALID EdgeIterator. */
    EdgeIterator(const Graph<V>* gp, size_type i) {
      graphp = gp;
      index = i;
    }
    //--END

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference edge iterator and produce Edge.
     * @post: does not modify data.
     */
    Edge operator*() const {
      return Edge(const_cast<Graph*>(graphp), graphp->edges_i[index]);
    }
    /** Increment EdgeIterator.
     * @pre: Does not try to increment past graph.edge_end()
     */
    EdgeIterator& operator++() {
      assert(index < graphp->number_edges);
      index += 1;
      return *this;
    }
    /** Tests equality of this edge iterator and another.
     */
    bool operator==(const EdgeIterator& ei) const {
      return ((graphp == ei.graphp) && (index == ei.index)); 
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Edge iterator consists of pointer to graph and edge index
    const Graph<V>* graphp;
    size_type index;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Returns iterator to "first" edge in graph. Does not alter data.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  /** Returns iterator to "last" edge in graph. Does not alter data.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, number_edges);
  }

  std::string print() {
    return "(edges: " + std::to_string(number_edges) + ", nodes: " + std::to_string(size()) + ")";
  }

 private:

  // HW0: YOUR CODE HERE
  // The node data is stored in 2 containers:
  // 1) A vector of unique identifiers (uids) (index in vector is index)
  // 2) An unordered hashmap from uids -> struct(Point, index)
  // Together, these two containers allow fast access to graph
  // data from indices OR from uids AS WELL as fast Node.index() calls
  std::vector<size_type> nodes; // vector of uids
  struct internal_node {
      Point loc; // location of node
      size_type index; // index of node in vector of uids (for fast .index())
      node_value_type val; // value of node
  };
  std::unordered_map<size_type, internal_node> node_data; // uids -> internal nodes
  size_type next_node_uid;
  //--design_0
  //--This is good. The outer std::unordered_map allows O(1) access to the vector
  //--of neighboring nodes. But then checking for has_edge(a, b) is O(degree(a)).
  //--Does the collection of neighbors of a node need to be a vector (ordered)?
  //--Can it be something with even better performance?
  //--START
  std::unordered_map<size_type, std::vector<size_type>> edges; // uids -> adjacency list of incident uids
  //--END
  std::vector<std::tuple<size_type, size_type>> edges_i; // (uid, uid) -> list of edges for fast iteration
  size_type number_edges;

  /** Return Point object corresponding to uid. 
   * @pre: uid is valid for graph
   */
  Point& get_point_by_uid(size_type uid) {
    return node_data[uid].loc;
  }

  /** Return value by uid.
   * @pre: uid is valid for graph
  */
  node_value_type& get_value_by_uid(size_type uid) {
    return node_data[uid].val;
  }

  /** Return index of internal node corresponding to uid. 
   * @pre: uid is valid for graph
  */
  size_type get_index_by_uid(size_type uid) {
    return node_data[uid].index;
  }

  bool has_uid(size_type node_uid) const {
    // O(1) with a hashmap:
    return !(node_data.find(node_uid) == node_data.end());
  }

};

//--functionality_0
//--I'm a tad confused---your Graph class seems great! But subgraph.cpp was not
//--attempted and shortest_path.cpp didn't compile. If this was simply a result
//--of having too many things to do and needing to let something drop, i
//--completely understand. If you turned this in expecting it to be fine, then
//--i'd recommend reaching out or dropping by office hours.
//--END

// for debug
template <typename V>
std::ostream& operator<< (std::ostream &out, typename Graph<V>::Edge data) {
    out << "(" << data.node1().index() << ',' << data.node2().index() << ")";
    // and so on... 
    return out;
    }





#endif // CME212_GRAPH_HPP
