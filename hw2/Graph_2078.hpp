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
#include <math.h>

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
template <typename V, typename E>
class Graph {
  typedef V node_value_type;
  typedef E edge_value_type;

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
  //using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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

  using incident_edge_iterator = typename std::unordered_map<size_type, edge_value_type>::iterator;

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

    /** Return this node's position as a non-const Point. */
    Point& position() const {
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
    */
    const node_value_type& value() const {
      return graph->get_value_by_uid(uid);
    }
    
    /**
     * Return the number of incident edges.
     * @post: Node is unchanged.
     */
    size_type degree() const {
      return graph->edge_data[uid].size(); 
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
        return IncidentIterator(this, 0); // returns pointer to this node and index
      } else {//degree == 0 
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

    const size_type get_uid() const {
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
    edge_data[uid]; 
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

    edge_value_type& value() {
      return graph->edge_data[std::get<0>(uids)][std::get<1>(uids)];
    }

    const edge_value_type& value() const {
      return graph->edge_data[std::get<0>(uids)][std::get<1>(uids)];
    }

    /* Compute the length of this edge */
    double length() const{
      return norm(node1().position() - node2().position());
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
  
    // constructor of valid Edges from uids
    Edge(Graph* g, std::tuple<size_type, size_type> ids){
      graph = g;
      uids = ids;
    }

    // constructor of valid Edges from two nodes
    Edge(Graph* g, Node& a, Node& b) {
      graph = g;
      uids = std::make_tuple(a.uid, b.uid);
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
   * Achieve complexity: O(num_edges())
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
      // use edge iterator
      EdgeIterator ei = edge_begin();
      size_type j = 0;
      while (j < i) {
        ei++;
        j++;
      }
      return (*ei);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // use hashmap of hashmap of edge data
    // O(1)
    return !( (edge_data.at(a.uid)).find(b.uid) == (edge_data.at(a.uid)).end() );
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
  Edge add_edge(const Node& a, const Node& b, edge_value_type val = edge_value_type()) {
    // If no value is passed, assume default value of edge_value_type
    // HW2
    if (!has_edge(a, b)) {
      // add new edge to adjacency list in 2 places
      edge_data[a.uid][b.uid] = val; 
      edge_data[b.uid][a.uid] = val;
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
    //edges.clear();
    edge_data.clear();
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
     * @post: Returns equal only if graph, index, and uid are the same.
    */
    bool operator==(const NodeIterator& ni) const {
      return ((graph == ni.graph) && (index == ni.index) && (graph->nodes[index] == ni.graph->nodes[ni.index]));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph<V, E>* graph; // pointer to a constant graph
    size_type index; // index into vector of nodes inside graph

    // private constructor for a valid NodeIterator
    NodeIterator(const Graph<V, E>* g, size_type i) {
      graph = g;
      index = i;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return a NodeIterator pointing to "beginning" of the nodes.
   * @post: Nothing about graph has changed.
   */
  node_iterator node_begin() const {
    NodeIterator ni = NodeIterator(this, 0);
    return ni;
  }
  /** Return a NodeIterator pointing to "end" of the nodes.
   * @post: Nothing about graph has changed.
   */
  node_iterator node_end() const {
    NodeIterator ni = NodeIterator(this, num_nodes());
    return ni;
  }

  // return an interator to "middle" of nodes (for use with adjacency list)
  // an alternate "node_end()"
  node_iterator node_middle() const {
    return NodeIterator(this, ceil(num_nodes() / 2));
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
      //iterator=nullptr;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference an IncidentIterator to get an edge incident to this Node.
     * @pre: This is a valid IncidentIterator.   
     * @post: IncidentIterator is unchanged.
     */
    Edge operator*() const {
      size_type node2uid = iterator->first;
      return Edge(node->graph, std::make_tuple(node->uid, node2uid));
    }
    /** Increment this incident iterator.
     * @pre: Incrementing will not cause the incident iterator to go out of bounds.
     */
    IncidentIterator& operator++() {
      assert(iterator != node->graph->edge_data[node->uid].end());
      ++iterator;
      return *this;
    }
    /** Compare another inicident iterator to this one. Return true only if the node
     * and iterator are the same.
     * @post: Nothing has changed about this or the other IncidentIterator.
     */
    bool operator==(const IncidentIterator& ii) const {
      // check if node is the same and if the iterator are pointing to the same edge for this node
      // the first element is the node_uid at the far end of the edge
      return ((node->uid == ii.node->uid) && (iterator->first == ii.iterator->first)); 
    }

    // erase incident edge pointed to by this IncidentIterator
    // returns valid iterator to next unvisited element in unordered_set
    IncidentIterator erase() {
      // erase other part of adjacency list (2->1), by uid
      size_type node2uid = iterator->first;
      node->graph->edge_data[node2uid].erase(node->uid);
      // then erase (1->2) using iterator.erase() functions
      iterator = node->graph->edge_data[node->uid].erase(iterator);
      return (*this);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // things for accessinig graph.edges adjacency list
    const Node* node; // has uid and graph pointer (mostly for safekeeping)
    incident_edge_iterator iterator; // iterator into hashmap of edges for this node

    /** Construct a valid incident iterator. */
    IncidentIterator(const Node* n, size_type i) {
      node=n;
      assert(i <= n->graph->edge_data[n->uid].size());
      iterator=incident_edge_iterator();
      if (i==0) {
        iterator = n->graph->edge_data[n->uid].begin(); // iterator to start of hashmap
      }
      else if(i == n->graph->edge_data[n->uid].size()) {
        iterator = n->graph->edge_data[n->uid].end(); // iterator to end
      }
      else {
        incident_edge_iterator it= n->graph->edge_data[n->uid].begin(); // create iterator to somewhere in the middle
        size_type j = 0;
        while (j < i) {
          it++;
          j++;
        }
        iterator = it;
      }
          }
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference edge iterator and produce Edge.
     * @post: does not modify data.
     */
    Edge operator*() const {
      return *edge_i;
    }
    /** Increment EdgeIterator.
     * @pre: Does not try to increment past graph.edge_end()
     */
    EdgeIterator& operator++() {
      if (!(at_edge_end())) {
        // first try and increment edge_i
        if (edge_i != (*node_i).edge_end()) {
          ++edge_i;
          if (has_been_visited((*edge_i))) {
            ++(*this); // recursive call
          }
        }
        // if cannot, increment node_i and get new edge_i
        if (edge_i == (*node_i).edge_end()) {
            ++node_i;
            edge_i = (*node_i).edge_begin();

            if (has_been_visited((*edge_i)) && (!at_edge_end())){
            ++(*this); // recursive call
            }
        }

        // only mark visitation for middle node of graph, when overlap may start to occur
        if (at_edge_endm1()) {
          visited[(*node_i).uid][(*edge_i).node2().uid] = true;
          visited[(*edge_i).node2().uid][(*node_i).uid] = true;
        }
      }
      return *this;
    }

    bool has_been_visited(Edge e) {
      size_type uid1 = e.node1().uid;
      size_type uid2 = e.node2().uid;
      bool e12 = visited[uid1].find(uid2) != visited[uid1].end(); // true if can find key
      bool e21 = visited[uid2].find(uid1) != visited[uid2].end(); // true if can find key
      return (e12 || e21); // if either have been visited, this edge has been visited
    }

    bool at_edge_endm1() {
      size_type last_node_index = ceil(graphp->num_nodes() / 2);  // bc adjacency list
      return (*node_i).index() == last_node_index - 1;
    }

    bool at_edge_end() {
      size_type last_node_index = ceil(graphp->num_nodes() / 2);  // bc adjacency list
      return (*node_i).index() == last_node_index;
    }

    /** Tests equality of this edge iterator and another.
     */
    bool operator==(const EdgeIterator& ei) const {
      return ((graphp == ei.graphp) && (node_i == ei.node_i) && (edge_i == ei.edge_i)); 
    }

    // Remove the edge pointed to by this EdgeIterator and return an iterator to the next unvisited Edge
    // precondition: never call this function with a node.index() >= node.middle().index()
    EdgeIterator erase() {
      edge_i = edge_i.erase(); // returns new, valid incident iterator
      // if we've iterated through all of this node's incident edges, go to next node
      if (edge_i == (*node_i).edge_end()) {
        if (!at_edge_end()) {
          // if we're not at the end of the nodes
          ++node_i;
          edge_i = (*node_i).edge_begin();
        }
      }
      // check if this edge has been visited yet
      if (has_been_visited((*edge_i)) && (!at_edge_end())) {
        ++(*this);
      }
      return (*this);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Edge iterator consists of pointer to graph, node iterator, and incident edge iterator
    const Graph<V,E>* graphp;
    NodeIterator node_i;
    IncidentIterator edge_i;
    std::unordered_map<size_type, std::unordered_map<size_type, bool>> visited;

    /** Construct a VALID EdgeIterator. 
     * @pre: always passed 0 <= i <= graph.number_edges
    */
    EdgeIterator(const Graph<V, E>* gp, size_type i) {
      graphp = gp;
      edge_i = IncidentIterator();
      node_i = NodeIterator();
      if (i == 0) { //start of edges
        NodeIterator ni = gp->node_begin();
        edge_i = (*ni).edge_begin();
        node_i = ni;
      } else if(i == gp->number_edges) { // edge of  edges
        NodeIterator ni = gp->node_middle();
        edge_i = (*ni).edge_begin(); // chosen custom
        node_i = ni;
      } else { // somewhere in the middle
        throw std::invalid_argument(" cannot ask for EdgeIterator to somewhere in the middle.");
      }
    }
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

  // remove node from the  graph's data structures
  // precondition: this is a valid node in the graph
  // returns: 1 indicating success/that 1 node was removed
  size_type remove_node(const Node& n) {
    size_type uid_to_delete = nodes[n.index()];
    // delete from vector of uids by swapping with end and deleting end
    nodes[n.index()]  = nodes[size()-1];
    nodes.pop_back();
    // delete entry from node_data, and alter index of moved node in node data
    node_data.erase(uid_to_delete);
    node_data[n.index()].index = n.index();
    // delete this node from edge data structures:
    // first delete from node2 positions O(degree)
    // then delete from node1 position
    // ^ both compactly done using IncidentIterator.erase()
    for (IncidentIterator i = n.edge_begin(); i != n.edge_end(); ++i) {
      i = i.erase();
    }
    // THEN delete whole entry from hashmap
    edge_data.erase(uid_to_delete);
    return 1;
  }

  // same as remove_node(Node) above, but takes iterator and returns
  // iterator to next valid Node
  // precondition: node iterator is valid
  node_iterator remove_node(node_iterator n_it) {
    remove_node((*n_it));
    return n_it; // because remove_node(Node) is done by swapping indices, index in this NodeIterator is now valid and unvisited
  }

  // samee as other remove_edge functions
  size_type remove_edge(const Node& a, const Node& b) {
    Edge e = Edge(this, a, b);
    return remove_edge(e);
  }

  // remove edge
  // reurn number of edges removed (1)
  size_type remove_edge(const Edge& e) {
    edge_data[e.node1().uid].erase(e.node2().uid);
    edge_data[e.node2().uid].erase(e.node1().uid);
    number_edges -= 1;
    return 1;
  }

  // remove edge
  // returns valid ierator to next unvisited edge
  edge_iterator remove_edge(edge_iterator e_it){
    return e_it.erase();
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
  size_type number_edges;
  std::unordered_map<size_type, std::unordered_map<size_type, edge_value_type>> edge_data; // map from uid -> (map uid -> data)


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

// for debug
template <typename V, typename E>
std::ostream& operator<< (std::ostream &out, typename Graph<V,E>::Edge data) {
    out << "(" << data.node1().index() << ',' << data.node2().index() << ")";
    // and so on... 
    return out;
    }


#endif // CME212_GRAPH_HPP
