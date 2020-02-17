#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <vector>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //


  using node_value_type = V;
  using edge_value_type = E;
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
  Graph() : points(), values(), edg_by_ord(), edg_by_idx(), node_edges(), 
      edge_value_vec(), id_compl_to_exist(), id_exist_to_compl(), nr_nodes(0) , 
      nr_edges(0), nr_nodes_created(0) {
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
    Node() {}

    /** standard Node constructor
     * param[in] g, pointer to a Graph
     *           i, index for which Node was initialised
     */
    Node(const Graph* g, size_type i) {
      idx = i;
      graph = g;
    }

    /** Return this node's position. */
    const Point& position() const { 
      // position saved inside the Graph's points vector
      return (*graph).points[this->idx];
    }

    /** Return this node's position (Position is modifiable) */
    Point& position() { 
      // position saved inside the Graph's points vector
      Graph* gr= const_cast<Graph*>(graph);
      size_type idx_exist = (*gr).id_compl_to_exist[this->idx];
      return (*gr).points[idx_exist];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      size_type idx_exist = (*graph).id_compl_to_exist[this->idx];
      return idx_exist;
    }

    // YOUR HW1 #1 CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    // return a reference to the value of the node (value can be changed)
    node_value_type& value() {
      Graph* gr= const_cast<Graph*>(graph);
      size_type idx_exist = (*gr).id_compl_to_exist[this->idx];
      return (*gr).values[idx_exist];
    }

    // return a reference to the constant value of the node
    const node_value_type& value() const {
      size_type idx_exist = (*graph).id_compl_to_exist[this->idx];
      return (*graph).values[idx_exist];
    }

    //the number of adjacent nodes
    size_type degree() const {
      return (*graph).node_edges[this->idx].size();
    }

    // initialization of the IncidentIterator
    incident_iterator edge_begin() {
      return (IncidentIterator(const_cast<Graph*>(graph), idx));
    }

    // end of the IncidentIterator
    incident_iterator edge_end() const {
      return IncidentIterator(nullptr, 0);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (this->idx != n.index()) {
        return false;
      }
      if (this->graph != n.graph) {
        return false;
      }
      return true;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * 
     * Nodes are order by their index and by their
     * Graph pointer if the indices are the same
     */
    bool operator<(const Node& n) const {
      // calculate 2-norms
      //double norm1 = norm(position());
      //double norm2 = norm(n.position());

      // determine lower index
      if (this->idx < n.idx) {
        return true;
      }
      else if (this->idx > n.idx) {
        return false;
      }
      // indices equal: compare Graphs
      else {
        return (std::less<const Graph*>{}(graph, n.graph));
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    // stored index is the id of the complete order (ignoring removals)
    size_type idx;
    const Graph* graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nr_nodes;
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // add this point to the saved points
    points.push_back(const_cast<Point&>(position));
    // add the value to the saved values
    values.push_back(const_cast<node_value_type&>(value));
    // add empty vector to the node_edges vector
    node_edges.push_back(std::vector<size_type> {});
    // add empty map to edge_value_vec
    edge_value_vec.push_back(std::map<size_type, edge_value_type> {});
    // add to id_compl_to_exist
    id_compl_to_exist.push_back(nr_nodes);
    // add to id_exist_to_compl
    if (id_exist_to_compl.size() <= nr_nodes) {
      id_exist_to_compl.push_back(nr_nodes_created);
    }
    else {
      id_exist_to_compl[nr_nodes] = nr_nodes_created;
    }

    Node ret_node = node(nr_nodes_created);
    // and increment the nodes counter
    nr_nodes += 1;
    nr_nodes_created += 1;
    return ret_node;        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // check if index is inside the vector's range
    if (n.index() >= nr_nodes) {
      return false;
    }
    // check if node gives the correct Graph
    else if (this != n.graph) {
      return false;
    }
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const  {
    // give out a new Node with the known properties
    return Node(this, id_exist_to_compl[i]);
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
    Edge() {}

    /** Normal Edge constructor
     * param[in] g, a pinter to the Graph
     *           i, index of the first Node
     *           j, index of the second node
     */
    Edge(const Graph* g, size_type i, size_type j) {
      idx_1 = i;
      idx_2 = j;
      graph = g;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // generate new Node with knwon properties
      return (*graph).node(idx_1);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // generate new Node with knwon properties
      return (*graph).node(idx_2); 
    }

    // return a reference to the value of the node
    edge_value_type& value() {
      Graph* gr= const_cast<Graph*>(graph);
      size_type first_idx = std::min(idx_1, idx_2);
      size_type second_idx = std::max(idx_1, idx_2);
      return (*gr).edge_value_vec[first_idx][second_idx];
    }

    // return a reference to the constant value of the node
    const edge_value_type& value() const {
      size_type first_idx = std::min(idx_1, idx_2);
      size_type second_idx = std::max(idx_1, idx_2);
      return (*graph).edge_value_vec[first_idx][second_idx];
    }

    double length() const {
      Point diff = node1().position() - node2().position();
      return norm(diff);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // 2 equality cases: 
      // * node 1 of this edge equals node 1 of e
      // * node 2 of this graph equals node 1 of e
      if (node1() == e.node1()) {
        if (node2() == e.node2()) {
          return true;
        }
        else {
          return false;
        }
      }
      else {
        if (node2() == e.node1()) {
          if (node1() == e.node2()) {
            return true;
          }
        }
        else {
          return false;
        }
      }
      // to surpress compiler warnings
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * 
     * First compares by the Node's indices, then by their Graph pointers
     */
    bool operator<(const Edge& e) const {
      
      // check which edge has the smallest index
      // if equal smallest indices, check for second smallest index

      size_type help_11 = this->idx_1;
      size_type help_12 = this->idx_2;
      if (this->idx_2 < this->idx_1) {
        help_11 = this->idx_2;
        help_12 = this->idx_1;
      }

      size_type help_21 = e.idx_1;
      size_type help_22 = e.idx_2;
      if (e.idx_2 < e.idx_1) {
        help_21 = e.idx_2;
        help_22 = e.idx_1;
      }

      if (help_11 < help_21) {
        return true;
      }
      else if (help_11 > help_21) {
        return false;
      }
      else if (help_12 < help_22) {
        return true;
      }
      else if (help_12 > help_22) {
        return false;
      }
      // equal dindices: compare the Graph pointers
      else {
        return (std::less<const Graph*>{}(graph, e.graph));
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    size_type idx_1;
    size_type idx_2;
    const Graph* graph;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * number of edges saved as incremented variable -> complexity O(1)
   */
  size_type num_edges() const {
    return nr_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Separate data structure for node as grrouped by order of creation
   * Call in O(1)
   */
  Edge edge(size_type i) const {
    // call indices from edg_by_ord vector
    // and generate the 2 Nodes with known information
    size_type idx_1 = edg_by_ord[i][0];
    size_type idx_2 = edg_by_ord[i][1];
    return Edge(this, idx_1, idx_2);    
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // if no edges: can not have this edge
    if (nr_edges == 0) {
      return false;
    }
    // order the 2 indices of theb 2 Nodes by their size
    size_type idx_1 = a.index();
    size_type idx_2 = b.index();

    size_type idx1_exist = id_exist_to_compl[idx_1];
    size_type idx2_exist = id_exist_to_compl[idx_2];

    size_type idx_a;
    size_type idx_b;
    if (idx1_exist > idx2_exist) {
      idx_a = idx2_exist;
      idx_b = idx1_exist;
    }
    else {
      idx_a = idx1_exist;
      idx_b = idx2_exist;
    }

    // try to find lower index 
    // (which would be a key of the map)
    auto idx_check = edg_by_idx.find(idx_a);
    if (idx_check != edg_by_idx.end()) {
      // key exists
      // check if the higher index is in the corresponding value set
      if ((idx_check->second).find(idx_b) != (idx_check->second).end()) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& e_value = edge_value_type()) {
    // get the indices for the 2 Nodes
    size_type idx_1 = a.index();
    size_type idx_2 = b.index();

    size_type idx1_compl = id_exist_to_compl[idx_1];
    size_type idx2_compl = id_exist_to_compl[idx_2];
    
    if (!has_edge(a, b)) {
      // Node not already in the Graph
      // increment edge counter
      nr_edges += 1;
      // append Node to the vector storing the indices by order
      // of creation
      std::vector<size_type> indices = {idx1_compl, idx2_compl};
      edg_by_ord.push_back(indices);

      // sort the 2 indices
      size_type idx_a;
      size_type idx_b;
      if (idx1_compl > idx2_compl) {
        idx_a = idx2_compl;
        idx_b = idx1_compl;
      }
      else {
        idx_a = idx1_compl;
        idx_b = idx2_compl;
      }

      // check whether lower index is already a key in the map
      auto idx_check = edg_by_idx.find(idx_a);
      if (idx_check != edg_by_idx.end()) {
        // is key
        // append new corresponding index to the value set
        std::set<size_type> to_add = idx_check->second;
        to_add.insert(idx_b);
        // replace old value set with new one for this key
        edg_by_idx[idx_a] = to_add;
      }
      else {
        // lower index not a key yet
        // generate new value set and append to map
        std::set<size_type> to_add = {idx_b };
        edg_by_idx[idx_a] = to_add;
      }

      // add vectors to vector - vector storage
      std::vector<size_type> vec_help = node_edges[idx_a];
      vec_help.push_back(idx_b);
      node_edges[idx_a] = vec_help;

      vec_help = node_edges[idx_b];
      vec_help.push_back(idx_a);
      node_edges[idx_b] = vec_help;

      // add default value to edge vector
      edge_value_vec[idx_a][idx_b] = const_cast<edge_value_type&>(e_value);
      
    }
    return Edge(this, idx1_compl, idx2_compl);
  }

  /** Delete an edge of this graph
   * @param[in]     a         a boundary node of the edge
   * @param[in]     b         the other boundary node of the edge
   * @return        size_type the index of the replaced edge
   *                or 0 if the edge is not on the graph
   *
   * @tparam @a and @b are both of type Node
   * @pre @a and @b are both Nodes belonging to this graph, that is,
   * a.index() < num_nodes, b.index() < num_nodes(). Also, to delete an
   * existing node, @a and @b should have an edge between them, such that
   * there exists an edge e s.t. e.node1() is in {a, b} and e.node2() is the 
   * other node
   * @post There no longer exists said this edge e
   *        num_edges() is reduced by 1.
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (has_edge(a, b)) {
      // delete the corresponding row in edg_by_ord
      // find Node index combination
      size_type idx1_compl = id_exist_to_compl[a.index()];
      size_type idx2_compl = id_exist_to_compl[b.index()];

      std::vector<size_type> inds1 = {idx1_compl, idx2_compl};
      std::vector<size_type> inds2 = {idx2_compl, idx1_compl};

      auto it1 = std::find(edg_by_ord.begin(), edg_by_ord.end(), inds1);
      auto it2 = std::find(edg_by_ord.begin(), edg_by_ord.end(), inds2); 
      size_type ind;
      
      // find edge index
      if (it1 != edg_by_ord.end()) {
        ind = std::distance(edg_by_ord.begin(), it1);
      }
      else {
        ind = std::distance(edg_by_ord.begin(), it2);
      }
      // swap with last element, delete last element
      if (edg_by_ord.size() > 1) {
        unsigned last_ind = (edg_by_ord.size()-1);
        iter_swap(edg_by_ord.begin()+ind, edg_by_ord.begin()+last_ind);
      }
      edg_by_ord.pop_back();

      // delete from edg_by_idx
      unsigned min_ind = std::min(idx1_compl, idx2_compl);
      unsigned max_ind = std::max(idx1_compl, idx2_compl);

      auto idx_check = edg_by_idx.find(min_ind);
      if (idx_check != edg_by_idx.end()) {
        // get set of adjacent indices, delete max_ind
        std::set<size_type> to_add = idx_check->second;
        to_add.erase(max_ind);
        edg_by_idx[min_ind] = to_add;
      }

      // delete both combinations from node_edges
      std::vector<size_type> a_adj = node_edges[idx1_compl];
      auto it_adj = std::find(a_adj.begin(), a_adj.end(), idx2_compl);
      a_adj.erase(it_adj);
      node_edges[idx1_compl] = a_adj;

      std::vector<size_type> b_adj = node_edges[idx2_compl];
      it_adj = std::find(b_adj.begin(), b_adj.end(), idx1_compl);
      b_adj.erase(it_adj);
      node_edges[idx2_compl] = b_adj;

      // delete from edge_value_vec
      std::map<size_type, edge_value_type> ma = edge_value_vec[min_ind];
      //auto it_mapf = ma.find(max_ind);
      ma.erase(max_ind);
      edge_value_vec[min_ind] = ma;

      // decrease number of edges
      nr_edges -= 1;

      return ind;
    }
    else {
      return 0;
    }
  }

  /** Delete an edge of this graph
   * @param[in]     e         an edge from this graph
   * @return        size_type the index of the replaced edge
   *                or 0 if the edge is not on the graph
   *
   * @tparam @e is of type Edge
   * @pre @e is an Edge belonging to this graph, that is,
   * e.node1().index() < num_nodes, e.node2().index() < num_nodes() and 
   * has_edge(e.node1(), e.node2()) = true
   * @post There no longer exists said this edge e
   *        num_edges() is reduced by 1.
   */
  size_type remove_edge(const Edge& e) {
    Node n1 = e.node1();
    Node n2 = e.node2();
    return remove_edge(n1, n2);
  }

  /** Delete an edge of this graph
   * @param[in]     e_it         an edge iterator from this graph
   * @return        edge_iterator, an iterator from the graph without this node
   *                with distance(edge_begin(), return) 
   *                        = distance(edge_begin(), e_it) or the next value
   *                if such that (*return).node1() <= (*return).node2()
   *
   * @tparam @e_it is of type EdgeIterator
   * @pre @e_it is an EdgeIterator belonging to this graph, that is,
   * (*e_it).node1().index() < num_nodes, (*e_it).node2().index() < num_nodes() and 
   * has_edge((*e_it).node1(), (*e_it).node2()) = true
   * @post There no longer exists said this edge (*e_it)
   *       num_edges() is reduced by 1.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type ind = remove_edge(*e_it);
    if (nr_edges == 0) {
      return edge_end();
    }
    else {
      if ((*e_it).node1().index() > (*e_it).node2().index()) {
        ++e_it;
      }
      return e_it;
    }
  }

  /** Delete a node from this graph
   * @param[in]     n         a node
   * @return        size_type the index of the replaced edge (n.index())
   *                or 0 if the node is not on the graph
   *
   * @tparam @n is of type Node
   * @pre @n is a Node belonging to this graph, that is,
   * n.index() < num_nodes. 
   * @post There no longer exists said this Node n
   *        However, it will be replaced by the last node according to the
   *        current order, that is, if m.index() = num_nodes() - 1 before,
   *        then m.index() = a if a was n.index() before the removal. Also,
   *        all edges e with e.node1().index() = a or e.node2().index() = a 
   *        are removed from the graph.
   *        Also, if num_nodes() = c before removal, num_nodes() = c-1 now and 
   *         num_edges() is rediced by n.degree()
   */
  size_type remove_node(const Node& n) {
    size_type ind = n.index();
    if (ind < nr_nodes) {
      // find the edges to be removed
      size_type id_compl = id_exist_to_compl[ind];
      std::vector<size_type> adj_ids = node_edges[id_compl];
      for (unsigned int i_adj = 0; i_adj < adj_ids.size(); i_adj++) {
        size_type id_other = id_compl_to_exist[adj_ids[i_adj]];
        // remove them (ignoring which node is node 1 or node 2)
        remove_edge(node(ind), node(id_other));
        remove_edge(node(adj_ids[i_adj]), node(id_compl));
      }

      if (nr_nodes > 1) {
        // swap index
        iter_swap(points.begin()+ind, points.begin()+(nr_nodes-1));
        iter_swap(values.begin()+ind, values.begin()+(nr_nodes-1));
      }
      points.pop_back();
      values.pop_back();

      nr_nodes -= 1;

      //size_type ind_compl_curr = id_exist_to_compl[ind];
      size_type id_exist_end = nr_nodes;
      size_type id_compl_end = id_exist_to_compl[id_exist_end];

      id_exist_to_compl[ind] = id_compl_end;
      id_compl_to_exist[id_compl_end] = ind;

      return ind;
    }
    else {
      return 0;
    }
  }

  /** Delete a node from this graph
   * @param[in]     n_it        a node iterator frorm this graph with
   *                            (*n_it).graph = this
   * @return        node_iterator the index of the replaced edge (n.index())
   *                or 0 if the node is not on the graph
   *
   * @tparam @n is of type Node
   * @pre @n is a Node belonging to this graph, that is,
   * (*n_it).index() < num_nodes. 
   * @post There no longer exists said this Node (*n_it)
   *        However, it will be replaced by the last node according to the
   *        current order, that is, if m.index() = num_nodes() - 1 before,
   *        then m.index() = a if a was n.index() before the removal 
   *        and (*return) = m.
   *        Also, all edges e with e.node1().index() = a or 
   *        e.node2().index() = a are removed from the graph.
   *        Also, if num_nodes() = c before removal, num_nodes() = c-1 now and 
   *         num_edges() is rediced by n.degree()
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    if (nr_nodes == 0) {
      return node_end();
    }
    else {
      return n_it;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // set all variables back to default
    points = {};
    values = {};
    edg_by_ord = {};
    edg_by_idx = {};
    node_edges = {};
    edge_value_vec = {};
    id_compl_to_exist = {};
    id_exist_to_compl = {};
    nr_nodes = 0;
    nr_edges = 0;
    nr_nodes_created = 0;
  }

  void check() {
    std::cout << 1;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Construct valid Node Iterator */
    NodeIterator(Graph* g) : graph(g){
      // check if graph is nullptr
      if (g == nullptr) {
        idx = 0;
      }
      // intialized by node_begin() function
      else {
        idx = 0;
        // check if we have nodes
        // otherwise, this should be equal to end()
        if ((*graph).nr_nodes == 0) {
          graph = nullptr;
        }
      }
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return the node object */
    Node operator*() const {
      return (graph->node(idx));
    }

    /** Iterate to the next Node */
    NodeIterator& operator++() {
      // next node exists
      if ((idx + 1) < (*graph).nr_nodes) {
        idx += 1;
      }
      // no next node
      // go to default
      else {
        idx = 0;
        graph = nullptr;
      }
      // return the incremented iterator
      return (*this);
    }

    // iterator is equal if it has same graph pointer, 
    //  same index (and same number of nodes)
    bool operator==(const NodeIterator& n) const {
      bool check1 = (idx == n.idx);
      bool check2 = (graph == n.graph);

      return (check1 & check2);
    }


   private:
    friend class Graph;
    Graph* graph;
    size_type idx;
    

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Initialize NodeIterator at the first instance */
  NodeIterator node_begin() {
    return NodeIterator(this);
  }

  /** Initialize NodeIterator's end */
  NodeIterator node_end() {
    return NodeIterator(nullptr);
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
    IncidentIterator() {}

    /** Construct a valid IncidentIterator. */
    IncidentIterator(Graph* g, size_type id) : graph(g), this_idx(id) {
      // end of INcidentIterator is initialized
      if (g == nullptr) {
        other_idx = 0;
      }
      // initialization by edge_begin() function
      else {
        other_idx = 0;

        // if there are no adjacent ondes, iterator is set to end
        if ((*graph).node_edges[this_idx].size() == 0) {
          graph = nullptr;
          this_idx  = 0;
        }
      }
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return the Node element */
    Edge operator*() const {
      return (Edge(graph, this_idx, (*graph).node_edges[this_idx][other_idx]));
    }

    /** Iterate to next element */
    IncidentIterator& operator++() {
      if ((other_idx + 1) < (*graph).node_edges[this_idx].size()) {
        other_idx += 1;
      }
      // reached the end
      else {
        graph = nullptr;
        this_idx = 0;
        other_idx = 0;
      }

      return (*this);
    }

    /** check for equality
     * equality holds if all IncidentIterator variables are equal 
     */
    bool operator==(const IncidentIterator& inc) const {
      bool cond1 = (graph == inc.graph);
      bool cond2 = (this_idx == inc.this_idx);
      bool cond3 = (other_idx == inc.other_idx);

      return (cond1 & cond2 & cond3);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph;
    size_type this_idx;
    size_type other_idx;

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
    EdgeIterator() {}

    /** Construct a valid EdgeIterator. */
    EdgeIterator(graph_type* g, bool is_end) : 
        graph(g), node_iter((*graph).node_end()), 
        inc_iter(Node(graph, 0).edge_end()) {

      // we initialize an end iterator   
      if (is_end) {
        // construct end of NodeIterator
        node_type help_node = Node();
        node_iter = (*graph).node_end();
        // construct end of EdgeIterator
        inc_iter = help_node.edge_end();
        graph = nullptr;
      }
      else {
        if (((*graph).nr_nodes > 0) & ((*graph).nr_edges > 0)) {
          // start node iterator
          node_iter = (*graph).node_begin();
          // iterate unitil we found a node with edges
          bool searching = true;
          while (searching) {
            if ((*node_iter).degree() > 0) {
              // found such a node
              // initialize EdgeIterator
              searching = false;
              inc_iter = (*node_iter).edge_begin();
            }
            else {
              ++node_iter;
            }
          }
        }
        // no nodes or no edges in graphin graph
        // thus initialize an end iterator
        else {
          node_type help_node = Node();
          node_iter = (*graph).node_end();
          inc_iter = help_node.edge_end();
          graph = nullptr;
        }
      }
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Return the currrent edge */
    Edge operator*() const {
      return (*inc_iter);
    }

    /** Increment the EdgeIterator */
    EdgeIterator& operator++() {
      bool new_node = false;

      // increment the IncidentIterator until we found a new node
      ++inc_iter;
      for (; node_iter != (*graph).node_end(); ++node_iter) {
        if (new_node) {
          // start a new IncidentIterator for the new Node
          inc_iter = (*node_iter).edge_begin();
        }
        // iterate through all new incident edges
        for (; inc_iter != (*node_iter).edge_end(); ++inc_iter) {
          // only select edges where the smaller node index occurs first
          if ((*inc_iter).node1().index() <= (*inc_iter).node2().index()) {
            return *this;
          }
        }
        new_node = true;
      }
      graph = nullptr;
      return *this;
    }

    /** Check for equality of EdgeIterators
     * Equality of same graph pointer, same NodeIterator
     * and same IncidentIterator
     */
    bool operator==(const EdgeIterator& e) const {
      bool cond1 = (graph == e.graph);
      bool cond2 = (node_iter == e.node_iter);
      bool cond3 = (inc_iter == e.inc_iter);

      return (cond1 & cond2 & cond3);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph;
    NodeIterator node_iter;
    IncidentIterator inc_iter;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** start a new EdgeIterator */
  edge_iterator edge_begin() {
    return EdgeIterator(this, false);
  }

  /** initialize the end of an EdgeIterator */
  edge_iterator edge_end() {
    return EdgeIterator(this, true);
  }

 private:

  std::vector<Point> points;
  std::vector<node_value_type> values;
  std::vector<std::vector<size_type>> edg_by_ord;
  std::map<size_type, std::set<size_type>> edg_by_idx;
  std::vector<std::vector<size_type>>  node_edges;
  std::vector<std::map<size_type, edge_value_type>> edge_value_vec;
  std::vector<size_type> id_compl_to_exist;
  std::vector<size_type> id_exist_to_compl;
  size_type nr_nodes;
  size_type nr_edges;
  size_type nr_nodes_created;

};

#endif // CME212_GRAPH_HPP