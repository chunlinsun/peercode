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
template <typename V>
class Graph {
 private:

  struct internalNode;

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
  /** Synonym for Node value type. */
  using node_value_type = V;

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
  Graph(): node_vec_(), neighbors_() {}

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

    /** Return this node's position. */
    const Point& position() const {
      return graph_->node_vec_[idx_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        return idx_;
    }

    // HW1: YOUR CODE HERE
    /** @brief Access the value of a node
     * @return the node's value as a reference. (can be changed)
     */
    node_value_type& value() {
        return graph_->node_vec_[idx_].val_;
    }

    /** @brief Access the value of a node
     * @return the node's value as a reference. (cannot be changed)
     */
    const node_value_type& value() const {
        return graph_->node_vec_[idx_].val_;
    }

    /** @brief Get the number of edges incident to this node
     * @return the number of incident edges of this node
     */
    size_type degree() const {
        return graph_->neighbors_[idx_].size();
    }

   /** @brief Start of the incident iterator
    * @return the incident iterator corresponding to the first edge visited
    */
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_, idx_, 0);
    }

    /** @brief End of the incident iterator
     * @return the incident iterator where the edge traversal stops
     */
    incident_iterator edge_end() const{
        return IncidentIterator(graph_, idx_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        if(graph_ == n.graph_ and idx_ == n.idx_){
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
        if(graph_ == n.graph_ and idx_ < n.idx_) {
            return true;
        }
        else if (graph_ < n.graph_) {
            return true;
        }
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // the graph the node belongs to
    Graph* graph_;
    // node's idx
    size_type idx_;

    Node(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {

    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_vec_.size();
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
      // new node's idx
      size_type result_node_idx = node_vec_.size();
      // construct a new node and add to node_vec_
      internalNode resultInternalNode = internalNode(position, result_node_idx, value);
      node_vec_.push_back(resultInternalNode);
      neighbors_.push_back(std::vector<size_type>());

      return Node(this, result_node_idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      if (n.graph_ == this and n.idx_ <= size()){
          return true;
      }
      return false;
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, idx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, idx2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        // first check if the graph equals
        bool graph_equal = graph_ == e.graph_;

        // then check if the two nodes equals
        bool nodes_equal = (((e.idx1_ == idx1_) && (e.idx2_ == idx2_)) || ((e.idx1_ == idx2_) && (e.idx2_ == idx1_)));

        return (graph_equal && nodes_equal);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        size_type this_min = std::min(idx1_, idx2_);
        size_type this_max = std::max(idx1_, idx2_);
        size_type e_min = std::min(e.idx1_, e.idx2_);
        size_type e_max = std::max(e.idx1_, e.idx2_);

        // two edges are in different graph
        if (graph_ < e.graph_) {
            return true;
        }
        else if (graph_ > e.graph_) {
            return false;
        }

        // two edges are in the same graph
        else {
            if (this_min < e_min) {
                return true;
            }
            else if (this_min == e_min) {
                return this_max < e_max;
            } else {
                return false;
            }
        }
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type idx1_;
    size_type idx2_;

    Edge(const Graph* graph, size_type idx1, size_type idx2):
        graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2) {

    };
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
      size_type total_edges = 0;
      for (auto n : neighbors_) {
          total_edges += n.size();
      }
      total_edges /= 2;
      return total_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
      size_type edge_idx = -1;
      for (size_type n1 = 0; n1 < neighbors_.size(); ++n1) {
          for (size_type n2 = 0; n2 < neighbors_[n1].size(); ++n2) {
              if (n1 < neighbors_[n1][n2]) {
                  ++edge_idx;
              }
              if (edge_idx == i) {
                  return Edge(this, n1, neighbors_[n1][n2]);
              }
          }
      }
  }



  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check the neighbors of node a
    for (auto n : neighbors_[a.idx_]) {
        if (n == b.idx_) {
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
    // If edge exists, return the current edge
    if (has_edge(a, b)) {
        return Edge(this, a.idx_, b.idx_);
    } else {
        // If the edge doesn't exist, create a new one, and add to graph and return
        neighbors_[a.idx_].push_back(b.idx_);
        neighbors_[b.idx_].push_back(a.idx_);

        return Edge(this, a.idx_, b.idx_);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_vec_.clear();
    neighbors_.clear();
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    /** @brief Operator for de-referencing a NodeIterator
     * @return The node the NodeIterator traverses at
     */
    Node operator*() const {
        return graph_->node(node_iter_idx_);
    }

    /** @brief Operator for increasing a NodeIterator
     * @return A NodeIterator pointing to the next node
     * @post   @node_iter_idx_ of this NodeIterator increases by 1
     */
    NodeIterator& operator++() {
        node_iter_idx_++;
        return *this;
     }

    /** @brief Operator for testing the equality between this NodeIterator and another NodeIterator @node_iter
     * @param[in] @node_iter  NodeIterator to compare with
     * @return    true iff this NodeIterator and @node_iter are in the same graph and have the same index
     * @pre       @node_iter is a valid NodeIterator
     */
    bool operator==(const NodeIterator& node_iter) const {
        bool graph_equal = graph_ == node_iter.graph_;
        bool idx_equal = node_iter_idx_ == node_iter.node_iter_idx_;
        return (graph_equal && idx_equal);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_iter_idx_;

    /** @brief Valid constructor for NodeIterator.
      * @param[in] graph              Pointer to iterator's Graph
      * @param[in] nodeiterator_idx   Starting iterator index (default is 0)
      */
    NodeIterator(const Graph* graph, const size_type node_iter_idx = 0):
        graph_(const_cast<Graph*>(graph)), node_iter_idx_(node_iter_idx) {

    }
  };

  // HW1 #2: YOUR CODE HERE
  /** @brief Start of the NodeIterator
   * @return the NodeIterator corresponding to the first node in the graph visited
   */
  NodeIterator node_begin() const {
      return NodeIterator(this, 0);
  }

  /** @brief End of the NodeIterator
   * @return the NodeIterator where the node traversal stops
   */
  NodeIterator node_end() const {
      return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>  {
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
    /** @brief Operator for de-referencing the IncidentIterator
     * @return The edge the IncidentIterator points to
     */
    Edge operator*() const {
        size_type neighbor_id = graph_->neighbors_[target_node_id_][in_edges_visited_];
        return Edge(graph_, target_node_id_, neighbor_id);
    }

    /** @brief Operator that increases this IncidentIterator
     * @return the IncidentIterator pointing to the next edge incident to the target node
     *          or the IncidentIterator with @in_edges_visited_ == degree(target node)
     * @post @in_edges_visited_ of this IncidentIterator increases by 1
     */
    IncidentIterator& operator++() {
        in_edges_visited_++;
        return *this;
    }
    /** @brief Operator for testing the equality between this IncidentIterator and another IncidentIterator @in_iter
     * @param[in] @in_iter  IncidentIterator to compare with
     * @return    true iff this IncidentIterator and @in_iter are in the same graph,
     *            have the same target node and have the same number of incident edges visited
     * @pre @in_iter is a valid IncidentIterator
     */
    bool operator==(const IncidentIterator& in_iter) const {
        bool graph_equal = graph_ == in_iter.graph_;
        bool target_node_equal = target_node_id_ == in_iter.target_node_id_;
        bool in_edges_visited_equal = in_edges_visited_ == in_iter.in_edges_visited_;

        return (graph_equal && target_node_equal && in_edges_visited_equal);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type target_node_id_;
    size_type in_edges_visited_;

    IncidentIterator(const Graph* graph, const size_type target_node_id, const size_type in_edges_visited = 0):
        graph_(const_cast<Graph*>(graph)), target_node_id_(target_node_id), in_edges_visited_(in_edges_visited) {
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
    /** @brief Operator for de-referencing an EdgeIterator
     * @return The edge currently indexed by EdgeIterator
     */
    Edge operator*() const {
        size_type neighbor_id = graph_->neighbors_[current_target_idx_][neighbors_visited_];
        return Edge(graph_, current_target_idx_, neighbor_id);

    }

    /** @brief Operator for increasing an EdgeIterator
     * @return An EdgeIterator pointing to the next EdgeIterator
     * @post   (@neighbors_visited_ of EdgeIterator increases by 1 if @current_target_idx_ != number of total nodes in this graph
     *         and @neighbors_visited_ != total neighbors of @current_target_idx_), otherwise, @neighbors_visited_ resets to 0
     */
    EdgeIterator& operator++() {
        neighbors_visited_++;
        size_type current_target_degree;
        size_type neighbor_id;
        for(; current_target_idx_ != graph_->num_nodes(); ++current_target_idx_) {
            current_target_degree = graph_->node(current_target_idx_).degree();
            // Iterate over incident edges of current target node
            for (; neighbors_visited_ != current_target_degree; ++neighbors_visited_) {
                neighbor_id = graph_->neighbors_[current_target_idx_][neighbors_visited_];

                // ensure each edge is only traversed once due to undirected graph
                if (current_target_idx_ < neighbor_id){
                    return *this;
                }
            }

            // all neighbors of current target node visited;
            // move to the next target node and reset the number of neighbors visited to be 0
            neighbors_visited_ = 0;
        }

        // all target nodes visited in the graph
        neighbors_visited_ = 0;
        return *this;
    }

   /** @brief Operator for testing the equality between this EdgeIterator and another EdgeIterator @edge_iter
    * @param[in] @edge_iter  EdgeIterator to compare with
    * @return    true iff this EdgeIterator and @edge_iter are in the same graph, have the same target node
    *            and have the same number of visited neighbors
    * @pre       @edge_iter is a valid EdgeIterator
    */
    bool operator==(const EdgeIterator& edge_iter) const {
        bool graph_equal = graph_ == edge_iter.graph_;
        bool cur_target_equal = current_target_idx_ == edge_iter.current_target_idx_;
        bool neighbor_visited_equal = neighbors_visited_ == edge_iter.neighbors_visited_;
        return (graph_equal && cur_target_equal && neighbor_visited_equal);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type current_target_idx_;
    size_type neighbors_visited_;

    EdgeIterator(const Graph* graph, const size_type current_target_idx = 0, const size_type neighbors_visited = 0) :
        graph_(const_cast<Graph*>(graph)), current_target_idx_(current_target_idx), neighbors_visited_(neighbors_visited) {}

  };


  // HW1 #5: YOUR CODE HERE
  /** @brief Start of the EdgeIterator
   * @return the EdgeIterator corresponding to the first edge in the graph visited
   */
  EdgeIterator edge_begin() const {
      return EdgeIterator(this, 0, 0);
  }

  /** @brief End of the EdgeIterator
   * @return the EdgeIterator corresponding to where the edge traversal stops
   */
  EdgeIterator edge_end() const {
      return EdgeIterator(this, num_nodes(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internalNode {
      // node position
      Point position_;
      // node index
      size_type idx_;
      // Value of node
      node_value_type val_;
      // constructor for internalNode
      internalNode(const Point& position, const size_type idx, node_value_type val):
        position_(position), idx_(idx), val_(val) {

      }

  };

    // all nodes in the graph, stored in a vector
    std::vector<internalNode> node_vec_;
    // all nodes' neighbors
    std::vector<std::vector<size_type>> neighbors_;



};

//--functionality_0
//--Great job!
//--END
#endif // CME212_GRAPH_HPP
