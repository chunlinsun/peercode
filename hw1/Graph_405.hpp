#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <utility>
#include <cassert>
#include <map>
#include <unordered_map>

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
        // HW0: YOUR CODE HERE
        // Use this space for declarations of important internal types you need
        // later in the Graph's definition.
        // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
        // code here. Just use the space if you need it.)

    public:

        //
        // PUBLIC TYPE DEFINITIONS
        //

        /** Type of this graph. */
        using graph_type = Graph;

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

        // My defs:
        using um = std::unordered_map<size_type, size_type>;
        using mum = std::map<size_type, um>;
        using ii = std::pair<size_type, size_type>;
        using vii = std::vector<ii>;
        using vp = std::vector<Point>;

        //
        // CONSTRUCTORS AND DESTRUCTOR
        //

        /** Construct an empty graph. */
        Graph() {
            // HW0: YOUR CODE HERE
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
            private:


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
                Node() { }

                /** Return this node's position. */
                const Point& position() const {
                    return g_->points[idx];
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    return idx;
                }

                node_value_type& value() {
                    return const_cast<node_value_type&>(g_->vals.at(idx));
                }

                const node_value_type& value() const {
                    return g_->vals.at(idx);
                }
//--functionality_1
//--Compilation error: Node has no member adj_list. once you fix this, make sure
//--it works even if the node has degree 0.
//--START
                size_type degree() const {
                    return adj_list[idx].size();
                }
//--END

                IncidentIterator edge_begin() const {
                    return IncidentIterator(g_, g_->adj_list.at(idx).begin(), idx);
                }

                IncidentIterator edge_end() const {
                    return IncidentIterator(g_, g_->adj_list.at(idx).end(), idx);
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    return g_ == n.g_ and n.idx == idx;
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
                    if (g_ == n.g_)
                        return idx < n.idx;
                    return g_ < n.g_;
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;

                Node(Graph const *g, size_type i) {
                    this->g_ = g;
                    this->idx = i;
                }

                Graph const *g_;
                size_type idx;
        };

        /** Return the number of nodes in the graph.
         *
         * Complexity: O(1).
         */
        size_type size() const {
            return points.size();
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
        Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
            points.push_back(position);
            vals.push_back(val);
            return Node(this, points.size()-1);
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const {
            return this == n.g_;
        }

        /** Return the node with index @a i.
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         */
        Node node(size_type i) const {
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
        class Edge : private totally_ordered<Edge> {
            public:
                /** Construct an invalid Edge. */
                Edge() { }

                /** Return a node of this Edge */
                // Complexity: O(1)
                Node node1() const {
                    if (flip)
                        return Node(g_, g_->edges[idx].second); 
                    else
                        return Node(g_, g_->edges[idx].first); 
                }

                /** Return the other node of this Edge */
                // Complexity: O(1)
                Node node2() const {
                    if (flip)
                        return Node(g_, g_->edges[idx].first); 
                    else
                        return Node(g_, g_->edges[idx].second); 
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const {
                    return g_ == e.g_ and idx == e.idx;
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const {
                    if (g_ == e.g_) {
                        return idx < e.idx;
                    }
                    return g_ < e.g_;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;

                Edge(Graph const *g, size_type i, bool flip = false) {
                    this->g_ = g;
                    this->idx = i;
                    this->flip = flip;
                }

                Graph const *g_;
                size_type idx;
                bool flip;
        };

        /** Return the total number of edges in the graph.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        // Complexity: O(1)
        size_type num_edges() const {
            return edges.size();
        }

        /** Return the edge with index @a i.
         * @pre 0 <= @a i < num_edges()
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        // Complexity: O(1)
        Edge edge(size_type i) const {
            return Edge(this, i);
        }

        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        // Complexity: log(num_nodes())
//--functionality_0
//--Compilation error: lookup_a returns an iterator to adj_list but you're comparing
//--it to an iterator of edges. Once fixed, note that *lookup_a returns a pair, so
//--you need to grab the second element (the unordered map) to look for b.idx.
//--START
        bool has_edge(const Node& a, const Node& b) const {
            auto lookup_a = adj_list.find(a.idx);
            if (lookup_a != edges.end()) {
                return lookup_a->count(b.idx);
            } else 
                return false;
        }
//--END

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
        // Complexity: log(num_nodes())
        Edge add_edge(const Node& a, const Node& b) {
            auto lookup_a = adj_list.find(a.idx); // Check if we have encountered this node before
            if (lookup_a != adj_list.end()) {
                // In case we have seen a before, check if is connected to b
                auto lookup_b = lookup_a->second.find(b.idx);
                if (lookup_b != lookup_a->second.end()) {
                    return Edge(this, lookup_b->second); // Edge already exists, return it
                } else {
                    // Edge doesn't exist, add it, using our previous query for optimization
                    edges.push_back(make_edge_pair(a, b));
                    lookup_a->second.insert({b.idx, edges.size()-1});
                    adj_list[b.idx].insert({a.idx, edges.size()-1});
                }
            } else {
                // Edge doesn't exist, add it
                edges.push_back(make_edge_pair(a, b));
                adj_list[a.idx].insert({b.idx, edges.size()-1});
                adj_list[b.idx].insert({a.idx, edges.size()-1});
            }
            return Edge(this, edges.size()-1);
        }

        /** Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() {
             edges.clear();
             adj_list.clear();
             points.clear();
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
                NodeIterator() { }

//--documentation_2
//--No documentation for HW1 Graph methods
//--START
                Node operator*() const {
                    return Node(g_, idx);
                }

                NodeIterator& operator++() {
                    idx++;
                    return *this;
                }

                bool operator==(const NodeIterator& it) const {
                    return g_ == it.g_ and idx == it.idx;
                }

                bool operator!=(const NodeIterator& iter) const {
                    return not (*this == iter);
                }
//--END
            private:
                friend class Graph;

                NodeIterator(Graph const* g, size_type idx) {
                    this->g_ = g;
                    this->idx = idx;

                }

                Graph const *g_;
                size_type idx;

        };

        NodeIterator node_begin() const {
            return NodeIterator(this, 0);
        }

        NodeIterator node_end() const {
            return NodeIterator(this, points.size());
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
                using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid IncidentIterator. */
                IncidentIterator() { }

                Edge operator*() const {
                    if (from < it->first) // make sure the orientation is correct
                        return Edge(g_, it->second);
                    else
                        return Edge(g_, it->second, true);
                }

                IncidentIterator& operator++() {
                    it++;
                    return *this;
                }

                bool operator==(const IncidentIterator& inc_it) const {
                    return it == inc_it.it;
                }

                bool operator!=(const IncidentIterator& inc_it) const {
                    return it != inc_it.it;
                }

            private:
                friend class Graph;
                
                Graph const *g_;
                um::const_iterator it;
                size_type from;

                IncidentIterator(Graph const* g, um::const_iterator it, size_type from) {
                    this->g_ = g;
                    this->it = it;
                    this->from = from;
                }
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
                EdgeIterator() { }

                Edge operator*() const {
                    return Edge(g_, idx);
                }

                EdgeIterator& operator++() {
                    idx++;
                    return *this;
                }

                bool operator==(const EdgeIterator& it) const {
                    return g_ == it.g_ and idx == it.idx;
                }

                bool operator!=(const EdgeIterator& it) const {
                    return not (*this == it);
                }

            private:
                friend class Graph;

                Graph const* g_;
                size_type idx;

                EdgeIterator(Graph const* g, size_type idx) {
                    this->g_ = g;
                    this->idx = idx;
                }

        };

        EdgeIterator edge_begin() const {
            return EdgeIterator(this, 0);
        }
        EdgeIterator edge_end() const {
            return EdgeIterator(this, edges.size());
        }

    private:

        // Creates a pair of indices with the lesser one first
        ii make_edge_pair(const Node& a, const Node& b) const {
            return std::make_pair(std::min(a, b).idx, std::max(a, b).idx);
        }

        vp points; 
        vii edges;
        std::vector<node_value_type> vals;
        mum adj_list;
};

#endif // CME212_GRAPH_HPP
