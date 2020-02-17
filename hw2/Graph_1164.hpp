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
template <typename V, typename E>
class Graph {
    private:
        // HW0: YOUR CODE HERE
        // Use this space for declarations of important internal types you need
        // later in the Graph's definition.
        // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
        // code here. Just use the space if you need it.)
        
        typedef V node_value_type;
        typedef E edge_value_type;

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
        Graph() {}

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
                Node() { }

                /** Return this node's position. */
                Point& position() {
                    // return const_cast<Point&>(g_->points.at(idx));
                    return const_cast<Point&>(g_->get_point(uid));
                }

                /** Return this node's position. */
                const Point& position() const {
                    // return g_->points.at(idx);
                    return g_->get_point(uid);
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    // return idx;
                    return g_->n_uids[uid];
                }

                node_value_type& value() {
                    // return const_cast<node_value_type&>(g_->n_vals.at(idx));
                    return const_cast<node_value_type&>(g_->get_node_val(uid));
                }

                const node_value_type& value() const {
                    return g_->get_node_val(uid);
                }

                size_type degree() const {
                    try {
                        return adj_list[uid].size();
                    } catch(std::out_of_range e) {
                        return size_type{0};
                    }
                }

                IncidentIterator edge_begin() const {
                    return IncidentIterator(g_, g_->adj_list.at(uid).begin(), uid);
                }

                IncidentIterator edge_end() const {
                    return IncidentIterator(g_, g_->adj_list.at(uid).end(), uid);
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    return g_ == n.g_ and n.uid == uid;
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
                        return uid < n.uid;
                    return g_ < n.g_;
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;

                Node(Graph const *g, size_type i) {
                    this->g_ = g;
                    this->uid = i;
                }

                Graph const *g_;
                size_type uid;
        };
        
        /** Remove a node from the Graph
         *  @post g.has_node(n) -> false
         *  @post g.node_num() one fewer than @pre
         *  @post all edges that are connected are removed as well
         *  @post iterations past node_end()-1 are invalidated
         *
         * Complexity: n.num_adj() log(num_nodes())
         */
        size_type remove_node(const Node& n) {
            if (has_node(n)) {
                size_type idx = n.index(); // index of element to be removed
                auto adj_nodes = adj_list.find(n.uid);
                if (adj_nodes != adj_list.end()) { // Remove all edges connected to the node
                    std::vector<Edge> to_rm;
                    for (auto& val : adj_nodes->second) { 
                        to_rm.push_back(Edge(this, val.second));
                    }
                    while (not to_rm.empty()) {
                        remove_edge(to_rm.back());
                        to_rm.pop_back();
                    }
                    adj_list.erase(n.uid); // Remove the node entry itself
                }
                size_type uid_of_last_el = n_idxs.back();
                // update the index of the element about to be moved in its place 
                n_uids[uid_of_last_el] = idx;
                swap_n_pop(points, idx);
                swap_n_pop(n_vals, idx);
                swap_n_pop(n_idxs, idx);
                return idx;
            }
            return 0;
        }

        /** Remove a node by reference to its iterator
         *  @post g.has_node(n) -> false
         *  @post g.node_num() one fewer
         *  @post all edges that are connected are removed as well
         *  @post iterations past node_end()-1 are invalidated
         *
         * Complexity: n.num_adj() log(num_nodes())
         */
        node_iterator remove_node(node_iterator n_it) {
            remove_node(*n_it);
            return n_it;
        }

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
            n_vals.push_back(val);
            n_uids.push_back(points.size()-1);
            n_idxs.push_back(n_uids.size()-1);
            return Node(this, n_uids.size()-1);
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const {
            return this == n.g_ and  0 <= n.index() and n.index() < num_nodes();
        }

        /** Return the node with index @a i.
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         */
        Node node(size_type idx) const {
            assert(0<= idx and idx < num_nodes());
            return Node(this, n_idxs[idx]);
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

                edge_value_type& value() {
                    // return const_cast<edge_value_type&>(g_->e_vals.at(idx));
                    return const_cast<edge_value_type&>(g_->get_edge_val(uid));
                }

                const edge_value_type& value() const {
                    // return g_->e_vals.at(idx);
                    return g_->g_->get_edge_val(uid);
                }

                /** Return a node of this Edge */
                // Complexity: O(1)
                Node node1() const {
                    if (flip)
                        // return Node(g_, g_->edges[idx].second); 
                        return Node(g_, g_->get_pair(uid).second); 
                    else
                        // return Node(g_, g_->edges[idx].first); 
                        return Node(g_, g_->get_pair(uid).first); 
                }

                /** Return the other node of this Edge */
                // Complexity: O(1)
                Node node2() const {
                    if (flip)
                        // return Node(g_, g_->edges[idx].first); 
                        return Node(g_, g_->get_pair(uid).first); 
                    else
                        // return Node(g_, g_->edges[idx].second); 
                        return Node(g_, g_->get_pair(uid).second); 
                }

                double length() const {
                    return norm(node1().position() - node2().position());
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const {
                    return g_ == e.g_ and uid == e.uid;
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const {
                    if (g_ == e.g_) {
                        return uid < e.uid;
                    }
                    return g_ < e.g_;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;

                Edge(Graph const *g, size_type i, bool flip = false) {
                    this->g_ = g;
                    this->uid = i;
                    this->flip = flip;
                }

                Graph const *g_;
                size_type uid;
                bool flip;
        };

        /** Remove an edge between two nodes
         *
         * @post num_edge() one less than @pre
         * @post edge_end()-1 invalidated from @pre
         *
         * Complexity: log(num_nodes())
         */
        size_type remove_edge(const Node& a, const Node& b) {
            if (has_edge(a, b)) {
                size_type uid = adj_list[a.uid][b.uid];
                size_type idx = e_uids[uid];
                // Remove the adj list entry
                adj_list[a.uid].erase(b.uid);
                adj_list[b.uid].erase(a.uid);
                size_type uid_of_last_el = e_idxs.back();
                // update the index of the element about to be moved in its place 
                e_uids[uid_of_last_el] = idx;
                swap_n_pop(edges, idx);
                swap_n_pop(e_vals, idx);
                swap_n_pop(e_idxs, idx);
                return idx;
            }
            return 0;
        }

        /** Remove an edge
         *
         * Complexity: log(num_edges())
         */
        size_type remove_edge(const Edge& e) {
            return remove_edge(e.node1(), e.node2());
        }

        edge_iterator remove_edge(edge_iterator e_it) {
            remove_edge(*e_it);
            return e_it;
        }

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
            return Edge(this, e_idxs[i]);
        }

        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        // Complexity: log(num_nodes())
        bool has_edge(const Node& a, const Node& b) const {
            auto lookup_a = adj_list.find(a.uid);
            if (lookup_a != adj_list.end()) {
                return lookup_a->second.count(b.uid);
            } else 
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
        // Complexity: log(num_nodes())
        Edge add_edge(const Node& a, const Node& b) {
            auto lookup_a = adj_list.find(a.uid); // Check if we have encountered this node before
            if (lookup_a != adj_list.end()) {
                // In case we have seen a before, check if is connected to b
                auto lookup_b = lookup_a->second.find(b.uid);
                if (lookup_b != lookup_a->second.end()) {
                    return Edge(this, lookup_b->second); // Edge already exists, return it
                } else {
                    // Edge doesn't exist, add it, using our previous query for optimization
                    edges.push_back(make_edge_pair(a, b));
                    e_vals.push_back({});
                    e_uids.push_back(edges.size()-1);
                    e_idxs.push_back(e_uids.size()-1);
                    lookup_a->second.insert({b.uid, e_uids.size()-1});
                    adj_list[b.uid].insert({a.uid, e_uids.size()-1});
                }
            } else {
                // Edge doesn't exist, add it
                edges.push_back(make_edge_pair(a, b));
                e_vals.push_back({});
                e_uids.push_back(edges.size()-1);
                e_idxs.push_back(e_uids.size()-1);
                adj_list[a.uid].insert({b.uid, e_uids.size()-1});
                adj_list[b.uid].insert({a.uid, e_uids.size()-1});
            }
            return Edge(this, e_uids.size()-1);
        }

        /** Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() {
             edges.clear();
             n_idxs.clear();
             n_vals.clear();
             e_idxs.clear();
             e_vals.clear();
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

                // Dereference the iterator
                Node operator*() const {
                    return Node(g_, *it);
                }

                // Increment the iterator
                NodeIterator& operator++() {
                    ++it;
                    return *this;
                }

                bool operator==(const NodeIterator& n_it) const {
                    return g_ == n_it.g_ and it == n_it.it;
                }

                bool operator!=(const NodeIterator& n_it) const {
                    return not (*this == n_it);
                }

            private:
                friend class Graph;

                NodeIterator(Graph const* g, bool begin=true) {
                    g_ = g;
                    if (begin)
                        it = g_->n_idxs.begin();
                    else
                        it = g_->n_idxs.end();
                }

                Graph const *g_;
                std::vector<size_type>::const_iterator it;

        };

        // Returns an iterator over all the nodes, pointing to the first node
        // Iterates over the nodes in an undefined order
        NodeIterator node_begin() const {
            return NodeIterator(this);
        }

        // Retunrs an iterator over all the nodes, pointing after the last node
        // Iterates over the nodes in an undefined order
        NodeIterator node_end() const {
            return NodeIterator(this, false);
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
                    return Edge(g_, *it);
                }

                EdgeIterator& operator++() {
                    ++it;
                    return *this;
                }

                bool operator==(const EdgeIterator& e_it) const {
                    return g_ == e_it.g_ and it == e_it.it;
                }

                bool operator!=(const EdgeIterator& e_it) const {
                    return not (*this == e_it);
                }

            private:
                friend class Graph;

                Graph const* g_;
                std::vector<size_type>::const_iterator it;

                EdgeIterator(Graph const* g, bool begin=true) {
                    g_ = g;
                    if (begin)
                        it = g_->e_idxs.begin();
                    else
                        it = g_->e_idxs.end();
                }

        };

        // Returns an iterator over all the edges, pointing to the first edge
        // Iterates over the nodes in an undefined order
        EdgeIterator edge_begin() const {
            return EdgeIterator(this);
        }

        // Return an iterator over all the edges, pointing after the last edge
        // Iterates over the nodes in an undefined order
        EdgeIterator edge_end() const {
            return EdgeIterator(this, false);
        }

    private:

        // Creates a pair of indices with the lesser one first
        ii make_edge_pair(const Node& a, const Node& b) const {
            return std::make_pair(std::min(a, b).uid, std::max(a, b).uid);
        }

        template <typename D>
            void swap_n_pop(std::vector<D>& v, typename std::vector<D>::const_iterator i) {
                std::iter_swap(v.end()-1, i);
                v.pop_back();
            }

        template <typename D>
            void swap_n_pop(std::vector<D>& v, size_type i) {
                std::iter_swap(v.end()-1, v.begin() + i);
                v.pop_back();
            }

        ii& get_pair(size_type uid) {
            return const_cast<ii&>(edges.at(e_uids[uid]));
        }

        const ii& get_pair(size_type uid) const {
            return edges.at(e_uids[uid]);
        }

        edge_value_type& get_edge_val(size_type uid) {
            return const_cast<edge_value_type&>(e_vals.at(e_uids[uid]));
        }

        const edge_value_type& get_edge_val(size_type uid) const {
            return e_vals.at(e_uids[uid]);
        }

        Point& get_point(size_type uid) {
            return const_cast<Point&>(points.at(n_uids[uid]));
        }

        const Point& get_point(size_type uid) const {
            return points.at(n_uids[uid]);
        }

        node_value_type& get_node_val(size_type uid) {
            return const_cast<node_value_type&>(n_vals.at(n_uids[uid]));
        }

        const node_value_type& get_node_val(size_type uid) const {
            return n_vals.at(n_uids[uid]);
        }

        vp points; 
        vii edges;
        std::vector<size_type> e_uids; // uid -> idx
        std::vector<size_type> e_idxs; // idx -> uid
        std::vector<size_type> n_uids; // uid -> idx
        std::vector<size_type> n_idxs; // idx -> uid
        std::vector<node_value_type> n_vals;
        std::vector<edge_value_type> e_vals;
        mum adj_list;
};

#endif // CME212_GRAPH_HPP
