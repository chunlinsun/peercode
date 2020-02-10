#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/**
 * @file Graph.hpp
 *
 * @author Chih-Hsuan (Carolyn) Kao
 * Contact: chkao831@stanford.edu
 * Date: Jan 31, 2020
 *
 * @brief This file, Graph.hpp, contains an undirected graph type, class Graph.
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <vector>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/**
 * @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {

    private:
        struct NodeInternal;

    public:

        /** PUBLIC TYPE DEFINITIONS */

        /** Type of this graph. */
        using graph_type = Graph;

        /** Predeclaration of Node type. */
        class Node;
        /** Synonym for Node (following STL conventions). */
        using node_type = Node;
        /** Synonym for Node Value Type as Graph is a template */
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

        /** Invalid node index value */
        static size_type constexpr invalid_index = size_type(-1);

        /** CONSTRUCTORS AND DESTRUCTOR of Graph */

        /** Construct an empty graph. */
        Graph()
            : vec_of_points_(),
            num_of_nodes_(0),
            map_of_edges_(),
            num_of_edges_(0),
            smallernode_to_edgeidx_(),
            map_of_adjnodes_(),
            node_internals_(){
            }

        /** Default destructor */
        ~Graph() = default;

        //
        // NODE
        //
        /**
         * @class Graph::Node
         * @brief Class representing the graph's nodes.
         *
         * Node objects are used to access information about the Graph's nodes.
         */
        class Node : private totally_ordered<Node>{

            public:
                /*
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
                /** Construct an invalid node. */
                Node()
                    : ptr_to_graph_(nullptr),
                    node_idx_(invalid_index){
                    }

                /**
                 * @brief Return this node's position.
                 */
                const Point& position() const {
                    return ptr_to_graph_-> vec_of_points_[this->node_idx_];
                }

                /**
                 * @brief Return this node's index, a number in the range [0, graph_size).
                 */
                size_type index() const {
                    return node_idx_;
                }

                /**
                 * @brief Return ref of (or set) this node's value (from internal struct)
                 * @return the value of this node
                 */
                node_value_type& value() {
                    return ptr_to_graph_->node_internals_[node_idx_].val_;
                }

                /**
                 * @brief Return ref of this node's value from internal struct; value is unchangeable
                 * @return the value of this node
                 */
                const node_value_type& value() const {
                    return ptr_to_graph_->node_internals_[node_idx_].val_;
                }

                /**
                 * @brief Return the num of edges incident to this node
                 * @return degree of this node (connecting edges)
                 */
                size_type degree() const {
                    return ptr_to_graph_-> map_of_adjnodes_.at(node_idx_).size();
                }

                /**
                 * @brief For this node's edges, return first @a incident_iterator
                 * @return @a incident_iterator corresponding to first edge of this node
                 */
                IncidentIterator edge_begin() const {
                    return IncidentIterator(ptr_to_graph_, node_idx_);
                }

                /**
                 * @brief For this node's edges, return last @a incident_iterator
                 * @return @a incident_iterator corresponding to pass-the end edge of this node
                 */
                IncidentIterator edge_end() const{
                    return IncidentIterator(ptr_to_graph_, node_idx_, degree());
                }

                /**
                 * @brief Test whether this node and @a n are equal.
                 * @param[in] @a n  Node to test
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    return (n.ptr_to_graph_ == this->ptr_to_graph_) \
                                             && (n.node_idx_ == this->node_idx_);
                }

                /**
                 * @brief Test whether this node is less than @a n in a global order.
                 * @param[in] @a n  Node to test
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any geometric meaning.
                 *
                 * The node ordering relation must obey trichotomy: For any two nodes x
                 * and y, exactly one of x == y, x < y, and y < x is true.
                 */
                bool operator<(const Node& n) const {
                    if(this->node_idx_ == n.node_idx_ && this->ptr_to_graph_ != n.ptr_to_graph_){
                        return std::less<Graph*>{}(this->ptr_to_graph_, n.ptr_to_graph_);
                    } else {
                        return this->node_idx_ < n.node_idx_;
                    }
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;

                /** private attributes of Node class*/
                // pointer to Graph
                Graph* ptr_to_graph_;
                //the index of this node
                size_type node_idx_;

                /** Node's private constructor */
                Node(const Graph* graph, size_type nodeid)
                    : ptr_to_graph_(const_cast<Graph*>(graph)),
                    node_idx_(nodeid){
                    }
        }; //end Node class

        /** public methods of Graph related to Node */
        /**
         * @brief Return the number of nodes in the graph.
         * Complexity: O(1).
         */
        size_type size() const {
            return num_of_nodes_;
        }

        /** Synonym for size(). */
        size_type num_nodes() const {
            return size();
        }

        /**
         * @brief Add a node to the graph, returning the added node.
         * @param[in] position The new node's position
         * @param[in] value  Value of the adding node
         * @post new num_nodes() == old num_nodes() + 1
         * @post result_node.index() == old num_nodes()
         *
         * Complexity: O(1) amortized operations.
         */
        Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
            //the pre-number-of-nodes is the index of newly added node
            size_type new_index = num_of_nodes_ ;

            //add point to vector of points for Graph
            vec_of_points_.push_back(position);
            //key is node index; current value is empty vector (later pushing its neighbor nodes)
            map_of_adjnodes_.insert(std::pair<size_type,std::vector<size_type> >(new_index,
                        std::vector<size_type>()));
            //For internal struct, create a node internal with corresponding node index, position, and value type
            NodeInternal newInternal = NodeInternal(new_index, position, value);
            node_internals_.push_back(newInternal);

            //increment number of nodes
            num_of_nodes_+=1;

            return Node(this,new_index);
        }

        /**
         * @brief Determine if a Node belongs to this Graph
         * @param[in] @a n  Node to test
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const {
            return (n.node_idx_ < num_of_nodes_ && n.ptr_to_graph_== this);
        }

        /**
         * @brief Return the node with index @a i.
         * @param[in] @a i  Index of the node
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         */
        Node node(size_type i) const {
            //precondition: the passed-in index is within range
            assert((i>=0) && (i<num_nodes()));
            return Node(this,i);
        }

        //
        // EDGE
        //
        /**
         * @class Graph::Edge
         * @brief Class representing the graph's edges.
         *
         * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
         * are considered equal if they connect the same nodes, in either order.
         */
        class Edge : private totally_ordered<Edge>{

            public:
                /** Construct an invalid Edge. */
                Edge()
                    : ptr_to_graph_(nullptr),
                    edge_idx_(invalid_index),
                    node1_idx_(invalid_index),
                    node2_idx_(invalid_index){
                    }

                /** Return a node of this Edge */
                Node node1() const {
                    return ptr_to_graph_ -> node(node1_idx_);
                }

                /** Return the other node of this Edge */
                Node node2() const {
                    return ptr_to_graph_ -> node(node2_idx_);
                }

                /**
                 * @brief Test whether this edge and @a e are equal.
                 * @param[in] @a e  Edge to test
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const {
                    bool same_nodes = (e.node1_idx_ == node1_idx_ && e.node2_idx_ == node2_idx_);
                    bool same_graphs = e.ptr_to_graph_ == ptr_to_graph_;
                    //return true if connecting same nodes and pointing to same graph
                    return (same_nodes && same_graphs);
                }

                /**
                 * @brief Test whether this edge is less than @a e in a global order.
                 * @param[in] @a e  Edge to test
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const {
                    if (this->edge_idx_ == e.edge_idx_ && this->ptr_to_graph_ != e.ptr_to_graph_) {
                        return std::less<Graph*>{}(this->ptr_to_graph_, e.ptr_to_graph_);
                    } else {
                        return this->edge_idx_ < e.edge_idx_;
                    }
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;

                /** private attributes for Edge */
                //pointer to graph
                Graph* ptr_to_graph_;
                //edge index
                size_type edge_idx_;
                //first node (smaller) index
                size_type node1_idx_;
                //second node (bigger) index
                size_type node2_idx_;

                /** private constructor for Edge */
                Edge(const Graph* graph,
                        size_type edgeindex,
                        size_type node1index,
                        size_type node2index)
                    : ptr_to_graph_(const_cast<Graph*>(graph)),
                    edge_idx_(edgeindex),
                    node1_idx_(node1index),
                    node2_idx_(node2index){
                    }
        };//end Edge

        /** public methods of Graph related to Edge */
        /**
         * @brief Return the total number of edges in the graph.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        size_type num_edges() const {
            return num_of_edges_;
        }

        /**
         * @brief Return the edge with index @a i.
         * @param[in] @a i  Index of this edge
         * @pre 0 <= @a i < num_edges()
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        Edge edge(size_type i) const {
            //precondition: passed-in index within range
            assert((i>=0) && (i<num_of_edges_));
            //checked within range (edge is valid), then check content
            //get index of smaller and bigger nodes
            size_type index1 = map_of_edges_.at(i)[0];
            size_type index2 = map_of_edges_.at(i)[1];

            return Edge(this,i,index1,index2);
        }

        /**
         * @brief Test whether two nodes are connected by an edge.
         * @param[in] @a a  First node of this edge
         * @param[in] @a a  Second node of this edge
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * In my design, the nodes are checked in ascending order (smaller to bigger)
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        bool has_edge(const Node& a, const Node& b) const {
            //precondition: nodes are valid
            assert(has_node(a));
            assert(has_node(b));

            //two cases: if the passed-in arguments are in ascending order
            if(a.index() < b.index()){
                //check map of map, from smaller node to bigger node and to edge
                if(smallernode_to_edgeidx_.count(a.index())){
                    return smallernode_to_edgeidx_.at(a.index()).count(b.index());
                } else {
                    return false;
                }
                //if the passed-in arguments are in descending order
            } else if(b.index() < a.index()){
                //check map of map, from smaller node to bigger node and to edge
                if(smallernode_to_edgeidx_.count(b.index())){
                    return smallernode_to_edgeidx_.at(b.index()).count(a.index());
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }

        /**
         * @brief Add an edge to the graph, or return the current edge if it already exists.
         * @param[in] @a a  First node of this edge
         * @param[in] @a a  Second node of this edge
         * @pre @a a and @a b are distinct valid nodes of this graph
         * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
         * @post has_edge(@a a, @a b) == true
         * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
         *       Else,                        new num_edges() == old num_edges() + 1.
         *
         * Can invalidate edge indexes -- in other words, old edge(@a i) might not
         * equal new edge(@a i). Must not invalidate outstanding Edge objects.
         *
         * In my design, the passed-in nodes are added to structures in ascending order.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        Edge add_edge(const Node& a, const Node& b) {

            //precondition: nodes are valid
            assert(has_node(a));
            assert(has_node(b));
            //precondition: nodes are distinct
            assert(!(a==b));

            //nodes passed-in in ascending order, no re-ordering
            if(a.index() < b.index()){
                //edge already added, then simply return this edge
                if(has_edge(a,b)){
                    return Edge(this,
                            smallernode_to_edgeidx_.at(a.index()).at(b.index()),
                            a.index(),
                            b.index());
                }

                //otherwise, capture new index
                size_type new_index = num_of_edges_;

                //add to map (key is edge; val is pair (vec) of node indices)
                map_of_edges_[new_index] = {a.index(), b.index()};
                //add to map of map (key is smaller node index; val is bigger node index;
                //                   inner val is edge index)
                smallernode_to_edgeidx_[a.index()][b.index()] = num_of_edges_;

                //add to map of adjacency in both direction for two nodes
                map_of_adjnodes_[a.index()].push_back(b.index());
                map_of_adjnodes_[b.index()].push_back(a.index());

                //increment number of edges
                num_of_edges_ += 1;

                return Edge(this,num_of_edges_-1,a.index(),b.index());

            } else { //nodes passed-in in descending order, re-ordering needed
                //edge already added, then simply return this edge
                if(has_edge(b,a)){
                    return Edge(this,
                            smallernode_to_edgeidx_.at(b.index()).at(a.index()),
                            b.index(),
                            a.index());
                }

                //otherwise, capture new index
                size_type new_index = num_of_edges_;

                //add to map (key is edge; val is pair (vec) of node indices)
                map_of_edges_[new_index] = {b.index(), a.index()};
                //add to map of map (key is smaller node index; val is bigger node index;
                //                   inner val is edge index)
                smallernode_to_edgeidx_[b.index()][a.index()] = num_of_edges_;

                //add to map of adjacency in both direction for two nodes
                map_of_adjnodes_[b.index()].push_back(a.index());
                map_of_adjnodes_[a.index()].push_back(b.index());

                //increment number of edges
                num_of_edges_ += 1;

                return Edge(this,num_of_edges_-1,b.index(),a.index());
            }
        }

        /**
         * @brief Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() {
            num_of_nodes_ = 0;
            num_of_edges_ = 0;
            vec_of_points_.clear();
            map_of_edges_.clear();
            smallernode_to_edgeidx_.clear();
            map_of_adjnodes_.clear();
        }

        //
        // Node Iterator
        //
        /**
         * @class Graph::NodeIterator
         * @brief Iterator class for nodes. A forward iterator.
         */
        class NodeIterator: private totally_ordered<NodeIterator>{

            public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Node;                     // Element type
                using pointer           = Node*;                    // Pointers to elements
                using reference         = Node&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid NodeIterator. */
                NodeIterator()
                    : ptr_to_graph_(nullptr),
                    nodeiter_currentidx_(invalid_index){
                    }

                /**
                 * @brief A dereferencing operator for @a node_iterator
                 * @return the value of underlying node of this iterator
                 */
                Node operator*() const {
                    //using node function under Graph (returning a node)
                    return ptr_to_graph_->node(nodeiter_currentidx_);
                }

                /**
                 * @brief An incrementing operator for @node_iterator to the next position
                 * @return this iterator that wraps to the next underlying element
                 */
                NodeIterator& operator++() {
                    nodeiter_currentidx_++;
                    return *this;
                }

                /**
                 * @brief An equality operator for @a node_iterator to test if two iters are equal
                 * @param @a nodeiter node_iterator to be compared with this
                 * @return true if two are pointing to the same graph and have same current index
                 */
                bool operator==(const NodeIterator& nodeiter) const {
                    bool same_graph = (nodeiter.ptr_to_graph_ == ptr_to_graph_);
                    bool same_currentindex = (nodeiter.nodeiter_currentidx_ == nodeiter_currentidx_);
                    return (same_graph && same_currentindex);
                }

            private:
                friend class Graph;

                Graph* ptr_to_graph_;
                size_type nodeiter_currentidx_;

                /** Nodeiterator  private constructor */
                NodeIterator(const Graph* graph, const size_type nodeiterator_idx = 0)
                    : ptr_to_graph_(const_cast<Graph*>(graph)),
                    nodeiter_currentidx_(nodeiterator_idx) {
                    }
        }; //end class NodeIterator

        /** public methods of Graph related to node iterator */
        /** @brief Return the node iterator at the first starting position */
        NodeIterator node_begin() const {
            return NodeIterator(this);
        }

        /** @brief Return the node iterator at the pass-the-end position */
        NodeIterator node_end() const {
            return NodeIterator(this, num_nodes());
        }

        //
        // Incident Iterator
        //
        /**
         * @class Graph::IncidentIterator
         * @brief Iterator class for edges incident to a node. A forward iterator.
         */
        class IncidentIterator: private totally_ordered<IncidentIterator>{

            public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Edge;                     // Element type
                using pointer           = Edge*;                    // Pointers to elements
                using reference         = Edge&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid IncidentIterator. */
                IncidentIterator()
                    : ptr_to_graph_(nullptr),
                    root_node_idx_(invalid_index),
                    num_visited_edges_(invalid_index){
                    }

                /**
                 * @brief A dereferencing operator for @a incident_iterator
                 * @return the value of underlying edge of this iterator
                 */
                Edge operator*() const {
                    //firstly capture the adjacent node relative to the root node
                    size_type adjnode_idx = ptr_to_graph_->map_of_adjnodes_.at(root_node_idx_)[num_visited_edges_];
                    //declare a variable that captures the connecting edge index of root & adjacent node
                    size_type connecting_edge_index;
                    //if rootnode < adjacentnode
                    if(root_node_idx_ < adjnode_idx){
                        //get edge value from the map of map (which has ascending order)
                        connecting_edge_index = ptr_to_graph_->smallernode_to_edgeidx_.at(root_node_idx_).at(adjnode_idx);
                    } else { //otherwise, adjacentnode < rootnode
                        //get edge value from the map of map (which has ascending order)
                        connecting_edge_index = ptr_to_graph_->smallernode_to_edgeidx_.at(adjnode_idx).at(root_node_idx_);
                    }
                    return Edge(ptr_to_graph_, connecting_edge_index, root_node_idx_, adjnode_idx);
                }

                /**
                 * @brief An incrementing operator for @incident_iterator to the next neighbor node
                 * @return this iterator that wraps to the next underlying element
                 */
                IncidentIterator& operator++() {
                    num_visited_edges_++;
                    return *this;
                }

                /**
                 * @brief An equality operator for @a incident_iterator to test if two iters are equal
                 * @param @a inc_iter incident_iterator to be compared with this
                 * @return true if two are pointing to the same graph, have same root note, and visited same neighbors
                 */
                bool operator==(const IncidentIterator& inc_iter) const {
                    bool same_graph = inc_iter.ptr_to_graph_ == ptr_to_graph_;
                    bool same_rootnode = inc_iter.root_node_idx_ == root_node_idx_;
                    bool same_visited_edges = inc_iter.num_visited_edges_ == num_visited_edges_;
                    return (same_graph && same_rootnode && same_visited_edges);
                }

            private:
                friend class Graph;

                Graph* ptr_to_graph_;
                size_type root_node_idx_;
                size_type num_visited_edges_;

                /** Incident Iterator private constructor */
                IncidentIterator(const Graph* graph,
                        const size_type root_node_id,
                        const size_type num_visited_edges = 0)
                    : ptr_to_graph_(const_cast<Graph*>(graph)),
                    root_node_idx_(root_node_id),
                    num_visited_edges_(num_visited_edges) {
                    }
        };//end IncidentIterator class

        //
        // Edge Iterator
        //
        /**
         * @class Graph::EdgeIterator
         * @brief Iterator class for edges. A forward iterator.
         */
        class EdgeIterator : private totally_ordered<EdgeIterator> {

            public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Edge;                     // Element type
                using pointer           = Edge*;                    // Pointers to elements
                using reference         = Edge&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid EdgeIterator. */
                EdgeIterator()
                    : ptr_to_graph_(nullptr),
                    current_edge_idx_(invalid_index){
                    }

                /**
                 * @brief A dereferencing operator for @a edge_iterator
                 * @return the value of underlying edge of this iterator
                 */
                Edge operator*() const{
                    size_type minornode = ptr_to_graph_->map_of_edges_.at(current_edge_idx_)[0];
                    size_type majornode = ptr_to_graph_->map_of_edges_.at(current_edge_idx_)[1];
                    return Edge(ptr_to_graph_,current_edge_idx_,minornode,majornode);
                }

                /**
                 * @brief An incrementing operator for @a edge_iterator to the next position
                 * @return the value of underlying edge of this iterator
                 */
                EdgeIterator& operator++() {
                    current_edge_idx_++;
                    return *this;
                }

                /**
                 * @brief An equality operator for @a edge_iterator to test if two iters are equal
                 * @param @a edge_iter edge_iterator to be compared with this
                 * @return true if two are pointing to the same graph and at same edge
                 */
                bool operator==(const EdgeIterator& edge_iter) const {
                    bool same_graphs = (edge_iter.ptr_to_graph_ == ptr_to_graph_);
                    bool same_current_edge = (edge_iter.current_edge_idx_ == current_edge_idx_);
                    return same_graphs && same_current_edge;
                }

            private:
                friend class Graph;

                Graph* ptr_to_graph_;
                size_type current_edge_idx_;

                /** Edge Iterator private constructor */
                EdgeIterator(const Graph* graph,
                        const size_type current_edge = 0)
                    : ptr_to_graph_(const_cast<Graph*>(graph)),
                    current_edge_idx_(current_edge){
                    }

        };//end EdgeIterator class

        /** public methods of Graph related to edge iterator */
        /** @brief Return the edge iterator at the first starting position */
        EdgeIterator edge_begin() const {
            return EdgeIterator(this);
        }
        /** @brief Return the edge iterator at the pass-the-end position */
        EdgeIterator edge_end() const {
            return EdgeIterator(this, num_edges());
        }

    /** PRIVATE BLOCK FOR GRAPH */
    private:

        /** Internal struct for node */
        struct NodeInternal {
            size_type node_idx_;
            Point points_;
            node_value_type val_;

            //private constructor
            NodeInternal(const size_type node_id,
                    const Point& position,
                    node_value_type val)
                : node_idx_(node_id),
                points_(position),
                val_(val){
                }
        };//end NodeInternal

        /** Other private attributes for Graph */
        //vector of points
        std::vector<Point> vec_of_points_;
        //number of nodes
        size_type num_of_nodes_;
        //keys are edge indices, value is a vector of 2 elements (smaller and bigger node indices)
        std::unordered_map< size_type, std::vector<size_type> > map_of_edges_;
        //number of edges
        size_type num_of_edges_;
        //map of map: smaller node index is key; bigger node index is value; edge index is inner value
        std::unordered_map<size_type, std::unordered_map<size_type, size_type>> smallernode_to_edgeidx_;
        //get all neighbor nodes by root node index
        std::unordered_map<size_type, std::vector<size_type> > map_of_adjnodes_;

        // Vector containing the graph nodes
        std::vector<NodeInternal> node_internals_;

};//end Graph class

//--functionality_0
//--Great job!
//--END
#endif // CME212_GRAPH_HPP


