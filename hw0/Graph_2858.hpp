#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP
/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>


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
    // To solve the problem of re-index, we give those nodes two kind of
    // indentification, nodeid_(the position in the internalnode) and index.
    // Each time we add or del a point, we only need to add or del the point in
    // the internalnode and, then, adjust map bewteen index and nodeid.
    class internalnode;
    class internaledges;



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
            // Initilization the graph pointer
            graph_ = nullptr;
        }

        /** Return this node's position. */
        const Point& position() const {
            // HW0: YOUR CODE HERE
            if (graph_==nullptr || nodeid_> graph_->size()) {
                std::cout << nodeid_ << std::endl;
                throw std::invalid_argument( "invalid node to return position" );

            }
            return(graph_->InternalNodes[nodeid_].nodeposition_);
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            if (graph_==nullptr || nodeid_> graph_->size()) {
                std::cout << nodeid_ << std::endl;
                std::cout << "invalid node to return index" << std::endl;
                return -1;
            }
            return(graph_->InternalNodes[nodeid_].nodeindex_);
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
            (void) n;          // Quiet compiler warning
            if (n.graph_ == this->graph_  && n.nodeid_ == nodeid_){
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
            // HW0: YOUR CODE HERE
            (void) n;           // Quiet compiler warning
            if (n.graph_ == this->graph_  && graph_->InternalNodes[n.nodeid_].nodeindex_ > graph_->InternalNodes[nodeid_].nodeindex_){
                return true;
            }
            return false;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects
        // Pointer back to the Graph container
        graph_type* graph_;
        // This nodet's unique identification number
        size_type nodeid_;
        /** Private Constructor */
        Node(const graph_type * graph, size_type nodeid): graph_(const_cast<Graph*>(graph)) , nodeid_(nodeid) {
        }


    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        // HW0: YOUR CODE HERE
        return InternalNodes.size();
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
        (void) position;      // Quiet compiler warning

        NodesID.push_back(size());
        NodesIndex.push_back(size());
        InternalNodes.push_back( internalnode(position, size()));
        return Node(this,size()-1);        // Invalid node
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        (void) n;            // Quiet compiler warning
        if (n.graph_ == this  && n.nodeid_ < size()){
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
        // HW0: YOUR CODE HERE
        (void) i;             // Quiet compiler warning
        if(i<size())
        {
            return Node(this,NodesID[i]);
        }
        // Invalid node
        return Node();
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
            graph_ = nullptr;
            // HW0: YOUR CODE HERE
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            if (graph_ == nullptr)
            {
                return Node();      // Invalid Node

            }
            return (Node(graph_,graph_->InternalEdges[edgeid_].firstnodeid_));
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            if (graph_ == nullptr)
            {
                return Node();      // Invalid Node

            }
            return (Node(graph_,graph_->InternalEdges[edgeid_].secondnodeid_));
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            (void) e;           // Quiet compiler warning
            if (graph_ == e.graph_ && ( (e.node1() == node1() && e.node2()==node2()) || (e.node1() == node2() && e.node2()==node1()) ) )
            {
                return true;
            }
            return false;
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            (void) e;           // Quiet compiler warning
            if (e.graph_ == this->graph_  && graph_->InternalEdges[e.edgeid_].edgeindex_ > graph_->InternalEdges[edgeid_].edgeindex_){
                return true;
            }
            return false;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects
        graph_type* graph_;
        size_type edgeid_;
        /** Private Constructor */
        Edge(const graph_type * graph, size_type edgeid): graph_(const_cast<Graph*>(graph)) , edgeid_(edgeid) {
        }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return InternalEdges.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        (void) i;             // Quiet compiler warning
        if (i< num_edges())
        {
            return Edge(this,EdgesID[i]);
        }
        return Edge();        // Invalid Edge
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE
        (void) a; (void) b;   // Quiet compiler warning
        for(size_type i=0;i<num_edges();i++)
        {
            if( (a.nodeid_==InternalEdges[i].firstnodeid_ && b.nodeid_==InternalEdges[i].secondnodeid_)||  (b.nodeid_==InternalEdges[i].firstnodeid_ && a.nodeid_==InternalEdges[i].secondnodeid_))
            {
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
        (void) a, (void) b;   // Quiet compiler warning
        if (a==b)
        {
            return Edge();        // Invalid Edge
        }
        for(size_type  i=0;i<num_edges();i++)
        {
            if( (a.nodeid_==InternalEdges[i].firstnodeid_ && b.nodeid_==InternalEdges[i].secondnodeid_)||  (b.nodeid_==InternalEdges[i].firstnodeid_ && a.nodeid_==InternalEdges[i].secondnodeid_))
            {
                return Edge(this,EdgesID[i]);
            }
        }
        EdgesID.push_back(num_edges());
        EdgesIndex.push_back(num_edges());
        InternalEdges.push_back(internaledge(a.nodeid_,b.nodeid_,num_edges()));
        return (Edge(this,num_edges()-1));

    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        InternalNodes.clear();
        NodesIndex.clear();
        NodesID.clear();
        InternalEdges.clear();
        EdgesIndex.clear();
        EdgesID.clear();
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
    class internalnode
    {
        // variables
        Point nodeposition_;
        size_type nodeindex_;

        // constructors
        internalnode() {}
        internalnode(Point nodeposition,size_type nodeindex)
                : nodeposition_(nodeposition), nodeindex_(nodeindex) {}

        friend class Graph;
    };
    class internaledge
    {
        // variables
        size_type firstnodeid_;
        size_type secondnodeid_;
        size_type edgeindex_;

        // constructors
        internaledge() {}
        internaledge(size_type firstnodeid, size_type secondnodeid, size_type edgeindex)
                : firstnodeid_(firstnodeid), secondnodeid_(secondnodeid), edgeindex_(edgeindex) {}

        friend class Graph;
    };

    std::vector<internalnode> InternalNodes;
    std::vector<size_type> NodesIndex;
    std::vector<size_type> NodesID;

    std::vector<internaledge> InternalEdges;
    std::vector<size_type> EdgesIndex;
    std::vector<size_type> EdgesID;

};


#endif // CME212_GRAPH_HPPff00f51dff
