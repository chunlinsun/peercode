#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

/* Uncomment when compile*/
//#include "CME212/Util.hpp"
//#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

//--functionality_1
//--No shortest_path/subgraph
//--END

//--style_0
//--variable names are very descriptive but perhaps longer than they need to be
//--END

// HW1:
template <typename V>
class Graph {
    private:
        // HW0:
        // Use this space for declarations of important internal types you need
        // later in the Graph's definition.
        // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
        // code here. Just use the space if you need it.)

        //Predeclare the internal structs
        struct internal_node;
        struct internal_edge; //Re-think this
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
        // HW1:
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

        // HW0:
        /**
         * @brief Construct an empty graph.
         *
         * By default, the graph has zero edges; graph by definition can have zero edges.
         */
        Graph():
            globalNodeSet(),
            allNodeAdjacencies(),
            globalEdgeSet() {
        }

        /** Default destructor
         *
         * @brief Deletes graph object.
         *
         */
        // !!Re-examine the complexity of this
        ~Graph() = default;
        /* ~Graph(){
         *      delete[] globalNodeset;
         * }
         */


        //
        // NODES
        //

        /** @class Graph::Node
         * @brief Class implementation of a _Node_ within the _Graph_ class.
         *
         * Node objects are used to access information about the Graph's nodes.
         */
        // HW1:
        class Node : private totally_ordered<Node>{
            public:
                /** Construct an invalid node.
                 *
                 * Valid nodes are obtained from the _Graph_ class, but it
                 * is occasionally useful to declare an invalid _node_, and assign a
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
                 *
                 * We let this declaration remain an invalid node. _Node()_ has a private definition
                 * for work specifically in this class.
                 */
                Node() {
                // HW0:
                }

                /** Return this node's position.
                 *
                 * The _globalNodeSet_ is a vector of _interal_node_ structs.
                 * A _size_type_ value is assigned as the _nodeUID_ attribute of each _internal_node_.
                 * _Position_ is an attribute of _internal_node_, inheirited from _Point.hpp_.
                 *
                 */
                const Point& position() const {
                // HW0:
                    return graph->globalNodeSet[nodeUID].position;
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                /**
                 * @brief _Index()_ is built into the _Node()_ struct as the _nodeUID_ attribute.
                 * @return _NodeUID_ is an attribute of our _Node()_ class and based on range _[0, graph_size)_.
                 */
                size_type index() const {
                // HW0:
                    return nodeUID;
                }
                // HW1: YOUR CODE HERE
                // Supply definitions AND SPECIFICATIONS for:
                // node_value_type& value();
                // const node_value_type& value() const;
                // size_type degree() const;
                // incident_iterator edge_begin() const;
                // incident_iterator edge_end() const;
                
                // Return date for user-specified value type
                node_value_type& value() {
                      return graph->globalNodeSet[nodeUID].value;
                }
                const node_value_type& value() const {
                      return graph->globalNodeSet[nodeUID].value;
                }
                // return the number of incident edges
                size_type degree() const {
                      return graph->allNodeAdjacencies[nodeUID].size();
                }
                // start the incident iterator
                incident_iterator edge_begin() const {
                      return IncidentIterator(graph, nodeUID, 0);
                }
                // End of incident iterator
                incident_iterator edge_end() const{
                      return IncidentIterator(graph, nodeUID, degree());
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                // HW0:
                    return (n.graph == graph) && (n.nodeUID == nodeUID);
                }

                /** Test whether this node is less than @a n in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * _std::map<>_. It need not have any geometric meaning.
                 *
                 * The node ordering relation must obey trichotomy: For any two nodes _x_
                 * and _y_, exactly one of x == y, x < y, and y < x is true.
                 *
                 * Because a node is assigned an ID that matches the _graph.size()_
                 * at the time the node is added, each node within a given graph
                 * has a unique identity. This should satisfy the trichotomy requirement.
                 */
                bool operator<(const Node& n) const {
                  // HW0:
                    if (n.graph != graph) return false;
                    return nodeUID < n.nodeUID;
                }
            // !!End of Node Class Public Members

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;
                // HW0:
                // Use this space to declare private data members and methods for Node
                // that will not be visible to users, but may be useful within Graph.
                // i.e. Graph needs a way to construct valid Node objects
                /**
                 * @brief Private:Node() Our _Node()_ has _graph_ and _nodeUID_ attributes.
                 * We assign a _pointer_ attribute to associate the _node_ to a specific graph.
                 * We assign a unique ID to the node of _size_type_ value.
                 */
                Graph* graph;
                size_type nodeUID;
                Node(const Graph* g, size_type nID):
                    graph(const_cast<Graph*>(g)),
                    nodeUID(nID){
                }
            // !!End of Node Class Private Members
        };
        // !!End of Node Class

        /** Return the number of nodes in the graph.
         *
         * Complexity: O(1).
         *
         */
        size_type size() const {
            // HW0:
            /**
             * @brief Returns number of nodes in the graph.
             *
             * We meet the _size()_ complexity requirement of O(1)
             * since we use vector to define our set of nodes.
             * _Vector.size()_ has constant complexity.
             */
            return globalNodeSet.size();
        }

        // !!Possibly delete unless required for future assignments
        /** Synonym for size(). */
        size_type num_nodes() const {
            return size();
        }

        // HW1:
        /** Add a node to the graph, returning the added node.
         * @param[in] position The new node's position
         * @post new num_nodes() == old num_nodes() + 1
         * @post result_node.index() == old num_nodes()
         *
         * Complexity: O(1) amortized operations
         *
         * A node is added to the graph object using _.push_back_.
         *
         * There is no need to re-write a _size_ attribute,
         * as the above post-condition comment suggests, because we access
         * the graph's size using a member function, which is based on _vector.size()_.
         *
         * We update the vector of vectors that describes's each node's adjacencies;
         * We use _.push_back_ to add an empty vector here.
         */
        Node add_node(const Point& position,
                      const node_value_type& value = node_value_type()) {
            // HW0:
            /**
             * @brief newNodeUID The new node's ID is based on _size()_.
             *
             * Assigning the new node's unique ID as the _graph.size()_
             * works fine because our _globalNodeSet_ vector starts at index 0
             * and _size()_ starts counting from 1.
             */
            size_type newNodeUID = globalNodeSet.size();
            internal_node newNode = internal_node(position, newNodeUID, value);
            globalNodeSet.push_back(newNode);
            allNodeAdjacencies.push_back(std::vector<size_type>());

            return Node(this, newNodeUID);
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const {
            // HW0:
            return n.graph == this;
        }

        /** Return the node with index @a i.
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         *
         * Perforn some basic user filtering before carrying out return.
         */
        Node node(size_type i) const {
            // HW0:
            assert((i >= 0) && (i < num_nodes()));
            return Node(this, i);
        }

        //
        // EDGES
        //

        // HW1:
        /** @class Graph::Edge
         * @brief Implementation of _edges_ for a graph.
         *
         * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
         * are considered equal if they connect the same nodes, in either order.
         */
        class Edge : private totally_ordered<Edge>{
            public:
                /** Construct an invalid Edge. */
                Edge() {
                // HW0:
                }

                /** Return a node of this Edge */
                Node node1() const {
                    // HW0:
                    return Node(graph, node1ID);
                }

                /** Return the other node of this Edge */
                Node node2() const {
                    // HW0:
                    return Node(graph, node2ID);
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 *
                 * We first assess that the edges belong to the same graph.
                 *
                 */
                bool operator==(const Edge& e) const {
                    // HW0:
                    if (e.graph != graph) return false;
                    return (e.node1ID == node2ID && e.node2ID == node1ID)
                           ||
                           (e.node1ID == node1ID && e.node2ID == node2ID);
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 *
                 * This member function was difficult to implement in a meaningful way.
                 * I initially wanted to compare the edges by their index in the _globalEdgeSet_
                 * vector. But Edge _e_ may come from a different graph with a separate _globalEdgeSet_.
                 *
                 * I deferred to comparing edges by the index of their first node.
                 * And if the first nodes are equal for the two edges, then the second nodes are
                 * compared.
                 */
                bool operator<(const Edge& e) const {
                  // HW0:
                    // !!Re-examine this
                    if (e.node1ID == node1ID){
                        return node2ID  < e.node2ID;
                    }
                    return node1ID < e.node1ID;
                }
            // !!End of Edge Class Public Members

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;
                // HW0:
                // Use this space to declare private data members and methods for Edge
                // that will not be visible to users, but may be useful within Graph.
                // i.e. Graph needs a way to construct valid Edge objects

                /**
                 * @brief Internal attributes for defining the _edge_ class.
                 *
                 * The _Edge_ class consists of attributes:
                 * _graph_ _pointer_, indicating which graph the edge belongs to,
                 * and _node1ID_ and _node2ID_ of _size_type_ value to indicate the edge nodes.
                 */
                Graph* graph;
                size_type node1ID;
                size_type node2ID;
                Edge(const Graph* g, size_type n1ID, size_type n2ID):
                    graph(const_cast<Graph*>(g)),
                    node1ID(n1ID),
                    node2ID(n2ID){
                    }
            // !!End of Edge Class Private Members
        };
        // !!End of Edge Class

        /** Return the total number of edges in the graph.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        size_type num_edges() const {
            // HW0:
            return globalEdgeSet.size();
        }

        /** Return the edge with index @a i.
         * @pre 0 <= @a i < num_edges()
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        Edge edge(size_type i) const {
            // HW0:
            // !! Re-examine this
            internal_edge nodesOfEdge = globalEdgeSet[i];
            return Edge(this, nodesOfEdge.node1, nodesOfEdge.node2);
        }

        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        bool has_edge(const Node& a, const Node& b) const {
            // HW0:
            std::vector<size_type> nodeEdgeSet = allNodeAdjacencies[a.nodeUID];
            for(size_type i = 0; i < nodeEdgeSet.size(); i++){
                if (nodeEdgeSet[i] == b.nodeUID) return true;
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
         *
         * We explicitly collect the pairings of nodes by creating a
         * _internal_edge_ struct, which is then added to the _globalEdgeSet_ vector attribute.
         *
         * The node with a lower _size_type_ value is listed as the first node in the pairing.
         */
        Edge add_edge(const Node& a, const Node& b) {
        // HW0:
        // !! Re-examine this
            if (has_edge(a, b) == true) return Edge(this, a.nodeUID, b.nodeUID);

            allNodeAdjacencies[a.nodeUID].push_back(b.nodeUID);
            allNodeAdjacencies[b.nodeUID].push_back(a.nodeUID);

            internal_edge nodeIDs = internal_edge(a.nodeUID, b.nodeUID);
            if (b.nodeUID < a.nodeUID){
                nodeIDs = internal_edge(b.nodeUID, a.nodeUID);
            }

            globalEdgeSet.push_back(nodeIDs);
            return Edge(this, a.nodeUID, b.nodeUID);
        }

        /** Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() {
        // HW0:
            allNodeAdjacencies.clear();
            globalEdgeSet.clear();
            globalNodeSet.clear();
        }

        //
        // Node Iterator
        //

        // HW1:
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
                // Supply definitions AND SPECIFICATIONS for:
                // Node operator*() const
                // NodeIterator& operator++()
                // bool operator==(const NodeIterator&) const
//--documentation_1
//--docs should be in doxygen style
//--START
                /** Dereference a node through the iterator pointer */
                Node operator*() const{
                    return graph->node(nodeIterIndx);
                }
                /** Increments to the iterator pointer to the next node */
                NodeIterator& operator++(){
                    nodeIterIndx += 1;
                    return (*this);
                }
                /** Assesses equality between two iterators */
                bool operator==(const NodeIterator& otherNodeIterator) const {
                    if (otherNodeIterator.graph != graph) return false;
                    return otherNodeIterator.nodeIterIndx == nodeIterIndx;
                }
//--END
            private:
                friend class Graph;
                // HW1 #2: YOUR CODE HERE
                Graph* graph;
                size_type nodeIterIndx;
                NodeIterator(const Graph* graph, const size_type startingIndex = 0):
                    graph(const_cast<Graph*>(graph)),
                    nodeIterIndx(startingIndex){
                }
        };


        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // node_iterator node_begin() const
        // node_iterator node_end() const

        /** Our index is based on size_type, so we can begin and end at
         * at length of our _globalNodeSet_ vector.
         */
        node_iterator node_begin() const{
            return NodeIterator(this, 0);
        }
        node_iterator node_end() const{
            return NodeIterator(this, size());
        }

        //
        // Incident Iterator
        //
        // HW1:
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
                }

                // HW1 #3: YOUR CODE HERE
                // Supply definitions AND SPECIFICATIONS for:
                // Edge operator*() const
                // IncidentIterator& operator++()
                // bool operator==(const IncidentIterator&) const

                // Dereference edge for current index
                Edge operator*() const{
                    size_type nextEdge = graph->allNodeAdjacencies[thisNodeID][localEdgesIndex];
                    return Edge(graph, thisNodeID, nextEdge);
                }
                // Increment our iterator to the next incident edge
                IncidentIterator& operator++(){
                    localEdgesIndex += 1;
                    return (*this);
                }
                // Assess equality between two incident iterators
                bool operator==(const IncidentIterator& otherIncidentIterator) const{
                    if (otherIncidentIterator.graph != graph) return false;
                    if (otherIncidentIterator.thisNodeID != thisNodeID) return false;
                    return otherIncidentIterator.localEdgesIndex == localEdgesIndex;
                }
            private:
                friend class Graph;
                // HW1 #3: YOUR CODE HERE
                Graph* graph;
                size_type thisNodeID;
                size_type localEdgesIndex;
                IncidentIterator(const Graph* g,
                                 const size_type nID,
                                 const size_type startingIndex = 0):
                    graph(const_cast<Graph*>(g)),
                    thisNodeID(nID),
                    localEdgesIndex(startingIndex){
                }
        };

        //
        // Edge Iterator
        //

        // HW1:
//--functionality_2
//--since you have your edges in a vector, the hint is not necessary - you can implement
//--EdgeIterator in the same way as NodeIterator. Using the hint is perfectly fine too;
//--it would probably be easiest to keep a NodeIterator and IncidentIterator as data
//--members and step through those, while skipping (i,j) if i>j (and also being careful
//--at the end).
//--START
        // Couldn't quite figure out how to combine NodeIterator and IncidentIterator to create this iterator
        /** @class Graph::EdgeIterator
         * @brief Iterator class for edges. A forward iterator. */
        class EdgeIterator : private totally_ordered<EdgeIterator>{
//--END
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

                /**
                 * This dereferencing assumes some order based on how these edges were built as opposed
                 * to an order based on their nodes' properties such as their node_type_value
                 */
                Edge operator*() const{
                    // what are the fringe cases for this? re-check how you're building allNodeAdjacencies
                    size_type adjacent = graph->allNodeAdjacencies[currentNodeID][currentPathDegree];
                    return Edge(graph, currentNodeID, adjacent);
                }

                /*** !!This needs work; needs to account for non-direction; be alert on complexity O(n)
                // Increment to the next edge
                EdgeIterator& operator++(){
                    currentPathDegree += 1;
                    for(int i = 0; i < graph->globalNodeSet.size(); i++){
                        size_type rootDegree =  graph->globalNodeSet[i];
                        for ( ; currentPathDegree != rootDegree; currentPathDegree++){
                            size_type nextNode = graph->allNodeAdjacencies[currentNodeID][currentPathDegree];
                            if (currentNodeID < nextNode) {return *this;}
                        }
                        currentPathDegree = 0;
                    }
                    currentPathDegree = 0;
                    return (*this);
                }
                **/
                // Assess for equality of edge iterators. I'm unsure where this would come up (??)
                bool operator==(const EdgeIterator& otherEdgeIterator) const{
                    if (graph != otherEdgeIterator.graph) return false;
                    if (currentNodeID != otherEdgeIterator.currentNodeID) return false;
                    return currentPathDegree == otherEdgeIterator.currentPathDegree;
                }
            private:
                friend class Graph;
                // HW1 #5: YOUR CODE HERE
                /**
                 * We require these three attributes because we need to identify that two edges are
                 * of the (1) same graph for _==_ operator, (2) that an edge has a given first point,
                 * i.e. currentNodeID, and
                 *
                 */
                Graph *graph;
                size_type currentNodeID;
                size_type currentPathDegree;
                EdgeIterator(const Graph* g,
                             const size_type nID = 0,
                             const size_type startingIndex = 0):
                    graph(const_cast<Graph*>(g)),
                    currentNodeID(nID),
                    currentPathDegree(startingIndex){
                }
        };

        // HW1 #5: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // edge_iterator edge_begin() const
        // edge_iterator edge_end() const
        /** We base our edge iterator on the size of our _globalEdgeSet_,
         * so the iterator can span from 0 to this size.
         */
        edge_iterator edge_begin() const{
            return EdgeIterator(this, 0, 0);
        }
        edge_iterator edge_end() const{
            return EdgeIterator(this, size(),0);
        }
        // !!End of Graph Class public members

    private:
        // HW0:
        // HW1:
        // Use this space for your Graph class's internals:
        //   helper functions, data members, and so forth.
        /**
         * @brief The internal structs are based on the provided struct in proxy_example.cpp.
         *
         * _Point_ is inheirited from _Point.hpp_.
         * A _size_type_ value provides the unique identification of each _internal_node_ (_nodeUID_).
         * A _node_value_type_ provides user-specified info for the node.
         */
        struct internal_node{
            Point position;
            size_type nodeUID;
            node_value_type value;  // HW1:
            internal_node (const Point& pos,
                           const size_type nUID,
                           node_value_type v):
                position(pos),
                nodeUID(nUID),
                value(v){
            }
        };

        struct internal_edge{
            size_type node1;
            size_type node2;
            internal_edge (const size_type n1, const size_type n2):
                node1(n1),
                node2(n2){
            }
        };

        std::vector<internal_node> globalNodeSet;
        std::vector<std::vector<size_type>> allNodeAdjacencies;
        std::vector<internal_edge> globalEdgeSet;
    // !!End of Graph Class Private members
};
// !!End of Graph Class

#endif // CME212_GRAPH_HPP
