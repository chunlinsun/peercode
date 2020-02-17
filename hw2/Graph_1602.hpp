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
template <typename V, typename E>

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

        Point& position() {
            // HW0: YOUR CODE HERE
            if (graph_==nullptr || nodeid_> graph_->size()) {
                std::cout << nodeid_ << std::endl;
                throw std::invalid_argument( "invalid node to return a modifiable position" );
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

        /** Find the node's value.
         * @return A reference to the value of this node.
         */
        node_value_type& value() {
            if(graph_== nullptr || nodeid_> graph_->size())
            {
                throw std::invalid_argument( "invalid node to return value" );
            }
            return graph_->InternalNodes[nodeid_].nodevalue;
        }


        /** Find the node's value.
         * @return A constant reference to the value of this node.
         */
        const node_value_type& value() const {
            if(graph_== nullptr || nodeid_> graph_->size())
            {
                throw std::invalid_argument( "invalid node to return value" );
            }
            return graph_->InternalNodes[nodeid_].nodevalue;
        }

        /** Find the node's degree
        * @return degree of a node.
        *
        * Complexity: O(1).
        */

        size_type degree() const
        {
            return(graph_->IncidentEdgestoNodes[graph_->NodesIndex[nodeid_]].size());
        }

        /** Return the iterator pointing to the first element
        * @return incident_iterator last incident edge of a node.
        *
        * Complexity: O(1).
        */
        incident_iterator edge_begin() const
        {
            return IncidentIterator(graph_,graph_->NodesIndex[nodeid_],
                graph_->IncidentEdgestoNodes[graph_->NodesIndex[nodeid_]].begin());
        }

        /** Return the iterator pointing to the end
        * @return incident_iterator an end iterator
        *
        * Complexity: O(1).
        */
        incident_iterator edge_end() const
        {
            return IncidentIterator(graph_,graph_->NodesIndex[nodeid_],
                graph_->IncidentEdgestoNodes[graph_->NodesIndex[nodeid_]].end());
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            // HW0: YOUR CODE HERE

            return (n.graph_ == this->graph_  && n.nodeid_ == nodeid_);
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
            if(n.graph_ != this->graph_)
            {
                throw std::invalid_argument( "invalid node to compare, not in the same graph!" );
            }
            return ( graph_->InternalNodes[n.nodeid_].nodeindex_ > graph_->InternalNodes[nodeid_].nodeindex_);
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
    Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
        // HW0: YOUR CODE HERE
        // Here, although NodesID and NodesIndex are same, we can solve the reindex porblem by finding
        // proper mapping from NodesID and NodesIndex. For example, if we need more functions, such
        // as delete and insert, we delete a node in the middle of the vector and adjust the map between
        // NodesID and nodesIndex can fix the problem.
        NodesID.push_back(size());
        NodesIndex.push_back(size());
        InternalNodes.push_back( internalnode(position, value, size()));
        std::vector<size_type> new_node_incidentedgestonode;
        IncidentEdgestoNodes.push_back(new_node_incidentedgestonode);

        return Node(this,size()-1);        // Invalid node
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        return n.graph_ == this  && n.nodeid_ < size();
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
        if(i<size())
        {
            return Node(this,NodesID[i]);
        }
        // Invalid node
        return Node();
    }

    // HW2 section 7: nodes related removal
    /** Delete a node from the graph, and return an integer indicating whether the deletion is successful.
     * @pre @a n is a valid or invalid node of this graph.
     * @return 1 if the node is sucessfully deleted.
     *         0 if the node is invalid.
     * @post has_node(@a n) == false
     * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
     *       Else,                  new num_nodes() == old num_nodes().
     *
     * Complexity: O(num_nodes()) if we can assume the number of adjacent nodes of the node can be bounded
     * by a constant.
     */
    size_type remove_node(const Node& n){
        size_type nodeid = n.nodeid_;
        if(n.nodeid_>=num_nodes())
        {
            std::cout<< "do not have this node" <<std::endl;
            return 0;
        }
        size_type num_adj = IncidentEdgestoNodes[nodeid].size();
        for(size_type it=0;it<num_adj;)
        {
            Edge e = edge(IncidentEdgestoNodes[nodeid][it]);
            remove_edge(e);
            num_adj--;
        }
        for(size_type it=nodeid+1;it<num_nodes();it++)
        {
           InternalNodes[it].nodeindex_ = InternalNodes[it].nodeindex_-1;
        }
        for(size_type it=0;it<num_edges();++it)
        {
            if(InternalEdges[it].firstnodeid_>nodeid)
                InternalEdges[it].firstnodeid_--;
            if(InternalEdges[it].secondnodeid_>nodeid)
                InternalEdges[it].secondnodeid_--;

        }
        IncidentEdgestoNodes.erase(IncidentEdgestoNodes.begin()+nodeid);
        InternalNodes.erase(InternalNodes.begin()+nodeid);
        NodesID.pop_back();
        NodesIndex.pop_back();
        return 1;
    }
    /** Delete a node from the graph, and return an integer indicating whether the deletion is successful.
     * @pre @a n_it is a vlid node iterator of this graph.
     * @return the iterator pointing to the next node.
     * @post has_node(@a (*n_it)) == false
     * @post If old has_node(@a (*n_it)), new num_nodes() == old num_nodes() - 1.
     *       Else,                        new num_nodes() == old num_nodes().
     *
     * Complexity: O(num_nodes()) if we can assume the number of adjacent nodes of the node can be bounded
     * by a constant.
     */
    node_iterator remove_node(node_iterator n_it)
    {
        Node n = InternalNodes[(*n_it)];
        unsigned result = remove_node(n);
        if (result)
        {
            return n_it;
        }else
        {
            return NodeIterator();
        }

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
            graph_ = nullptr;
            // HW0: YOUR CODE HERE
        }

        // /** Return the length of the edge**/
        //double length() const{
        //    return graph_->InternalEdges[edgeid_].edgevalue_;
        //}

        /** Find the edge's value.
         * @return A reference to the value of this edge.
         */
        edge_value_type& value() {
            if(graph_== nullptr || edgeid_> graph_->num_edges())
            {
                throw std::invalid_argument( "invalid edge to return value" );
            }
            return graph_->InternalEdges[edgeid_].edgevalue_;
        }


        /** Find the edge's value.
         * @return A constant reference to the value of this edge.
         */
        const edge_value_type& value() const {
            if(graph_== nullptr || edgeid_> graph_->num_edges())
            {
                throw std::invalid_argument( "invalid edge to return value" );
            }
            return graph_->InternalEdges[edgeid_].edgevalue_;
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            if (graph_ == nullptr)
            {
                return Node();      // Invalid Node

            }
            return (graph_->node(node1_index));
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            if (graph_ == nullptr)
            {
                return Node();      // Invalid Node

            }
            return (graph_->node(node2_index));
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            return (graph_ == e.graph_ && ( (e.node1() == node1() &&
                e.node2()==node2()) || (e.node1() == node2() && e.node2()==node1()) ));
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            if (e.graph_ != this->graph_)
            {
                throw std::invalid_argument( "invalid edge to compare, not in the same graph!" );
            }

            return (graph_->InternalEdges[e.edgeid_].edgeindex_ > graph_->InternalEdges[edgeid_].edgeindex_);
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        //friend class IncidentIterator;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects
        graph_type* graph_;
        size_type edgeid_;
        size_type node1_index;
        size_type node2_index;
        /** Private Constructor */
        Edge(const graph_type * graph, size_type edgeid, size_type n1_index, size_type n2_index):
        graph_(const_cast<Graph*>(graph)) , edgeid_(edgeid), node1_index(n1_index), node2_index(n2_index) {
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
     * Actrully, the complexity is O(1)
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        if (i< num_edges())
        {
            unsigned edgeid = EdgesID[i];
            return Edge(this,edgeid,NodesIndex[this->InternalEdges[edgeid].firstnodeid_],
                        NodesIndex[this->InternalEdges[edgeid].secondnodeid_]);
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
        // HW0: YOUR CODE HERE// Quiet compiler warning
        for(size_type i=0;i<num_edges();i++)
        {
            if( (a.nodeid_==InternalEdges[i].firstnodeid_ && b.nodeid_==InternalEdges[i].secondnodeid_)||
                (b.nodeid_==InternalEdges[i].firstnodeid_ && a.nodeid_==InternalEdges[i].secondnodeid_))
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
        if (a==b)
        {
            return Edge();        // Invalid Edge
        }
        for(size_type  i=0;i<num_edges();i++)
        {
            if( (a.nodeid_==InternalEdges[i].firstnodeid_ && b.nodeid_==InternalEdges[i].secondnodeid_)||
                (b.nodeid_==InternalEdges[i].firstnodeid_ && a.nodeid_==InternalEdges[i].secondnodeid_))
            {
                return Edge(this,EdgesID[i],NodesIndex[a.nodeid_],NodesIndex[b.nodeid_]);
            }
        }

        EdgesID.push_back(num_edges());
        EdgesIndex.push_back(num_edges());
        InternalEdges.push_back(internaledge(a.nodeid_,b.nodeid_,num_edges()));
        // Next line is used for hw2 section3
        //InternalEdges[num_edges()-1].edgevalue_ = norm(a.position()-b.position());


        IncidentEdgestoNodes[NodesIndex[a.nodeid_]].push_back(num_edges()-1);
        IncidentEdgestoNodes[NodesIndex[b.nodeid_]].push_back(num_edges()-1);

        return Edge(this,num_edges()-1,NodesIndex[a.nodeid_],NodesIndex[b.nodeid_]);

    }

    Edge add_edge(const Node& a, const Node& b,edge_value_type edgedata) {
        // HW0: YOUR CODE HERE
        if (a==b)
        {
            return Edge();        // Invalid Edge
        }
        for(size_type  i=0;i<num_edges();i++)
        {
            if( (a.nodeid_==InternalEdges[i].firstnodeid_ && b.nodeid_==InternalEdges[i].secondnodeid_)||
                (b.nodeid_==InternalEdges[i].firstnodeid_ && a.nodeid_==InternalEdges[i].secondnodeid_))
            {
                return Edge(this,EdgesID[i],NodesIndex[a.nodeid_],NodesIndex[b.nodeid_]);
            }
        }

        EdgesID.push_back(num_edges());
        EdgesIndex.push_back(num_edges());
        InternalEdges.push_back(internaledge(a.nodeid_,b.nodeid_,num_edges(),edgedata));
        // Next line is used for hw2 section3
        //InternalEdges[num_edges()-1].edgevalue_ = norm(a.position()-b.position());


        IncidentEdgestoNodes[NodesIndex[a.nodeid_]].push_back(num_edges()-1);
        IncidentEdgestoNodes[NodesIndex[b.nodeid_]].push_back(num_edges()-1);

        return Edge(this,num_edges()-1,NodesIndex[a.nodeid_],NodesIndex[b.nodeid_]);

    }

    // HW2 section 7: edges related removal
    /** Delete an edge from the graph, and return an integer indicating whether the deletion is successful.
     * @pre @a e is a valid or invalid edge of this graph.
     * @return 1 if the edge is sucessfully deleted.
     *         0 if the edge is invalid.
     * @post has_edge(@a e.node1(), @a e.node2()) == false
     * @post has_edge(@a e.node1(), @a e.node2()) , new num_edges() == old num_edges() - 1.
     *       Else,                                  new num_edges() == old num_edges().
     *
     *
     * Complexity: O( num_edges()+num_nodes() )
     */
    size_type remove_edge(const Edge& e) {
        size_type edgeid = e.edgeid_;
        size_type number_nodes = num_nodes();
        if(edgeid>=num_edges())
        {
            std::cout << "do not have this edge" << std::endl;
            return 0;
        }
        for(size_type  it=0;it<number_nodes;++it)
        {
            size_type inci_size = IncidentEdgestoNodes[it].size();
            for(size_type in_it=0; in_it<inci_size;)
            {
                if(IncidentEdgestoNodes[it][in_it]==edgeid)
                {
                    IncidentEdgestoNodes[it].erase(IncidentEdgestoNodes[it].begin()+in_it);
                    inci_size--;
                    in_it--;
                }else if(IncidentEdgestoNodes[it][in_it]>edgeid)
                {
                    IncidentEdgestoNodes[it][in_it]--;
                }
                ++in_it;
            }
        }
        InternalEdges.erase(InternalEdges.begin()+edgeid);
        EdgesIndex.pop_back();
        EdgesID.pop_back();
        return 1;
    }
    /** Delete an edge from the graph, and return an integer indicating whether the deletion is successful.
     * @pre @a a and @a b are two valid nodes of the graph.
     * @return 1 if the edge is sucessfully deleted.
     *         0 if the edge is invalid or does not exist.
     * @post has_edge(@a a, @a b) == false
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
     *       Else,                        new num_edges() == old num_edges().
     *
     *
     * Complexity: Complexity: O( num_edges()+num_nodes() )
     */
    size_type remove_edge(const Node& a, const Node& b){
        size_type has_edge = 0;
        size_type edgeid = 0;
        for(size_type it=0;it<num_edges();it++)
        {
            if( (a.nodeid_==InternalEdges[it].firstnodeid_ && b.nodeid_==InternalEdges[it].secondnodeid_)||
                (b.nodeid_==InternalEdges[it].firstnodeid_ && a.nodeid_==InternalEdges[it].secondnodeid_))
            {
                edgeid = it;
                has_edge = 1;
                break;
            }
        }
        if (has_edge)
        {
            return remove_edge(edge(edgeid));
        }else
        {
            std::cout <<"do not have this edge" <<std::endl;
            return 0;
        }
    }

    /** Delete an edge from the graph, and return an integer indicating whether the deletion is successful.
     * @pre @a e_it is an edge iterator of this graph.
     * @return the iterator pointing to the next edge.
     * @post If old has_edge(@a (*e_it).node1(), @a (*e_it).node2()), new num_edges() == old num_edges() - 1.
     *       Else,                                                    new num_edges() == old num_edges().
     *
     *
     * Complexity: Complexity: O( num_edges()+num_nodes() )
     */
    edge_iterator remove_edge(edge_iterator e_it){
        size_type edgeid = *e_it;
        if(edgeid>=num_edges())
        {
            std::cout << "do not have this edge" << std::endl;
            return EdgeIterator();
        }else
        {
            remove_edge(edge(edgeid));
            return e_it;
        }
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
        IncidentEdgestoNodes.clear();
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
            nodeiterator = nullptr;
            graph_ = nullptr;
        }

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Node operator*() const
        // NodeIterator& operator++()
        // bool operator==(const NodeIterator&) const

        /** This function returns an Node by getting the internal Node
         * with specific index.
         * @return node for the specific edge iterator.
         *
         * Complexity: O(1).
         */
        Node operator*() const
        {
            size_type nodeid = *nodeiterator;
            return graph_->node(graph_->NodesIndex[nodeid]);
        }

        /** This function returns a reference to the Node iterator that was incremented by 1.
         * @return Edgeiterator incremented by 1.
         *
         * Complexity: O(1).
         */

        NodeIterator& operator++()
        {
            nodeiterator++;
            return *this;
        }

        /** This function checks whether two iterators are the same.
         * @param[in] new_NodeIterator This is another incident iterator that will be
         * compared to the one defined in the private attributes of the class NodeIterator.
         * @return true if the iterators are the same and false if they aren't.
         *
         * Complexity: O(1).
         */

        bool operator==(const NodeIterator& new_NodeIterator) const
        {
            return (nodeiterator == new_NodeIterator.nodeiterator);
        }



    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE


        const Graph* graph_;
        typename std::vector<size_type>::const_iterator nodeiterator;
        // Constructor
        NodeIterator(const Graph* graph, typename std::vector<size_type>::const_iterator iter)
                : graph_(const_cast<Graph*>(graph)), nodeiterator(iter) {
        }
    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_iterator node_begin() const
    // node_iterator node_end() const
    /** Return the iterator pointing to the first node in the graph
    * @return _node_iterator_ first node in the graph.
    *
    * Complexity: O(1).
    */
    node_iterator node_begin() const
    {
        return(NodeIterator(this,this->NodesID.begin()));
    }

    /** Return the iterator pointing to the end iterator in the graph
    * @return _node_iterator_ end in the graph.
    *
    * Complexity: O(1).
    */
    node_iterator node_end() const
    {
        return(NodeIterator(this,this->NodesID.end()));
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
        }

        // HW1 #3: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
        // IncidentIterator& operator++()
        // bool operator==(const IncidentIterator&) const

        /** This function returns an Edge by getting the internal Edge
         * at the specific iD.
         * @return Edge for the specific incident edge iterator.
        *
        * Complexity: O(1).
        */
        Edge operator*() const
        {
            size_type edgeindex = *incidentiterator;
            Edge currentedge = graph_->edge(edgeindex);
            unsigned aindex,bindex;
            aindex = currentedge.node1().index();
            bindex = currentedge.node2().index();
            if(aindex == nodeindex_)
            {
                return Graph::Edge(graph_,graph_->EdgesID[edgeindex],nodeindex_,bindex);
            }
            else
            {
                return Graph::Edge(graph_,graph_->EdgesID[edgeindex],nodeindex_,aindex);
            }
        }

        /** This function returns a reference to the Edge iterator that was incremented by 1.
         * @return Edgeiterator incremented by 1.
         *
         * Complexity: O(1).
         */

        IncidentIterator& operator++()
        {
            incidentiterator++;
            /** This function checks whether two iterators are the same.
             * @param[in] new_EdgeIterator This is another incident iterator that will be
             * compared to the one defined in the private attributes of the class EdgeIterator.
             * @return true if the iterators are the same and false if they aren't.
             *
             * Complexity: O(1).
             */
            return *this;
        }

        /** This function checks whether two iterators are the same.
         * @param[in] new_IncidentIterator This is another incident iterator that will be
         * compared to the one defined in the private attributes of the class EdgeIterator.
         * @return true if the iterators are the same and false if they aren't.
         *
         * Complexity: O(1).
         */
        bool operator==(const IncidentIterator& new_IncidentIterator) const
        {
            bool same_graph,same_node;
            same_graph = (graph_==new_IncidentIterator.graph_) ;
            same_node = (nodeindex_ ==new_IncidentIterator.nodeindex_);
            if (same_graph && same_node)
            {
                return (incidentiterator == new_IncidentIterator.incidentiterator);
            }
            return false;
        }

    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE

        friend class Node;
        const Graph* graph_;
        size_type nodeindex_;
        typename std::vector<size_type>::const_iterator incidentiterator;
        // Constructor
        IncidentIterator(const Graph* graph, size_type nodeindex, typename std::vector<size_type>::const_iterator iter)
                : graph_(const_cast<Graph*>(graph)), nodeindex_(nodeindex), incidentiterator(iter) {
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
            edgeiterator = nullptr;
            graph_ = nullptr;
        }

        // HW1 #5: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
        // EdgeIterator& operator++()
        // bool operator==(const EdgeIterator&) const
        /** This function returns an Edge by getting the internal Edge
         * at the specific iD.
         * @return Edge for the specific edge iterator.
         *
         * Complexity: O(1).
         */
        Edge operator*() const
        {
            size_type edgeid = *edgeiterator;
            return graph_->edge(graph_->EdgesIndex[edgeid]);
        }

        /** This function returns a reference to the Edge iterator that was incremented by 1.
         * @return Edgeiterator incremented by 1.
         *
         * Complexity: O(1).
         */

        EdgeIterator& operator++()
        {
            edgeiterator++;
            return *this;
        }

        /** This function checks whether two iterators are the same.
         * @param[in] new_EdgeIterator This is another incident iterator that will be
         * compared to the one defined in the private attributes of the class EdgeIterator.
         * @return true if the iterators are the same and false if they aren't.
         *
         * Complexity: O(1).
         */

        bool operator==(const EdgeIterator& new_EdgeIterator) const
        {
            return (edgeiterator == new_EdgeIterator.edgeiterator);
        }

    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
        const Graph* graph_;
        typename std::vector<size_type>::const_iterator edgeiterator;
        // Constructor
        EdgeIterator(const Graph* graph, typename std::vector<size_type>::const_iterator iter)
                : graph_(const_cast<Graph*>(graph)), edgeiterator(iter) {
        }
    };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // edge_iterator edge_begin() const
    // edge_iterator edge_end() const

    edge_iterator edge_begin() const
    {
        return(EdgeIterator(this,this->EdgesID.begin()));
    }
    edge_iterator edge_end() const
    {
        return(EdgeIterator(this,this->EdgesID.end()));
    }

private:


    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.
    class internalnode
    {
        // variables
        Point nodeposition_;
        node_value_type nodevalue;
        size_type nodeindex_;

        // constructors
        internalnode() {}
        internalnode(Point nodeposition, node_value_type value, size_type nodeindex)
                : nodeposition_(nodeposition), nodevalue(value), nodeindex_(nodeindex) {}

        friend class Graph;
    };
    class internaledge
    {
        // variables
        size_type firstnodeid_;
        size_type secondnodeid_;
        size_type edgeindex_;
        edge_value_type edgevalue_;

        // constructors
        internaledge() {}
        internaledge(size_type firstnodeid, size_type secondnodeid, size_type edgeindex)
                : firstnodeid_(firstnodeid), secondnodeid_(secondnodeid), edgeindex_(edgeindex){}
        internaledge(size_type firstnodeid, size_type secondnodeid, size_type edgeindex, edge_value_type edgevalue)
                : firstnodeid_(firstnodeid), secondnodeid_(secondnodeid), edgeindex_(edgeindex), edgevalue_(edgevalue){}


        friend class Graph;
    };
    // To solve the problem of re-index, we give those nodes two kind of
    // indentification, nodeid_(the position in the internalnode) and index.
    // Each time we add or del a point, we only need to add or del the point in
    // the internalnode and, then, adjust map bewteen index and nodeid.

    std::vector<internalnode> InternalNodes;
    // NodesIndex[i] is the index of the node in the i-th location in the InternalNodes, i.e. id i to inex
    std::vector<size_type> NodesIndex;
    // NodesID[i] is the  location in the InteralNodes of the nodex indexed by i, i.e. index i to id
    std::vector<size_type> NodesID;

    std::vector<internaledge> InternalEdges;
    std::vector<size_type> EdgesIndex;
    std::vector<size_type> EdgesID;

    std::vector<std::vector<size_type>> IncidentEdgestoNodes;

};


#endif // CME212_GRAPH_HPPff00f51dff
