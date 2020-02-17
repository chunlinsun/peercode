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
template <typename V, typename E>
class Graph {
private:
	struct PrvNode; // Forward declaration of structure
public:

	//
	// PUBLIC TYPE DEFINITIONS
	//

	/** Template type for node. */
	using node_value_type = V;

	/** Template type for edge. */
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
		this->nodes_ = std::vector<PrvNode>();
		this->size_ = 0;
		this->nodeSuccessors_ = std::vector<std::vector<size_type>>();
		this->edges_ = std::vector<std::pair<size_type, size_type>>();
		this->numEdges_ = 0;
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
			//No init needed
		}

		/**
		*@brief Return this node's position. 
		*/
		const Point& position() const {
			return((graph_->nodes_)[index_].position_);
		}

		/** 
		* @brief Return a referernce of this node's position. 
		*/
		Point& position() {
			return((graph_->nodes_)[index_].position_);
		}

		/** 
		 * @brief Return this node's index, a number in the range [0, graph_size). 
		*/
		size_type index() const {
			return(index_);
		}

		/** 
		 * @brief Returns this node's value, of type node_value_type. 
		*/
		node_value_type& value() {
			return((graph_->nodes_)[index_].val_);
		}

		/** 
		* @brief Returns this node's value, of type const node_value_type
		*/
		const node_value_type& value() const {
			return((graph_->nodes_)[index_].val_);
		}

		/** 
		 * @brief Returns the number of adjacent nodes of this node. 
		*/
		size_type degree() const {
			return((graph_->nodeSuccessors_)[index_].size());
		}

		/** 
		 * @brief Returns the first incident_iterator of this node. 
		*/
		incident_iterator edge_begin() const {
			return(IncidentIterator(graph_, index_, 0));
		}

		/** 
		 * @brief Returns the last incident_iterator of this node.
		*/
		incident_iterator edge_end() const {
			return(IncidentIterator(graph_, index_, degree()));
		}
		

		/** Test whether this node and @a n are equal.
		 *
		 * Equal nodes have the same graph and the same index.
		 */
		bool operator==(const Node& n) const {
			return (index_ == n.index() && graph_ == n.graph_);
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
			return(index_ < n.index());// && graph_ == n.graph_);
		}

	private:
		// Allow Graph to access Node's private member data and functions.
		friend class Graph;

		Graph* graph_; //Pointer to the node's graph
		size_type index_; //Index of the node in graph nodes

		/** Construct a valid node. */
		Node(const Graph* graph, size_type index) {
			this->graph_ = const_cast<Graph*>(graph);
			this->index_ = index;
		}
	};

	/** Return the number of nodes in the graph.
	 *
	 * Complexity: O(1).
	 */
	size_type size() const {
		return size_;
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
		PrvNode privateNode = PrvNode(size_, position);
		nodes_.push_back(privateNode);
		nodeSuccessors_.push_back(std::vector<size_type>());
		Node node = Node(this, size_);
		this->size_ += 1;
		return node;
	}

	/** Add a node to the graph, returning the added node.
	 * @param[in] position The new node's position
	 * @param[in] val The new node's value
	 * @post new num_nodes() == old num_nodes() + 1
	 * @post result_node.index() == old num_nodes()
	 *
	 * Complexity: O(1) amortized operations.
	 */
	Node add_node(const Point& position, node_value_type val) {
		Node node = add_node(position);
		node.value() = val;
		return node;
	}

	/** Determine if a Node belongs to this Graph
	 * @return True if @a n is currently a Node of this Graph
	 *
	 * Complexity: O(1).
	 */
	bool has_node(const Node& n) const {
		return(n.index() < size_ && n.graph_ == this);
	}

	
	/** Return the node with index @a i.
	 * @pre 0 <= @a i < num_nodes()
	 * @post result_node.index() == i
	 *
	 * Complexity: O(1).
	 */
	Node node(size_type i) const {
		Node node = Node();
		//If index not in graph we return invalid node
		if (i < this->size_) {
			node = Node(this, i);
		}
		return node;
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
			//No init needed
		}
	
		/** Return a node of this Edge */
		Node node1() const {
			return Node(graph_, index1_);
		}

		/** Return the other node of this Edge */
		Node node2() const {
			return Node(graph_, index2_);
		}

		/** 
		* @brief Return a reference of this edge's value .
		*/
		edge_value_type& value() {
			size_type edgeIdx = find_index();
			return(graph_->edgeValues_[edgeIdx]);
		}

		/**
		* @brief Return a this edge's value .
		*/
		const edge_value_type& value() const{
			size_type edgeIdx = find_index();
			return(graph_->edgeValues_[edgeIdx]);
		}

		/**
		* @brief Finds global index of this edge.
		*/
		size_type find_index() {
			return(graph_->find_index(*this));
		}

		/**
		* @brief Returns the distance between the two nodes of this edge.
		*/
		double length() const{
			return norm(node2().position() - node1().position());
		}

		/** @brief Test whether this edge and @a e are equal.
		 *
		 * Equal edges represent the same undirected edge between two nodes.
		 */
		bool operator==(const Edge& e) const {
			//As we force the indexe of node1 to be inferior or equal to the one of node2,
			//We can restrict ourselves to this simple check
			bool sameGraph = (graph_ == e.graph_);
			bool sameNodes = (index1_ == e.index1_ && index2_ == e.index2_) ||
				(index1_ == e.index2_ && index2_ == e.index1_);
			return(sameGraph && sameNodes);
		}

		/** @brief Test whether this edge is less than @a e in a global order.
		 *
		 * This ordering function is useful for STL containers such as
		 * std::map<>. It need not have any interpretive meaning.
		 */
		bool operator<(const Edge& e) const {
			//We choose the lexicographic order on (node1, node2) as global order
			if (graph_ == e.graph_){
				if (index1_ == e.index1_)
					return(index2_ < e.index2_);
				else
					return(index1_ < e.index1_);
			}
			else{
				return (graph_ < e.graph_);
			}
		}

	private:
		// Allow Graph to access Edge's private member data and functions.
		friend class Graph;

		Graph* graph_;
		size_type index1_;
		size_type index2_;

		/** @brief Return a valid edge based on its nodes. */
		Edge(const Graph* graph, size_type index1, size_type index2) {
			this->graph_ = const_cast<Graph*>(graph);
			this->index1_ = index1;
			this->index2_ = index2;
		}

		/** @brief Return a valid edge base on its edge index. */
		Edge(const Graph* graph, size_type edgeIdx) {
			this->graph_ = const_cast<Graph*>(graph);
			std::pair<size_type, size_type> edge = graph_->edges_[edgeIdx];
			this->index1_ = edge.first;
			this->index2_ = edge.second;
		}

	};


	/** @brief Return the total number of edges in the graph.
	 *
	 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	 */
	size_type num_edges() const {
		return numEdges_;
	}

	/** @brief Return the edge with index @a i.
	 * @pre 0 <= @a i < num_edges()
	 *
	 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	 */
	Edge edge(size_type i) const {
		Edge edge = Edge();
		//If edge not in graph, we return invalid edge
		if (i < numEdges_) {
			edge = Edge(this, i);
		}
		return edge;
	}

	/**
	* @brief Finds global index of @a edge in edges_.
	*/
	size_type find_index(Edge edge) {
		return(find_index(edge.node1(), edge.node2()));
	}

	/**
	* @brief Finds global index of @a edge in deges_.
	*/
	size_type find_index(const Node& a, const Node& b) {
		size_type mn = a.index() <= b.index() ? a.index() : b.index(); //smaller index
		size_type mx = a.index() < b.index() ? b.index() : a.index(); //bigger index
		for (size_type k = 0; k < numEdges_; k++) {
			if (edges_[k].first == mn && edges_[k].second == mx)
				return k;
		}
		return(numEdges_);
	}




	/** @brief Test whether two nodes are connected by an edge.
	 * @pre @a a and @a b are valid nodes of this graph
	 * @return True if for some @a i, edge(@a i) connects @a a and @a b.
	 *
	 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	 */
	bool has_edge(const Node& a, const Node& b) const {
		if (!has_node(a) || !has_node(b)){
			return false;
		}
		for(auto iter = a.edge_begin(); iter != a.edge_end(); ++iter){
			if((*iter).node2() == b)
				return(true);
		}
		return false;
	}

	/** @brief Add an edge to the graph, or return the current edge if it already exists.
	 * @pre @a a and @a b are distinct valid nodes of this graph
	 * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
	 * @post has_edge( @a a, @a b) == true
	 * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
	 *       Else,                        new num_edges() == old num_edges() + 1.
	 *
	 * Can invalidate edge indexes -- in other words, old edge(@a i) might not
	 * equal new edge(@a i). Must not invalidate outstanding Edge objects.
	 *
	 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	 */
	Edge add_edge(const Node& a, const Node& b) {
		Edge edge = Edge(this, a.index_, b.index_);
		//Adding Edge to graph
		if (!this->has_edge(a, b)) {
			size_type mn = a.index_ <= b.index_ ? a.index_ : b.index_; //smaller index
			size_type mx = a.index_ < b.index_ ? b.index_ : a.index_; //bigger index
			std::pair<size_type, size_type> indexPair = std::make_pair(mn, mx);
			edges_.push_back(indexPair);
			edgeValues_.push_back(edge_value_type());
			nodeSuccessors_[mn].push_back(mx);
			nodeSuccessors_[mx].push_back(mn);
			numEdges_ += 1;
		}
		return edge;
	}
	/** @brief Remove all nodes and edges from this graph.
	 * @post num_nodes() == 0 && num_edges() == 0
	 *
	 * Invalidates all outstanding Node and Edge objects.
	 */
	void clear() {
		nodes_.clear();
		size_ = 0;
		nodeSuccessors_.clear();
		edges_.clear();
		numEdges_ = 0;
	}

	//
	// Node Iterator
	//

	/** @class Graph::NodeIterator
	 * @brief Iterator class for nodes. A forward iterator. */
	class NodeIterator : private totally_ordered<NodeIterator> {
	public:
		// These type definitions let us use STL's iterator_traits.
		using value_type = Node;                     // Element type
		using pointer = Node*;                    // Pointers to elements
		using reference = Node&;                    // Reference to elements
		using difference_type = std::ptrdiff_t;           // Signed difference
		using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

		/** @brief Construct an invalid NodeIterator. */
		NodeIterator() {
			//No init needed
		}

		/** @brief Return the underlying node of the iterator. Dereferrencing. */
		Node operator*() const {
			return(Node(graph_, index_));
		}

		/** @brief Return the iterator incremented. */
		NodeIterator& operator++() {
			index_++;
			return(*this);
		}

		/** @brief Return the iterator decremented. */
		NodeIterator& operator--() {
			index_--;
			return(*this);
		}

		/** @brief Test whether this NodeIterator and @NodeIterator node_iter are equal.
		 *
		 * Equal NodeIterators represent the same underlying node in the same graph.
		 */
		bool operator==(const NodeIterator& node_iter) const {
			return(index_ == node_iter.index_ && graph_ == node_iter.graph_);
		}

	private:
		friend class Graph;
		Graph* graph_;
		size_type index_;

		/** @brief Return a valide NodeIterator. */
		NodeIterator(const Graph* graph, size_type index) {
			this->graph_ = const_cast<Graph*>(graph);
			this->index_ = index;
		}
	};

	/** @brief Return the first NodeIterator of the graph. */
	node_iterator node_begin() const {
		return(NodeIterator(this, 0));
	}

	/** @brief Return the last NodeIterator of the graph. */
	node_iterator node_end() const {
		return(NodeIterator(this, size_));
	}


	//
	// Incident Iterator
	//

	/** @class Graph::IncidentIterator
	 * @brief Iterator class for edges incident to a node. A forward iterator. */
	class IncidentIterator : private totally_ordered<IncidentIterator> {
	public:
		// These type definitions let us use STL's iterator_traits.
		using value_type = Edge;                     // Element type
		using pointer = Edge*;                    // Pointers to elements
		using reference = Edge&;                    // Reference to elements
		using difference_type = std::ptrdiff_t;           // Signed difference
		using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

		/** Construct an invalid IncidentIterator. */
		IncidentIterator() {
			// No init needed
		}

		/** @brief Return the underlying Edge of the iterator. Dereferrencing. */
		Edge operator*() const {
			size_type nodeSuccessorIdx = (graph_->nodeSuccessors_)[rootIdx_][edgeIdx_];
			return(Edge(graph_, rootIdx_, nodeSuccessorIdx));
		}

		/** @brief Return the iterator incremented. */
		IncidentIterator& operator++() {
			edgeIdx_++;
			return(*this);
		}

		/** @brief Return the iterator decremented. */
		IncidentIterator& operator--() {
			edgeIdx_--;
			return(*this);
		}

		/** @brief Test whether this IncidentIterator and @IncidentIterator incident_iter are equal.
		 *
		 * Equal IncidentIterators represent the same underlying edge 
		 * starting from the same root node in the same graph.
		 */
		bool operator==(const IncidentIterator& incident_iter) const {
			return(graph_ == incident_iter.graph_ &&
				rootIdx_ == incident_iter.rootIdx_ &&
				edgeIdx_ == incident_iter.edgeIdx_);
		}


	private:
		friend class Graph;
		Graph* graph_;
		size_type rootIdx_;
		size_type edgeIdx_; //Index of the edge within incident edges

		/** Return a valid IncidentIterator. */
		IncidentIterator(const Graph* graph, size_type rootIdx, size_type edgeIdx) {
			this->graph_ = const_cast<Graph*>(graph);
			this->rootIdx_ = rootIdx;
			this->edgeIdx_ = edgeIdx;
		}

	};

	//
	// Edge Iterator
	//

	/** @class Graph::EdgeIterator
	 * @brief Iterator class for edges. A forward iterator. */
	class EdgeIterator : private totally_ordered<EdgeIterator> {
	public:
		// These type definitions let us use STL's iterator_traits.
		using value_type = Edge;                     // Element type
		using pointer = Edge*;                    // Pointers to elements
		using reference = Edge&;                    // Reference to elements
		using difference_type = std::ptrdiff_t;           // Signed difference
		using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

		/** @brief Construct an invalid EdgeIterator. */
		EdgeIterator() {
			//No init needed.		
		}

		/** @brief Return the underlying Edge of the iterator. Dereferrencing. */
		Edge operator*() const {
			return(Edge(graph_, index_));
		}

		/** @brief Return the iterator incremented. */
		EdgeIterator& operator++() {
			index_++;
			return(*this);
		}

		/** @brief Test whether this EdgeIterator and @EdgeIterator edge_iter are equal.
		 *
		 * Equal EdgeIterators represent the same underlying edge in the same graph.
		 */
		bool operator==(const EdgeIterator& edge_iter) const {
			return(index_ == edge_iter.index_ && graph_ == edge_iter.graph_);
		}


	private:
		friend class Graph;
		Graph* graph_;
		size_type index_;

		/** Return a valid iterator. */
		EdgeIterator(const Graph* graph, size_type index) {
			this->graph_ = const_cast<Graph*>(graph);
			this->index_ = index;
		}
	};
	 /** Return the first edge_iterator in the graph. */
	edge_iterator edge_begin() const {
		return(EdgeIterator(this, 0));
	}

	/** Return the last edge_iterator in the graph. */
	edge_iterator edge_end() const {
		return(EdgeIterator(this, numEdges_));
	}

	//
	// NODE & EDGE REMOVER
	//
	

	 /** @brief Remove @a node from the graph (if it exists in the
     * graph), and return 1 if node deleted OR 0 otherwise.
     *
     * @param[in]  a    Node to delete
     *
	 * @return 			number of nodes removed
	 * 
	 * @post            IF has_node( @a node) == true,
     *                         size_ <-   size_ - 1
     *                  ELSE
     *                         size_ <-   size_
	 * 
     * Complexity:  - has_node is O(1).
	 * 				- removing all edges connected to node is O(numEdges_).
	 * 				- changing the swaping node's index in all structures is O(numEdges_).
	 * 				- swaping node is O(1).
     */
	size_type remove_node (const Node &node){
				
		if(!has_node(node))
			return false;

		//Remove all edges connected to node
		std::vector<size_type> successors = nodeSuccessors_[node.index()];
		for (auto iter = successors.begin(); iter != successors.end(); ++iter){
			remove_edge(node, Node(this, *iter));	
		}
		nodeSuccessors_[node.index()].clear();

		//Swap node
		if(node.index() != size_-1){
			
			Node swapNode = Node(this,size_-1);

			//Look for all the edges connected to swap node and change its index
			std::vector<size_type>& swapSuccessors = nodeSuccessors_[size_-1];
			for( auto iter = swapSuccessors.begin(); iter != swapSuccessors.end(); ++iter){
				size_type edgeIdx = find_index(swapNode,Node(this,*iter));
				if (edges_[edgeIdx].first == swapNode.index())
					 edges_[edgeIdx].first = node.index();
				else
					edges_[edgeIdx].second = node.index();

				if (edges_[edgeIdx].first > edges_[edgeIdx].second){
					size_type temp = edges_[edgeIdx].second;
					edges_[edgeIdx].second = edges_[edgeIdx].first;
					edges_[edgeIdx].first = temp;
				}
								
				std::vector<size_type>& newSuccessors = nodeSuccessors_[*iter];
				for (auto secondIter = newSuccessors.begin(); secondIter != newSuccessors.end(); ++secondIter){
					if(*secondIter == size_-1){
						*secondIter = node.index();
						break;
					}
				}
			}
			nodes_[size_-1].index_ = node.index();
			nodes_[node.index()] = nodes_[size_-1];
			nodeSuccessors_[node.index()] = nodeSuccessors_[size_-1];
		}

		nodeSuccessors_.erase(nodeSuccessors_.end()-1);
		nodes_.erase(nodes_.end()-1);
		size_--;
		return(1);
	}


	/** @brief Remove node referenced by iterator @a n_it from the graph (if it exists in the
     * graph), and return 1 if node deleted OR 0 otherwise.
     *
     * @param[in]  a    iterator referencing the node to delete
     *
	 * @return 			iterator with a new reference.
	 * 
	 * @post            IF has_node( @a node) == true,
     *                         size_ <-   size_ - 1
     *                  ELSE
     *                         size_ <-   size_
	 * 
     * Complexity:  - has_node is O(1).
	 * 				- removing all edges connected to node is O(numEdges_).
	 * 				- changing the swaping node's index in all structures is O(numEdges_).
	 * 				- swaping node is O(1).
     */
	node_iterator remove_node ( node_iterator n_it ){
		remove_node(*n_it);
		return(n_it);
	}

	 /** @brief Remove the edge between two nodes @a a and @a b from the graph (if the
     * edge exists), returning the number of edges successfully deleted (0 or 1).
     *
     * @param[in,out] a   Node of the edge
     * @param[in,out] b   Node of the edge
	 * 
	 * @returns         Number of edges deleted (0 or 1)
     *
     * @post            IF has_edge( @a a, @a b) == true,
     *                         num_edges_ <-   num_edges_ - 1
     *                  ELSE
     *                         num_edges_ <-   num_edges_
	 * 
	 * Runtime: - Finding edge index takes O(numEdges_).
	 * 			- Removing edge from edges_ and edgeValues_ takes O(1) on average.
	 * 			- Removing edge from nodeSuccessors takes O(1) on average.
	 */	
	size_type remove_edge (const Node &a, const Node &b) {
		
		if(!has_edge(a,b))
			return(0);

		size_type edgeIdx = find_index(a,b);
		
		// Remove from edges_ and edgeValues_
		if(edgeIdx != numEdges_){
			edges_[edgeIdx] = edges_[numEdges_ - 1];
			edgeValues_[edgeIdx] = edgeValues_[numEdges_ -1];
		}
		edges_.erase(edges_.end()-1);
		edgeValues_.erase(edgeValues_.end()-1);

		// Remove from nodeSuccessors
		std::vector<size_type>& aSuccessors = nodeSuccessors_[a.index()];
		for(auto iter = aSuccessors.begin(); iter != aSuccessors.end(); ++iter){
			if (*iter==b.index()){
				if(iter != aSuccessors.end()-1){
					*iter = *(aSuccessors.end()-1);
				}
				aSuccessors.erase(aSuccessors.end()-1);
				break;
			}
		}
		std::vector<size_type>& bSuccessors = nodeSuccessors_[b.index()];
		for(auto iter = bSuccessors.begin(); iter != bSuccessors.end(); ++iter){
			if (*iter==a.index()){
				if(iter != bSuccessors.end()-1){
					*iter = *(bSuccessors.end()-1);
				}
				bSuccessors.erase(bSuccessors.end()-1);
				break;
			}
		}
		numEdges_--;
		return(1);
	}

	/** @brief Remove the edge @a edge from the graph (if the
     * edge exists), returning the number of edges successfully deleted (0 or 1).
     *
     * @param[in,out] edge	Edge to delete
	 * 
	 * @returns         Number of edges deleted (0 or 1)
     *
     * @post            IF has_edge( @a a, @a b) == true,
     *                         num_edges_ <-   num_edges_ - 1
     *                  ELSE
     *                         num_edges_ <-   num_edges_
	 * 
	 * Runtime: - Finding edge index takes O(numEdges_).
	 * 			- Removing edge from edges_ and edgeValues_ takes O(1) on average.
	 * 			- Removing edge from nodeSuccessors takes O(1) on average.
	 */	
	size_type remove_edge (const Edge &edge){
		return(remove_edge(edge.node1(), edge.node2()));
	}

	/** @brief Remove the edge  referenced by the edge iterator @a e_it from the graph (if the
     * edge exists), returning the number of edges successfully deleted (0 or 1).
     *
     * @param[in,out] e_it	edge iterator referencing the edge to delete
	 * 
	 * @returns         e_it with new reference.
     *
     * @post            IF has_edge( @a a, @a b) == true,
     *                         num_edges_ <-   num_edges_ - 1
     *                  ELSE
     *                         num_edges_ <-   num_edges_
	 * 
	 * Runtime: - Finding edge index takes O(numEdges_).
	 * 			- Removing edge from edges_ and edgeValues_ takes O(1) on average.
	 * 			- Removing edge from nodeSuccessors takes O(1) on average.
	 */	
	edge_iterator remove_edge (edge_iterator e_it){
		remove_edge(*e_it);
		return(e_it);
	}




private:

	/** Structure that links all the information pertainning to a specific node. */
	struct PrvNode {

		size_type index_;
		Point position_;
		node_value_type val_;

		/** Construct a PrvNode based on index and position. */
		PrvNode(const size_type index, const Point& position) {
			this->index_ = index;
			this->position_ = position;
			this->val_ = node_value_type();
		}

		/** Construct a PrvNode based on index, position and value. */
		PrvNode(const size_type index, const Point& position, node_value_type val) {
			this->index_ = index;
			this->position_ = position;
			this->val_ = val;
		}
	};
	std::vector<PrvNode> nodes_; //Container of our datastructure.
	size_type size_; //Number of nodes in graph.
	std::vector<std::vector<size_type>> nodeSuccessors_; //Container of the all the nodes' successors.
	std::vector<std::pair<size_type, size_type>> edges_; //Container of the edges.
	std::vector<edge_value_type> edgeValues_; //Container of the edges values
	size_type numEdges_; //Number of edges in graph.
};

#endif // CME212_GRAPH_HPP
