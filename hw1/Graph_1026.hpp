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
	struct PrvNode; // Forward declaration of structure
public:

	//
	// PUBLIC TYPE DEFINITIONS
	//

	/** Template type for node. */
	using node_value_type = V;

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

		/** Return this node's position. */
		const Point& position() const {
			return((graph_->nodes_)[index_].position_);
		}

		/** Return this node's index, a number in the range [0, graph_size). */
		size_type index() const {
			return(index_);
		}

		/** Returns this node's value, of type node_value_type. */
		node_value_type& value() {
			return((graph_->nodes_)[index_].val_);
		}

		/** Returns this node's value, of type const node_value_type. */
		const node_value_type& value() const {
			return((graph_->nodes_)[index_].val_);
		}

		/** Returns the number of adjacent nodes of this node. */
		size_type degree() const {
			return((graph_->nodeSuccessors_)[index_].size());
		}

		/** Returns the first incident_iterator of this node. */
		incident_iterator edge_begin() const {
			return(IncidentIterator(graph_, index_, 0));
		}

		/** Returns the last incident_iterator of this node. */
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
			return(index_ < n.index() && graph_ == n.graph_);
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
//--functionality_1
//--no add_node function that lets you add node with given value
//--START
	Node add_node(const Point& position) {
		PrvNode privateNode = PrvNode(size_, position);
		nodes_.push_back(privateNode);
		nodeSuccessors_.push_back(std::vector<size_type>());
		Node node = Node(this, size_);
		this->size_ += 1;
		return node;
	}
//--END

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

		/** Test whether this edge and @a e are equal.
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

		/** Test whether this edge is less than @a e in a global order.
		 *
		 * This ordering function is useful for STL containers such as
		 * std::map<>. It need not have any interpretive meaning.
		 */
		bool operator<(const Edge& e) const {
			//We choose the lexicographic order on (node1, node2) as global order
			if (graph_ != e.graph_)
				return false;
			if (index1_ == e.index1_)
				return(index2_ < e.index2_);
			else
				return(index1_ < e.index1_);
		}

	private:
		// Allow Graph to access Edge's private member data and functions.
		friend class Graph;

		Graph* graph_;
		size_type index1_;
		size_type index2_;

		/** Return a valid edge based on its nodes. */
		Edge(const Graph* graph, size_type index1, size_type index2) {
			this->graph_ = const_cast<Graph*>(graph);
			this->index1_ = index1;
			this->index2_ = index2;
		}

		/** Return a valid edge base on its edge index. */
		Edge(const Graph* graph, size_type edgeIdx) {
			this->graph_ = const_cast<Graph*>(graph);
			std::pair<size_type, size_type> edge = graph_->edges_[edgeIdx];
			this->index1_ = edge.first;
			this->index2_ = edge.second;
		}

	};


	/** Return the total number of edges in the graph.
	 *
	 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	 */
	size_type num_edges() const {
		return numEdges_;
	}

	/** Return the edge with index @a i.
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

	/** Test whether two nodes are connected by an edge.
	 * @pre @a a and @a b are valid nodes of this graph
	 * @return True if for some @a i, edge(@a i) connects @a a and @a b.
	 *
	 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	 */
//--design_0
//--you should speed this up now that you have nodeSuccessors_
//--START
	bool has_edge(const Node& a, const Node& b) const {
		assert(has_node(a));
		assert(has_node(b));
		size_type mn = a.index_ <= b.index_ ? a.index_ : b.index_; //smaller index
		size_type mx = a.index_ < b.index_ ? b.index_ : a.index_; //bigger index
		for (size_type k = 0; k < this->numEdges_; k++) {
			if (edges_[k].first == mn && this->edges_[k].second == mx)
				return true;
		}
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
	Edge add_edge(const Node& a, const Node& b) {
		Edge edge = Edge(this, a.index_, b.index_);
		//Adding Edge to graph
		if (!this->has_edge(a, b)) {
			size_type mn = a.index_ <= b.index_ ? a.index_ : b.index_; //smaller index
			size_type mx = a.index_ < b.index_ ? b.index_ : a.index_; //bigger index
			std::pair<size_type, size_type> indexPair = std::make_pair(mn, mx);
			edges_.push_back(indexPair);
			nodeSuccessors_[mn].push_back(mx);
			nodeSuccessors_[mx].push_back(mn);
			numEdges_ += 1;
		}
		return edge;
	}
	/** Remove all nodes and edges from this graph.
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

		/** Construct an invalid NodeIterator. */
		NodeIterator() {
			//No init needed
		}

//--documentation_1
//--docs should be in doxygen style
//--START
		/** Return the underlying node of the iterator. Dereferrencing. */
		Node operator*() const {
			return(Node(graph_, index_));
		}

		/** Return the iterator incremented. */
		NodeIterator& operator++() {
			index_++;
			return(*this);
		}

		/** Test whether this NodeIterator and @NodeIterator node_iter are equal.
		 *
		 * Equal NodeIterators represent the same underlying node in the same graph.
		 */
		bool operator==(const NodeIterator& node_iter) const {
			return(index_ == node_iter.index_ && graph_ == node_iter.graph_);
		}
//--END

	private:
		friend class Graph;
		Graph* graph_;
		size_type index_;

		/** Return a valide NodeIterator. */
		NodeIterator(const Graph* graph, size_type index) {
			this->graph_ = const_cast<Graph*>(graph);
			this->index_ = index;
		}
	};

	/** Return the first NodeIterator of the graph. */
	node_iterator node_begin() const {
		return(NodeIterator(this, 0));
	}
//--functionality_0
//--incorrect end iterator for node and edge; should be one past the last element
//--START
	/** Return the last NodeIterator of the graph. */
	node_iterator node_end() const {
		return(NodeIterator(this, size_ - 1));
	}
//--END


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

		/** Return the underlying Edge of the iterator. Dereferrencing. */
		Edge operator*() const {
			size_type nodeSuccessorIdx = (graph_->nodeSuccessors_)[rootIdx_][edgeIdx_];
			return(Edge(graph_, rootIdx_, nodeSuccessorIdx));
		}

		/** Return the iterator incremented. */
		IncidentIterator& operator++() {
			edgeIdx_++;
			return(*this);
		}

		/** Test whether this IncidentIterator and @IncidentIterator incident_iter are equal.
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
		// HW1 #3: YOUR CODE HERE
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

		/** Construct an invalid EdgeIterator. */
		EdgeIterator() {
			//No init needed.		
		}

		/** Return the underlying Edge of the iterator. Dereferrencing. */
		Edge operator*() const {
			return(Edge(graph_, index_));
		}

		/** Return the iterator incremented. */
		EdgeIterator& operator++() {
			index_++;
			return(*this);
		}

		/** Test whether this EdgeIterator and @EdgeIterator edge_iter are equal.
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

//--functionality_0
//--also needs to be one past the last
//--START
	/** Return the last edge_iterator in the graph. */
	edge_iterator edge_end() const {
		return(EdgeIterator(this, numEdges_ - 1));
	}
//--END

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
	size_type numEdges_; //Number of edges in graph.
};

#endif // CME212_GRAPH_HPP
