#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>
#include <iostream>
#include <typeinfo>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

template <typename V, typename E>
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and s. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph 
{
	public:

		//
		// PUBLIC TYPE DEFINITIONS
		//
		using size_type = unsigned;
		using node_value_type = V;
		using edge_value_type = E;
		/** Type of this graph. */
		using graph_type = Graph;

		/** Predeclaration of Node type. */
		class Node;
		/** Synonym for Node (following STL conventions). */
		using node_type = Node;
		/** making a tuple for endpoints of an edge. */
		typedef std::tuple<size_type, size_type> node_tuple;

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

		/** Construct an empty graph. */
		Graph() 
		{}
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
		class Node : private totally_ordered<Node>
		{
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

				Point& position() 
				{
					return graph_ptr->node_positions.at(node_index);
				}

				const Point& position()  const
				{
					return const_cast<Point&>(graph_ptr->node_positions.at(node_index));
				}

				node_value_type& value()
				{
					return graph_ptr->node_values.at(node_index);
				}

				// return value corresponding to the node
				const node_value_type& value() const
				{ 
					return const_cast<node_value_type&> (graph_ptr->node_values.at(node_index));
				}

				// TEST
				size_type degree() const
				{
					if(graph_ptr->adj_list.find(node_index) == graph_ptr->adj_list.end())
						return 0;
					return graph_ptr->adj_list.at(node_index).size();
				}

				// return uid of the node
				size_type index()
				{
					return graph_ptr->node_uid_map.at(node_index);
				}

				/** Test whether this node and @a n are equal.
				 *
				 * Equal nodes have the same graph and the same index.
				 */
				bool operator==(const Node& n) const 
				{
					if(graph_ptr != n.graph_ptr)
						return false;
					return n.node_index == node_index;
				}

				/** Test whether this node is less than @a n in a global order.
				 *
				 * This ordering function is useful for STL containers such as
				 * std::map<>. It need not have any geometric meaning.
				 *
				 * The node ordering relation must obey trichotomy: For any two nodes x
				 * and y, exactly one of x == y, x < y, and y < x is true.
				 */
				bool operator<(const Node& n) const 
				{
					// graph_1 smaller than graph_2
					if(std::less<const Graph*>{} (graph_ptr, n.graph_ptr))  
						return true;
					// graph_1 is larger than graph_2
					if(std::greater<const Graph*>{} (graph_ptr, n.graph_ptr))  
						return false;
					// now we are comparing nodes of same graph
					return node_index < n.node_index;
				}

				// create incident iterator pointing at first edge incident on node
				incident_iterator edge_begin() const
				{
					return IncidentIterator(node_index, 0, graph_ptr);
				}

				// create incident iterator pointing at one past the last edge incident on node
				incident_iterator edge_end() const
				{
					return IncidentIterator(node_index, degree(), graph_ptr);
				}

			 private:
				size_type node_index;
				graph_type* graph_ptr;
				// private constructor
				Node(size_type ind, const graph_type* gp)
				{
					node_index = ind;
					graph_ptr = const_cast<graph_type*> (gp);
				}
				
				friend class Graph;
		};

		class Edge : private totally_ordered<Edge>
		{
			public:
				/** Construct an invalid Edge. */
				Edge() 
				{}

				/** Return a node of this Edge */
				Node node1() const 
				{
					size_type uid = graph_ptr->node_uid_map.at(node_1_index);
					return graph_ptr->nodes.at(uid);      
				}

				/** Return the other node of this Edge */
				Node node2() const 
				{
					size_type uid = graph_ptr->node_uid_map.at(node_2_index);
					return graph_ptr->nodes.at(uid);      
				}

				/** Test whether this edge and @a e are equal.
				 *
				 * Equal edges represent the same undirected edge between two nodes.
				 */
				bool operator==(const Edge& e) const 
				{
					if(graph_ptr != e.graph_ptr)
						return false;
					// if(node_1_index == e.node_1_index and node_2_index == e.node_2_index)
					// 	return true;
					// if(node_1_index == e.node_2_index and node_2_index == e.node_1_index)
					// 	return true;
					return edge_index == e.edge_index;
				}

				/** Test whether this edge is less than @a e in a global order.
				 *
				 * This ordering function is useful for STL containers such as
				 * std::map<>. It need not have any interpretive meaning.
				 */
				bool operator<(const Edge& e) const 
				{
					// graph_1 smaller than graph_2
					if(std::less<const Graph*>{} (graph_ptr, e.graph_ptr))  
						return true;
					// graph_1 is larger than graph_2
					if(std::greater<const Graph*>{} (graph_ptr, e.graph_ptr))  
						return false;
					return edge_index < e.edge_index;
				}

				// return value associated with the edge
				const edge_value_type& value() const
				{
					return const_cast<edge_value_type&>(graph_ptr->edge_values.at(edge_index));
				}

				edge_value_type& value() 
				{
					return graph_ptr->edge_values.at(edge_index);
				}

			private:
				// node_1_index stores the smaller index
				size_type node_1_index;
				// node_2_index stores the larger index
				size_type node_2_index;
				// edge_index stores unique id for the edge
				size_type edge_index;
				// store the graph pointer
				Graph* graph_ptr;

				// private constructor
				Edge(size_type a, size_type b, size_type ind, const Graph* gp)
				{
					// TODO : Check for self loops
					// the smaller node goes first
					node_1_index = std::min(a, b);
					node_2_index = std::max(a, b);
					edge_index = ind;
					graph_ptr = const_cast<Graph*> (gp);
				}
				friend class Graph;	
		};

		//
		// Node Iterator
		//
		/** @class Graph::NodeIterator
		 * @brief Iterator class for nodes. A forward iterator. */
		class NodeIterator : private totally_ordered<NodeIterator> 
		{
			public:
				// These type definitions let us use STL's iterator_traits.
				using value_type        = Node;                     // Element type
				using pointer           = Node*;                    // Pointers to elements
				using reference         = Node&;                    // Reference to elements
				using difference_type   = std::ptrdiff_t;           // Signed difference
				using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy


				/** Construct an invalid NodeIterator. */
				NodeIterator() {}

				// dereference operator
				Node operator*() const
				{
					return graph_ptr->node(current_index);
				}

				NodeIterator& operator++() 
				{
					if(current_index < graph_ptr->nodes.size())
						current_index++;
					return *this;
				}

				bool operator==(const NodeIterator& it) const
				{
					return current_index == it.current_index;
				}

			private:
					size_type current_index;
					graph_type* graph_ptr;

				// private constructor
				NodeIterator(size_type i, const graph_type* g) 
				{
					current_index = i;
					graph_ptr = const_cast<graph_type*> (g);

				}
				friend class Graph;

		};

		//
		// Incident Iterator
		//
		/** @class Graph::IncidentIterator
		 * @brief Iterator class for edges incident to a node. A forward iterator. */
		class IncidentIterator : private totally_ordered<IncidentIterator>
		{
			public:
				// These type definitions let us use STL's iterator_traits.
				using value_type        = Edge;                     // Element type
				using pointer           = Edge*;                    // Pointers to elements
				using reference         = Edge&;                    // Reference to elements
				using difference_type   = std::ptrdiff_t;           // Signed difference
				using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

				/** Construct an invalid IncidentIterator. */
				IncidentIterator() 
				{}

				// dereference operator
				Edge operator*() const
				{
					return graph_ptr->adj_list.at(node_index).at(current_index);	
				}

				IncidentIterator& operator++()
				{
					if(current_index < graph_ptr->adj_list.at(node_index).size())
						current_index++;
					return *this;
				}

				// equality operator
				bool operator==(const IncidentIterator& it) const
				{
					return node_index == it.node_index and current_index == it.current_index;
				}

			private:
				size_type node_index;
				size_type current_index;
				graph_type* graph_ptr;

				IncidentIterator(size_type index, size_type i, const graph_type* g)
				{
					node_index = index;
					current_index = i;
					graph_ptr = const_cast<graph_type*> (g);
				}

				friend class Graph;

		};

//
	// Edge Iterator
	//
	/** @class Graph::EdgeIterator
	 * @brief Iterator class for edges. A forward iterator. */
	class EdgeIterator : private totally_ordered<EdgeIterator>
	{
		public:
			// These type definitions let us use STL's iterator_traits.
			using value_type        = Edge;                     // Element type
			using pointer           = Edge*;                    // Pointers to elements
			using reference         = Edge&;                    // Reference to elements
			using difference_type   = std::ptrdiff_t;           // Signed difference
			using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

			/** Construct an invalid EdgeIterator. */
			EdgeIterator() {}

			// dereference operator
			Edge operator*() const
			{ 
				return graph_ptr->edges.at(current_index);
			}
			
			// pre increment operator
			EdgeIterator& operator++()
			{
				if(current_index < graph_ptr->num_edges())
					current_index++;
				return *this;
			}
			
			// equality operator
			bool operator==(const EdgeIterator& iter) const
			{ 
				return  current_index == iter.current_index;
			}

		private:
			size_type current_index;
			graph_type* graph_ptr;

			// private constructor
			EdgeIterator(size_type i, const graph_type* g)
			{
				current_index = i;
				graph_ptr = const_cast<graph_type*>(g);
			}
			friend class Graph;
	};


		/** Return the number of nodes in the graph.
		 *
		 * Complexity: O(1).
		 */
		size_type size() const 
		{
			return nodes.size();
		}

		/** Synonym for size(). */
		size_type num_nodes() const 
		{
			return size();
		}

		/** Determine if a Node belongs to this Graph
		 * @return True if @a n is currently a Node of this Graph
		 *
		 * Complexity: O(1).
		 */
		bool has_node(const Node& n) const 
		{
			return adj_list.find(n.node_index) != adj_list.end();
		}

		/** Add a node to the graph, returning the added node.
		 * @param[in] position The new node's position
		 * @post new num_nodes() == old num_nodes() + 1
		 * @post result_node.index() == old num_nodes()
		 *
		 * Complexity: O(1) amortized operations.
		 */
		Node add_node(const Point& position, const node_value_type& val = node_value_type()) 
		{
			// create a node
			size_type node_index = num_inserted_nodes;
			node_type new_node = Node(node_index, this);
			
			node_uid_map[node_index] = num_nodes();
			nodes.push_back(new_node);
			node_positions[node_index] = position;
			node_values[node_index] = val;
			
			adj_list[node_index] = std::vector<edge_type>();
			
			num_inserted_nodes++;
			
			return new_node;       
		}

		/** Return the node with index @a i.
		 * @pre 0 <= @a i < num_nodes()
		 * @post result_node.index() == i
		 *
		 * Complexity: O(1).
		 */
		Node node(size_type i) const 
		{
			return nodes.at(i);
		}

		/** Return the total number of edges in the graph.
		 *
		 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
		 */
		size_type num_edges() const 
		{
			return edges.size();
		}

		/** Return the edge with index @a i.
		 * @pre 0 <= @a i < num_edges()
		 *
		 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
		 */
		Edge edge(size_type i) const 
		{
			return edges.at(i);
		}

		/** Test whether two nodes are connected by an edge.
		 * @pre @a a and @a b are valid nodes of this graph
		 * @return True if for some @a i, edge(@a i) connects @a a and @a b.
		 *
		 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
		 */
		bool has_edge(const Node& a, const Node& b) const 
		{
			size_type node_1_index = std::min(a.node_index, b.node_index);
			size_type node_2_index = std::max(a.node_index, b.node_index);
			node_tuple tup = std::make_tuple(node_1_index, node_2_index);

			return edge_map.find(tup) != edge_map.end();
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
		Edge add_edge(const Node& a, const Node& b, edge_value_type val = edge_value_type()) 
		{
			size_type node_1_index = std::min(a.node_index, b.node_index);
			size_type node_2_index = std::max(a.node_index, b.node_index);
			node_tuple tup = std::make_tuple(node_1_index, node_2_index);

			if (has_edge(a,b))
				return Edge(node_1_index, node_2_index, edge_map.at(tup), this);

			size_type edge_index = num_inserted_edges;
			edge_type new_edge = Edge(node_1_index, node_2_index, edge_index, this);
			edge_uid_map[edge_index] = num_edges();
			edge_values[edge_index] = val;
			edge_map[tup] = edge_index;
			adj_list.at(node_1_index).push_back(new_edge);
			adj_list.at(node_2_index).push_back(new_edge);
			edges.push_back(new_edge);

			num_inserted_edges++;
			return new_edge;
		}

		// begin() method for edge_iterator
		edge_iterator edge_begin() const
		{
			return EdgeIterator(0, this);
		}

		// end() method for edge_iterator
		edge_iterator edge_end() const
		{
			return EdgeIterator(num_edges(), this);
		}

		// begin() method for edge_iterator
		node_iterator node_begin() const
		{
			return NodeIterator(0, this);
		}

		// end() method for edge_iterator
		node_iterator node_end() const
		{
			return NodeIterator(num_nodes(), this);
		}


		edge_iterator remove_edge(edge_iterator e_it)
		{
			std::vector<size_type> indices {};
			indices.push_back((*e_it).node1().node_index);
			indices.push_back((*e_it).node2().node_index);

			size_type first = (*e_it).node1().node_index;
			size_type second = (*e_it).node2().node_index;

			size_type edge_idx = (*e_it).edge_index;
			size_type uid_current = edge_uid_map.at(edge_idx);

			// step 1 : update adjacency lists
			for(auto& x : indices)
			{
				unsigned n = adj_list.at(x).size();
				for(unsigned i = 0; i < n; i++)
				{
					if(adj_list.at(x).at(i) == *e_it)
					{
						std::swap(adj_list.at(x).at(i), adj_list.at(x).at(n-1));
						adj_list.at(x).pop_back();
						break;
					}
				}
			}


			// step 2 : update edge_values
			edge_values.erase(edge_idx);
			// step 3 : update edge_map
			node_tuple tup = std::make_tuple(first, second);
			edge_map.erase(tup);;

			// step 4 : update edge_uid_map and edges
			edge_uid_map.erase(edge_idx);
			// *e_it = edges.back();
			std::swap(edges[uid_current], edges[num_edges() -1]);
			edges.pop_back();

			// if there are no edges, don't update uid
			if(num_edges() > 0 and uid_current < num_edges())
			{	
				edge_uid_map.at(edges[uid_current].edge_index) = uid_current;
			}

			if(uid_current == num_edges())
				return EdgeIterator(num_edges()-1, this);

			return  EdgeIterator(uid_current, this);
		}


		size_type remove_edge(const edge_type& e)
		{
			if(has_edge(e.node1(), e.node2()))
			{
				size_type uid = edge_uid_map.at(e.edge_index);
				remove_edge(EdgeIterator(uid, this));
			}
			return 0;
		}


		size_type remove_edge(const node_type& a, const node_type& b)
		{
			if (has_edge(a,b))
				remove_edge(add_edge(a,b));
 			return 0;
		}

		size_type remove_node(const node_type& n)
		{
			if(has_node(n))
			{
				size_type uid = node_uid_map.at(n.node_index);
				remove_node(NodeIterator(uid, this));
			}
			return 0;
		}

		node_iterator remove_node(node_iterator n_it)
		{
			node_type n = *n_it;
			size_type idx = n.node_index;
			size_type uid = node_uid_map.at(idx);

			std::vector<edge_type> edges_for_removal {};
			for(unsigned i = 0; i < adj_list.at(idx).size(); i++)
				edges_for_removal.push_back(adj_list.at(idx).at(i));

			// step 1 : remove all incident edges
			for(size_type i = 0; i < edges_for_removal.size(); i++)
			{
				edge_type e = edges_for_removal[i];
				remove_edge(e);
			}
		
			// step 2 : remove node from adjacency list
			adj_list.erase(idx);

			// step 3 : update values and positions
			node_values.erase(idx);
			node_positions.erase(idx);


			// step 4 : update edge_uid_map and edges
			node_uid_map.erase(idx);
			std::swap(nodes[uid], nodes[num_nodes() -1]);
			nodes.pop_back();

			// if there are no edges, don't update uid
			if(num_nodes() > 0 and uid < num_nodes())
			{	
				node_uid_map.at(nodes[uid].node_index) = uid;
			}

			if(uid == num_nodes())
			{
				return NodeIterator(uid-1, this);
			}
			return NodeIterator(uid, this);

		}


		void clear()
		{
			num_inserted_edges = 0;
			num_inserted_nodes = 0;

			nodes.clear();
			node_values.clear();
			node_positions.clear();
			node_uid_map.clear();

			adj_list.clear();

			edges.clear();
			edge_values.clear();
			edge_uid_map.clear();
			edge_map.clear();
		}

	 private:

		// stores valid nodes in the graph
		std::vector< node_type > nodes {};
		// stores values associated with the nodes
		std::unordered_map< size_type, node_value_type > node_values;
		// stores 3D position of the node
		std::unordered_map< size_type, Point > node_positions;
		// stores mapping from node_index to node_uid
		std::unordered_map< size_type, size_type > node_uid_map;

		// stores valid edges of the graph
		std::unordered_map< size_type, std::vector < edge_type > > adj_list;

		// stores a vector of edges
		std::vector<edge_type> edges {};
		// stores values associated with edges
		std::unordered_map< size_type, edge_value_type > edge_values;
		// stores mapping from edge_index to <node_1_index, node_2_index>
		std::map< node_tuple, size_type > edge_map;
		// stores mapping from edge_index to edge_uid
		std::unordered_map<size_type, size_type> edge_uid_map;
		
		// total number of nodes added from creation of the graph
		size_type num_inserted_nodes {};
		// total number of edges added from creation of the graph
		size_type num_inserted_edges {}; 

};
#endif // CME212_GRAPH_HPP


