


#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP





/** @file Graph.hpp
 * @brief An undirected graph typ
 */



#include <algorithm>
#include <vector>
#include <cassert>
#include <map>


#include "CME212/Util.hpp"
#include "CME212/Point.hpp"



using namespace std;



/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template  <typename V, typename E>
class Graph {

	private:

		struct internal_node;
		struct internal_edge;


		// HW0: YOUR CODE HERE
		// Use this space for declarations of important internal types you need
		// later in the Graph's definition.
		// (As with all the "YOUR CODE HERE" markings, you may not actually NEED
		// code here. Just use the space if you need it.

	public:

		//
		// PUBLIC TYPE DEFINITIONS
		//

		/** Type of this graph. */
		using graph_type = Graph;

		/** Type of node value to support a user specified value  */
		using  node_value_type = V;

                /** Type of edge value to support a user specified value  */
                using  edge_value_type = E;

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
		Graph() : vect_nodes_(), vect_edges_(), map_adj_(), i2u_(), edge_i2u_() {

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

		class Node : private totally_ordered<Node>  {
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
				 *   x= some other node using a complicated calculation
				 * do_something(x);
				 * @endcode
				 */

				Node() {

				}



				/** Return this node's position. */
				const Point& position() const {

					return graph_->vect_nodes_[node_id_].position_;

				}

                                /** Return the reference of the node's position. */
                                Point& position() {

                                        return graph_->vect_nodes_[node_id_].position_;

                                }



				/** Return this node's index, a number in the range [0, graph_size). */
				size_type index() const {

					return graph_->vect_nodes_[node_id_].idx_;

				}

				// HW1: YOUR CODE HERE
				// Supply definitions AND SPECIFICATIONS for:
				// node_value_type& value();
				// const node_value_type& value() const;
				// size_type degree() const;
				// incident_iterator edge_begin() const;
				// incident_iterator edge_end() const;

				// @brief this method returns the reference of the node value type  and it can be changed
				node_value_type& value(){

					return graph_->vect_nodes_[node_id_].value_;

				}

				// @brief this method returns the reference of the node value type and it cannot be changed
				const  node_value_type& value() const  {

					return graph_->vect_nodes_[node_id_].value_;

				}

				//  @brief this method returns the number of incident edges
				size_type  degree() const{

					return graph_->map_adj_[node_id_].size();

				}

				// @brief  Start  of the  incident  iterator.
				incident_iterator  edge_begin() const{

					return IncidentIterator(graph_,node_id_,0);

				}

				// @brief  End of  incident  iterator
				incident_iterator  edge_end() const{

					return IncidentIterator(graph_,node_id_,graph_->map_adj_[node_id_].size());

				}




				/** Test whether this node and @a n are equal
				 *
				 * Equal nodes have the same graph and the same index.
				 */

				bool operator==(const Node& n) const {

					return  (graph_==n.graph_ && node_id_==n.node_id_);

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

					return  (node_id_<n.node_id_);

				}



			private:

				// Allow Graph to access Node's private member data and functions.

				friend class Graph;
				friend class NodeIterator;

				// Use this space to declare private data members and methods for Node
				// that will not be visible to users, but may be useful within Graph.
				// i.e. Graph needs a way to construct valid Node objects

				graph_type* graph_;
				size_type node_id_;



				/** Private Constructor */

				Node(const Graph* graph, size_type node_idx)

					: graph_(const_cast<graph_type*>(graph)), node_id_(node_idx) {

					}


		};

		/** Return the number of nodes in the graph.
		 *
		 * Complexity: O(1).
		 */

		size_type size() const {
			return i2u_.size() ;

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

			vect_nodes_.push_back(internal_node(num_nodes(),position,value));
			//vect_adj_.push_back(std::vector<internal_edge>());
                        i2u_.push_back(vect_nodes_.size()-1);
			return Node(this,num_nodes()-1);

		}


		/** Determine if a Node belongs to this Graph
		 * @return True if @a n is currently a Node of this Graph
		 *
		 * Complexity: O(1).
		 */

		bool has_node(const Node& n) const {

			if (this==n.graph_){
				return true;
			}
			else{
				return false;
			}
		}


		/** Return the node with index @a i.
		 * @pre 0 <= @a i < num_nodes()
		 * @post result_node.index() == i
		 *
		 * Complexity: O(1).
		 */

		Node node(size_type i) const {

			// HW0: YOUR CODE HERE

			assert(0<=i && i<num_nodes());

			return Node(this,i2u_[i]);

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

		class Edge : private totally_ordered<Edge>  {
			public:

				/** Construct an invalid Edge. */
				Edge() {

				}

				/** Return a node of this Edge */
				Node node1() const {

					return Node(graph_,node1_idx_);

				}



				/** Return the other node of this Edge */
				Node node2() const {

					return Node(graph_,node2_idx_);

				}



				/** Test whether this edge and @a e are equal.
				 *
				 * Equal edges represent the same undirected edge between two nodes.
				 */

				bool operator==(const Edge& e) const  {
					// Here we check if this edge  and @ e  have the same nodes 
					// instead of comparing only the edges' indices , using the operator 
					//  == of class Node

					if (node1() == e.node1()  && node2() == e.node2() ){
						return true;
					}
					else if (node1() == e.node2()  && node2() == e.node1() ){
						return true;
					}
					else {
						return false;
					}

				}

				/** Test whether this edge is less than @a e in a global order.
				 *
				 * This ordering function is useful for STL containers such as
				 * std::map<>. It need not have any interpretive meaning.
				 */

				bool operator<(const Edge& e) const {

					if (length()<e.length()){
						return true;
					}
					else {
						return false;
					}
				}


                               edge_value_type& value(){
                                   return graph_->vect_edges_[edge_id_].value_;
                              }

                               const  edge_value_type& value() const  {
                                   return graph_->vect_edges_[edge_id_].value_;
                              }

                              double length() const {
                                    return norm(graph_->vect_nodes_[node1_idx_].position_-graph_->vect_nodes_[node2_idx_].position_);
                              }

			private:

				// Allow Graph to access Edge's private member data and functions.

				friend class Graph;

				// Use this space to declare private data members and methods for Edge
				// that will not be visible to users, but may be useful within Graph.
				// i.e. Graph needs a way to construct valid Edge objects

				graph_type* graph_;
				size_type edge_id_;
				size_type node1_idx_;
				size_type node2_idx_;

				/** Private Constructor */

				Edge(const Graph* graph, size_type edge_idx, size_type node1_idx, size_type node2_idx)

					: graph_(const_cast<graph_type*>(graph)), edge_id_(edge_idx), node1_idx_(node1_idx), node2_idx_(node2_idx) {

					}



		};



		/** Return the total number of edges in the graph.
		 *
		 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
		 */

		size_type num_edges() const {

			return edge_i2u_.size();

		}



		/** Return the edge with index @a i.
		 * @pre 0 <= @a i < num_edges()
		 *
		 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
		 */

		Edge edge(size_type i) const {

			assert(0<=i && i<num_edges());
                        size_type id = edge_i2u_[i];
			return Edge(this,id,vect_edges_[id].node1_idx_,vect_edges_[id].node2_idx_);

		}



		/** Test whether two nodes are connected by an edge.
		 * @pre @a a and @a b are valid nodes of this graph
		 * @return True if for some @a i, edge(@a i) connects @a a and @a b.
		 *
		 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
		 */

		bool has_edge(const Node& a, const Node& b) const {

			// HW0: YOUR CODE HERE
                        assert(vect_edges_.size()>=num_edges());
			size_type  i=0;
			while (i<num_edges()){
                                //std::cout<<"1";
                                size_type id = edge_i2u_[i];
				if (vect_edges_[id].node1_idx_ == a.node_id_ && vect_edges_[id].node2_idx_ == b.node_id_){
					return true;
				}
				else if (vect_edges_[id].node1_idx_ == b.node_id_ && vect_edges_[id].node2_idx_ == a.node_id_){
					return true;
				}
				i++;
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

		Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {

			if (has_edge(a,b) == true){

				// Here, we look for the index of this edge, in order to return 
				// edge(this,index) 
				size_type i =0;
				while (i<num_edges()){
                                        size_type id = edge_i2u_[i];
					size_type node1_idx = vect_edges_[id].node1_idx_;
					size_type node2_idx = vect_edges_[id].node2_idx_;

					if (node1_idx == a.node_id_ && node2_idx == b.node_id_){
						break;
					}
					else if (node1_idx == b.node_id_ && node2_idx == a.node_id_){
						break;
					}
					i++;
				}
				return Edge(this,edge_i2u_[i],a.node_id_,b.node_id_);
			}

			else{ // Here the index of the added edge will be the old  num_edges()  and the num_edges()
				//  will increase by 1
				size_type NumEdges = edge_i2u_.size();
				vect_edges_.push_back(internal_edge(NumEdges,a.node_id_,b.node_id_,value));
                                // Here we add the internal edge to map_adj_ such that the root is node1
				map_adj_[a.node_id_].push_back(internal_edge(NumEdges,a.node_id_,b.node_id_,value));
				map_adj_[b.node_id_].push_back(internal_edge(NumEdges,b.node_id_,a.node_id_,value));
                                edge_i2u_.push_back(vect_edges_.size()-1);
                               // std::cout<<std::string(edge_i2u_[vect_edges_.size()-1])<<std::endl;
 				return Edge(this,vect_edges_.size()-1,a.node_id_,b.node_id_);

			}

		}


		/** Remove all nodes and edges from this graph.
		 * @post num_nodes() == 0 && num_edges() == 0
		 *
		 * Invalidates all outstanding Node and Edge objects.
		 */

		void clear() {

			vect_nodes_.clear();
			vect_edges_.clear();
			map_adj_.clear();
                        i2u_.clear();
                        edge_i2u_.clear();
		}



		//
		// Node Iterator
		//



		/** @class Graph::NodeIterator
		 * @brief Iterator class for nodes. A forward iterator. */

		class NodeIterator : private totally_ordered<NodeIterator>  {
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


				Node operator*() const {

					return graph_->node(node_iter_);

				}

				NodeIterator& operator++(){

					node_iter_++;
					return *this;

				}

				bool operator==(const NodeIterator& Nodeiterator) const{

					return (Nodeiterator.graph_==graph_ && Nodeiterator.node_iter_==node_iter_);

				}

			private:

				friend class Graph;

				// HW1 #2: YOUR CODE HERE
				graph_type* graph_;
				size_type node_iter_ ;

				NodeIterator(const graph_type* graph, size_type node_iter=0 ) : graph_(const_cast<graph_type*>(graph)), node_iter_(node_iter) {};


		};



		// HW1 #2: YOUR CODE HERE
		// Supply definitions AND SPECIFICATIONS for:
		// node_iterator node_begin() const
		// node_iterator node_end() const

		node_iterator node_begin() const  {

			return NodeIterator(this, 0);

		}

		node_iterator node_end()  const  {

			return NodeIterator(this, num_nodes());

		}



		//
		// Incident Iterator
		//



		/** @class Graph::IncidentIterator
		 * @brief Iterator class for edges incident to a node. A forward iterator. */

     		class IncidentIterator : private totally_ordered<IncidentIterator> {
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


				Edge operator*() const {
                                        size_type edge_idx = graph_->map_adj_[mother_node_id_][num_edge_].edge_idx_;
					return Edge(graph_, graph_->edge_i2u_[ edge_idx ] ,
                                                            mother_node_id_ ,
                                                            graph_->map_adj_[mother_node_id_][num_edge_].node2_idx_);

				}


				IncidentIterator& operator++(){

					++num_edge_;
					return *this;

				}

				bool operator==(const IncidentIterator& inc_iter)  const {

					return (inc_iter.graph_==graph_ && inc_iter.mother_node_id_==mother_node_id_ && inc_iter.num_edge_==num_edge_);

				}



			private:

				friend class Graph;

				// HW1 #3: YOUR CODE HERE

				graph_type* graph_;
				size_type mother_node_id_;
				size_type num_edge_;

				IncidentIterator(const graph_type* graph, size_type mother_node_idx, size_type num_edge) :
					graph_(const_cast<graph_type*>(graph)) , mother_node_id_(mother_node_idx)  , num_edge_(num_edge) {};


		};



		//
		// Edge Iterator
		//



		/** @class Graph::EdgeIterator
		 * @brief Iterator class for edges. A forward iterator. */

		class EdgeIterator : private totally_ordered<EdgeIterator> {
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

				Edge operator*() const {

					return graph_->edge(edge_iter_);

				}

				EdgeIterator& operator++(){

					++edge_iter_;
					return *this;

				}

				bool operator==(const EdgeIterator Edgeiterator) const {

					return (Edgeiterator.graph_==graph_ && Edgeiterator.edge_iter_==edge_iter_);

				}

			private:

				friend class Graph;

				// HW1 #5: YOUR CODE HERE
				graph_type* graph_;
				size_type edge_iter_;


				EdgeIterator(const graph_type* graph, size_type edge_iter):
					graph_(const_cast<graph_type*>(graph)), edge_iter_(edge_iter) {}





		};



		// HW1 #5: YOUR CODE HERE
		// Supply definitions AND SPECIFICATIONS for:
		// edge_iterator edge_begin() const
		// edge_iterator edge_end() const

		edge_iterator edge_begin() const {

			return EdgeIterator(this,0);

		}


		edge_iterator edge_end() const {

			return EdgeIterator(this,num_edges());

		}



  	/** @brief This method removes an edge from the graph,
  	*  and return 1 if one of the edges is removed,
  	*  0 if no edge is removed
  	* @pre @a n1 @a n2 are valid nodes of this graph
  	* @return 1 if @pre has_edge(n1,n2), 0 otherwise
  	* @post has_edge(n1,n2) == false
  	* @post IF the edge is removed new edge_i2u_.size() = old edge_i2u_.size() - 1
  	*       ELSE,                 new edge_i2u_.size() = old edge_i2u_.size()
  	* Complexity: No more than O(num_nodes()+num_edges())
  	*/
               size_type  remove_edge(const  Node& n1 , const  Node & n2){

                   if (has_edge(n1,n2)){
                     size_type edge_id;
                     // THis is the index of the edge in edge_i2u_
                     size_type edge_idx;

                     for (auto e_it = map_adj_[n1.node_id_].begin(); e_it!= map_adj_[n1.node_id_].end(); ++e_it){

                             size_type neighbor_node_id = (*e_it).node2_idx_;
                             if (neighbor_node_id == n2.node_id_ ){

                                 edge_id = (*e_it).edge_idx_;
                                 edge_idx = vect_edges_[edge_id].edge_idx_;
                                 // we remove the edge from the map_adj_[node1_idx]
                                 map_adj_[n1.node_id_].erase(e_it);
                                 break;
                             }
                     }

                     for (auto e_it = map_adj_[n2.node_id_].begin(); e_it!= map_adj_[n2.node_id_].end(); ++e_it){
                             size_type neighbor_node_id = (*e_it).node2_idx_;
                             if (neighbor_node_id == n1.node_id_ ){
                                 // we remove the edge from the map_adj_[node2_idx]
                                 map_adj_[n2.node_id_].erase(e_it);
                                 break;
                             }
                     }


                    assert(edge_i2u_.size()>=edge_idx);
                    // we remove the edge_id from active edges
                    edge_i2u_.erase(edge_i2u_.begin()+edge_idx);
                    // we update the indices of edges in edge_i2u that are big than edge_idx
                    for (size_type idx = edge_idx; idx < edge_i2u_.size(); ++idx ){
                         vect_edges_[edge_i2u_[idx]].edge_idx_ = idx;
                     }

                   return 1;
                   }

              return 0;
              };

  	/** @brief This method removes an edge from the graph,
  	* for an edge @a given as argument
  	*  and return 1 if one of the edges is removed,
  	*  0 if no edge is removed
  	* @pre @a edge is a valid edge of this graph
        * @return 1 if @pre has_edge(@a.node1(),@a.node2()), 0 otherwise
  	* @post has_edge(@a.node1(),@a.node2()) == false
  	* @post IF the edge is removed new edge_i2u_.size() = old edge_i2u_.size() - 1
  	*       ELSE,                 new edge_i2u_.size() = old edge_i2u_.size()
  	* Complexity: No more than O(num_nodes()+num_edges())
  	*/


               size_type remove_edge(const Edge& edge){
                   return remove_edge(edge.node1(),edge.node2());
               }

  	/** @brief This method removes an edge from the graph,
  	* for an edge_iterator @a given as argument
  	*  and return 1 if one of the edges is removed,
  	*  0 if no edge is removed
  	* @pre @a e_it an edge_iterator
  	* 
  	* @post e_it is not a valid edge iterator
  	* @post IF the edge (*e_it) is removed new edge_i2u_.size() = old edge_i2u_.size() - 1
  	*       ELSE                 new edge_i2u_.size() = old edge_i2u_.size()
  	* Complexity: No more than O(num_nodes()+num_edges())
  	*/


              edge_iterator remove_edge(edge_iterator e_it){
                   remove_edge(*e_it);
                   return e_it;
              };

  	/** @brief This function removes a @a node given as argument
   	* from the graph (if it exists in the graph), and
   	* return  1 if @a is removed  OR 0.
   	*
   	* @pre @a n is a node of this graph
   	* @return 1 if n is a valid node of this graph, 0 if n is not a valid node of this graph
   	* @post IF @a has_node(a)  == true  @returns   1
   	*       ELSE                        @returns   0
   	*
   	* Complexity:No more than  O(num_nodes()).
	   */



             size_type remove_node(const Node& node){
               if (has_node(node)) {
                   while (node.degree() > 0){
                       remove_edge(*(node.edge_begin()));
                   }
                   // we remove the node from th vector of active nodes indices
                   i2u_.erase(i2u_.begin() + node.index());
                   // we update the active nodes indices for nodes with index bigger than node.index()
                   for (size_type idx = node.index(); idx!= num_nodes(); ++idx){
                      vect_nodes_[i2u_[idx]].idx_ = idx;
                  }

                   return 1;
               }
               return 0;
             };


  	/** Remove an node from the graph, and return the new location followed by the erased element
  	* @pre @a n_it is a node iterator
  	* @return the new iterator followed by the erase element
  	* @post n_it is not a valid node iterator
        * @post If @a n is removed, new num_nodes() == old num_nodes()-1.
        *       Else,               new num_nodes() == old num_nodes()
        * Complexity: No more than O(num_nodes())
        */
            node_iterator  remove_node(node_iterator  n_it){
                remove_node(*n_it);
                return n_it;
           };












	private:

		// Struct internal-node

		struct internal_node{

			size_type idx_;
			Point position_;
			node_value_type value_;

			internal_node(size_type idx,const  Point& position, node_value_type value): idx_(idx) , position_(position), value_(value)  {}

		};

		// Struct internal_edge

		struct internal_edge{

			size_type edge_idx_;
			size_type node1_idx_;
			size_type node2_idx_;
                        edge_value_type value_;

			internal_edge(size_type edge_idx, size_type node1_idx, size_type node2_idx, edge_value_type value): \
                                              edge_idx_(edge_idx) , node1_idx_(node1_idx), node2_idx_(node2_idx), value_(value) {}

		}; 


		// vector of nodes
		std::vector<internal_node> vect_nodes_;
		// vector of edges
		std::vector<internal_edge> vect_edges_;
		// map of adjencies
		std::map<size_type, std::vector<internal_edge> > map_adj_;
                // vector indexed by node_id_
                std::vector<size_type> i2u_;
                // vector indexed by edge_idx_
                std::vector<size_type> edge_i2u_;
};



#endif // CME212_GRAPH_HPP
