#ifndef __GRAPH_HPP
#define __GRAPH_HPP

#include <exception>
#include <stdexcept>
#include <algorithm> // std::find, std::remove, std::remove_if
#include <vector>
#include <queue>
#include <string>
#include <utility> // std::pair
#include "constants.hpp"

/* <-------- Exceptions --------> */

class GraphException : public std::logic_error {
protected:
	GraphException(const std::string& message) : std::logic_error(message) {}
public:
	virtual bool has_trigger() const { return false; }
};

template<typename T>
class HasTriggerGraphException : public GraphException {
	T _trigger;
protected:
	HasTriggerGraphException(const T& t, const std::string& message) : 
		GraphException(error_message::format(message, t)), _trigger(t) {}
public:
	virtual bool has_trigger() const { return true; }
	virtual T trigger() const { return _trigger; }
};

template<typename T>
class VertexNotFoundException : public HasTriggerGraphException<T> {
public:
	VertexNotFoundException(const T& t, const std::string& message) : 
		HasTriggerGraphException<T>(t, "VertexNotFoundException: " + message) {}
};

template<typename T>
class EdgeNotFoundException : public HasTriggerGraphException<T> {
public:
	EdgeNotFoundException(const T& t, const std::string& message) : 
		HasTriggerGraphException<T>(t, "EdgeNotFoundException: " + message) {}
};

template<typename T>
class VertexExistsException : public HasTriggerGraphException<T> {
public:
	VertexExistsException(const T& t, const std::string& message) : 
		HasTriggerGraphException<T>(t, "VertexExistsException: " + message) {}
};

template<typename T>
class EdgeExistsException : public HasTriggerGraphException<T> {
public:
	EdgeExistsException(const T& t, const std::string& message) : 
		HasTriggerGraphException<T>(t, "EdgeExistsException: " + message) {}
};

template<typename T>
class IncidentVertexNotFoundException : public HasTriggerGraphException<T> {
public:
	IncidentVertexNotFoundException(const T& t, const std::string& message) :
		HasTriggerGraphException<T>(t, "IncidentVertexNotFoundException: " + message) {}
};

/* <----------------------------> */

template<typename VertexElementBase>
class Edge{
	const VertexElementBase * const _from;
	const VertexElementBase * const _to;
	std::string _label;
	
	double* _weight;

public:
	explicit Edge(const VertexElementBase& from, const VertexElementBase& to, const double& weight) :
		_from(&from), _to(&to), _label(""), _weight(new double(weight)) {}

	explicit Edge(const VertexElementBase& from, const VertexElementBase& to, const std::string& label = "") :
		_from(&from), _to(&to), _label(label), _weight(nullptr) {}

	explicit Edge(const VertexElementBase& from, const VertexElementBase& to, const std::string& label, const double& weight) :
		_from(&from), _to(&to), _label(label), _weight(new double(weight)) {}


	Edge(const Edge<VertexElementBase>& other) : _from(other._from), _to(other._to), _label(other._label), _weight(nullptr) {
		if(other._weight != nullptr) {
			_weight = new double(*other._weight);
		}
	}

	Edge(Edge<VertexElementBase>&& other) : _from(other._from), _to(other._to), _label(other._label), _weight(other._weight) {
		other._weight = nullptr;
	}

	Edge<VertexElementBase>& operator=(const Edge<VertexElementBase>& other){
		if(this != &other){
			if(_weight != nullptr){
				delete _weight;
				_weight = nullptr;
			}
			_from = other._from; _to = other._to;
			_label = other._label;
			if(other._weight != nullptr){
				_weight = new double(*other._weight);
			}
		}
		return *this;
	}

	Edge<VertexElementBase>& operator=(Edge<VertexElementBase>&& other){
		if(this != other){
			if(_weight != nullptr){
				delete _weight;
			}
			_from = other._from; _to = other._to;
			_label = other._label; _weight = other._weight;
			other._weight = nullptr;
		}
		return *this;
	}

	inline bool operator==(const Edge& other) const {
		return (*_from == *(other.from())) && (*_to == *(other.to()));
	}

	bool incidents(const VertexElementBase& vertex) const {
		return incidents_from(vertex) || incidents_to(vertex);
	}

	bool incidents_from(const VertexElementBase& vertex) const {
		return *_from == vertex;
	}

	bool incidents_to(const VertexElementBase& vertex) const{
		return *_to == vertex;
	}

	const VertexElementBase* from() const { return _from; }
	const VertexElementBase* to() const { return _to; }
	std::string label() const { return _label; }
	double* weight() const { 
		return _weight;
	}

	void set_label(const std::string& label) { _label = label; }
	void set_weight(const double& weight) { 
		if(_weight == nullptr){
			_weight = new double;
		}
		*(_weight) = weight;
	}

	virtual ~Edge(){
		if(_weight != nullptr){
			delete _weight;
		}
	}
};

template<typename VertexElementBase>
std::ostream& operator<<(std::ostream& out, const Edge<VertexElementBase>& edge){
	out << *(edge.from()) << " -> " << *(edge.to());
	return out;
}

template<typename VertexElementBase>
class Graph{
	std::vector<VertexElementBase*> _vertices;
	std::vector<Edge<VertexElementBase>*> _edges;

	typename std::vector<VertexElementBase*>::const_iterator _find_vertex(const VertexElementBase& vertex) const {
		return std::find_if(_vertices.begin(), _vertices.end(),
							[&vertex](const VertexElementBase* maybe_vertex){
								return *maybe_vertex == vertex;
							});
	}

	typename std::vector<Edge<VertexElementBase>*>::const_iterator _find_edge(const Edge<VertexElementBase>& edge) const {
		return std::find_if(_edges.begin(), _edges.end(),
							[&edge](const Edge<VertexElementBase>* maybe_edge){
								return *maybe_edge == edge;
							});
	}

	bool _adjacent(const VertexElementBase& first, const VertexElementBase& second) const {
		return has_edge(Edge<VertexElementBase>(first, second)) || has_edge(Edge<VertexElementBase>(second, first));
	}

	std::vector<VertexElementBase*> _all_vertices() const {
		return _vertices;
	}

	std::vector<Edge<VertexElementBase>*> _all_edges() const {
		return _edges;
	}

	/* Adds a vertex if it is not contained by the graph,
	else throws an exception. */
	template<typename VertexElementDerived>
	void _add_vertex(const VertexElementDerived& vertex){
		if(! has_vertex(vertex)){
			_vertices.reserve(_vertices.size() + 1); // Avoid memory leak if push_back throws an exception
			_vertices.push_back(new VertexElementDerived(vertex));
		}
		
		else{
			throw VertexExistsException<VertexElementDerived>(vertex, error_message::kAddedVertexExists);
		}
	}

	/* Removes a vertex if it is contained by the graph,
	else throws an exception. */
	void _remove_vertex(const VertexElementBase& vertex){
		if(has_vertex(vertex)){
			_remove_all_edges(vertex);
			_vertices.erase(std::remove_if(_vertices.begin(), _vertices.end(), 
						[&vertex](const VertexElementBase* maybe_vertex){
							if(vertex == *maybe_vertex){
								delete maybe_vertex;
								return true;
							}
							return false;
						}), _vertices.end());
		}
		else{
			throw VertexNotFoundException<VertexElementBase>(vertex, error_message::kRemoveVertexNotFound);
		}
	}

	/* Removes the edges that are incident to the given vertex. */
	void _remove_all_edges(const VertexElementBase& vertex){
		_edges.erase(std::remove_if(_edges.begin(), _edges.end(),
									[&vertex](const Edge<VertexElementBase>* edge){
										if(edge->incidents(vertex)){
											delete edge;
											return true;
										}
										return false;
									}), _edges.end());	
	}

	/* Removes an edge if it is contained by the graph,
	else throws an exception. */
	void _remove_edge(const Edge<VertexElementBase>& edge){
		if(has_edge(edge)){
			_edges.erase(std::remove_if(_edges.begin(), _edges.end(), 
						[&edge](const Edge<VertexElementBase>* maybe_edge){
							if(edge == *maybe_edge){
								delete maybe_edge;
								return true;
							}
							return false;
						}), _edges.end());
		}
		else{
			throw EdgeNotFoundException<Edge<VertexElementBase>>(edge, error_message::kRemoveEdgeNotFound);
		}
	}

	/* Adds an edge if it is not contained by the graph,
	else throws an exception. 
	If an incident vertex is not contained by the graph, 
	throws en exception. */
	void _add_edge(const Edge<VertexElementBase>& edge){
		typename std::vector<VertexElementBase*>::const_iterator it_from;
		typename std::vector<VertexElementBase*>::const_iterator it_to;
		if((it_from = _find_vertex(*(edge.from()))) == _vertices.end()) throw IncidentVertexNotFoundException<VertexElementBase>(*(edge.from()), error_message::kIncidentVertexNotFound);
		if((it_to = _find_vertex(*(edge.to()))) == _vertices.end()) throw IncidentVertexNotFoundException<VertexElementBase>(*(edge.to()), error_message::kIncidentVertexNotFound);
		
		if(! has_edge(edge)){
			_edges.reserve(_edges.size() + 1);
			(edge.weight() != nullptr) ? _edges.push_back(new Edge<VertexElementBase>(**it_from, **it_to, *edge.weight())) : 
										 _edges.push_back(new Edge<VertexElementBase>(**it_from, **it_to));
		}
		else{
			throw EdgeExistsException<Edge<VertexElementBase>>(edge, error_message::kAddedEdgeExists);
		}
	}

	/* Returns the successors of the given vertex. Throws an exception if the 
	given vertex is not contained by the graph. */
	std::vector<VertexElementBase*> _get_out_vertices(const VertexElementBase& vertex) const {
		if(has_vertex(vertex)){
			std::vector<VertexElementBase*> out_vertices;
			std::for_each(_edges.begin(), _edges.end(),
							[&vertex, &out_vertices](const Edge<VertexElementBase>* edge){
								if(edge->incidents_from(vertex)){
									out_vertices.push_back(const_cast<VertexElementBase*>(edge->to()));
								}
							});
			return out_vertices;	
		}
		else{
			throw VertexNotFoundException<VertexElementBase>(vertex, error_message::kGetVertexNotFound);
		}
	}

	/* Returns the predecessors of the given vertex. Throws an exception if the 
	given vertex is not contained by the graph. */
	std::vector<VertexElementBase*> _get_in_vertices(const VertexElementBase& vertex) const {
		if(has_vertex(vertex)){
			std::vector<VertexElementBase*> in_vertices;
			std::for_each(_edges.begin(), _edges.end(),
							[&vertex, &in_vertices](const Edge<VertexElementBase>* edge){
								if(edge->incidents_to(vertex)){
									in_vertices.push_back(const_cast<VertexElementBase*>(edge->from()));
								}
							});
			return in_vertices;
		}
		else{
			throw VertexNotFoundException<VertexElementBase>(vertex, error_message::kGetVertexNotFound);
		}
	}

	/* Returns the neighbours of the given vertex. Duplicate neighbours are ignored.
	Throws an exception if the given vertex is not contained by the graph. */
	std::vector<VertexElementBase*> _get_neighbours(const VertexElementBase& vertex) const {
		std::vector<VertexElementBase*> out = _get_out_vertices(vertex);
		std::vector<VertexElementBase*> in = _get_in_vertices(vertex);
		out.reserve(out.size() + in.size());
		std::for_each(in.begin(), in.end(),
						[&out](VertexElementBase* in_vertex){
							if(std::find_if(out.begin(), out.end(),
								[in_vertex](VertexElementBase* out_vertex) {
									return *in_vertex == *out_vertex;
								}) == out.end()) {
									out.push_back(in_vertex);
							}
						});
		out.shrink_to_fit();
		return out;
	}

	std::vector<Edge<VertexElementBase>*> _get_out_edges(const VertexElementBase& vertex) const {
		if(has_vertex(vertex)){
			std::vector<Edge<VertexElementBase>*> out_edges;
			std::for_each(_edges.begin(), _edges.end(),
							[&vertex, &out_edges](Edge<VertexElementBase>* edge){
								if(edge->incidents_from(vertex)){
									out_edges.push_back(edge);
								}
							});
			return out_edges;
		}
		else{
			throw VertexNotFoundException<VertexElementBase>(vertex, error_message::kGetOutEdgesVertexNotFound);
		}
	}

	std::vector<Edge<VertexElementBase>*> _get_in_edges(const VertexElementBase& vertex) const {
		if(has_vertex(vertex)){
			std::vector<Edge<VertexElementBase>*> in_edges;
			std::for_each(_edges.begin(), _edges.end(),
							[&vertex, &in_edges](Edge<VertexElementBase>* edge){
								if(edge->incidents_to(vertex)){
									in_edges.push_back(edge);
								}
							});
			return in_edges;
		}
		else{
			throw VertexNotFoundException<VertexElementBase>(vertex, error_message::kGetInEdgesVertexNotFound);
		}
	}

	Graph<VertexElementBase> _sub_graph(const std::vector<VertexElementBase>& vertices) const {
		Graph<VertexElementBase> sub;
		/* Add vertices. */
		for(const VertexElementBase& vertex : vertices){
			sub.add_vertex(vertex);
		}
		/* Add edges. */
		for(const VertexElementBase& vertex_from : vertices){
			for(const VertexElementBase& vertex_to : vertices){
				typename std::vector<Edge<VertexElementBase>*>::const_iterator it;
				if((it = _find_edge(Edge<VertexElementBase>(vertex_from, vertex_to))) != _edges.end()) {
					sub._add_edge(**it);
				}
			}
		}
		return sub;
	}

	void _topological_sort() {
		std::map<VertexElementBase*, std::size_t> pred;
		std::vector<VertexElementBase*> L;
		std::queue<VertexElementBase*> Q;
		for(VertexElementBase* vertex : _vertices){
			std::size_t num_preds = get_in_edges(*vertex).size(); 
			pred[vertex] = num_preds;
			if(num_preds == 0) Q.push(vertex);
		}
		while(! Q.empty()){
			VertexElementBase* vertex = Q.front();
			Q.pop();
			L.push_back(vertex);
			for(Edge<VertexElementBase>* edge : get_out_edges(*vertex)){
				VertexElementBase* to_vertex = const_cast<VertexElementBase*>(edge->to()); 
				std::size_t& num_preds = pred[to_vertex];
				--num_preds;

				if(num_preds == 0) Q.push(to_vertex);
			}
		}
		if(_vertices.size() != L.size()) throw std::logic_error("fail toposort");
		else _vertices = L;
	}

	/* Removes all the edges. Deletes the dynamically allocated edges before
	calling the standard clear on the vector. */
	void _clear_all_edges() {
		std::for_each(_edges.begin(), _edges.end(), 
						[](const Edge<VertexElementBase>* edge){
							delete edge;
						});
		_edges.clear();
	}

	/* Removes all the vertices. Deletes the dynamically allocated vertices before
	calling the standard clear on the vector. */
	void _clear_all_vertices() {
		_clear_all_edges();
		std::for_each(_vertices.begin(), _vertices.end(), 
						[](const VertexElementBase* vertex){
							delete vertex;
						});
		_vertices.clear();
	}

public:
	/* Graph interface */
	Graph() : _vertices(), _edges() {}

	Graph(const Graph<VertexElementBase>& other) : _vertices(), _edges() {
		_vertices.reserve(other._vertices.size());
		for(std::size_t i = 0; i < other._vertices.size(); ++i){
			_vertices[i] = new VertexElementBase(*(other._vertices[i]));
		}
		_edges.reserve(other._edges.size());
		for(std::size_t i = 0; i < other._edges.size(); ++i){
			_edges[i] = new Edge<VertexElementBase>(*(other._edges[i]));
		}
	}

	Graph(Graph<VertexElementBase>&& other) :
		_vertices(std::move(other._vertices)), _edges(std::move(other._edges)) {}

	std::size_t num_vertices() const { return _vertices.size(); }
	std::size_t num_edges() const { return _edges.size(); }

	std::vector<VertexElementBase*> get_vertices() const { 
		return _all_vertices();
	}

	std::vector<Edge<VertexElementBase>*> get_edges() const {
		return _all_edges();
	}

	bool has_vertex(const VertexElementBase& vertex) const {
		return _find_vertex(vertex) != _vertices.end();
	}

	bool has_edge(const Edge<VertexElementBase>& edge) const {
		return _find_edge(edge) != _edges.end();
	}

	bool has_edge(const VertexElementBase& from_vertex, const VertexElementBase& to_vertex) const {
		return has_edge(Edge<VertexElementBase>(from_vertex, to_vertex));
	}

	bool adjacent(const VertexElementBase& first, const VertexElementBase& second) const {
		return _adjacent(first, second);
	}

	VertexElementBase* get_vertex(const VertexElementBase& vertex) const {
		typename std::vector<VertexElementBase*>::const_iterator it;
		if((it = _find_vertex(vertex)) != _vertices.end()) {
			return *it;
		}
		else{
			throw VertexNotFoundException<VertexElementBase>(vertex, error_message::kGetVertexNotFound);
		}
	}

	void add_vertex(const VertexElementBase& vertex){
		_add_vertex(vertex);
	}

	void remove_vertex(const VertexElementBase& vertex){
		_remove_vertex(vertex);
	}

	void add_edge(const VertexElementBase& from, const VertexElementBase& to){
		_add_edge(Edge<VertexElementBase>(from, to));
	}

	void add_edge(const VertexElementBase& from, const VertexElementBase& to, const std::string& label){
		_add_edge(Edge<VertexElementBase>(from, to, label));
	}

	void add_edge(const VertexElementBase& from, const VertexElementBase& to, const double& weight){
		_add_edge(Edge<VertexElementBase>(from, to, weight));
	}

	void add_edge(const VertexElementBase& from, const VertexElementBase& to, const std::string& label, const double& weight){
		_add_edge(Edge<VertexElementBase>(from, to, label, weight));
	}

	void remove_edge(const VertexElementBase& from, const VertexElementBase& to){
		_remove_edge(Edge<VertexElementBase>(from, to));
	}

	std::vector<VertexElementBase*> get_out_vertices(const VertexElementBase& vertex) const {
		return _get_out_vertices(vertex);
	}	

	std::vector<VertexElementBase*> get_in_vertices(const VertexElementBase& vertex) const {
		return _get_in_vertices(vertex);
	}

	std::vector<VertexElementBase*> get_neighbours(const VertexElementBase& vertex) const {
		return _get_neighbours(vertex);
	}

	std::vector<Edge<VertexElementBase>*> get_in_edges(const VertexElementBase& vertex) const {
		return _get_in_edges(vertex);
	}

	std::vector<Edge<VertexElementBase>*> get_out_edges(const VertexElementBase& vertex) const {
		return _get_out_edges(vertex);
	}

	Graph<VertexElementBase> sub_graph(const std::vector<VertexElementBase>& vertices) const {
		return _sub_graph(vertices);
	}

	void topological_sort() {
		_topological_sort();
	}

	void clear_all_edges() {
		_clear_all_edges();
	}

	void clear_all_vertices() {
		_clear_all_vertices();
	}
	
	virtual ~Graph() { 
		_clear_all_edges();
		_clear_all_vertices();
	}
};

#endif