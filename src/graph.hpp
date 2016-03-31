#ifndef __GRAPH_HPP
#define __GRAPH_HPP

#include <exception>
#include <stdexcept>
#include <algorithm> // std::find, std::remove, std::remove_if
#include <vector>
#include <string>
#include <utility> // std::pair
#include "constants.hpp"

/* <-------- Exceptions --------> */

class GraphException : public std::logic_error {
protected:
	GraphException(const std::string& message) :
		std::logic_error(message) {}
};

class VertexNotFoundException : public GraphException {
public:
	template<typename T>
	VertexNotFoundException(const T& t) : 
		GraphException(error_message::format("VertexNotFoundException: " + error_message::kVertexNotFound, t)) {}
};

class EdgeNotFoundException : public GraphException {
public:
	template<typename T>
	EdgeNotFoundException(const T& t) : 
		GraphException(error_message::format("EdgeNotFoundException: " + error_message::kEdgeNotFound, t)) {}
};

class VertexExistsException : public GraphException {
public:
	template<typename T>
	VertexExistsException(const T& t) : 
		GraphException(error_message::format("VertexExistsException: " + error_message::kAddedVertexExists, t)) {}
};

class EdgeExistsException : public GraphException {
public:
	template<typename T>
	EdgeExistsException(const T& t) : 
		GraphException(error_message::format("EdgeExistsException: " + error_message::kAddedEdgeExists, t)) {}
};

class IncidentVertexNotFoundException : public GraphException {
public:
	template<typename T>
	IncidentVertexNotFoundException(const T& t) :
		GraphException(error_message::format("IncidentVertexNotFoundException: " + error_message::kIncidentVertexNotFound, t)) {}
};

/* <----------------------------> */

template<typename VertexElementBase>
class Edge{
	const VertexElementBase* _from;
	const VertexElementBase* _to;
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
		return contains(Edge<VertexElementBase>(first, second)) || contains(Edge<VertexElementBase>(second, first));
	}

	std::vector<VertexElementBase*>& _all_vertices() const {
		return _vertices;
	}

	std::vector<Edge<VertexElementBase>*>& _all_edges() const {
		return _edges;
	}

	/* Adds a vertex if it is not contained by the graph,
	else throws an exception. */
	template<typename VertexElementDerived>
	void _add_vertex(const VertexElementDerived& vertex){
		if(! contains(vertex)){
			_vertices.reserve(_vertices.size() + 1); // Avoid memory leak if push_back throws an exception
			_vertices.push_back(new VertexElementDerived(vertex));
		}
		
		else{
			throw VertexExistsException(vertex);
		}
	}

	/* Removes a vertex if it is contained by the graph,
	else throws an exception. */
	void _remove_vertex(const VertexElementBase& vertex){
		if(contains(vertex)){
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
			throw VertexNotFoundException(vertex);
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
		if(contains(edge)){
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
			throw EdgeNotFoundException(edge);
		}
	}

	/* Adds an edge if it is not contained by the graph,
	else throws an exception. 
	If an incident vertex is not contained by the graph, 
	throws en exception. */
	void _add_edge(const Edge<VertexElementBase>& edge){
		if(! contains(*(edge.from()))) throw IncidentVertexNotFoundException(*(edge.from()));
		if(! contains(*(edge.to()))) throw IncidentVertexNotFoundException(*(edge.to()));
		
		if(! contains(edge)){
			_edges.reserve(_edges.size() + 1);
			_edges.push_back(new Edge<VertexElementBase>(edge));
		}
		else{
			throw EdgeExistsException(edge);
		}
	}

	/* Returns the successors of the given vertex. Throws an exception if the 
	given vertex is not contained by the graph. */
	std::vector<VertexElementBase*> _get_out_vertices(const VertexElementBase& vertex) const {
		if(contains(vertex)){
			std::vector<VertexElementBase*> out_vertices;
			std::for_each(_edges.begin(), _edges.end(),
							[&vertex, &out_vertices](const Edge<VertexElementBase>* edge){
								if(edge->incidents_from(vertex)){
									out_vertices.push_back(edge->to());
								}
							});
			return out_vertices;	
		}
		else{
			throw VertexNotFoundException(vertex);
		}
	}

	/* Returns the predecessors of the given vertex. Throws an exception if the 
	given vertex is not contained by the graph. */
	std::vector<VertexElementBase*> _get_in_vertices(const VertexElementBase& vertex) const {
		if(contains(vertex)){
			std::vector<VertexElementBase*> in_vertices;
			std::for_each(_edges.begin(), _edges.end(),
							[&vertex, &in_vertices](const Edge<VertexElementBase>* edge){
								if(edge->incidents_to(vertex)){
									in_vertices.push_back(edge->from());
								}
							});
			return in_vertices;
		}
		else{
			throw VertexNotFoundException(vertex);
		}
	}

	/* Returns the neighbours of the given vertex. Duplicate neighbours are ignored.
	Throws an exception if the given vertex is not contained by the graph. */
	std::vector<VertexElementBase*> _get_neighbours(const VertexElementBase& vertex) const {
		std::vector<VertexElementBase*> out = _get_out_vertices(vertex);
		std::vector<VertexElementBase*> in = _get_in_vertices(vertex);
		out.reserve(out.size() + in.size());
		std::for_each(in.begin(), in.end(),
						[&out](const VertexElementBase* in_vertex){
							if(std::find_if(out.begin(), out.end(),
								[in_vertex](const VertexElementBase* out_vertex) {
									return *in_vertex == *out_vertex;
								}) == out.end()) {
									out.push_back(in_vertex);
							}
						});
		out.shrink_to_fit();
		return out;
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
		std::for_each(_vertices.begin(), _vertices.end(), 
						[](const VertexElementBase* vertex){
							delete vertex;
						});
		_vertices.clear();
	}

public:
	/* Graph interface */
	Graph() : _vertices(), _edges() {}

	std::size_t num_vertices() const { return _vertices.size(); }
	std::size_t num_edges() const { return _edges.size(); }

	std::vector<VertexElementBase*>& get_vertices() const { 
		return _all_vertices();
	}

	std::vector<Edge<VertexElementBase>*>& get_edges() const {
		return _all_edges();
	}

	bool contains(const VertexElementBase& vertex) const {
		return _find_vertex(vertex) != _vertices.end();
	}

	bool contains(const Edge<VertexElementBase>& edge) const {
		return _find_edge(edge) != _edges.end();
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
			throw VertexNotFoundException(vertex);
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