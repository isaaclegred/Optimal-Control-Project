#include "Graph.hpp"

#include <functional>
#include <cmath>
#include <vector>


Graph::Graph(std::vector<std::vector<bool>> Adj, std::function<double(std::vector<double>)> Cost,
             std::vector<double> pos, std::function<bool(std::vector<double>)> permitted_direction):
  N(Adj), K(Cost), x(pos), A(permitted_direction) {
  u = std::vector<double>(x.size(), pow(10,9));
}
std::vector<double> Graph::Solve_minmax(){
  std::vector<double> u_vals(x.size(), pow(10,9));
  return u_vals;
};
std::vector<std::vector<bool>> Graph::get_adjacency(size_t num_nodes, std::string type,
                                                           size_t dim){
  std::vector<std::vector<bool>> adj;
  if (type == "cartesian" or type == "Cartesian"){
    double l = pow(num_nodes, 1/dim);
    if (l - floor(l)  > pow(10,-8)){
      throw("Dimension Mismatch in adjacency generation");
    }
    auto  L = static_cast<size_t>(l);
    std::vector<std::vector<bool>> local_adj;
    for(size_t i = 1; i < num_nodes - 1; i++){
      std::vector<bool> current_node(num_nodes);
      for(size_t d = 0; d < dim; d++){
        current_node[pow(L, d) + i] = 1;
        current_node[-pow(L, d) + i] = 1;
      }
      local_adj.push_back(current_node);
    }
    return local_adj;
  }
}
