#include <string>
#include <functional>
#include <vector>
// The idea is to have a graph which implements a "solve minmax problem function", hopefully the
// anisotropy can be handeled purely in the adjacency matrix part stage
struct Graph{
  // Create a graph using adjency matrix Adj (which has a 1 in Adj[i][j] if there is an
  // edge from node i to node j),
  Graph(std::vector<std::vector<bool>> adjacency, std::function<double(std::vector<double>)> cost,
        std::vector<double> positions, std::function<bool(std::vector<double>)> permitted_direction );

  std::vector<double> Solve_minmax();
  static std::vector<std::vector<bool>> get_adjacency(size_t num_nodes, std::string type, size_t dim);
  std::vector<std::vector<bool>> N;
  std::function<double(std::vector<double>)> K;
  std::vector<double> x;
  std::function<bool(std::vector<double>)> A;
  std::vector<double> u;
};
