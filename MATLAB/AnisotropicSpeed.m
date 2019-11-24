% This will hopefully contain the code which attempts to use Dijkstra's
% method on a non-isotropic grid, and it will compare to an iterative
% method on the same grid.  
[x, y] = meshgrid(linspace(-2,2,100), linspace(-2,2,100));
coords = [x(:) y(:)];
Adj = load('G.mat', 'G');
Adj = Adj.G;
values = Value(coords, @Cost, 4951, Adj)
function C = Cost(X, ~)
     C =  1/10*sum(X.^2,2).*cos(.112*sum(X, 2)).^2.*exp(-.5*(X(:,1).^2 + X(:,2).^2));
end
function V = Value(positions, K, exit, adjacency)
    num_coords = size(positions);
    nodes  = 1:1:num_coords(1);
    Values = 10^9*ones(size(nodes));
    Values(exit) = 0;
    size(adjacency)
    exit_neighbors = nonzeros(adjacency(exit,:).*nodes);
    % Check if there is a valid path from exit_neighbors to exit, 
    % if so, then proceed with the rest of the process
    % Initialize the value of the neighbors of the exit node
    Values(exit_neighbors) = K(positions(exit_neighbors,:), ...
     repmat(positions(exit,:), size(exit_neighbors))-positions(exit_neighbors,:));
    Tentative = transpose(exit_neighbors);
    [~,I] = sort(Values(Tentative));
    Tentative = unique(Tentative(I));
    counter = 0;
    while size(Tentative) ~= 0 
         current_node = Tentative(1)
         Tentative(1) = [];
         % Find which nodes have a connection to current_node
         relevant_nodes = adjacency(current_node,:).*nodes;
         neighbors = transpose(nonzeros(relevant_nodes));
         % Only consider those neighbors which can benefit from this node
         % and have a valid path to it which obeys the dynamics.
         pick = nonzeros((Values(neighbors) > Values(current_node)).*neighbors);
         use_neighbors = transpose(relevant_nodes(pick));
         size(use_neighbors)
         size(positions(use_neighbors,:))

         Values(use_neighbors) = max(Values(current_node), K(positions(use_neighbors,:), ...
              repmat(positions(current_node,:), size(use_neighbors))-positions(use_neighbors,:)));
         Tentative = [Tentative, transpose(use_neighbors)];
         
         Tentative = unique(Tentative);
         [~,I] = sort(Values(Tentative));
         Tentative = Tentative(I);
         counter = counter + 1; 
    end
    counter
    V = Values;
end
function G = make_graph(coords, h)
nodes = 1:size(coords ,1);
G = zeros(size(nodes, 1));
for i = nodes
    for j = 1:i-1
        if norm(coords(i,:)- coords(j,:)) < h
            G(i,j)=1;
            G(j,i)=1;
        end
    end
end
end
% Unless the cost is anisotropic this function cannot be useful given
% vurrent methods
function [direction, mincost] = SemiLagrange(tentative, accepted, settled, K)
    % make a line connecting accepted and settled
    x1 = coords(accepted, :);
    x2 = coords(settled, :);
    x = coords(tentative, :);
    f(x)
    % compute minimum allowed cost on that edge,
    % in the isotropic case this is just find the closest you
    % can get to the smaller node.  
    
end