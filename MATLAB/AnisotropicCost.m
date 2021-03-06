[x, y] = meshgrid(linspace(-2,2,100), linspace(-2,2,100));
coords = [x(:) y(:)];
% h = .77;
% Adj = make_graph(coords, h);
% gplot(Adj , coords);
% %Adj = diag(ones(24,1),1) + diag(ones(24,1),-1)+ diag(ones(20,1), 5) + diag(ones(20,1), -5);
% psize = size(coords);
% values = Value(coords, @Cost, 1, Adj)
% %contourf(linspace(-1,1,20), linspace(-1,1,20), reshape(values, [20,20]) )
Adj = load('G.mat', 'G');
Adj = Adj.G;
values = Value(coords, @App_Cost, 4951, Adj)


function C = Cost(current, targets)
    vdt = targets-current;
    v_0 = [.0 , 05 ];
    C =  sum(current.^2,2)*.4*((vdt(2)-v_0(2)*cos(current(2))).^2 +  (vdt(1) - v_0(1)*cos(current(1))).^2);
end
function C_1  = App_Cost(current, targets)
    %v = [cos(current(1)); sin(current(2))];
    v = [1;0];
    fdt = targets - current;
    norms  = sqrt(sum(fdt.^2,2));
    a = fdt./([norms, norms]);
    
    C_1 = norm(v)-  (fdt*v > 0).*a*v;
  
    
end
% target should be a list of nodes where exit is possible
function V = Value(positions, K, target, adjacency)
h = norm(positions(2,:) - positions(1,:));
num_coords = size(positions);
nodes  = 1:1:num_coords(1);
% Initialize all values to +infinity
Values = 10^9*ones(size(nodes));
% On the target, the value is zero
Values(target) = 0;
counter = 0;
% An array stores which nodes have been reevalyated in a certain sweep
total_counter = 0;
while(counter < 2)
    Values(target) = 0;
    %if(counter ~= 0)
    %    target = randi(num_coords(1));
    %end
    target_neighbors = nonzeros(adjacency(target,:).*nodes);
    % Initialize the value of the neighbors of the exit node
    Values(target_neighbors) = K(positions(target_neighbors,:), ...
    repmat(positions(target,:), size(target_neighbors))-positions(target_neighbors,:));
    To_Reevaluate = transpose(target_neighbors);
    [~,I] = sort(Values(To_Reevaluate));
    To_Reevaluate = unique(To_Reevaluate(I));

    reevaluated = zeros(num_coords(1), 1);
    reevaluated(target) =  1;
    while(size(To_Reevaluate) ~= 0)
        To_Reevaluate = unique(To_Reevaluate);
        current_node = To_Reevaluate(1);
        reevaluated(current_node) = 1;
        To_Reevaluate(1) = [];
        current_min = Values(current_node);
        relevant_nodes = adjacency(current_node,:).*nodes;
        neighbors = transpose(nonzeros(relevant_nodes));
        % If this nodes neighbors haven't been reevaluated, reevaluate them
        simplices = [];
        for neighbor = neighbors
            if norm(positions(neighbor,:) - positions(current_node, :))>1.8*h 
                break
            end
            if ~reevaluated(neighbor)
                To_Reevaluate = [To_Reevaluate, neighbor];
            end
            % Get Simplices for the `current_node`
            neighbor_relevant_nodes = adjacency(neighbor,:).*nodes;
            neighbor_neighbors = transpose(nonzeros(neighbor_relevant_nodes));
            for node = neighbor_neighbors
                if ismember(node, neighbors) && node < neighbor
                    simplices = [simplices; [neighbor, node]];
                end
            end
            
            current_min = 10^9;
            % The direction is given by the [simplex, xi_value]
            % [0,0] would mean don't go anywhere
            current_direction = [0, 0];
            num_simps = size(simplices);
            for i  = 1:num_simps(1)
                simplex = simplices(i,:);
                [try_direction, try_val]= SemiLagrange(positions, current_node, simplex(1), simplex(2), K, Values, 0);
                if try_val < current_min
                    current_min = try_val(1);
                    current_direction = [i, try_direction(1)];
                end
            end
            
        end
        Values(current_node) = current_min;
        total_counter  = total_counter+1;
        
    end
    counter = counter + 1;
end

%     This loop is probably not useful anymore, need to sweep instead of
%     march
%     while size(Tentative) ~= 0
%          current_node = Tentative(1)
%          Tentative(1) = [];
%          % Find which nodes have a connection to current_node
%          relevant_nodes = adjacency(current_node,:).*nodes;
%          neighbors = transpose(nonzeros(relevant_nodes));
%          % Only consider those neighbors which can benefit from this node
%          pick = nonzeros((Values(neighbors) > Values(current_node)).*neighbors);
%          use_neighbors = transpose(relevant_nodes(pick));
%          size(use_neighbors)
%          size(positions(use_neighbors,:))
%
%          Values(use_neighbors) = max(Values(current_node), K(positions(use_neighbors,:), ...
%               repmat(positions(current_node,:), size(use_neighbors))-positions(use_neighbors,:)));
%          Tentative = [Tentative, transpose(use_neighbors)];
%
%          Tentative = unique(Tentative);
%          [~,I] = sort(Values(Tentative));
%          Tentative = Tentative(I);
%          counter = counter + 1;
%     end
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
function [direction, min_cost] = SemiLagrange(coords, tentative, accepted, settled, K, values, minimization_routine)
% make a line connecting accepted and settled
x1 = coords(accepted, :);
u1 = values(accepted);
x2 = coords(settled, :);
u2 = values(accepted);
x0 = coords(tentative, :);
xi=0;
% if an analytical minimization routine is provided use it, otherwise
% interpolate
if minimization_routine
    [direction, min_cost] = minimization_routine(u1, u2, x0, x1, x2);
else
    interp_points = transpose([linspace(x1(1), x2(1), 20);linspace(x1(2), x2(2), 20)]);
    costs  = transpose(K(x0, interp_points));
    interp_values  = u1*linspace(0,1,20) + u2*(1-linspace(0,1,20));
    possible_values  = max(costs, interp_values);
    [test_val, test_place] = min(possible_values);
    direction = test_place/20;
    min_cost = test_val;
end
% compute minimum allowed cost on that edge,
% in the isotropic case this is just find the closest you
% can get to the smaller node.

end

