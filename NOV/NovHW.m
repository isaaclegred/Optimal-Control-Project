close all

global inner_contacts_up
global inner_contacts_down
global inner_contacts_left
global inner_contacts_right

global outer_contacts_up 
global outer_contacts_down
global outer_contacts_left
global outer_contacts_right

inner_contacts_up = [];
inner_contacts_down = [];
inner_contacts_left = [];
inner_contacts_right = [];

outer_contacts_up = [];
outer_contacts_down = [];
outer_contacts_left = [];
outer_contacts_right = [];
N = 640;
fig =figure()
for i = 1:1:N
    angle = 1/N * 2*i * pi;
   [times, trajectory, hits] = get_ray(angle);
end

fig = figure();
hist(outer_contacts_up, 20)
title("Illuminance of y_h =1")
saveas(fig, "out_up.png")
fig = figure();
hist(outer_contacts_down, 20)
title("Illuminance of y_h = -1")
saveas(fig, "out_down.png")
fig = figure();
hist(inner_contacts_up, 20)
title("Illuminance of y_h = .5")
saveas(fig, "in_up.png")
fig = figure();
hist(inner_contacts_down, 20)
title("Illuminance of y_h = -.5")
saveas(fig, "in_down.png")
function [times, trajectory, hits] = get_ray(theta)
  

    p = []; % parameters

    t = 0:0.01:5; % timespan

    Z0 = [0,0, cos(theta), sin(theta)]; % Initial condition
    
    tol.RelTol=1e-10; % Set ODE45 options
    tol.AbsTol=1e-10;
    % Stop when the trajctory hits the bouyndary of the domain
    function [value, isterminal, direction] = events(t, state)
        value = [ceil(min(.25-state(1)^2, .25-state(2)^2)),ceil(min(1-state(1)^2, 1-state(2)^2))];
        isterminal = [0, 1];
        direction = [0,0];
        
        
    end
 
    options  = odeset('Events', @events );
    f = @(t,Z)  dynamics(t,Z,p); % anonymous function handle
    [times, zOut,te,ye, ie] = ode45(f, t, Z0, options); % call ODE45
    global inner_contacts_up
    global inner_contacts_down
    global inner_contacts_left
    global inner_contacts_right

    global outer_contacts_up 
    global outer_contacts_down
    global outer_contacts_left
    global outer_contacts_right
    if (ye(1,2)>= .5)
        inner_contacts_up = [inner_contacts_up, ye(1,1)];
    elseif (ye(1,2) <= -.5 )
        inner_contacts_down = [inner_contacts_down, ye(1,1)];
    elseif (-.5 >= ye(1,1) )  
        inner_contacts_left = [inner_contacts_left, ye(1,2)];
    elseif (.5 <= ye(1,1))
        inner_contacts_right = [inner_contacts_right, ye(1,2)];
    end 
    if (ye(2,2)>= 1)
        outer_contacts_up = [outer_contacts_up, ye(2,1)];
    elseif (ye(2,2) <= -1 )
        outer_contacts_down = [outer_contacts_down, ye(2,1)];
    elseif (-1 >= ye(2,1)  )  
        outer_contacts_left = [outer_contacts_left, ye(2,2)];
    elseif (1 <= ye(2,1) )
        outer_contacts_right = [outer_contacts_right, ye(2,2)];
    end 
    hits = [];
    trajectory = zOut;
    % plot results
    hold on
    plot(zOut(:,1), zOut(:,2))
    axis([-1 1 -1 1])


    function dZ = dynamics(t,Z,p) % RHS
        dZ = zeros(4,1);
        epsilon  = .001;
        eta = @(x,y) 3/2*exp((-(x-.5)^2 - (y-.5)^2)/2/.1^2)+ 3/2*exp((-(x+.5)^2 - (y-.5)^2)/2/.1^2)...+
            -1/2*exp((-(x)^2 - (y-.5)^2)/2/.2^2) + 1;
        norm_p = sqrt(Z(3).^2 + Z(4).^2);
        x = Z(1)/norm_p;
        y = Z(2)/norm_p;
        dZ(1) = Z(3);
        dZ(2) = Z(4);
        dZ(3) = (eta(x + epsilon, y) - eta (x, y))/epsilon;
        dZ(4) = (eta(x, y + epsilon) - eta (x, y))/epsilon;
    end
     
end