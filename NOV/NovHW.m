close all

for angle = 0:.05:6.25
   [times, trajectory] = get_ray(angle);
end
function [times, trajectory] = get_ray(theta)
  

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
    ye(2,:)
    global inner_contacts
    global outer_contacts
    inner_contacts = [inner_contacts, ye(1,:)];
    outer_contacts = [outer_contacts, ye(1,:)];
    trajectory = zOut;
    % plot results
    hold on
    plot(zOut(:,1), zOut(:,2))
    axis([-1 1 -1 1])


    function dZ = dynamics(t,Z,p) % RHS
        dZ = zeros(4,1);
        epsilon  = .001;
        eta = @(x,y) 1;
        norm_p = sqrt(Z(3).^2 + Z(4).^2);
        x = Z(1)/norm_p;
        y = Z(2)/norm_p;
        dZ(1) = Z(3);
        dZ(2) = Z(4);
        dZ(3) = (eta(x + epsilon, y) - eta (x, y))/epsilon;
        dZ(4) = (eta(x, y + epsilon) - eta (x, y))/epsilon;
    end
end