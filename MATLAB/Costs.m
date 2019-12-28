[x, y] = meshgrid(linspace(-2,2,100), linspace(-2,2,100));
coords = [x(:) y(:)];
omega = .098;
Frames  = [];
while(omega < .099)
     contourf(linspace(-1,1,100), linspace(-1,1,100), reshape(Cost(coords, omega), [100,100]) );
     title("omega = " + omega)
     frame = getframe;
     omega = omega + .0005;
     Frames = [Frames, frame];
end
function C = Cost(X, omega)
         C =  1/10*sum(X.^2,2).*cos(omega*sum(X, 2)).^2.*exp(-.5*(X(:,1).^2 + X(:,2).^2));
end