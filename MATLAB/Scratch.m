a  = linspace(0,1,20);
x = 1;
vals = func(x,a)
function f = func(x,a)
    f = x.^2 + a.^2
end