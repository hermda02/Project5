function [ u ] = analytical( x,dt,steps )

N = 100;
u = zeros([1 length(x)]);
t = steps * dt;
v = 0;


for i = 1:length(x)
    for k = 1:N
        A = (2*(-1)^k/(k*pi));
        v = v + (A .* sin(pi .* k .* x(i) .* exp(-(k.*pi)^2.*t)));
        v = v + x(i);
    end
    u(i) = v;
end


end
