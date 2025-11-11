function sdot=slab_diff2(t,v,D,h)
% ODE file for diffusion out of a slab.
% D - diffusivity h - half thickness
% Method lines for the PDE dc/dt = D d^x/dx^2 
% over 0 <= x <= h and BC. dc/dx(0,t)=0 at x=0
% and c(h,t)=0. 
% v is a vector of length N+1 
% v(i) is c((i-1)h/N,t) i=1,..,N
% v(N+1) is Q(t)
% vdot(i) = D(N/h)^2 [v(i-1) + v(i+1) - 2 v(i)] i=2,,.N
% vdot(1) = D(N/h)^2 [2 v(2) - 2 v(1)]

N = length(v);
Dh=D/(h^2);
s = v(1:N); s(N+1) =0; 
sdot = zeros(N,1);
sdot(1) = N^2*Dh*(2*s(2)-2*s(1));
for i=2:N
    sdot(i)=N^2*Dh*(s(i+1)+s(i-1)-2*s(i));
end
end

