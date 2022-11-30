function j = jacob(r,x)
% j = jacob(r,x)
%
% Calculates the jacobian matrix of the function r at x.

r_vector = r(x);
lx = length(x);
rows = length(r_vector);
j = zeros(rows,lx);
for k = 1:rows
    for i = 1:lx
       xplus = x;
       xminus = x;
       xplus(i) = x(i) + 1.e-8;
       xminus(i) = x(i) - 1.e-8;
       upper = r(xplus);
       lower = r(xminus);
       j(k,i) = ( upper(k) - lower(k) )/2.e-8;
    end
end
