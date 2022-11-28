function [lambdastar]=Golden(F,a,b,n)
% For example F = @(x) 1-x*exp(-x)
% [a,b] interval
% tol = max ratio of final to initial interval lengths
% X output matrix containing final a, b and b-a from every iteration
% n = number of function evaluations
alpha = (sqrt(5)-1)/2; %golden ratio
initialL = b-a; %initial interval length
if(initialL<=0) %Makes sure intervall exists
   disp('Intervall has to be bigger than zero')
   return
end
N = 1;
X(1,1)=0;
X(1,2)=a; %lower bound (a) for every iteration
X(1,3)=b; % upper bound (b)
X(1,4)=b-a; %intervall length
for i = 1:n
   lambda = a + (1-alpha)*(b-a);
   my = a + alpha*(b-a);
   if F(lambda) <= F(my)
       b = my; %if F(lambda) is smaller, then new upper limit equals my
   else
       a = lambda; %if F(my) is smaller, then new lower limit equals
                   %lambda
   end
   N = N + 1;
   X(N,1)=N-1;
   X(N,2)=a;
   X(N,3)=b;
   X(N,4)=b-a;
end
%colNames = {'k', 'a', 'b', 'L'};
%sTable = array2table(X, 'VariableNames',colNames)
%N-1
lambdastar=(a+b)/2;
%{
Task 2.1
Function from task 2.1 F = @(x) 1-x*exp(-x)
N = 4 gives tol = alpha^(N-1) = 0.23607
Golden(F,0,2,0.236) gives L0 = 2, L1 = 1.2361, L2 = 0.7639, L3 = 0.4721,
L4 = 0.2918. Final a = 0.9443, final b = 1.2361
%}
%{
Task 2.3
Function from task 2.3 G = @(x) exp(-x) + x^2
N = 5 gives tol = alpha^(N-1) = 0.14590
Golden(G,-1,1,0.146) gives L0 = 2, L1 = 1.2361, L2 = 0.7639, L3 = 0.4721,
L4 = 0.2918, L5 = 0.1803. Final a = 0.2361, final b = 0.4164
%}

