function [fval] = fvalevaluator(phi,t,y,x)
%Evaluates the objective functions value.
numofpoints=length(t);
if(length(t)~=length(y))
   print("The amount of points t and y must be equal!")
   return
end
fval=0;
for n=1:numofpoints
   fval=fval+(phi(x,t(n))-y(n))^2;
end

