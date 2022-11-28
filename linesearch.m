function [lambda] = linesearch(phi,t,y,x,d,numberofiterations,searchmethod)
objectivefunction = @(lambda) 0;
lambda=0;
numofpoints=length(t);
for n=1:numofpoints
           % objectivefunction = @(lambda) objectivefunction(lambda) +  (phi([x(1)+lambda*d(1),x(2)+lambda*d(2)],t(n))-y(n))^2;
         
   objectivefunction = @(lambda) objectivefunction(lambda) +  (phi(x+lambda*d,t(n))-y(n))^2;
end
if(searchmethod=="golden")
   lambda=Golden(objectivefunction,0,100,numberofiterations);
end
if(searchmethod=="newton")
   lambda=Newton(objectivefunction,0.5,0.000001,0.000001);
end

