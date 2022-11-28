function [x,no_iter] = gaussnewton(phi,t,y,x0,tol,printout,plotout)
%Nu provar jag först med en kod jag vet är sjukt ineffektiv bara för att se
%så övrig logik osv funkisch
%This first section is to determine a search direction
x=x0; %Initial X-values
numofpoints=length(t); %Will be used multiple times, thus saved as variable
numofparameters=length(x);%--||--
Jacobian=zeros(numofpoints,numofparameters); %Creates empty matrix
change=1; %Used for initial stop criteria
numberofiterations=0; %A counter
while change>tol && numberofiterations<275%A rather primitive stop criteria
   lastf=fvalevaluator(phi,t,y,x); %saving functionvalue after last iteration
   for n=1:numofpoints
     r=@(x) phi(x,t(n))-y(n); %this is non-square for every point. Used for the gradient 
   
     gradient=grad(r,x); %calculating the gradient
   
        for m=1:numofparameters
         Jacobian(n,m)=gradient(m); %Saving in the jacobian
  
        end
   end
   %Jacobian is now determined
   rows=zeros(numofpoints,1); %creating empty vector to save each row function is
   for n=1:numofpoints
       rows(n)=phi(x,t(n))-y(n); %Adding all the row functions to a vector
   end
   d=-(Jacobian'*rows)\(Jacobian'*Jacobian)  % search direction, using 3.6
%Search direction funnen
   lambda=linesearch(phi,t,y,x,d,50,"golden") %using golden search. The linesearchfunction essentially takes our objectivefunction and converts it into a function of just lambda
   x=x+lambda*d; %new x
   functionvalue=fvalevaluator(phi,t,y,x) %lets see how big the squared error is
   change=lastf-functionvalue; %used for stop criteria
   numberofiterations=numberofiterations+1; %counting
end
numberofiterations

