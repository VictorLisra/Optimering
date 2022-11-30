function [X,N,optimal_x]=Golden(F,a,b,tol)
arguments
    F
    a {mustBeNumeric}
    b {mustBeNumeric}
    tol {mustBeNumeric}
end 
% For example F = @(x) 1-x*exp(-x)
% [a,b] interval
% X output matrix containing final a, b and b-a from every iteration
% N = number of function evaluations
% optimal_x is the optimal x

initialL = b-a; %initial interval length
if(initialL<=0) %Makes sure intervall exists
    error('Intervall has to be bigger than zero.')
end

alpha = (sqrt(5)-1)/2; %golden ratio
N = 1;
X(1,1)=0;
X(1,2)=a; %lower bound (a) for every iteration
X(1,3)=b; % upper bound (b)
X(1,4)=b-a; %intervall length

while (b-a) > tol
    lambda = a + (1-alpha)*(b-a);
    my = a + alpha*(b-a);
    FLambda = F(lambda);
    FMy = F(my);

    if isnan(FMy) || isnan(FLambda)
        error(['Golden: Function values are NaN. ' ...
            'Try another starting x.'])
    end

    

    if FLambda <= FMy
        b = my; %if F(lambda) is smaller or equal, then new upper limit 
                % equals my. This includes the case where both FLambda and
                %FMy are Inf
    elseif FLambda > FMy
        a = lambda; %if F(my) is smaller, then new lower limit equals 
                    %lambda  
    elseif abs(FLambda-FMy) < tol
        optimal_x = a;
        return
    end

    N = N + 1;
    X(N,1)=N-1;
    X(N,2)=a;
    X(N,3)=b;
    X(N,4)=b-a;
end
optimal_x = (a+b)/2;


if isnan(F(optimal_x)) || F(optimal_x)>F(0)
    error('Bad job of the line search!')
end

%colNames = {'k', 'a', 'b', 'L'};
%sTable = array2table(X, 'VariableNames',colNames)
end