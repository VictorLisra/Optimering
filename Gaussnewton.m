function [x,no_iter] = Gaussnewton(phi,t,y,x0,tol,printout,plotout)

%initiazing variables and functions
tic;

no_iter = 0;
r=@(x) phi(x,t)-y;
f=@(x) r(x)'*r(x);
x = x0;
jacobian = zeros(length(t),length(x)); %pre-allocating to save time

%allocating array sizes (to save time)
allocated_size = 1000;

x_matrix = zeros(allocated_size,length(x));
d_array = zeros(allocated_size);
maxr_array = zeros(allocated_size);
ls_iter_array = zeros(allocated_size);
lambda_array = zeros(allocated_size);
f_array = zeros(allocated_size);
f_norm_array = zeros(allocated_size);

%adding row 1 values to table
for i = 1:length(x)
    x_matrix(i,1) = x(i); %every column is an x
end
d_array(1) = 0;
maxr_array(1) = max(abs(r(x)));
ls_iter_array(1) = 0;
lambda_array(1) = 0;
f_array(1) = f(x);
f_norm_array(1) = norm(grad(f,x));

change = 1; %in order to enter the while loop

while change > tol

    jacobian = jacob(r,x);
    A = jacobian'*jacobian;
    B = jacobian'*r(x);
    d = -A\B;
    gradient = grad(f,x);

    %Levenberg-Marquardt algorithm is used if d contains NaN (A is
    % singular)
    LevenbergMarquardt = false;
    while sum(isnan(d)) ~= 0
        A = A + 0.001*eye(length(A),length(A));
        d = -A\B;
       LevenbergMarquardt = true;
    end

    if LevenbergMarquardt
         disp(['Levenberg-Marquardt has been used.'])
    end


    %Investigating if d is a descent direction
    if transpose(gradient)*d >= 0
        error(['d is not a descent direction. Gauss-Newton terminates.' ...
            ' Try another starting point.'])
    end

    %Line search
    F = @(lambda) f(x+lambda*d);

    %Armijo's rule to estimate upper bound
    upperBound = 1;
    while F(0) >= F(upperBound)
        upperBound = upperBound*2;
    end

    [~,no_line_search,lambda] = Golden(F,0,upperBound,tol);
   
    %Preparing next iteration
    x0 = x;
    x = x + lambda*d;
    no_iter = no_iter + 1;
    change = f(x0)-f(x);

    %adding values to table
    for i = 1:length(x)
        x_matrix(i,no_iter+1) = x(i);
    end    
    d_array(no_iter+1) = norm(d);
    maxr_array(no_iter+1) = max(abs(r(x)));
    ls_iter_array(no_iter+1) = no_line_search;
    lambda_array(no_iter+1) = lambda;
    f_array(no_iter) = f(x);
    f_norm_array(no_iter) = norm(gradient);

end %End of while loop
toc %(Added here since plotting should not count)
if printout == 1

    fprintf(['%4.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s ' ...
        '%12.12s\n'],"iter","x","norm(d)","max(abs(r))","ls iter", ...
        "lambda","f", "norm(f)")    
    
    for n = 1:no_iter
        
        fprintf('%4.0f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
            n-1,x_matrix(1,n),d_array(n), maxr_array(n), ls_iter_array(n), ...
            lambda_array(n),f_array(n), f_norm_array(n));
        
        for k = 2:length(x)
           
            fprintf('%4.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
                '',x_matrix(k,n),'','','','','');
        end
    end
     temp = sum(ls_iter_array);
     total_iterations = temp(1)
end

if plotout == 1
    lin = linspace(min(t),max(t));
    phivalue = phi(x,lin);
    plot(lin,phivalue)
    hold on
    plot(t,y,'*')
    xlim([min(t) max(t)])
    ylim([min(y) max(y)])
end
end