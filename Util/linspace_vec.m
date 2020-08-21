function y = linspace_vec(x1, x2, n)
% linspace_vec linspace with vector input
% generates size(x1) * n matrix with each row i
% y(i,:) = linspace(x1(i), x2(i), n);
    dx = (x2-x1)/(n-1); 
    y = repmat(dx,1,n); 
    y(:,1) = x1; 
    y = cumsum(y,2);
end