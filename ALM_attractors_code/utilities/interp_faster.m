function yu = interp_faster(x,y,u)

yu = zeros(length(u),size(y,2));
indl = u<x(1); 
yu(indl,:) = 0;
indp = u>x(end); 
yu(indp,:) = u(indp);
ind = ~indl & ~indp;

u=u(ind);
[~,k] = histc(u,x);
n = length(x);
k(k == n) = n - 1;
t = (u - x(k))./(x(k+1) - x(k));
yu(ind,:) = (1-t).*y(k,:) + t.*y(k+1,:);

