function [Cov gamma] = shcov(X)

N = size(X,2);
T = size(X,1);
tic;
% Classical covariance estimate
C = cov(X);
mu = mean(X);
% evalutation of shrinkage parameter
z = zeros(N,N,T);
for i=1:N
    for j=1:N
        for k=1:T
            z(i,j,k) = (X(k,i)-mu(i))*(X(k,j)-mu(j));
        end
    end
end

vz = var(z,1,3);
v = trace(C)/N;
den = sum(sum(vz));
num1 = sum(sum((C-diag(diag(C))).^2));
num2 = sum((diag(C)-v).^2);

gamma = (T/(T-1)^2)*( den/(num1+num2));
toc;
% shrinkage
Cov = (1-gamma)*C + gamma*v*eye(N);
