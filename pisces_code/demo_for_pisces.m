%one easy example for PisCES and PisCES_cv

%simulate adjacency matrices
N = 200;
T = 6;
K = 4;
P = eye(K)*0.3 + ones(K)*0.1;
A = zeros(N,N,T);
z = [repmat([1:K],1, floor(N/K))];
for t = 1:T;
    %z = randsample(z0,N);
    a = double(rand(N) <= P(z,z));
    A(:,:,t) = triu(a,1) + triu(a,1)';
end

%run PisCES with alpha = 0.1
[Z] = PisCES(A,'T', 0.1*ones(T,2)); 

%run PisCES_cv
alphalist = [0.05 0.1 0.15];
[modu like] = pisces_cv(A, alphalist);
