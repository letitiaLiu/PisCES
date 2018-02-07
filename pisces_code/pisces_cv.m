function [modu like] = pisces_cv(A, alphalist, degree, Kmax, T_n, Kfold)
%This is a function for cross validation of PisCES method

%%Input
%A: adjacency matrices with dimention N*N*T, N is the number of nodes and T
%%%is the number of networks
%degree: 'True' for degree correction
%func: choose from 'PCM' and 'PCQ'. 'PCM' is PisCES while 'PCQ' is kernel
%%%additive method
%alpha: tuning parameter, default as 0.1
%Kmax: maximum number of communities, degault as N/10
%T_n: number of iteration of pisces, default is 10
%Kfole: number of folds in cross validation
%alphalist: possible alpha to choose from

%%Output
%modu: modularity of different alpha
%like: log likelihood of different alpha

if nargin > 6
    error('myfuns:somefun2:TooManyInputs', ...
        'requires at most 6 inputs');
end

S = size(A);
T = S(3);
N = S(2);

% Fill in unset optional values.
switch nargin
    case 1
        alphalist = [0.05 0.1];
        degree = 'T';
        Kmax = floor(N/10);
        T_n = 50;
        Kfold = 5;
        
    case 2
        degree = 'T';
        Kmax = floor(N/10);
        T_n = 50;
        Kfold = 5;
    case 3
        Kmax = floor(N/10);
        T_n = 50;
        Kfold = 5;
        
    case 4
       
        T_n = 50;
        Kfold = 5;
        
    case 5
        
        Kfold = 5;
       
       
end

n_idx = randsample(N*N,N*N);
Atrain = zeros(N,N,T,Kfold);
Atrain2 = Atrain;
for t = 1:T
    
    for k = 1:Kfold
        at = A(:,:,t);
        cvidxtemp = zeros(N,N);
        test = n_idx((k-1)*floor(N*N/Kfold)+1: floor(N*N/Kfold)*k);
        at(test) = 0;
        cvidxtemp(test) = 1;
        cvidxtemp = triu(cvidxtemp) + triu(cvidxtemp)';
        cvidx(:,:,k) = cvidxtemp;
        at = triu(at) + triu(at)';
        Atrain(:,:,t,k) = at;
        Atrain2(:,:,t,k) = eigen_complet(at, cvidxtemp, 10, 10);
    end
end

l = length(alphalist);
modu = zeros(l,1);
like = zeros(l,1);

for a = 1:l
    alpha = ones(T,2)*alphalist(a);
    Z = zeros(T,N,Kfold);
    for k = 1:Kfold
        [Z(:,:,k)] = PisCES(Atrain2(:,:,:,k), degree, alpha, Kmax, T_n);
        for t = 1:T
        modu(a) = modu(a) + wLoss(A(:,:,t), Atrain(:,:,t,k), Z(t,:,k), 1, cvidx(:,:,k));
        like(a) = like(a) + wLoss(A(:,:,t), Atrain(:,:,t,k), Z(t,:,k), 2, cvidx(:,:,k));
        end
    end
    fprintf('modularity for alpha = %d is %d \n',alphalist(a),modu(a))
    fprintf('loglikelihood for alpha = %d is %d \n',alphalist(a),like(a))
end





