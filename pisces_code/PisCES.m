function [Z obj obj2] = PisCES(A, degree, alpha, Kmax, T_n)
%A: adjacency matrices with dimention N*N*T, N is the number of nodes and T
%%%is the number of networks
%degree: 'True' for degree correction
%func: choose from 'PCM' and 'PCQ'. 'PCM' is PisCES while 'PCQ' is kernel
%%%additive method
%alpha: tuning parameter, default as 0.1
%Kmax: maximum number of communities, degault as N/10
%T_n: number of iteration of pisces, default is 50

if nargin > 5
    error('myfuns:somefun2:TooManyInputs', ...
        'requires at most 6 inputs');
end

S = size(A);
T = S(3);
N = S(2);

% Fill in unset optional values.
switch nargin
    case 1
        degree = 'T';
        
        alpha = 0.1*ones(T,2);
        Kmax = floor(N/10);
        T_n = 50;
    case 2
       
        alpha = 0.1*ones(T,2);
        Kmax = floor(N/10);
        T_n = 50;
    case 3
        
        Kmax = floor(N/10);
        T_n = 50;
    case 4
       
        T_n = 50;
   
end

Z = zeros(T,N);
k = zeros(T,1) + Kmax;
obj = zeros(T_n,1);
obj2 = zeros(T_n,1);
V = zeros(N,Kmax,T);
R = zeros(N,Kmax,T);
D = A;
Aori = A;

%degree correction
if degree == 'T'
    for t = 1:T
        a = A(:,:,t);
        d = diag(sum(a) + 10^-10);
        A(:,:,t) = sqrt(inv(d))*a*sqrt(inv(d));
        D(:,:,t) = d;
    end
else
    for t = 1:T
        D(:,:,t) = eye(N);
    end
end    
 

% initial V and R
for t = 1:T
      a = A(:,:,t);
      k(t) = choosek(a, D(:,:,t), Kmax,0.01,a);
      [V(:,1:k(t),t),e] = eigs(a,k(t));
end

for trial = 1:T_n
    R_old = R;
    V_old = V;
for t = 1:T
    if t == 1        
        X = V_old(:,1:k(t+1),t+1);
        At = A(:,:,t);
        S = A(:,:,t) + alpha(t,2)*X*X';
        k(t) = choosek(S, D(:,:,t), Kmax,0.01,At); 
        [V(:,1:k(t),t),e] = eigs(At,k(t));
        obj(trial) = obj(trial) + sum(abs(eig(V(:,1:k(t),t)'*V_old(:,1:k(t),t))));        
    elseif t == T
        
        X = V_old(:,1:k(t-1),t-1);
       
        At = A(:,:,t);
        S = At + alpha(t,1)*X*X';
        Ss = At*At' + 2*alpha(t,1)*X*X';
        k(t) = choosek(S,D(:,:,t), Kmax,0.01,At);
        [V(:,1:k(t),t),e] = eigs(S,k(t));
        obj(trial) = obj(trial) + sum(abs(eig(V(:,1:k(t),t)'*V_old(:,1:k(t),t))));
    else
       
        X1 = V_old(:,1:k(t-1),t-1);
        X2 = V_old(:,1:k(t+1),t+1);       
        At = A(:,:,t);       
        S = At + alpha(t,1)*X1*X1' + alpha(t,2)*X2*X2';
        Ss = At*At' + alpha(t,1)*X1*X1' + alpha(t,2)*X2*X2';
        k(t) = choosek(S,D(:,:,t), Kmax,0.01,At);       
        [V(:,1:k(t),t),e] = eigs(S,k(t));
        obj(trial) = obj(trial) + sum(abs(eig(V(:,1:k(t),t)'*V_old(:,1:k(t),t))));
    end
end
%disp(obj)

%stop 
if trial >1
if abs(obj(trial) - obj(trial-1)) < 0.001
    break
end
end

end

if (obj(trial) - obj(trial - 1) >= 0.001)
    error('method does not converge for alpha = %d, please try smaller alpha', alpha(1,1));
end

Z(1,:) = kmeans(V(:,1:k(1),1),k(1));
for i = 2:T-1
    Z(i,:) = kmeans(V(:,1:k(i),i),k(i));
end
Z(T,:) = kmeans(V(:,1:k(T),T),k(T));
end   
    

  
