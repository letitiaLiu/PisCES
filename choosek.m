function [K] = choosek(A, D, Kmax, d, At)

%method of choosing number of modules K

%%Input
%A: smoothed adjacency matrix from iteration of PisCES with dimension N*N
%D: degree matrix with dimention N*N
%Kmax: then maximum number of modules
%At: original adjacency matrix with dimension N*N

%%Output
%K: number of modules

S = size(A);
n = S(1);
if(D == eye(n))
    ad = diag(sum(A)+10^-10) ;
    erw = (eig(D - (sqrt(inv(ad))*A*sqrt(inv(ad)))));
else
erw = (eig(eye(n) - A));
end
erw = sort(erw);
gaprw = erw(2:Kmax+1) - erw(1:Kmax);
[M I] = sort(gaprw(2:Kmax));
pin = sum(sum(sqrt(D)*At*sqrt(D)))/nchoosek(n,2)/2;

threshold = 3.5/pin^(0.58)/n^(1.15);

i = find(gaprw > threshold);

K = max(i);
if length(i) == 0
    K = 1;
end
    


