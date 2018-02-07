function [A2] = eigen_complet(A, cvidx, epsilon,k)
   M = zeros(size(A));
while (1)
    A2 = A.*(1-cvidx) + M.*(cvidx);
    [U S V] = svd(A2);
    s = diag(S);
    s(k+1:length(s)) = 0;
    S = diag(s);
    M2 = U*S*V';
    M2 = min(1, max(0,M2));
    %disp(sqrt(sum(sum((M - M2).^2))))
    if sqrt(sum(sum((M - M2).^2))) < epsilon
        break
    end
    M = M2;
end

end

