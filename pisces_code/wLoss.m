function [loss] = wLoss(Adjtest, Adjtrain, zhat, type, cvidx)
%Calculate the modularity or loglikelihood

%%Input
%Adjtest: test matrix with dimention N*N; training edges with value 0;
%Adjtrain: training matrix with dimention N*N; test edges with value 0;
%zhat: estimated community assignment with dimension 1*N;
%type: modulatiry for 1 and loglikelihood for 2
%cvidx: N*N marix indicates the index of test edges: 1 for test and 0 for
%training

%%Output
%loss: modulatiry or loglikelihood on test data

k = max(zhat);
loss = 0;
M = mean(mean(Adjtrain));
minvalue = 0;
maxvalue =  max(max(max(Adjtrain)));

n = length(zhat);
W = sum(sum(Adjtest));

     
if type == 1
    Abin = Adjtrain;
    Kts = sum(Abin);
    W = sum(Kts);
    
    [row col] = find(cvidx == 0);
    Num = hist(col, n);
    NE = sum(Num);
    
    for k1 = 1:n
            for k2 = 1:n
                if cvidx(k1,k2) > 0 && zhat(k1) == zhat(k2) && k1 ~= k2
                    loss = loss + (Adjtest(k1,k2) - (Kts(k1)/Num(k1)*Kts(k2)/Num(k2))/W*NE);
                end
            end
    end

else 
    Abin = Adjtrain;
    Kts = sum(Abin);
    W = sum(Kts);
    hatO = zeros(k,k);
    theta = zeros(n,1);
    for i = 1:k
            for j = 1:k
                Ak = Abin(zhat == i, zhat == j);
                Aj = Adjtrain(zhat == i, zhat == j);
                ajx = cvidx(zhat == i, zhat == j);
                hatO(i,j) = sum(sum(Ak));
                
            end
    end
    for k1 = 1:n
        kk = zhat(k1);
        theta(k1) = Kts(k1)/sum(hatO(kk,:));
    end
    for k1 = 1:n
           for k2 = 1:n
                if cvidx(k1, k2) > 0
                %prob = Kts(k1)*Kts(k2)*hatB(zhat(k1),zhat(k2));
                prob = theta(k1)*theta(k2)*hatO(zhat(k1),zhat(k2))/4*5;
                
                %disp(prob)
                if isnan(prob) == 1 | prob == 0
                     prob = 10^(-5);
                end
                if prob >= 1
                     prob = 1 - 10^(-5);
                end
                loss = loss - log(prob)*(Adjtest(k1, k2) > 0.7^6) -  log(1-prob) * (Adjtest(k1, k2) <= 0.7^6);
                end
           end
    end
    loss = loss/W;
end


