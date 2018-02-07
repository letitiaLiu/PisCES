%
% function: layout_timeline
% #########################
% This function relabels the Z's to heuristically minimize the number of 
% crossing flows in a sankey plot
%

% David Choi
% 10/01/2017

function [Z] = layout_timeline(A_rec, Z, param)

T = size(Z,2);
n = size(Z,1);

min_flow_size = param.min_flow_size;
frac_min_flow_size = param.frac_min_flow_size;

% initialize the ordering of the first time step
Z(:,1) = init_time_ordering(Z(:,1),A_rec{1});

% use barycenter method to adapt the ordering the clusters in each time
% step, going forwards in time and then backwards

power = 3; % parameter of barycenter method
MAX_ITER = 50;
for j=1:MAX_ITER
    for i=1:T-1
        Z(:,i+1) = adapt_layout(Z(:,i), Z(:,i+1), power, min_flow_size, frac_min_flow_size);
    end
    for i=(T-1):-1:1
        Z(:,i) = adapt_layout(Z(:,i+1), Z(:,i), power, min_flow_size, frac_min_flow_size);
    end
end

end



% function: adapt layout
% #####################
% changes the order of the clusters in z2 so that they are
% more aligned with the clusters in z1

function [z2] = adapt_layout(z1, z2, power, min_flow_size, frac_min_flow_size)

[z1, class_size1, K1] = clean_z(z1);
[z2, class_size2, K2] = clean_z(z2);

z12 = accumarray([z1 z2], 1, [K1 K2]);

% sets each cluster position at time 2 to be the weighted centroid of the
% clusters in time 1 that it flows from 
% weights are proportional to the flow sizes (raised to some power), with small flows thresholded
% out 
[z2, no_change2] = find_barycenter(class_size1, K1, K2, z12, power, z2, ...
    min_flow_size, frac_min_flow_size);

end


% function: clean_z
% ####################
% just relabels the clusters to remove any empty classes that might exist
% and returns auxilliary information (# of nodes in each class, and K)

function [z, class_size, K] = clean_z(z)
class_map = [];
class_map(unique(z)) = 1:length(unique(z));
z = class_map(z)';
class_size = accumarray(z,1);
K = length(class_size);
end


% function: init_time_ordering
% ##############################
% before using the barycenter method to adapt the class orderings, we 
% have to initialize the ordering at time 1 somehow.

function [z1] = init_time_ordering(z1, A1)

% create a distance matrix for hierarchical clustering
% distance is based on laplacianized K x K blockmodel parameter matrix

class_size1 = accumarray(z1,1);
K1 = length(class_size1);

[i1 i2] = ndgrid(z1, z1);
i1 = reshape(i1,[],1);
i2 = reshape(i2,[],1);
edge_counts1 = accumarray([i1 i2], A1(:), [K1 K1]);
P1 = edge_counts1 ./ (class_size1*class_size1');
L1 = P1 - diag(diag(P1));
deg_L = sum(L1,2);
deg_A = sum(edge_counts1,2)./class_size1;
deg_A = deg_A ./ max(deg_A);
L1 = diag(deg_L.^(-1/2))*L1*diag(deg_L.^(-1/2));

distance_method = 'sqrt';
linkage_method = 'complete';

if distance_method == 'sqrt'
    dist_fun = @(i1, i2) 1-sqrt(L1(i1,i2));
    dist_vec = pdist((1:size(L1,1))', dist_fun);
end
clust = linkage(dist_vec, linkage_method);
leafOrder = optimalleaforder(clust, dist_vec);
h = dendrogram(clust,0, 'Reorder', leafOrder);
labels = str2num(get(gca, 'XTickLabel'));
z_order = [];
z_order(labels) = 1:K1;
z1 = z_order(z1)';

end


% function: find_barycenter
% #########################
% this function tries to orders the class labels at time 2 so that they are lined up
% with the class labels at time 1. It does this by assigning each cluster
% at time 2 to a y-coordinate that is the weighted centroid of the classes in
% time 1 that it connects to (after removing small flows), and then
% ordering the cluster labels in order of their y-coordinates.

function [z2, noChange] = find_barycenter(class_size1, K1, K2, z12, power, z2, ...
    min_flow_size, frac_min_flow_size)

vec = cumsum(class_size1);
vec = [1; vec];
midpoint1 = (vec(1:end-1) + vec(2:end))./2;
frac_z12 = z12./(sum(z12,2)*ones(1,K2));
frac_z21 = z12./(ones(K1,1)*sum(z12,1));

frac_flow = min(frac_z12, frac_z21);

P = z12.^power;
P(z12 < min_flow_size & frac_flow < frac_min_flow_size) = 0;
P = P./ (ones(K1,1)*sum(P,1));
midpoint2 = midpoint1'*P;
[val, ind] = sort(midpoint2);
z_order = [];
z_order(ind) = (1:K2);

old_z2 = z2;
z2 = z_order(z2)';
if sum(z2 ~= old_z2)==0
    noChange=1;
else
    noChange=0;
end
end

