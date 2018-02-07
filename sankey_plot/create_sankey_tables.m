% function: create_sankey_tables
% ######################
% wrapper function to create the table of x-y coordinates for the flows and
% for the clusters respectively

% David Choi
% 10/01/2017

function [flow_rec, cluster_rec] = create_sankey_tables(Z, A_rec)

[n T] = size(Z);

% create flow table for each (t, t+1)
for i=1:T-1
    z1 = Z(:,i);
    z2 = Z(:,i+1);
    class_size1 = accumarray(z1,1);
    class_size2 = accumarray(z2,1);
    K1 = length(class_size1);
    K2 = length(class_size2);
    z12 = accumarray([z1 z2], 1, [K1 K2]);
    color_mat = create_colormat(z1, z2, Z(:,1));
    flow_rec{i} = create_flow_table(class_size1, class_size2, K1, K2, z12, ...
        color_mat);
end


% create cluster table for each t
for i=1:T
    class_size = accumarray(Z(:,i),1);
    K = length(class_size);
    cluster_rec{i} = create_cluster_table(class_size);
end

end

% function: create_cluster_table
% #############################
% creates the table of x-y values needed to draw rectangles for each
% cluster. Each cluster is a rectangle.
%
% clust.y_start: the y-coordinate of the top side of the cluster (note: later we
% may edit this quantity to insert spaces between the clusters)
% clust.class_size: the number of nodes in the cluster (the height of
% the rectangle)
% clust.density: the within class density suitably normalized (the fill
% color of the rectangle)

function [c] = create_cluster_table(class_size)

K = length(class_size);
y_start = [0; cumsum(class_size(1:end-1))];

%c = [y_start class_size];
c.y_start = y_start;
c.class_size = class_size;

end


% function: create_colormat
% ###########################
% this function assigns to each flow a color
% the color of a flow is the mixture of the colors assigned to its constituent nodes
% (which are assigned by the hsv palette, in order of their classes at time 1
% (classes at time 1 = init_z). 

function [color_mat] = create_colormat(z1, z2, init_z)

initK = max(init_z);
init_color_mat = hsv(initK)';
for i=1:max(z1)
    for j=1:max(z2)
        ind = find(z1 == i & z2 == j);
        weights = accumarray(init_z(ind), 1, [initK 1]);
        weights = weights./sum(weights);
        color_vec = init_color_mat*weights;
        color_array(i,j,:) = color_vec;
    end
end
color_mat1 = reshape(color_array(:,:,1),[],1);
color_mat2 = reshape(color_array(:,:,2),[],1);
color_mat3 = reshape(color_array(:,:,3),[],1);
color_mat = [color_mat1 color_mat2 color_mat3];
color_mat = max(min(color_mat,1),0);
end

% function: create_flow_table
% #############################
% creates the table of x-y values needed to draw the flows at (t,t+1). Each
% flow is a parallelogram going from a cluster at time 1 to a cluster at
% time 2.
%
% flow.y1_start: the y-coord of the top of the cluster at time 1
% flow.y1_end: the y_coord of the bottom of the cluster at time 1
% flow.y2_start: "              " top of the cluster at time 2
% flow.y2_end:   "              " bottomw of the cluster at time 2
% flow.width: the number of nodes in the cluster (or height of
% parallelogram)
% flow.slope: the change in y coordinates/the change in x coord of the
% parallelogram 
% flow.class1: which cluster that the flow belongs to at time 1
% flow.class2: same but at time 2
% flow.frac_width: the number of nodes in the flow, as a fraction of the
% number of nodes in the larger of the two clusters
% flow.color_mat: each row is a 3-dim vector giving the fill color of the flow

function [flow_table] = create_flow_table(class_size1, class_size2, K1, K2, z12, ...
    color_mat)

z1start = [0; cumsum(class_size1(1:end-1))];
z1end = z1start;
for i=1:K2
    z1end(:,i) = z1start(:,i) + z12(:,i);
    z1start(:,i+1) = z1end(:,i);
end

z2start = [0; cumsum(class_size2(1:end-1))]';
z2end = z2start;
for i=1:K1
    z2end(i,:) = z2start(i,:) + z12(i,:);
    z2start(i+1,:) = z2end(i,:);
end

frac_z12 = z12./(sum(z12,2)*ones(1,K2));
frac_z21 = z12./(ones(K1,1)*sum(z12,1));
min_frac_flow = min(frac_z12,frac_z21);

[kk1 kk2] = ndgrid(1:K1, 1:K2);
kk1 = reshape(kk1, [], 1);
kk2 = reshape(kk2, [], 1);
z1start = reshape(z1start(1:K1,1:K2),[],1);
z2start = reshape(z2start(1:K1,1:K2),[],1);
z1end = reshape(z1end(1:K1,1:K2), [], 1);
z2end = reshape(z2end(1:K1,1:K2), [], 1);

min_frac_flow = reshape(min_frac_flow,[],1);
width = z1end - z1start;
slope = z1start - z2start;

keep_ind = find(width > .5);

flow_table.y1_start = z1start(keep_ind);
flow_table.y1_end = z1end(keep_ind);
flow_table.y2_start = z2start(keep_ind);
flow_table.y2_end = z2end(keep_ind);
flow_table.width = width(keep_ind);
flow_table.slope = slope(keep_ind);
flow_table.class1 = kk1(keep_ind);
flow_table.class2 = kk2(keep_ind);
flow_table.frac_width = min_frac_flow(keep_ind);
flow_table.color_mat = color_mat(keep_ind,:);

end


