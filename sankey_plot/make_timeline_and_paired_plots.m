% function: make_timeline_and_paired_plots
% #######################################
% draws a sankey plot, and/or a series of plots (one for each time step)
% showing the adjacency matrices connected to each flow

% David Choi
% 10/01/2017

function [param] = make_timeline_and_paired_plots(Z, A_rec, flow_rec, clust_rec, param)

param.n = size(Z,1);
param.T = size(Z,2);

clust_rec = find_normalized_cluster_densities(Z, A_rec, clust_rec);

if param.draw_whole_timeline == 1
    param = init_timeline_param(param);
    draw_timeline(flow_rec, clust_rec, param);
    % clean up the axes and stuff like that
    format_figure([10 10]);
end

% for each time step, show adjacency matrices for times t and t+1, and the
% corresponding flow
if param.draw_paired_plots == 1
    param = init_adj_pair_param(param);
    if isempty(param.which_paired_plots)
        plot_list = 1:(param.T-1);
    else
        plot_list = param.which_paired_plots;
    end
    for i=1:length(plot_list)
        ind = plot_list(i);
        figure;
        param = draw_adj_pair(flow_rec{ind}, clust_rec{ind}, clust_rec{ind+1}, Z(:,ind), Z(:, ind+1), A_rec{ind}, A_rec{ind+1}, ...
            param);
        % clean up axes and stuff like that
        format_figure([16 10]);
    end
end
end

% function: draw_adj_pair
% #######################
% draws a plot with the following components
% A1: the sparsity pattern of the adjacency matrix at time 1
% stripe_1: a stripe used to label the clusters of each flow at time 1
% middle: flows between the clusters are drawn here
% stripe_2: a stripe used to label the clusters of each flow at time 2
% A2: the sparsity pattern of the adj matrix at time 2

function [p] = draw_adj_pair(s, c1, c2, z1, z2, A1, A2, p)

% the adjacency matrices are permuted in order of the flows.
% Flows are ordered by (start cluster, end cluster), so the flow
% from class 1 to class 1 is on top and the flow from class K1 to class K2
% is on the bottom

[A1perm, A2perm, A1ind, A2ind] = permute_adj_matrices(s, z1, z2, A1, A2);
p.adj_indices{length(p.adj_indices)+1} = [A1ind' A2ind'];

% flows are drawn thickest first
s = set_drawing_order(s, p.min_flow_size, p.frac_min_flow_size);

% make a large matrix with A1, A2, and space for middle section
bigA = [A1perm zeros(p.n, p.middle_width+2*(p.stripe_width+p.stripe_spacer)) A2perm];

% okay now let's draw the figure:
subplot(1,1,1);
hold off
% draw sparsity pattern of bigA:
spy(bigA);
% draw each flow in order appearing in s:
draw_flows(s, p.middle_start, p.middle_end, p.min_flow_size, p.frac_min_flow_size, 0);

% draw a stripe for each cluster
draw_clusters(p.stripe_1_start, p.stripe_width, c1, 0, p);
draw_clusters(p.stripe_2_start, p.stripe_width, c2, 0, p);

% add text for the class labels in the adj and in the stripe
draw_class_labels( z1, p.A1_start, p.stripe_1_start, ...
    p.stripe_1_start, p.stripe_width, p);
draw_class_labels( z2, p.A2_start, p.stripe_2_start, ...
    p.stripe_2_start, p.stripe_width, p);

end

% function: draw_class_labels
% #########################
% writes the class label over each diagonal block in the sparsity patterns
% and also in the stripe of cluster rectangles

function [] = draw_class_labels(z1, A1_start, middle_start, ...
    stripe_1_start, stripe_width, param)

class_size1 = accumarray(z1,1);
K1 = max(z1);
n = length(z1);

vec = cumsum(class_size1);
vec = [1; vec];
midpoint = (vec(1:end-1) + vec(2:end))./2;
boundary_y = vec;
hold on
for j=1:length(midpoint)
    h = rectangle('Position', ...
        [A1_start+vec(j) vec(j) vec(j+1)-vec(j)*[1 1]]);
    set(h, 'EdgeColor', [1 1 1]*.5);
    set(h, 'Linewidth', .1);
    
    if param.add_class_labels == 1
        h = [];
        h(1) = text(A1_start+midpoint(j), midpoint(j), num2str(j));
        h(2) = text(stripe_1_start+stripe_width/2, midpoint(j), num2str(j));
        set(h, 'color', [.4 0 0]);
        set(h(1), 'HorizontalAlignment', 'center');
        set(h(2), 'HorizontalAlignment', 'center');
        set(h(1), 'FontSize', 10*(1 + max(0, (class_size1(j)./(n/K1)-1)/3)));
        set(h(2), 'FontSize', 10);
    end
end
hold off

end

% function: draw_timeline
% #######################
% this function is a wrapper for drawing the timeline
% it returns some parameters giving the x-coordinates of the different
% communities and flows so that we can draw additional things on top of the
% picture later if we want

function [] = draw_timeline(flow_rec, c_rec, p)

T = length(flow_rec)+1;

figure(1);
delete(gca);
temp = zeros(2);
spy(temp);

for i=1:T-1
    s = set_drawing_order(flow_rec{i}, p.min_flow_size, p.frac_min_flow_size);
    draw_flows(s, p.flow_start(i), p.flow_end(i), p.min_flow_size, p.frac_min_flow_size, p.cluster_spacer_factor);
end

for i=1:T
    draw_clusters(p.cluster_start(i), p.cluster_width, c_rec{i}, p.cluster_spacer_factor, p);
end
axis([-inf inf -inf inf]);
end


% function: draw_clusters
% #####################
% this function draws a rectangle for each cluster in the table c.
% the clusters can be spaced out or touching, depending on
% cluster_spacer_factor (>0 or =0, respectively)

function [] = draw_clusters(x_start, cluster_width, c, cluster_spacer_factor, param)

K = size(c.y_start,1);
n = c.y_start(end) + c.class_size(end);

% add space between the y_start coord if cluster_spacer_factor > 0
c.y_start = c.y_start + (0:K-1)'*n/K*cluster_spacer_factor;

hold on
for i=1:K
    h = rectangle('Position', [x_start c.y_start(i) cluster_width c.class_size(i)]);
    set(h, 'EdgeColor', 0*[1 1 1]);
    if param.show_density_by_greyscale == 1
        set(h, 'FaceColor', (.7 - .4*c.density(i))*[1 1 1]);
    else
        set(h, 'FaceColor', .7*[1 1 1]);
    end
    set(h, 'LineWidth', .1);
end
hold off

end


% function: draw_flows
% ####################
% this function draws a parallelogram for each flow

function [] = draw_flows(flow_table, middle_start, middle_end, ...
    min_flow_size, frac_min_flow_size, cluster_spacer_factor)

% if the flow is too small, we might not draw a black outline, or not even
% draw it at all
min_width_for_boundary = min_flow_size;
min_frac_width_for_boundary = frac_min_flow_size;
min_width_to_draw = min_flow_size;
min_frac_width_to_draw = frac_min_flow_size;

K1 = max(flow_table.class1);
K2 = max(flow_table.class2);
n= max(flow_table.y1_end);

% if space has been inserted between the clusters vertically, then we also
% need to space the flows so that they line up with the clusters
flow_table.y1_start = flow_table.y1_start + (flow_table.class1-1)*n./K1.*cluster_spacer_factor;
flow_table.y1_end = flow_table.y1_end + (flow_table.class1-1)*n./K1.*cluster_spacer_factor;
flow_table.y2_start = flow_table.y2_start + (flow_table.class2-1)*n./K2.*cluster_spacer_factor;
flow_table.y2_end = flow_table.y2_end + (flow_table.class2-1)*n./K2.*cluster_spacer_factor;

for i = 1:length(flow_table.y1_start)
    current_width = flow_table.width(i);
    frac_width = flow_table.frac_width(i);
    current_color = flow_table.color_mat(i,:);
    hold on
    if current_width >= min_width_to_draw | frac_width >= min_frac_width_to_draw
        % draw a parallelogram and fill in with color
        h = fill([middle_start middle_end middle_end middle_start], ...
            [flow_table.y1_start(i) flow_table.y2_start(i) flow_table.y2_end(i) flow_table.y1_end(i)], ...
            current_color);
        set(h, 'EdgeColor', current_color);
        set(h, 'LineWidth', .1);
        if current_width >= min_width_for_boundary  | frac_width >= min_frac_width_for_boundary
            % outline the parallelogram
            set(h, 'EdgeColor', .5*[1 1 1]);
%             h = plot([middle_start*[1 1]' middle_end*[1 1]']', ...
%                 [flow_table.y1_start(i) flow_table.y2_start(i); flow_table.y1_end(i) flow_table.y2_end(i)]', 'k');
%             set(h, 'LineWidth', .1);
%             set(h, 'color', .5*[1 1 1]);
        end
    else       
        % otherwise draw two thin dotted lines, to signify the flow at the
        % start and destination clusters 
        y1 = flow_table.y1_start(i);
        y2 = flow_table.y2_start(i);

        new_x2 = middle_start + (middle_end - middle_start)*.1;
        new_y2 = y1 + .1.*(y2 - y1);
        h = plot([middle_start*ones(size(y1)) new_x2*ones(size(y2))]', [y1 new_y2]', 'g--');
        set(h, 'color', current_color);

        new_x1 = new_x2;
        new_y1 = new_y2;
        new_x2 = middle_start + (middle_end - middle_start)*.2;
        new_y2 = y1 + .2.*(y2 - y1);
        h = plot([new_x1*ones(size(y1)) new_x2*ones(size(y2))]', [new_y1 new_y2]', 'g:');
        set(h, 'color', current_color);

        new_x1 = middle_start + (middle_end - middle_start)*.9;
        new_y1 = y1 + .9.*(y2 - y1);
        h = plot([new_x1*ones(size(y1)) middle_end*ones(size(y2))]', [new_y1 y2]', 'g--');        
        set(h, 'color', current_color);

        new_x2 = new_x1;
        new_y2 = new_y1;
        new_x1 = middle_start + (middle_end - middle_start)*.8;
        new_y1 = y1 + .8.*(y2 - y1);
        h = plot([new_x1*ones(size(y1)) new_x2*ones(size(y2))]', [new_y1 new_y2]', 'g:');        
        set(h, 'color', current_color);        
    end
    hold off
end
end

% function: find_normalized_cluster_densities
% ##########################################
% finds the within class density of each cluster at each time step
% then scales them between 0 and 1, where 1 is densest cluster which
% has more than n/K nodes (so not too small), and 0 is the least dense
% cluster that has more than n/K nodes. small clusters that exceed this
% are thresholded to 0 and 1.

function [clust_rec] = find_normalized_cluster_densities(Z, A_rec, clust_rec)

[n T] = size(Z);
for i=1:T
    class_size = accumarray(Z(:,i),1);
    K = length(class_size);
    [z1 z2] = ndgrid(Z(:,i));
    z1 = reshape(z1,[],1);
    z2 = reshape(z2,[],1);    
    edge_counts = accumarray([z1 z2], A_rec{i}(:), length(class_size)*[1 1]);
    cluster_densities(1:K,i) = diag(edge_counts) ./ class_size.^2;
    large_class_indicator(1:K,i) = class_size./n > 1/K;
end

min_density = 0;
max_density = max(cluster_densities(large_class_indicator == 1));
cluster_densities = cluster_densities - min_density;
cluster_densities = cluster_densities./(max_density - min_density);
cluster_densities = max(min(cluster_densities,1),0);

for i=1:T
    K = max(Z(:,i));
    clust_rec{i}.density = cluster_densities(1:K, i);
end

end


% function: flow_table_2_flow_mat
% ###############################
% sometimes it is convenient to convert the structure flow_table into a
% matrix, so that you can index all the columns at once

function [f_mat, table_names] = flow_table_2_flow_mat(flow_table)
table_names = fieldnames(flow_table);
f_cell = struct2cell(flow_table);
f_cell = f_cell';
f_mat = cell2mat(f_cell);
end

% function: flow_mat_2_flow_table
% ###############################
% after doing something to the flow matrix, we then convert it back to a
% structure

function [flow_table] = flow_mat_2_flow_table(f_mat, table_names)
f_cell = mat2cell(f_mat, size(f_mat,1), [1 1 1 1 1 1 1 1 1 3]);
flow_table = cell2struct(f_cell, table_names, 2);
end


% function: format_figure
% #######################
% this function sets the axes and stuff to make a figure prettier for
% printing or saving to PDF

function [] = format_figure(fig_size)

set(gca, 'XLabel', []);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', fig_size);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 fig_size]);
end


% function: init_timeline_param
% #######################

function [param] = init_timeline_param(param)

if ~isfield(param, 'cluster_spacer_factor')
    param.cluster_spacer_factor = 1;
end
if ~isfield(param, 'cluster_width')
    param.cluster_width = .1;
end
if ~isfield(param, 'cluster_start')
    param.cluster_start = 0;
end
if ~isfield(param, 'show_density_by_greyscale')
    param.show_density_by_greyscale = 1;
end

for i=1:param.T-1
    param.flow_start(i) = param.cluster_start(i) + param.cluster_width;
    param.flow_end(i) = param.flow_start(i) + 1;
    param.cluster_start(i+1) = param.flow_end(i);
end

end

% function: init_adj_pair_param
% ##########################

function [param] = init_adj_pair_param(param)

n = param.n;

if ~isfield(param, 'add_class_labels')
    param.add_class_labels = 1;
end
if ~isfield(param, 'show_density_by_greyscale')
    param.show_density_by_greyscale = 1;
end
if ~isfield(param, 'which_paired_plots')
    param.which_paired_plots = [];
end
if ~isfield(param, 'adj_indices')
    param.adj_indices = {};
end
if ~isfield(param, 'middle_width_mult')
    param.middle_width_mult = .5;
end
if ~isfield(param, 'middle_width')
    param.middle_width = floor(n*param.middle_width_mult);
end
if ~isfield(param, 'stripe_width')
    param.stripe_width = ceil(param.middle_width/10);
end
if ~isfield(param, 'stripe_spacer')
    param.stripe_spacer = ceil(n*.005);
end
if ~isfield(param, 'A1_start')
    param.A1_start = 1;
end
if ~isfield(param, 'A1_end')
    param.A1_end = n;
end
if ~isfield(param, 'stripe_1_start')
    param.stripe_1_start = n + param.stripe_spacer;
end
if ~isfield(param, 'stripe_1_end')
    param.stripe_1_end = param.stripe_1_start + param.stripe_width;
end
if ~isfield(param, 'middle_start')
    param.middle_start = param.stripe_1_end+1;
end
if ~isfield(param, 'middle_end')
    param.middle_end = param.middle_start + param.middle_width;
end
if ~isfield(param, 'stripe_2_start')
    param.stripe_2_start = param.middle_end+1;
end
if ~isfield(param, 'stripe_2_end')
    param.stripe_2_end = param.stripe_2_start + param.stripe_width;
end
if ~isfield(param, 'A2_start')
    param.A2_start = param.stripe_2_end + param.stripe_spacer;
end


end


% function: permute_adj_matrices
% ###############################
% the adjacency matrices are permuted in order of the flows.
% Flows are ordered by (start cluster, end cluster), so the flow
% from class 1 to class 1 is on top and the flow from class K1 to class K2
% is on the bottom

function [A1perm, A2perm, ind1, ind2] = permute_adj_matrices(flow_table, z1, z2, A1, A2)

K1 = max(z1);
K2 = max(z2);

ind1 = [];
for i=1:K1
    for j=1:K2
        ind_temp = find(z1 == i & z2 == j);
        ind1(end+1:(end+length(ind_temp))) = ind_temp;
    end
end
A1perm = A1(ind1,ind1);

ind2 = [];
for i=1:K2
    for j=1:K1
        ind_temp = find(z2 == i & z1 == j);
        ind2(end+1:(end+length(ind_temp))) = ind_temp;
    end
end
A2perm = A2(ind2,ind2);

end


% function: set_drawing_order
% ###########################
% reorders the rows of the flow data table so that they are drawn in a nice
% order. Currently set to draw the thickest flows first and thinner ones on
% top of them.
% Since flow_table is a structure, its easier to convert it to a matrix and
% then reorder the rows, instead of reordering each field in the structure
% separately. (maybe it's not easier, i dunno)

function [flow_table] = set_drawing_order(flow_table, min_flow_size, frac_min_flow_size)


ind0 = find(flow_table.width < min_flow_size & flow_table.frac_width < frac_min_flow_size);
ind = find(flow_table.width >= min_flow_size | flow_table.frac_width >= frac_min_flow_size);

[f_mat, table_names] = flow_table_2_flow_mat(flow_table);

f_clutter = f_mat(ind0,:);
f_keep = f_mat(ind,:);

% what order to draw the flows:
%[val, ind] = sort(abs(s(:,6))); % sort by absolute slope
[val, ind] = sort(f_keep(:,5), 'descend'); % sort by width (thickest first)
%[val, ind] = sort(s(:,5), 'ascend'); % sort by width (thinnest first)
%[val, ind] = sort((s(:,1) + s(:,3))/2); %midpoint in y coordinate

f_keep = f_keep(ind,:);
f_mat = [f_clutter; f_keep];
flow_table = flow_mat_2_flow_table(f_mat, table_names);
end


