% function: make_compressed_paired_plot
% ###################################
% removes nodes that aren't in large flows

% David Choi
% 10/01/2017

function [small_flow, small_cluster, flow_Z, original_K, new_param] = make_compressed_paired_plots(Z, ...
    flow_rec, A_rec, param, common_ordering)


if ~isempty(common_ordering)
    common_ordering = 0;
end

param.n = size(Z,1);
param.T = size(Z,2);

% create the params for the compressed plot
new_param = create_new_param(param, common_ordering);

% remove the nodes from small flows
ind = filter_small_flows(Z(:,1), Z(:,2), flow_rec{1}, param);
% create corresponding small Z and small adj matrices
smallZ = Z(ind,[1 2]);
smallA_rec{1} = A_rec{1}(ind, ind);
smallA_rec{2} = A_rec{2}(ind, ind);
% save the indices for later
filter_ind = ind;

% create the sankey tables for the compressed figure
[small_flow, small_cluster] = create_sankey_tables(smallZ, smallA_rec);
% change the colors in the new sankey table to match the old one
small_flow{1} = match_colors(small_flow{1}, flow_rec{1});
% if (common_ordering == 1) then you are making a diff plot and want all
% the flows to be horizontal
if common_ordering == 1
    % make the flows horizontal,and get rid of the greyscale (bug not
    % feature)
    [small_flow{1}, small_cluster, flow_Z, original_K] = make_common_ordering(small_flow{1}, small_cluster, smallZ);
    % make the plot
    new_param = make_timeline_and_paired_plots(flow_Z, smallA_rec, small_flow, small_cluster, new_param);
    % class labels are a little different for diff plots
    add_custom_class_labels(original_K, small_cluster, new_param);
else
    % otherwise just make the plot
    new_param = make_timeline_and_paired_plots(smallZ, smallA_rec, small_flow, small_cluster, new_param);
    flow_Z = [];
    original_K = [];
end
format_figure([16 10]);
% save the filtered indices so you know which genes were kept
new_param.filter_ind = filter_ind;
end




% function: add_custom_class_labels
% #################################
% the class labels are a little different, since the stripe is grouped by
% flows, but the labels should still be by cluster

function [] = add_custom_class_labels(labels, small_cluster, p)

for i=1:size(labels,1)
        h = text(p.stripe_1_start+p.stripe_width/2, small_cluster{1}.y_start(i)+small_cluster{1}.class_size(i)/2, ...
            num2str(labels(i,1)));
        set(h, 'color', [.4 0 0]);
        set(h, 'HorizontalAlignment', 'center');
        set(h, 'FontSize', 10);
    
        h = text(p.stripe_2_start+p.stripe_width/2, small_cluster{1}.y_start(i)+small_cluster{1}.class_size(i)/2, ...
            num2str(labels(i,2)));
        set(h, 'color', [.4 0 0]);
        set(h, 'HorizontalAlignment', 'center');
        set(h, 'FontSize', 10);
end

end


% function: create_new_param
% ##########################
% the compressed plot has slightly different plotting parameters

function [new_param] = create_new_param(param, common_ordering)


% create a new_param structure, some values are copied from param
new_param.min_flow_size = param.min_flow_size;
new_param.frac_min_flow_size = param.frac_min_flow_size;
new_param.show_density_by_greyscale = param.show_density_by_greyscale;
new_param.add_class_labels = param.add_class_labels;
% others are new
new_param.draw_whole_timeline = 0;
new_param.draw_paired_plots = 1;
new_param.middle_width_mult = .25;
if common_ordering == 1
    new_param.add_class_labels = 0;
end


end

% function: draw_class_labels
% ###########################
% draws the stripe with the cluster labels and adds them to the adjacency
% plot too

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

% function: make_common_ordering
% ##############################
% in a diff plot, the flows should all be horizontal. this function changes
% the sankey tables so that this is the case

function [small_flow_table, small_cluster, flow_Z, original_K] = make_common_ordering(small_flow_table, ...
    small_cluster, Z)

[original_K, ind] = sortrows([small_flow_table.class1 small_flow_table.class2]);
[flow_mat, table_names] = flow_table_2_flow_mat(small_flow_table);
[small_flow_table] = flow_mat_2_flow_table(flow_mat(ind,:), table_names);

for i=1:size(small_flow_table.y1_start,1)
    ind = find(Z(:,1) == small_flow_table.class1(i) & Z(:,2) == small_flow_table.class2(i));
    flow_Z(ind, 1:2) = i;    
end

small_flow_table.y2_start = small_flow_table.y1_start;
small_flow_table.y2_end = small_flow_table.y1_end;
small_flow_table.slope = small_flow_table.slope*0;
small_flow_table.class1 = (1:size(small_flow_table.y1_start,1))';
small_flow_table.class2 = (1:size(small_flow_table.y1_start,1))';
small_flow_table.frac_width = ones(size(small_flow_table.frac_width));

for i=1:2
    small_cluster{i}.y_start = small_flow_table.y1_start;
    small_cluster{i}.class_size = small_flow_table.width;
end

end

% function: match_colors
% #########################
% gets the colors from the old flow table and adds them to the new smaller
% table

function [small_table] = match_colors(small_table, flow_table)

for i=1:length(small_table.y1_start)
    ind = find(flow_table.class1 == small_table.class1(i) & flow_table.class2 == small_table.class2(i));
    small_table.color_mat(i,:) = flow_table.color_mat(ind(1),:);
end

end

% function: filter_small_flows
% ############################
% gets rid of all the nodes belonging to small flows

function [ind] = filter_small_flows(z1, z2, flow_table, p)

if isfield(p, 'max_clust1')
    max_clust1 = p.max_clust1;
else
    max_clust1 = Inf;
end
if isfield(p, 'max_clust2')
    max_clust1 = p.max_clust2;
else
    max_clust2 = Inf;
end

flow_ind = find((flow_table.width >= p.min_flow_size | ...
    flow_table.frac_width >= p.frac_min_flow_size) & ...
    (flow_table.class1 <= max_clust1 & flow_table.class2 <= max_clust2));

class1 = flow_table.class1(flow_ind);
class2 = flow_table.class2(flow_ind);

ind = find(ismember([z1 z2], [class1 class2], 'rows')==1);

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

