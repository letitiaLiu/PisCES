% example of functions to make sankey plots and adjacency matrix plots
% David Choi
% 10/01/2017

% make a small data example with two time steps

n = 300; % the number of nodes in network
K = 6; % the number of classes (can change from time to time if you want)
Z = ceil(rand(n,2)*K); % (the class labels)
for i=1:2
    theta = .05*ones(K) + diag(rand(K,1))*.5; % the blockmodel parameters
    P = theta(Z(:,i),Z(:,i)); % Expected value of data matrx
    A = rand(n) <= P; % the adjacency matrix
    A = triu(A,1); % make it symmetric
    A = A + A';
    A_rec{i} = A; % store in cell array
end

% 1) initialize some parameters
% (the layout algorithm needs to know these in order to find a good layout)
% ################
% to reduce clutter, only draw flows that have at least this many nodes
param.min_flow_size = n/K^2;
% or comprise at least this fraction of its start cluster and destiation cluster
param.frac_min_flow_size = .15;

% 2) get ready for plots
% ################
% find a good layout (relabeling the clusters) and assemble the tables
% needed to draw the plots
[newZ] = layout_timeline(A_rec, Z, param);
[flow_rec, cluster_rec] = create_sankey_tables(newZ, A_rec);

% 3) make some plots
% ################

% draw a sankey diagram (not useful if only two networks)
param.draw_whole_timeline = 1;
param.draw_paired_plots = 0;
param.show_density_by_greyscale = 1; 
[sankey_param] = make_timeline_and_paired_plots(newZ, A_rec, flow_rec, cluster_rec, param);

% draw the flows for two time steps and their adjacency matrices
param.add_class_labels = 1;
param.draw_whole_timeline = 0;
param.draw_paired_plots = 1;
param.which_paired_plots = 1; % just times 1 & 2, make [1 2] to add a plot for times 2 & 3, etc
[paired_param] = make_timeline_and_paired_plots(newZ, A_rec, flow_rec, cluster_rec, param);

% like before, but only include the nodes that are in the large flows
% note: output arguments are helpful if you want to draw more stuff on top
% of this plot later.
% sm_flow, sm_cluster, sm_param: analogous to previous tables and param
% sm_Z and sm_K are only used if set_common_ordering = 1
% sm_Z identifies which flow each node in the plot belongs to
% sm_K identifies which cluster each node in the plot belongs to
set_common_ordering = 0;
paired_param.show_density_by_greyscale = 0; 
[sm_flow, sm_cluster, sm_Z, sm_K, sm_param] = make_compressed_paired_plots(newZ(:,1:2), ...
    flow_rec(1), A_rec(1:2), paired_param, set_common_ordering);

% like previous, but make the adjacency matrices have a common ordering
% (may be easier to compare them visually)
set_common_ordering = 1;
paired_param.max_clust1 = 5; % if you want to only show first 5 clusters
[sm_flow, sm_cluster, sm_Z, sm_K, sm_param] = make_compressed_paired_plots(newZ(:,1:2), ...
    flow_rec(1), A_rec(1:2), paired_param, set_common_ordering);
