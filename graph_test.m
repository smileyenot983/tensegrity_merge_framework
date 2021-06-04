% edge_value - place of connection, defined in build_merging.m  [right,left,front,back,top,bottom]
% node_value - type of connected structure


matrix_path = 'merged_structs/merge_matrix.mat'
mmatrix = load(matrix_path);
mmatrix = mmatrix.merge_matrix;


% creating a graph which will describe merged structure

% adjacency matrix describes connections between nodes
% value of adjacency matrix = place of connection
% A = [0 1 0; 0 0 3; 0 0 0];

A = [0 1; 0 0 ]

% {1,2,3} - correspond to 'random_structure1.mat', 'random_structure2.mat' etc
% G = digraph(A,{'1','2','3'})
G = digraph(A,{'1','2'});

% 

structures = ['random_structure1.mat','random_structure2.mat','random_structure3.mat']

n_nodes = mmatrix.n_nodes;
% get shifts for each structure
% shift depends on sizes of all previous structures
shifts = [0];

for i=1:(size(n_nodes,2))
    prev_size=n_nodes(i);
    shifts(end+1) = shifts(end) + prev_size*2;

end

disp('shifts');
disp(shifts);

% creating matrices for merged structure
size_resulting = shifts(end);


R_resulting = zeros(size_resulting,size_resulting);
C_resulting = zeros(size_resulting,size_resulting);

% shift distance to avoid spawning all structures at the same place
DEFAULT_SHIFT_DISTANCE = 6;




%%

shifts_nodes = [0];
for i=1:size(n_nodes,2)
    shifts_nodes(end+1) = n_nodes(i)*i;
    
    
end

struct1 = str2num(G.Edges.EndNodes{1,1});
struct1_size = n_nodes(struct1);

% nodes_resulting = zeros(3,size_resulting);
nodes_resulting = zeros(3,struct1_size);
nodes_resulting(:,1:struct1_size) = mmatrix.nodes(:,shifts_nodes(struct1)+1:shifts_nodes(struct1)+n_nodes(struct1));

%%
% going through all connections
for i=1:size(G.Edges.EndNodes,1)
    
    
%GETTING INFO ABOUT CONNECTION BETWEEN STRUCTS     
%   take num of each node
    struct1 = str2num(G.Edges.EndNodes{i,1})
    struct2 = str2num(G.Edges.EndNodes{i,2})
    
%     getting their shapes
    struct1_size = n_nodes(struct1);
    struct2_size = n_nodes(struct2);
    
%     placement of connection(defined by edge weight)
%       [right,left,front,back,top,bottom]
    placement = G.Edges.Weight(i);
    
    shift_directions = ['x','x','y','y','z','z'];
    shift_idx = [1,1,2,2,3,3];
    shift_distances = [6,-6,6,-6,6,-6];
    
%     shifting nodes
    current_nodes = mmatrix.nodes(:,shifts_nodes(struct2)+1:shifts_nodes(struct2)+n_nodes(struct2))
    
    idx = shift_idx(placement);
    dist = shift_distances(placement);
    
    current_nodes(idx,:) = current_nodes(idx,:) + dist
        
    nodes_resulting(:,end+1:end+struct2_size) = current_nodes;
    
    
%     now get their connectivity matrix
    shift_x = shifts(struct2); 
    shift_y = shifts(struct1);
    C_ = mmatrix.C(shift_y+1:shift_y+struct1_size+struct2_size,shift_x+1:shift_x+struct1_size+struct2_size,placement);
    R_ = mmatrix.R(shift_y+1:shift_y+struct1_size+struct2_size,shift_x+1:shift_x+struct1_size+struct2_size,placement);
    
% FILLING RESULTING MATRIX WITH CONNECTIONS
    
    R_resulting(shifts(struct1)+1:shifts(struct1+1),shifts(struct2)+1:shifts(struct2+1)) = R_;
    C_resulting(shifts(struct1)+1:shifts(struct1+1),shifts(struct2)+1:shifts(struct2+1)) = C_;
    
end

size(G.Edges.EndNodes)

merge_sol.points = nodes_resulting;
merge_sol.R = R_resulting;
merge_sol.C = C_resulting;
filename = "graph merging test"
visualize_solution(merge_sol,3,filename,"separate");





