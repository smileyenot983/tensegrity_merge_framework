% TODO
% 
% 1. FIX DISTANCES(SHIFTING DEPENDING ON CENTER OF LEFT STRUCTURE)
% 2. NODE NAMES (FOR NOW NODE NAME = TYPE OF STRUCTURE)
% 3. CHECK CONNECTIVITY MATRIX(IN CASE OF 3 STRUCTURES, 3rd STRUCTURE SHOULD BE CONNECTED TO 2nd,
%                               BUT FOR SOME REASON IT IS CONNECTED TO 1st)




clc;
% close all;
clear all;
cvx_clear;
cvx_solver Gurobi_2;  


% TESTING ON CONNECTION OF 2 STRUCTURES

matrix_path = 'merged_structs/merge_matrix.mat'
mmatrix = load(matrix_path);
mmatrix = mmatrix.merge_matrix;


%% CREATING DIRECTED GRAPH BETWEEN NODES '1' AND '2' WITH WEIGHT 2
G = digraph
G = addnode(G,'1')
G = addnode(G,'2')
G = addnode(G,'3')


G = addedge(G,{'1'},{'2'},[1])
G = addedge(G,{'2'},{'3'},[1])

plot(G)
%%
n_nodes = mmatrix.n_nodes;

shifts = [0];

for i=1:(size(n_nodes,2))
    prev_size=n_nodes(i);
    shifts(end+1) = shifts(end) + prev_size;

end

size_resulting = shifts(end);

R_resulting = zeros(size_resulting,size_resulting);
C_resulting = zeros(size_resulting,size_resulting);

% shift distance to avoid spawning all structures at the same place
DEFAULT_SHIFT_DISTANCE = 6;


shifts_nodes = [0];
for i=1:size(n_nodes,2)
    shifts_nodes(end+1) = n_nodes(i)*i;

end

struct1 = str2num(G.Edges.EndNodes{1,1});
struct1_size = n_nodes(struct1);

% saving centroids of each new structure, to use it for shifting when
% merging structures
         
centroids = zeros(3,size(G.Edges.EndNodes,1))


% nodes_resulting = zeros(3,size_resulting);
nodes_resulting = zeros(3,struct1_size);
nodes_resulting(:,1:struct1_size) = mmatrix.nodes(:,shifts_nodes(struct1)+1:shifts_nodes(struct1)+n_nodes(struct1));

% nds = mmatrix.nodes(:,shifts_nodes(struct1)+1:shifts_nodes(struct1)+n_nodes(struct1))
% centroid = sum(nds,2)/(size(nds,2))

% centroids(end+1) = sum(mmatrix.nodes(:,shifts_nodes(struct1)+1:shifts_nodes(struct1)+n_nodes(struct1)));
%%
for i=1:size(G.Edges.EndNodes,1)
    
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
    shift_distances = [10,-10,10,-10,10,-10];
    
%     shifting nodes
    current_nodes = mmatrix.nodes(:,shifts_nodes(struct2)+1:shifts_nodes(struct2)+n_nodes(struct2))
    
    
    
    idx = shift_idx(placement);
    dist = shift_distances(placement);
    
    current_nodes(idx,:) = current_nodes(idx,:) + dist + centroids(idx,i)
        
%     s1 = size(centroids)
%     s2 = size(sum(current_nodes,2)/(size(current_nodes,2)))
    centroids(:,i+1) = sum(current_nodes,2)/(size(current_nodes,2))
    
    nodes_resulting(:,end+1:end+struct2_size) = current_nodes;
    
    
%     now get their connectivity matrix
    shift_x = shifts(struct2); 
    shift_y = shifts(struct1);
%     mmatrix.C()
    C_ = mmatrix.C(shift_y+1:shift_y+struct1_size,shift_x+1:shift_x+struct2_size,placement);
    R_ = mmatrix.R(shift_y+1:shift_y+struct1_size,shift_x+1:shift_x+struct2_size,placement);
    
% FILLING RESULTING MATRIX WITH CONNECTIONS
    
    R_resulting(shifts(struct1)+1:shifts(struct1+1),shifts(struct2)+1:shifts(struct2+1)) = R_;
    C_resulting(shifts(struct1)+1:shifts(struct1+1),shifts(struct2)+1:shifts(struct2+1)) = C_;
    
end


% ADDING INTERNAL CONNECTIONS
for i=1:size(G.Nodes.Name,1)

    path = strcat('merged_structs/','random_structure',num2str(i),'.mat');
    
    sol = load(path);
    C_ = sol.sol.C;
    R_ = sol.sol.R;
    
    R_resulting(shifts(i)+1:shifts(i)+n_nodes(i),shifts(i)+1:shifts(i)+n_nodes(i)) = R_;
    C_resulting(shifts(i)+1:shifts(i)+n_nodes(i),shifts(i)+1:shifts(i)+n_nodes(i)) = C_;
    
end


size(G.Edges.EndNodes)

merge_sol.points = nodes_resulting;
merge_sol.R = R_resulting;
merge_sol.C = C_resulting;
filename = "graph merging test"
visualize_solution(merge_sol,3,filename,"separate");
















