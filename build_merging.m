clc;
close all;
clear all;
cvx_clear;
cvx_solver Gurobi_2;    


% builds merging matrix of shape N x N x O
% where N - number of merged structures
% O - number of possible placemenets for
% merging(left,right,front,back,forward,backward) = 6


% get list of structures for which you want to build merge matrix
all_structures = dir('merged_structs');
% names of structures start from index 3

SEED=1;
N_ORIENTATIONS = 6;

shift_directions = ['x','x','y','y','z','z'];
shift_distances = [6,-6,6,-6,6,-6];

% shift_directions = ['x','x'];
% shift_distances = [6,-6];

% merge_solutions = [];



% n_structs = size(all_structures,1)-2;

% counting total number of nodes to build merging matrix
nodes_total=0;


sizes = [];
shifts = [0];
for i=3:size(all_structures,1)
    path1 = strcat('merged_structs/',all_structures(i).name);
    sol1 = load(path1);
    size1 = size(sol1.sol.C,1);
    nodes_total = nodes_total+size1;
    
%   also saving size of each node for indexing resulting matrix
    sizes(end+1) = size1;
    
    shifts(end+1) = shifts(end) + size1; 
end

points_total = zeros(3,nodes_total);
% getting info about nodes positions
size_prev = 0;
for i=3:size(all_structures,1)
    path1 = strcat('merged_structs/',all_structures(i).name);
    sol1 = load(path1);
    size1 = size(sol1.sol.C,1);
    points_total(:,size_prev*(i-3)+1:size_prev*(i-3)+size1) = sol1.sol.points;
    size_prev = size1;
    
end


C_total = zeros(nodes_total);
R_total = zeros(nodes_total);


%%


num_merges = 0;

last_result = 0;
n_sol = 0;
for k =1:size(shift_directions,2)
    shift_direction = shift_directions(k);
    shift_distance = shift_distances(k);

    
    shift_y=0;
    for i=3:size(all_structures,1)
        shift_x=0;
        for j=3:size(all_structures,1)
            
%            only filling the upper triangular part, lower triangle will be
%            same
            
                
            idx = i-2
            jdx = j-2

            path1 = strcat('merged_structs/',all_structures(i).name);
            path2 = strcat('merged_structs/',all_structures(j).name);

            merge_result = merge(SEED,path1,path2,shift_direction,shift_distance);

%                 extract slice of merged connectivity, responsible for
%                 external connections
            C_merge = merge_result.C(1:sizes(idx),sizes(jdx)+1:2*sizes(jdx));
            R_merge = merge_result.R(1:sizes(idx),sizes(jdx)+1:2*sizes(jdx));

            C_total(shifts(idx)+1:shifts(idx+1),shifts(jdx)+1:shifts(jdx+1),k) = C_merge;
            R_total(shifts(idx)+1:shifts(idx+1),shifts(jdx)+1:shifts(jdx+1),k) = R_merge;

%                 C_total(shift_y+1:shift_y+size(merge_result.C,1),shift_x+1:shift_x+size(merge_result.C,1),k) = merge_result.C;
%                 R_total(shift_y+1:shift_y+size(merge_result.C,1),shift_x+1:shift_x+size(merge_result.C,1),k) = merge_result.R;


%                 shift_x = shift_x + size(C_merge,1);

            num_merges=num_merges+1;
            
            
            disp(strcat("Made merges:",string(num_merges)))
        end
        
%         shift_y = shift_y + size(C_merge,1);
    
    end
end




merge_matrix.C = C_total;
merge_matrix.R = R_total;
merge_matrix.n_nodes = sizes;
merge_matrix.nodes = points_total;

filename = "merged_structs/merge_matrix";
save(filename,"merge_matrix");


% %%
% % Some tests to check if code is correct
% % test 1:
% 
% % checking structures 1 and 2
% 
% path1 = strcat('merged_structs/',all_structures(3).name);
% path2 = strcat('merged_structs/',all_structures(4).name);
% 
% sol1 = load(path1);
% sol2 = load(path2);
% 
% sol1 = sol1.sol;
% sol2 = sol2.sol;
% 
% shift=shift_directions(1);
% shift_dist=shift_distances(1);
% 
% switch shift
%     case 'x'
%         sol2.points(1,:) = sol2.points(1,:) + shift_dist;
%         disp("Shift in x");
%     case 'y'
%         sol2.points(2,:) = sol2.points(2,:) + + shift_dist;
%         disp("Shift in y");
%     case 'z'
%         sol2.points(3,:) = sol2.points(3,:) + + shift_dist;
%         disp("Shift in z");        
%     otherwise
%         disp("No shifting");
% end
% 
% p_1 = sol1.points;
% p_2 = sol2.points;
% 
% p_bar = horzcat(p_1,p_2);
% 
% C_bar = C_total(1:size(p_1,2)+size(p_2,2),1:size(p_1,2)+size(p_2,2));
% R_bar = R_total(1:size(p_1,2)+size(p_2,2),1:size(p_1,2)+size(p_2,2));
% 
% check_sol.points=p_bar;
% check_sol.R = R_bar;
% check_sol.C = C_bar;
% 
% filename="rand_filename";
% visualize_solution(check_sol,3,filename,"separate");
% 
% 
% %% test 1:
% 
% % checking structures 2 and 3
% 
% path1 = strcat('merged_structs/',all_structures(4).name);
% path2 = strcat('merged_structs/',all_structures(5).name);
% 
% sol1 = load(path1);
% sol2 = load(path2);
% 
% sol1 = sol1.sol;
% sol2 = sol2.sol;
% 
% shift=shift_directions(1);
% shift_dist=shift_distances(1);
% 
% switch shift
%     case 'x'
%         sol2.points(1,:) = sol2.points(1,:) + shift_dist;
%         disp("Shift in x");
%     case 'y'
%         sol2.points(2,:) = sol2.points(2,:) + + shift_dist;
%         disp("Shift in y");
%     case 'z'
%         sol2.points(3,:) = sol2.points(3,:) + + shift_dist;
%         disp("Shift in z");        
%     otherwise
%         disp("No shifting");
% end
% 
% p_1 = sol1.points;
% p_2 = sol2.points;
% 
% p_bar = horzcat(p_1,p_2);
% 
% C_bar = C_total(1:size(p_1,2)+size(p_2,2),1:size(p_1,2)+size(p_2,2));
% R_bar = R_total(1:size(p_1,2)+size(p_2,2),1:size(p_1,2)+size(p_2,2));
% 
% check_sol.points=p_bar;
% check_sol.R = R_bar;
% check_sol.C = C_bar;
% 
% filename="rand_filename";
% visualize_solution(check_sol,3,filename,"separate");
% 
% % shift_directions = ['x','x','y','y','z','z'];
% % shift_distances = [6,-6,6,-6,6,-6];
% 
% 
% % merge_sol.points = merge_result.points;
% % merge_sol.R = merge_result.R;
% % merge_sol.C = merge_result.C;
% % visualize_solution(merge_sol,3,filename,"separate");









