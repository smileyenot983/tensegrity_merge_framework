clc;
close all;
clear all;
cvx_clear;
cvx_solver Gurobi_2;


%  in case you don't have any pregenerated structures
                    
n_nodes = 16;
shift = 0;
sol = run_experiment(3,n_nodes,true,'cube stacks',false,0,false,0,false,0,0);

save("random_structure3","sol");
%                     
% sol = run_experiment(2,10,true,"cylinder",0,0,false,false,...
%                         false,0);
%                     
% save("random_structure2","sol");
% _____________________________

%  load structures
path1 = "forearm_generation.mat";
path2 = "hand_generation.mat";

sol1 = load("forearm_generation","sol");
sol2 = load("hand_generation","sol");


% merge
shift_direction = "x";
shift_distance = 9;

seed = 1;
merge_result = merge(seed,path1,path2,shift_direction,shift_distance);

filename = "forearm_hand_merging";
save(filename,"merge_result");

merge_sol.points = merge_result.points;
merge_sol.R = merge_result.R;
merge_sol.C = merge_result.C;

visualize_solution(merge_sol,3,filename,"separate");






