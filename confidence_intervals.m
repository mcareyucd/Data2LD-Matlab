function[CI] = confidence_intervals(Results_cell)

theta    = Results_cell{2,1}(:,1);
sdtheta  = Results_cell{2,1}(:,2);
df       = Results_cell{3,1}(3);
t_v      = tinv(0.975,size(Results_cell{4,1},1)-df);
CI       = num2cell([theta-t_v*sdtheta, theta, theta+t_v*sdtheta]); 
names    = {'LCI  theta ',' theta  ',' UCI  theta'};
Res      = cell2table(CI,'VariableNames',names');
Res