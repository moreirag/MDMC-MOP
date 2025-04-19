function ROIDWUmoea(Global)
% <algorithm> <D>
% Dominance-Weighted Uniformity multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% G. Moreira and L. Paquete, "Guiding under uniformity measure in the 
% decision space," 2019 IEEE Latin American Conference on Computational 
% Intelligence (LA-CCI), Guayaquil, Ecuador, 2019, pp. 1-6.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    %% Definição do cone
    axis = [1,1];
    theta = 0.3;
    %% Generate random population
    %[Z,Global.N] = UniformPoint(Global.N,Global.M);
    Population   = Global.Initialization();
    %Zmin         = min(Population(all(Population.cons<=0,2)).objs,[],1);

    [~,FrontNo,DWeight] = EnvironmentalSelection(Population,Global.N,axis,theta);

    %% Optimization
while Global.NotTermination(Population)
    MatingPool = TournamentSelection(2,Global.N,FrontNo,DWeight);
    Offspring  = GA(Population(MatingPool));
    [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,axis,theta);
end
end