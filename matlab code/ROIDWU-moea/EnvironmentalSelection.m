function [Population,FrontNo,DWeight] = EnvironmentalSelection(Population,N,axis,theta)
% The environmental selection of DWU-moea

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    angle = zeros(size(Population,2),1);
    %% Non-dominated sorting
    penalty = false(size(Population,2),1);
    [FrontNo,~] = NDSort(Population.objs,Population.cons,N);
    for i =1:size(Population,2)        
        angle(i) = acos(sum(Population(i).obj*axis')/(norm(Population(i).obj)*norm(axis)));
        if angle(i) > theta
            FrontNo(i) = FrontNo(i) + floor(exp(.3*(angle(i) - theta)));
            penalty(i) = true;
        end
    end

    
    %% Calculate the dominance information each solution
    DWeight = InfoDominance(Population.objs);
    
    
    % penalizacao tmb aqui
    %DWeight(penalty) = DWeight(penalty) + 30*floor(exp(.3*(angle - theta)));
    

    %% Environment Select
    Next = ReplacementUniformity(Population,N,FrontNo,DWeight,penalty,angle);
%     temp1 = length(Population);
%     Next = ismember(1:temp1,Next);
    
    %% Select the solutions in the last front
%     Last   = find(FrontNo==MaxFNo);
%     Choose = LastSelection(Population(Next).objs,Population(Last).objs,N-sum(Next),Z,Zmin);
%     Next(Last(Choose)) = true;
  
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    DWeight    = DWeight(Next);
end

function InfoD = InfoDominance(PopObj)
% Calculate the information dominance each solution

    N = size(PopObj,1);

    %% Dominance count each solution
    D = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
            if k == 1
                D(i,j) = true;
            elseif k == -1
                D(j,i) = true;
            end
        end
    end
    CountDominance = sum(D,2);
    %% Calculate information dominance each solution
    InfoD = D'*CountDominance;
end

