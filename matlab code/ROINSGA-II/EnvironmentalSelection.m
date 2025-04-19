function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N,axis,theta)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,~] = NDSort(Population.objs,Population.cons,N);
    
    for i =1:size(Population,2)        
        angle = acos(sum(Population(i).obj*axis')/(norm(Population(i).obj)*norm(axis)));
        if angle > theta
            FrontNo(i) = FrontNo(i) + floor(exp(.3*(angle - theta)));
        end
    end
    MaxFNo = max(FrontNo(~isinf(FrontNo)));
    Next = FrontNo < MaxFNo;
    
    temp = sum(Next);
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    if temp < N
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:N-sum(Next)))) = true;
    else
        Last     = find(FrontNo==(MaxFNo-1));
        [~,Rank] = sort(CrowdDis(Last));
        Next(Last(Rank(1:temp-N))) = false;
    end
        
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);   
    
 
end

