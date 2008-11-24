function [avgSize, totalMolecules, totalComplexes] = getAvgComplexSize(s)
%
%
%

%init the counters
totalMolecules = zeros(length(s),1);
totalComplexes = zeros(length(s),1);



%Calculates the avg size of each complex
for t=1:length(s)
    allComplexData= [];
    for i=1:length(s(t).data)
        data = s(t).data{i};
        allComplexData = [allComplexData; data(:,2)]; %#ok<AGROW>
    end
    
    totalMolecules(t)=length(allComplexData);
    totalComplexes(t)=length(unique(allComplexData));
end

avgSize = totalMolecules./totalComplexes;