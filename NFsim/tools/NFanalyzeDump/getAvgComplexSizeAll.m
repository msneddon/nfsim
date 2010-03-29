function [avgSize, totalMolecules, totalComplexes] = getAvgComplexSizeAll(s)
%   [avgSize, totalMolecules, totalComplexes] = getAvgComplexSizeAll(s)
%
%  Given the structure S created from reading an NFsim dump file either
%  using the function getTimeArray or getTimeArrayAtTime, this function 
%  calculates the average size of molecular complexes in the entire system.
%  This function also returns the total number of molecules in the system
%  and the total number of complexes.
%
%
%   Last Updated march, 2010
%   Michael Sneddon (michael.sneddon@yale.edu)
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