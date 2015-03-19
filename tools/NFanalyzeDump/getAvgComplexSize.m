function [avgSize, totalMolecules, totalComplexes] = getAvgComplexSize(s,molTypeName)
%  [avgSize, totalMolecules, totalComplexes] = getAvgComplexSize(s,molTypeName)
%
%  Given the structure S created from reading an NFsim dump file either
%  using the function getTimeArray or getTimeArrayAtTime, together with
%  the name of the molecule type to act on, this function calculates
%  the average size of molecular complexes that includes the moleculeType.
%  The function also returns arrays giving the total number of molecules 
%  at each time, the total as well as the total number of complexes.
%
%
%
%   Last Updated march, 2010
%   Michael Sneddon (michael.sneddon@yale.edu)
%



% first, determine which molecule type we are looking for
molTypeIndex = 1;
for i=1:length(s(1).molTypes)
    if strcmp(molTypeName,s(1).molTypes{i})
        molTypeIndex = i;
        fprintf(['molTypeIndex: ',num2str(i),'\n']);
        break;
    end
end



%init the counters
totalMolecules = zeros(length(s),1);
totalComplexes = zeros(length(s),1);



%Calculates the avg size of each complex
for t=1:length(s)
    data = s(t).data{molTypeIndex};
    allComplexData = data(:,2); 
    
    totalMolecules(t)=length(allComplexData);
    totalComplexes(t)=length(unique(allComplexData));
end

avgSize = totalMolecules./totalComplexes;




% 
% 
% %init the average size vector
% avgSize = zeros(length(s),1);
% 
% %Calculates the avg size of each complex
% for i=1:length(s)
%     data = s(i).data{molTypeIndex};
%     nComplexes = histc(data(:,2),(min(data(:,2))-0.5:1:max(data(:,2))+0.5));
%     avgSize(i)=mean(nComplexes(nComplexes~=0));
% end