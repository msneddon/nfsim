function [avgValue] = getAvgLocalFunctionValue(s,molTypeName, functionName)
%
%
%

% first, determine which molecule type we are looking for
molTypeIndex = 1;
for i=1:length(s(1).molTypes)
    if strcmp(molTypeName,s(1).molTypes{i})
        molTypeIndex = i;
        %fprintf(['molTypeIndex: ',num2str(i),'\n']);
        break;
    end
end




%init the counters
totalMolecules = zeros(length(s),1);
totalValueSum = zeros(length(s),1);



%Calculates the running sum
for t=1:length(s)
    data = s(t).data{molTypeIndex};
    allFunctionValues = data(:,end); 
    
    totalMolecules(t)=length(allFunctionValues);
    totalValueSum(t)=sum(allFunctionValues);
end

avgValue = totalValueSum./totalMolecules;


