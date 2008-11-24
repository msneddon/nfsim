function [avgSize] = getAvgComplexSize(s,molTypeName)
%
%
%

% first, determine which molecule type we are looking for
molTypeIndex = 1;
for i=1:length(s(1).molTypes)
    if molTypeName==s(1).molTypes{i}
        molTypeIndex = i;
    end
end


%init the average size vector
avgSize = zeros(length(s),1);

%Calculates the avg size of each complex
for i=1:length(s)
    data = s(i).data{molTypeIndex};
    nComplexes = histc(data(:,2),(min(data(:,2))-0.5:1:max(data(:,2))+0.5));
    avgSize(i)=mean(nComplexes(nComplexes~=0));
end