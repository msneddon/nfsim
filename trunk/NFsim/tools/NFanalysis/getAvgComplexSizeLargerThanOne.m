function [avgSize, totalMolecules, totalComplexes, sizeOneCount] = getAvgComplexSize(s,molTypeName)
%
%
%

% first, determine which molecule type we are looking for
molTypeIndex = 1;
for i=1:length(s(1).molTypes)
    if molTypeName==s(1).molTypes{i}
        molTypeIndex = i;
        fprintf(['molTypeIndex: ',num2str(i),'\n']);
        break;
    end
end



%init the counters
totalMolecules = zeros(length(s),1);
totalComplexes = zeros(length(s),1);
sizeOneCount = zeros(length(s),1);


%Calculates the avg size of each complex
for t=1:length(s)
    data = s(t).data{molTypeIndex};
    allComplexData = data(:,2); 
    
    soc = 0; lastValue = -1; counter = 0;
    allComplexData = sort(allComplexData);
    for i=1:length(allComplexData)
        if allComplexData(i)~=lastValue
            counter=counter+1;
        else
            if counter==1, soc=soc+1; end;
            counter = 0;
            lastValue = allComplexData(i);
        end
    end
    
    totalMolecules(t)=length(allComplexData);
    totalComplexes(t)=length(unique(allComplexData));
    sizeOneCount(t) = soc;
end

avgSize = totalMolecules./totalComplexes;



%x = length(unique(list))-sum(circshift(circshift(list,-1)-list,1)-(circshift(list,-1)-list)==-1);
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