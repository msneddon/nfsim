function [avgValue] = getAvgLocalFunctionValue(s,molTypeName,funcName)
%  [avgValue] = getAvgLocalFunctionValue(s,molTypeName, functionName)
%
%  Given the structure S created from reading an NFsim dump file either
%  the function name and moleculeTypeName, this function will return
%  the local function value averaged over each individual molecule for
%  each time point.  This is useful for extracting results when the result
%  is the averaged value of a local function.  The local function value is
%  the value of the function evaluated over the entire complex.  It has two
%  columns, for possibly in the future, outputting another value of the
%  local function evaluated over a different context.
%   
%
%   Last Updated march, 2010
%   Michael Sneddon (michael.sneddon@yale.edu)
%

% first, determine which molecule type we are looking for
molTypeIndex = -1;
for i=1:length(s(1).molTypes)
    if strcmp(molTypeName,s(1).molTypes{i})
        molTypeIndex = i;
        %fprintf(['molTypeIndex: ',num2str(i),'\n']);
        break;
    end
end
if molTypeIndex==-1
    error('Unable to find moleculeType with the given name.');
end;

% second, determine the index of the function, making sure it matches the
% moleculeType
funcIndex = -1;
for i=1:length(s(1).funcs)

    funcInfo = s(1).funcs{i};
    if ((molTypeIndex-1) == funcInfo{1}),
        if strcmp(funcName,funcInfo{3}),
            funcIndex = funcInfo{2};
        end;
    end;
end;


if funcIndex==-1
    error('Unable to find the function with the given name.');
end;

% figure out how many other columns are by taking into account the other
% components of the molecule
nComps = 0;
for i=1:length(s(1).comps)
    if ((molTypeIndex-1) == s(1).comps{i}{1}),
        nComps=nComps+1;
    end;
end;


%figure out what index into the data matrix the function values will be
indexIntoData = 2+nComps*2+(funcIndex*2)+1;


%init the counters
totalMolecules = zeros(length(s),1);
totalValueSum = zeros(length(s),2);


%Calculates the running sum
for t=1:length(s)
    data = s(t).data{molTypeIndex};
    allFunctionValues = data(:,indexIntoData:indexIntoData+1); 
    
    totalMolecules(t)=length(allFunctionValues);
    totalValueSum(t,:)=sum(allFunctionValues);
end

avgValue = [totalValueSum(:,1)./totalMolecules, totalValueSum(:,1)./totalMolecules];


