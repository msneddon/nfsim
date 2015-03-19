function [names] = getMolTypeNames(s,time)
%  [names] = getMolTypeNames(s,time)
%
%  Simple helper function that allows you to access the names of the
%  MoleculeTypes in the system at the specified time given in an array.
%
%   Last Updated march, 2010
%   Michael Sneddon (michael.sneddon@yale.edu)
%

%if time wasn't given, assume the first time is choosen
molTypes = [];
if ~exist('time')
    molTypes = s(1).molTypes;
else
    
    %Loop until the time we want comes up
    for i=1:length(s)
        if s(i).time > time
            molTypes = s(i).molTypes;
        end
    end
    
    if length(molTypes)==0
        fprintf('No data for sim time that high.\n');
        return
    end
    
end

names = molTypes;
% for n=1:length(molTypes)
%    if n==1, fprintf(['(',num2str(n),'):']);
%    else, fprintf([', (',num2str(n),'):']); end;
%    x = names{n};
%    fprintf(x);
%     
% end