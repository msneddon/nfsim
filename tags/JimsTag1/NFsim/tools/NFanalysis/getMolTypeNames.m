function [names] = getMolTypeNames(s,time)
%
%
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