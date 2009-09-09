function [time] = getTimeArray(s)
%
%
%

time = zeros(length(s),1);
for i=1:length(s)
    time(i)=s(i).time;
end