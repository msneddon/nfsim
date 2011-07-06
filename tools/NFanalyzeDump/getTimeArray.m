function [time] = getTimeArray(s)
%  [time] = getTimeArray(s)
%
%  Given the structure S created from reading an NFsim dump file either
%  using the function getTimeArray or getTimeArrayAtTime, this function
%  will return a column array of the times at each dump point.  This is
%  useful when you are plotting results from dump files.
%
%   Last Updated march, 2010
%   Michael Sneddon (michael.sneddon@yale.edu)
%

time = zeros(length(s),1);
for i=1:length(s)
    time(i)=s(i).time;
end