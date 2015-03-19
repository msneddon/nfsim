% RUNPARAMETERSCAN
%
%   This script, when run, will attempt to perform a parameter scan on an
%   NFsim model for the selected parameters.  To run this script, update
%   the scan options and setup below accordingly to point to the desired
%   model and scan the given parameters.  Output of this file is stored in
%   the trajectories that are saved to the given output directory.  For an
%   example for loop that reads all the output trajectories, see the very
%   end of this script.
%
%
%
%   created by Michael Sneddon, 8/19/2010



%  BASIC SETUP.  Change these options below so that they correctly point to
%  the correct directories, models, NFsim installation, and desired output
%  directory.
pathToModel =  'tlbrExample/';
bnglFileName = 'tlbr.bngl';
pathToNFsim =  '../../';
pathToOutput = 'tlbrExample/output/';



%  SCAN OPTIONS.  Add the name of the parameters you would like to scan,
%  togehter with the range and interval of values you would like to scan
%  below.

paramNames = {'Lig_conc',  'K1',  'K2'};

paramValues = { [1:1:10]; ...   % values to test for Lig_conc, range from 1 to 10 in steps of 1
                [0.5]; ...      % values to test for K1, note how it is valid to use only a single value
                [1:10:150]; };  % values to test for K2
            
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT THIS SCRIPT BELOW THIS LINE UNLESS YOU UNDERSTAND MATLAB PROGRAMMING

tic;
fprintf('\n+ + + + + + + + + + + + + + + + + + + + + + + + + +\n');
fprintf(['Running parameter scan on file: ',bnglFileName,'\n']);
totalRunCount = 1; currentRunNumber = 1;
for k=1:length(paramValues), totalRunCount = length(paramValues{k})*totalRunCount; end;
fprintf(['Total number of runs that will be generated: ', num2str(totalRunCount),'\n\n' ]);





currentPosition = 1;
currentIndex = ones(size(paramNames));

% run the first trial
fprintf(['\n------\n+ run (',num2str(currentRunNumber),' of ',num2str(totalRunCount),')\n']);
fprintf(['+ running with parameters:\n     ']);
for k=1:length(paramNames),
    fprintf(['   ',paramNames{k}]);
    fprintf(['=',num2str(paramValues{k}(currentIndex(k)))]);
end;
fprintf('\n');
fprintf(['+ trajectories will be saved to file with suffix: ',num2str(currentRunNumber),'.gdat\n']);
paramArray = size(paramNames);
for k=1:length(paramNames), paramArray(k) = paramValues{k}(currentIndex(k)); end;
console = runNFsimOnce(pathToModel,bnglFileName, pathToNFsim, pathToOutput, ...
    currentRunNumber, paramNames, paramArray);
thisRunTime = toc;
meanRunTime = thisRunTime/currentRunNumber;
fprintf(['+ done.  this run took ',num2str(thisRunTime),'s to complete.\n']);
fprintf(['+ estimated time to completion: ',num2str(meanRunTime.*(totalRunCount-currentRunNumber)./60),' min\n']);



while true,
    
     %reset the value at a particular position if we count too high
    allDone = 0;
    while currentIndex(currentPosition) == length(paramValues{currentPosition})
        currentIndex(currentPosition) = 1;
        currentPosition = currentPosition+1;
        %check if we are completely done (once the last digit roles over)
        if currentPosition>length(paramNames), allDone=1; break;  end;
        continue;
    end;
    if allDone, break; end;
    
     %if we got here, then we can increment the value at the current position
    currentIndex(currentPosition) = currentIndex(currentPosition)+1;
    currentPosition = 1;
   
    
    % run the actual simulation
    fprintf(['\n------\n+ run (',num2str(currentRunNumber+1),' of ',num2str(totalRunCount),')\n']);
    fprintf(['+ running with parameters:\n     ']);
    for k=1:length(paramNames), 
        fprintf(['   ',paramNames{k}]);
        fprintf(['=',num2str(paramValues{k}(currentIndex(k)))]);
    end;
    fprintf('\n');
    fprintf(['+ trajectories will be saved to file with suffix: ',num2str(currentRunNumber+1),'.gdat\n']);
    paramArray = size(paramNames);
    for k=1:length(paramNames), paramArray(k) = paramValues{k}(currentIndex(k)); end;
    console = runNFsimOnce(pathToModel,bnglFileName, pathToNFsim, pathToOutput, ...
        currentRunNumber+1, paramNames, paramArray);
    thisRunTime = toc;
    meanRunTime = thisRunTime/currentRunNumber;
    fprintf(['+ done.  this run took ',num2str(thisRunTime),'s to complete.\n']);
    fprintf(['+ estimated time to completion: ',num2str(meanRunTime.*(totalRunCount-currentRunNumber)./60),' min\n']);
    
    
    
    currentRunNumber = currentRunNumber +1;
    k=k+1;
end;
    
    

totalTime = toc;
fprintf(['\n\nParameter Scan Complete.  \nTotal Run time: ',num2str(totalTime./60),' minutes\n']);
fprintf(['Average Time Per Run: ',num2str(meanRunTime./60),' minutes\n']);
fprintf('+ + + + + + + + + + + + + + + + + + + + + + + + + +\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BELOW IS EXAMPLE CODE FOR READING IN THE RESULTS OF THE PARAMETER SCAN

% In general, you will have to modify this part of the script based on what
% you want to get out of the parameter scan.  This is problem specific, so
% cannot be generalized easily.  

% to run the example below, comment out this "return" line
return;


% this example reads in all the trajectories from the output files, saves
% the last value of the last observable (which is "xlinked" in the TLBR
% example that matches crosslinked receptors) into an array that you can
% then manipulate or plot.

% this is set to the number of runs, which was saved as the variable
% totalRunCount if the above script was used
numberOfRuns = totalRunCount;

% then create the output array of the same size as the number of runs that
% we performed during the parameter scan
xLinkedValues = zeros(numberOfRuns,1);

% this is the prefix in the names of the output files.  Remember to add the
% 'PS' before the filename that was added by this parameter scan script!
prefix = 'PS_tlbr_';

% this is the path to the output file (should be same as above)
pathToOutput = 'tlbrExample/output/';



% here we then loop over the runs
for k=1:numberOfRuns,
    % for each run, read in the file
    [data,headers] = tblread([pathToOutput,prefix,num2str(k),'.gdat']);
    
    % get the last output value from the 4th column, which is the
    % observable 'xlinked' in the TLBR example.  Here you could look at
    % more complicated things as well, such as the ratio between molecules,
    % averaged values over the simulation, etc.  We will save that value to
    % our output array
    xLinkedValues(k) = data(end,3);    
end;


% after this loop runs, you will have an array called xLinkedValues that
% you can plot or otherwise analyze using any of the other matlab tools.








