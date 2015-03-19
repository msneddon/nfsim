function [data, varNames, consoleOutput, figureHandles] = runbng(path,filename)
%  [data,varNames,consoleOutput,figureHandles] = bng(path,filename)
%
%  Runs BioNetGen on the given BNGL file located at the specified path. If
%  you are already in the correct path, then set path to [].  This method
%  now handles spaces in the path to BNG2.pl
%
%  The first time this function is called, it will ask you to identify the
%  path of the installation directory of BNG/NFsim so that it can find
%  the BNG executables.  It will also add this function to the path
%  so that you can call it from any directory.
%
%  You must also provide the path to the BNGL file so that the script 
%  can find the correct output directory.  This script assumes that 
%  the output GDAT files will be created in the same directory as the 
%  BNGL file.
%  
%  If there is only one GDAT output file created, it will return the data
%  and variable names from that file.  If there are multiple GDAT files,
%  data and varNames will be cell arrays where each cell contains the
%  output of a given GDAT file.  (note: you can use the iscell function to
%  determine if the returned variable is a cell array or not).
%
%  Last updated: 3/05/2009
%  michael.sneddon@yale.edu
%

%First, get the current time, so we can determine what files are generated
%by looking at the time stamps of files in the output directory
ctime=datestr(now);

%init the output variables
data =[]; varNames=[]; consoleOutput=[]; figureHandles=[];

%Try to identify the BNGNF path variable, which will identify
%where the BNG executable is
bngnf_path=getenv('BNGNF_PATH');
if isempty(bngnf_path)
    isValid=0; baseDir=0;
    while isValid~=2
        baseDir = uigetdir('','Identify the NFsim Installation Directory');
        if baseDir==0, return, end;
        if ispc, isValid=exist([baseDir,'\BNG2.pl'],'file');
        else isValid=exist([baseDir,'/BNG2.pl'],'file'); end;
        if isValid==2, break; end;
    
        %if we get here, then we could not find the BNG2.pl executable
        e=errordlg(['Could not find BNG2.pl in this directory! ', ...
             'Please try again.'],'Error finding BNG2.pl','modal');
        uiwait(e);
    end;
    
    %Good, we found it, so we can set the path and get on with the
    %simulations...
    addpath(baseDir);
    setenv('BNGNF_PATH',baseDir);
    bngnf_path=baseDir;
end;



%init the output variables and print a nice little message
fprintf('Running BioNetGen, please wait....');

%run BNG, and be sure to use the correct slash direction for the
%given operating system (annoying, isn't it!)

if(ispc)
    if ~isempty(path), path=[path,'\']; end;
    [status,consoleOutput]=system(['perl "',bngnf_path,'\BNG2.pl" "',path,filename,'"']);
else %if(isunix || ismac)
    if ~isempty(path), path=[path,'/']; end;
    [status,consoleOutput]=system(['perl "',bngnf_path,'/BNG2.pl" "',path,filename,'"']);
end

%Check the status of the run
if(status~=0)
    fprintf('\nBNG Error!\n');
    display(consoleOutput);
    return;
end;

fprintf('\n\n---\n');

%Go into the output directory, and see what files have been updated
if isempty(path), files=dir;
else files=dir(path); end;
for f=1:length(files)
    
    %first get rid of the folders '.' and '..'
    if(strlexcmp(files(f).name,'.')==0 || strlexcmp(files(f).name,'..')==0) 
        continue; 
    end; 
    
    % check if the file is new or updated
    if(strlexcmp(files(f).date,ctime)>=0)
        fprintf(['File: ',files(f).name,' is new or has been modified.\n']);
        
        %Check if it is a .gdat file
        loc=strfind(files(f).name, '.gdat');
        if isempty(loc), continue; end;
        if loc==(length(files(f).name)-4)
           
           %Since it is a .gdat file, plot it!
           fprintf('   -Processing .gdat file and making your plot...\n');
           [d,v]=tblread([path,files(f).name]);
           [figureHandle]=displayNicePlotOfResults(d,v,files(f).name);
           
           %Now check if we have past data or not! (this is so we can
           %return the output of multiple runs)
           if(isempty(data)) 
               %First result, so we can save it directly into the return values
               data=d; varNames=v;
               figureHandles=figureHandle;
           elseif(iscell(data))
               %Add to the cell array
               data{length(data)+1}=d; varNames{length(data)+1}=v;
               figureHandles=[figureHandles;figureHandle;]; %#ok<AGROW>
           else
               %Create a cell array to store all our results
               tempData=data; tempVarNames=varNames;
               data = cell(2,1); varNames=cell(2,1);
               data{1}=tempData;
               varNames{1}=tempVarNames;
               data{2}=d; varNames{2}=v;
               figureHandles=[figureHandles;figureHandle;]; %#ok<AGROW>
           end
           
        end;
    end;
end;


% and we're done....
fprintf('\n');
end



function [figureHandle] = displayNicePlotOfResults(data,varNames,fileName)
% Plots the results of a .gdat file nicely
%
    figureHandle=figure; hold on;
    plot(data(:,1),data(:,2:end));

    %Make sure to turn the interpreter off, so that we don't get unwanted
    %characters from underscores and such
    l=legend(varNames(2:end,:));
    set(l,'Interpreter','none');
    xlabel('Time');
    ylabel('Observable Concentration / Counts');
    title(['Output from: ',fileName],'Interpreter','none');

    %Make things pretty
    set(gcf, 'color', 'white'); box on;
    fontSize = 12;
    set(get(gca,'title'), 'fontSize', fontSize);
    set(get(gca, 'ylabel'), 'fontSize', fontSize);
    set(get(gca, 'xlabel'), 'fontSize', fontSize);

    hold off;
end


function tf = strlexcmp(a, b)
%STRLEXCMP Lexicographic comparison of two strings.
%
%   STRLEXCMP(A, B) returns -1, 0, or 1 depending on whether the left argument
%   is stringwise less than, equal to, or greater than the right argument.
%
%   This is a MATLAB version of the Perl `cmp' operator.
%
%   See also EQ, ISEQUAL.

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-09-22 19:49:47 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

% check arguments
error(nargchk(2, 2, nargin));
if ~ischar(a) || ~ischar(b)
    error('Both arguments must be char arrays (strings).');
end

% get lengths of strings
na = length(a);
nb = length(b);
n = min(na, nb);

% find characters that differ
k = find(a(1:n) ~= b(1:n));
if isempty(k)
    % all characters are identical -- compare lengths
    tf = sign(na - nb);
else
    % compare first character that is different
    k = k(1);
    tf = sign(a(k) - b(k));
end
end



