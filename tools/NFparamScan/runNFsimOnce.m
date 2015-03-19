function [ runOutput ] = runNFsimOnce(pathToBNGLFile,bnglFileName,pathToNFsim,outputDirectory,runNumber,paramNames,paramValues)
% RUNNFSIMONCE   Runs NFsim once on the given bngl file modified by the
%                given parameters and parameter values.  This allows you to
%                quickly incorporate NFsim runs into parameter scanning and
%                fitting scripts.  Note that this script will create a new
%                BNGL file named PS_[bnglfilename].bngl in the same
%                directory as the original BNGL file.  This is a temporary
%                file that can be deleted after running this function.
%
%   [ ] = runNFsimOnce(pathToBNGLFile,bnglFileName,pathToNFsim, ...
%               outputDirectory,runNumber,paramNames,paramValues)
%
%   pathToBNGLFile = the path to the BNGL file, including the last slash
%   bnglFileName   = the filename of the BNGL file to run
%   pathToNFsim    = the path to the base NFsim directory, including the last slash
%   outputDirectory = the path to the output directory
%   runNumber      = label of the output gdat file that will be appended to the filename
%   paramNames     = cellarray with the list of parameter names to change
%   paramValues    = for each parameter listed in paramNames, this gives
%                    the value that it should be changed to, this is in a
%                    vector format
%
%   Example usage:
%
%     
%     pathToModel =  '/home/me/Desktop/testModel/';
%     bnglFileName = 'tlbr.bngl';
%     pathToNFsim =  '/home/me/Desktop/NFsim_v1.052';
%     pathToOutput = '/home/me/Desktop/out/';
%
%     runNumber = 'example';  %note, this could also be an actual number
%    
%     paramNames = {'Lig_conc','K1','K2'};
%     paramValues = [0.01,1.2,0.8];
%
%     runNFsimOnce(pathToModel,bnglFileName, pathToNFsim, pathToOutput, ...
%            runNumber, paramNames, paramValues);
%
%
%
%   created by Michael Sneddon, 8/19/2010




% replace all parameters with the new values
if ~isempty(paramNames)
    if ~iscell(paramNames), error('paramNames in argument must be empty or a cell array'); end;
    if length(paramNames)~=length(paramValues), error('paramNames and paramValues must be the same length'); end;
    
    % read in the BNGL file and do the basic things
    baseModel = textread([pathToBNGLFile,bnglFileName],'%s','whitespace','','bufsize',8191*2);
    baseModel = baseModel{:};
    baseModel = strrep(baseModel,'%','%%');
    newOutput = baseModel;
    endlineChar = double(char(java.lang.System.getProperty('line.separator')));
    
    
    % first find where we start and end the parameter block, we should only
    % make changes within these bounds
    paramIndex = findstr(baseModel,'parameters');
    if length(paramIndex)~=2
        error(['could not find the parameters block!  make sure that you have one', ...
            ' and that you do not have a parameter named "parameters"!']);
    end
    
    
    % for each parameter value we want to replace
    for p=1:length(paramNames)
        
        %find the parameter
        %fprintf(['Looking for ',paramNames{p},'\n']);
        strIndex = findstr(baseModel,paramNames{p});
        
        % find the endlines
        
        endlines = findstr(baseModel,char(endlineChar));
       
        % for each occurence, we replace.
        found = 0;
        for k=1:length(strIndex)
            
            % check if the found string is within the parameters block
            if strIndex(k)<paramIndex(1) || strIndex(k)>paramIndex(2),
               continue;
            end
            
            % identify where that line ends by finding the endline
            startLineIndex = endlines(find(endlines<strIndex(k),1,'last'));
            endLineIndex = endlines(find(endlines>strIndex(k),1,'first'));
            
            % make sure it is not a comment, then we don't have to worry
            % about this find
            prefixOfLine = baseModel(startLineIndex:strIndex(k)-1);
            isComment = findstr(prefixOfLine,'#');
            if ~isempty(isComment), continue; end;
            
            % make sure this parameter is the one being defined, and not
            % just found because it is in the expression of another
            % parameter!
            isValid = 1;
            for r=1:length(prefixOfLine), 
                if(~isstrprop(prefixOfLine(r),'wspace')),
                    isValid = 0; break;
                end; 
            end
            if ~isValid, continue; end;
            
           
            
            % generate the replacement string
            replaceString = [paramNames{p},' ',num2str(paramValues(p))];
            
            % make sure the length is good and the new value will fit
            if length(replaceString)> (endLineIndex-startLineIndex-2), 
                error(['The parameter value to replace is longer than the ',...
                    'space available.  rerun and add extra characters to ',...
                    'the end of the line where the parameter ',paramNames{k},' was declared']);
            end;
            while length(replaceString)< (endLineIndex-startLineIndex-2)
                replaceString=[replaceString,' ']; %#ok<AGROW>
            end;
            
            
            % then replace that entire line in the BNGL file
            found = 1;
            newOutput(startLineIndex+1:endLineIndex-2) = replaceString;
            
        end;
        
         % if not found give warning and continue
        if found==0,
            fprintf(['WARNING: parameter named "',paramNames{p},'" was not found!\n']);
        end
        
    end;
    
    
    % save the file as temp, and remember to overwrite the original
    % bnglfile name so that we run the temp file
    outFile = fopen([pathToBNGLFile,'PS_',bnglFileName],'w');
    fprintf(outFile,newOutput);
    fclose(outFile);
    
    bnglFileName = ['PS_',bnglFileName];
end;



% run BioNetGen (which should run NFsim too, as that needs to be in there)
[status,runOutput]=system(['perl ',pathToNFsim,'/BNG/BNG2.pl "',pathToBNGLFile,bnglFileName,'"']);
if status~=0, error(runOutput); end


%Move the output file to the correct directory
strIndex = findstr(bnglFileName,'.');
[status,output]=system(['mv "',pathToBNGLFile,[bnglFileName(1:strIndex-1),'.gdat'],'" "', ...
    outputDirectory,[bnglFileName(1:strIndex-1),'_',num2str(runNumber),'.gdat'],'"']);
if status~=0, 
    error([output, 'make sure no prefix or suffix is added to the run command of the bngl file.']);
end


% return cause we're done here
return;