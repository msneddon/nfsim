function [results]=readNFdumpAtTime(baseDirectory,systemName,dumpTime)
%  [results]=readNFdumpAtTime(baseDirectory,systemName,dumpTime)
%
%  This function reads the dump output from an NFsim run of a system
%  named 'systemName' from the dump directory 'baseDirectory' at the
%  specified time.  It that time does not exist, this function gives an
%  error. The systemName is the generally set to the name of the bngl file 
%  you simulated, unless you explicilty changed the name. The result is
%  saved as a single struct.  For a description of the elements of the
%  struct, see the help for readNFdump(...).
%
%   Last Updated march, 2010
%   Michael Sneddon (michael.sneddon@yale.edu)
%


tic;
fprintf('\nRunning function readNFdumpAtTime ...\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First, get the list of times that we need to read in
[timeString] = getOutputTimes(baseDirectory);

%some variables to keep track of index values
N_MOLECULES = 1;
N_COMPONENTS = 2;
N_FUNCTIONS = 3;

results = [];
defResStruct = struct('time', [], 'molTypes', [], 'comps',[],'funcs',[],'data',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through the sets of output files
%for t=1:length(timeString)


timeDouble = str2double(timeString);
t=find(timeDouble==dumpTime);
if(prod(size(t))==0) 
    fprintf('Invalid dump time given.  Cannot find that dump.  Exiting\n');
    results=0;
    return;
end;

    %Initialize some arrays to store the MoleculeType info
    moleculeNames = {};
    moleculeTypeData = [];
    
    moleculeCompNames = {};
    moleculeFuncNames = {};

    %set up the filenames
    time=timeString{t};
    
 %   stopBar= progressbar(t./length(timeString),0);
 %   if (stopBar), break; end

    
   % fprintf(['<<<<<<<<<< Output Time: ',num2str(time), ...
   %     's (',num2str(t),' of ',num2str(length(timeString)),') >>>>>>>>>>\n']);
    resStruct=struct(defResStruct);
    resStruct.time=str2double(time);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % First, read in the header file for the set
    %fprintf(['Reading header file: ',systemName,'_nf.',time,'.dump.head\n']);
    headFileName = [baseDirectory,'/',systemName,'_nf.',time,'.dump.head'];

    % open the header file
    [fid, message] = fopen(headFileName,'r');
    if(fid==-1)
        fprintf(['Error when opening the header file named:\n\t', headFileName, '\n\n']);
        fprintf('Matlab says:\n');
        fprintf(['  ',message,'\n\n']);
        return;
     end

    %loop through the file getting the information we need
    readState = 'NONE';
    while ~feof(fid)
    
        %read the line
        line = fgetl(fid);
    
        %skip comment lines, as they are trash
        if strncmp(line,'##',2), continue; end;
    
        %switch to reading a particular type of information
        if strncmp(line,'>>',2),
            if ~isempty(findstr(line,'MoleculeTypeDef'))
                readState = 'TYPEDEF';
            elseif ~isempty(findstr(line,'MoleculeTypeComponents'))
                readState='COMP';
            elseif ~isempty(findstr(line,'MoleculeTypeFunctions'))
                readState='FUNC';
            elseif ~isempty(findstr(line,'EOF'))
                readState='NONE';
            end
            continue;
        end;
        
        %get rid of trailing / leading whitespace
        line = strtrim(line);
    
        %Read in information about the MoleculeTypes so we can parse the
        %other binary output files correctly
        if strcmp(readState,'TYPEDEF') && ~isempty(line)
            [idx,name,nMol,nComp,nFunc] = strread(line,'%d %s %d %d %d','delimiter','\t');
            moleculeNames{length(moleculeNames)+1} = name{1,:}; %#ok<AGROW>
            moleculeTypeData = [moleculeTypeData;[nMol,nComp,nFunc]]; %#ok<AGROW>
        end
    
        %Read in any other information here about component names and functions
        if strcmp(readState,'COMP') && ~isempty(line)
            [molTypeId,idx,name] = strread(line,'%d %d %s','delimiter','\t');
            moleculeCompNames{length(moleculeCompNames)+1} = {molTypeId,idx,name{1,:}}; %#ok<AGROW>
        end
        if strcmp(readState,'FUNC') && ~isempty(line)
            [molTypeId,idx,name] = strread(line,'%d %d %s','delimiter','\t');
            moleculeFuncNames{length(moleculeFuncNames)+1} = {molTypeId,idx,name{1,:}}; %#ok<AGROW>
        end
        
        resStruct.molTypes = moleculeNames;
        resStruct.comps = moleculeCompNames;
        resStruct.funcs = moleculeFuncNames;
    
    end
    fclose(fid);




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through the data files, and read in the information
    allRawData = cell(length(moleculeNames),1);

    for i=1:length(moleculeNames)
       % fprintf(['Reading data file: ',systemName,'_nf.',time,'.dump.',num2str(i-1),'\n']);
        dataFileName = [baseDirectory,'/',systemName,'_nf.',time,'.dump.',num2str(i-1)];

        %Number of columns = 1 for the molecule ID + 1 for the complex ID
        % +2 for each component + 1 for each local function
        columnCount = 1+1+2*moleculeTypeData(i,N_COMPONENTS)+moleculeTypeData(i,N_FUNCTIONS);
        rowCount = moleculeTypeData(i,N_MOLECULES);


        [fid, message] = fopen(dataFileName,'r'); %#ok<NASGU>
        if(fid==-1)
            fprintf('Could not read that file, skipping data for this MoleculeType.');
            allRawData{i}=[];
            %fprintf(['Error when opening the data file named:\n\t', dataFileName, '\n\n']);
            %fprintf('Matlab says:\n');
            %fprintf(['  ',message,'\n\n']);
            %return;
        end

        %If that was successful, then read in all the data with double
        %precision (which is the same way NFsim writes to binary files)
        data = fread(fid,[columnCount,rowCount],'double')';

        allRawData{i}=data;
        fclose(fid);

    end

    %Remember the data, and save it to our struct
    resStruct.data = allRawData;
    results = [results;resStruct]; %#ok<AGROW>
   % fprintf('\n');
   

%close(stopBar);
fprintf('done. '); toc;

end %%%% END OF FUNCTION



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeString] = getOutputTimes(pathToFolder)
%
%
    %first, extract out the actual files
    pathToFolder;
    folder = dir(pathToFolder);

    %set up the timeString to remember the results
    timeString={};

    % goal is to look at the header files, and extract out the
    % time string of each of the outputs.  We will use this to
    % read in all of the data....
    for i=1:length(folder)

        %skip directories
        if folder(i).isdir, continue; end;

        %grab the file name
        %fprintf(folder(i).name);
        %fprintf('\n');
        name = folder(i).name;

        %parse and split the string on the '.' character
        ind=findstr('.',name);
        if ~isempty(ind)
            if strcmp(name(ind(end)+1:end),'head')
                if length(ind)==3
                    timeString{length(timeString)+1} = ...
                        name((ind(1)+1):(ind(2)-1)); %#ok<AGROW>
                elseif length(ind)==4
                    timeString{length(timeString)+1} = ...
                        name((ind(1)+1):(ind(3)-1)); %#ok<AGROW>
                end
            end
        end
    end
    
    
    % Before we finish, we should sort these according to thier number 
    %values, not according to thier string values! so first sort the numbers...
    timeNumbers = str2double(timeString);
    timeNumbers = sort(timeNumbers);
    
    % Using the sorted numerical values of the time, put back in the
    % strings (this is needed in string form so as to extract the correct
    % file names
    timeStringSorted = cell(length(timeNumbers),1);
    for i=1:length(timeString)
        timeStringSorted(timeNumbers==str2double(timeString(i))) = timeString(i);
    end
    timeString = timeStringSorted;
end



