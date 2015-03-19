function [results]=readNFdump(baseDirectory,systemName)
%  [results]=readNFdump(baseDirectory,systemName)
%
%  This function reads all the dump output from an NFsim run of a system
%  named 'systemName' from the dump directory 'baseDirectory'.  The
%  systemName is the generally set to the name of the bngl file you
%  simulated, unless you explicilty changed the name.  When called, it
%  will read all the files in the directory.  Be careful, because
%  this may exhaust your system memory if you have too many output
%  files or too large of a system.  If this is the case, use the function
%  readNFdumpAtTime(...) instead.
%
%  When you call this function, it will save the results in a struct
%  array as:
%
%  struct array with fields:
%    time
%    molTypes
%    comps
%    funcs
%    data
%
%  You can access any of the elements in the array as:
%   results(index).time
%
%  Here, index will range from 1 to the size of the results struct array.
%  You can then access any of the fields by name.  Time gives the time of
%  the simulation at that index.  molTypes will give you a list of the
%  moleculeTypes that were saved.  The moleculeTypes are indexed from 0 to
%  n where n is the number of moleculeTypes.  That moleculeType index is
%  referenced in the comps field.  The comps field is a cell array of all
%  the components in the system.  In each position of the cell array, there
%  is another cell array of size three.  In this cell array, for each
%  component, you are given the moleculeType index in position 1, the
%  component index in position 2, and the component name in position 3.
%  This allows you to map names of components of moleculeTypes to thier
%  index which is used later.  The same syntax applies the functions block,
%  where the moleculeType index, followed by the function index, is joined
%  with the function name.  In this case, only local functions are
%  displayed (global functions can be output to gdat files using the -ogf 
%  flag, and are not generally a local property of individual molecules).
%  Local functions are associated with the molecule that requires thier
%  value to be evaluated.  The value of the local function is given for
%  the case where it is evaluated over the entire complex and over the
%  single molecule whether or not these values are used in rules.  Finally, 
%  all of the data is stored in the data field. The data field is a cell array,
%  containing the data for that time point for each moleculeType, numbered
%  according to the moleculeType index.  In each of these data matricies
%  for each molecule type, you are given the unique molecule Id in the
%  first column, the unique complex id in the second column (whehter or not
%  you used complex bookkeeping).  Then, for each component there are two
%  colums.  The first gives the state of the component, as an indexed
%  number value.  The indexed number value can be mapped to the component
%  state by looking at the bngl specification file.  The states are ordered
%  in the same way they were declared.  Next, for each component, you have
%  the unique id of the molecule that is bound at this site, or -1 if the
%  site is empty.  Finally, for each function that is local to the
%  molecule, you are provided with the value of the function evaluated
%  over the entire complex.
%
%
%   Last Updated march, 2010
%   Michael Sneddon (michael.sneddon@yale.edu)
%

tic;
fprintf('\nRunning function readNFdump...\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First, get the list of times that we need to read in
[timeString] = getOutputTimes(baseDirectory);

%some variables to keep track of index values
N_MOLECULES = 1;
N_COMPONENTS = 2;
N_FUNCTIONS = 3;

results = [];
defResStruct = struct('time', [], 'molTypes', [], 'comps',[],'funcs',[],'data',[]);

ticMarks = length(timeString);
if ticMarks >40, ticMarks=40; end;
fprintf('Progress:\n[');
for i=1:ticMarks, fprintf('-'); end;
fprintf(']\n[');

currentPos = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through the sets of output files
for t=1:length(timeString)

    if( t/length(timeString) > currentPos/ticMarks)
        fprintf('+');
        currentPos=currentPos+1;
    end;
    
    %Initialize some arrays to store the MoleculeType info
    moleculeNames = {};
    moleculeTypeData = [];
    
    moleculeCompNames = {};
    moleculeFuncNames = {};

    %set up the filenames
    time=timeString{t};
    
    
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
        columnCount = 1+1+2*moleculeTypeData(i,N_COMPONENTS)+2*moleculeTypeData(i,N_FUNCTIONS);
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
end;


fprintf('+]\n\n');
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



