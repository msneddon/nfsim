function [data, variableNames] = readNFsimBinary(dataFileName)
%  READNFSIMBINARY - Read in binary output of NFsim.  NFsim will provide two
%  output files if you choose to output in binary format.  The actual data
%  file which contains only data results and a header file (always ending
%  in the .head extension).  Give this matlab function the name of the data
%  file and it will automatically search for the header file as:
%  "[dataFileName].head".  The header file also gives the names of each
%  column.  This function will return a matrix of the results along with
%  the variableNames which is a vector of the names of each column.
%
%
%   [data, variableNames] = readNFsimBinary(dataFileName)
%
%   Created 9/10/08
%   Michael Sneddon (michael.sneddon@yale.edu)

%Declare the default output
data=[]; variableNames=[];

%get the name of the header file
headerFileName = [dataFileName,'.head'];


%First, make sure the files exist with some error checking...
f1=dir(dataFileName);
f2=dir(headerFileName);

if length(f2)~=1
    error('nfsim:readNFsimBinary:FileMissingError','   Could not find the binary data file.');
end
if length(f2)~=1
    error('nfsim:readNFsimBinary:FileMissingError','   Could not find the header file specified.');
end

%Get the size of the file we want to read (tells us how many rows
%there are in the file)
binaryFileSize = f1(1).bytes;

%Try to read the header file which tells us the names of each
%column (this also gives us the number of columns). Also here we'll
%do just a bit of error checking to make sure everything was formatted in
%the correct way
[rowCountInfo,variableNames] = tblread(headerFileName,'\t');
if(size(rowCountInfo)>1)
   error('nfsim:readNFsimBinary:RowCountError','   Error when reading header file: too many rows of information!');
end

%From the header file, extract out the information on column names, along
%with the number of columns and rows.  We need this info to read binary.
columnCount=size(variableNames,1);
rowCount=(binaryFileSize/8)/columnCount; %gives us number of elements in the
                                         %entire file (each of double size, 8
                                         %bytes) divided by the column
                                         %count thus giving the row count

%Try to open the binary data file in read only mode, with some error
%checking to make sure things worked
[fid, message] = fopen(dataFileName,'r');
if(fid==-1)
   error('nfsim:readNFsimBinary:FileCouldNotOpenError', ...
    ['Error when opening the .dat file named:\n\t', dataFileName, '\n\n', ...
   	'Matlab says:\n', ...
    '  ',message,'\n\n']);
   data = [];
   return;
end

%If that was successful, then read in all the data with double
%precision (which is the same way NFsim writes to binary files)
data = fread(fid,[columnCount,rowCount],'double')';

%Close up the file nicely
fclose(fid);
