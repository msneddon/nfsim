%  Validation script for NFsim
%
%  Runs all BioNetGen files in the models directory with names v1.bngl,
%  v2.bngl, v[...].bngl, ... v[N].bngl, using the run settings for NFsim from
%  the cooresponding r1.bngl, ... r[N].bngl files and flags any suspicious
%  looking runs based on the final observable count.  The

clear; clc; close all; drawnow;
fprintf('=======================\n');
fprintf('running nfsim validation...\n\n');

% Identify the path to NFsim and BioNetGen that you want to test
%bngPath = '/home/msneddon/bionetgen/bng2/BNG2.pl';  % <<< SET YOUR BIONETGEN PATH HERE <<<
bngPath = '/home/justin/bng-vers/BioNetGen_Testing/BNG2.pl';  % <<< SET YOUR BIONETGEN PATH HERE <<<
nfsimPath = '../bin/NFsim';                         % <<< SET YOUR NFSIM PATH HERE <<<
deleteTempFiles = 1; %set to 1 to delete temporary files, 0 otherwise
tol = 0.35; %this is the error tolerance when comparing nfsim's run to the ssa where 0.35 = 35%

viewAllPlots =1; % set to 1 to display, zero otherwise

% get the list of files
modelDirName = 'basicModels';
[files] = dir(modelDirName);




% loop over each of the files
testCounter = 1;
fprintf('status   model         Description\n');
for f=1:length(files)
    
    % make sure we are looking at a valid bngl file
    if files(f).isdir==1, continue; end;
    if length(files(f).name)<6, continue; end;
    if strcmp(files(f).name(end-4:end),'.bngl')~=1, continue; end;
    
    
    % set the filenames and find the associated run settings file
    bngFilename = files(f).name;
    runFilename = ['r',files(f).name(2:end-5),'.txt'];
    if exist([modelDirName,'/',runFilename],'file')~=2
        fprintf(2,['Cannot find associated run file for model: "',files(f).name,'".  skipping.\n']);
        continue;
    end;
    
    %get the run options and descriptions
    %give no delimeter so that we just always break at end of lines
    runFile=textread([modelDirName,'/',runFilename],'%s','delimiter','');
    description = runFile{1};
    runOptions = runFile{2};
    
    % run the test 3 times and average the results
    ssaDiff = []; nfDiff = []; maxRuns=10;
    for runs = 1:maxRuns
    
        %fprintf( 'perl %s -outdir %s %s/%s\n', bngPath, modelDirName, modelDirName, bngFilename );
        [status,bngTerminal]=system(['perl ',bngPath,' -outdir ',modelDirName,' ',modelDirName,'/',bngFilename]);
        [status,nfTerminal]=system([nfsimPath,' -xml ',modelDirName,'/',bngFilename(1:end-5),'.xml ', ...
            ' -o ',modelDirName,'/',bngFilename(1:end-5),'_nf.gdat ' runOptions]);

        %analyze the test
        [ode,headers] = tblread([modelDirName,'/',bngFilename(1:end-5),'_ode.gdat']);
        [ssa,headers] = tblread([modelDirName,'/',bngFilename(1:end-5),'_ssa.gdat']);
        [nf,headers] = tblread([modelDirName,'/',bngFilename(1:end-5),'_nf.gdat']);

        %to do our "naive" analysis, we simply compare at all time points and
        %all observables the difference between the ode and ssa, and the ode
        %and nfsim results.  If nfsim output is within the acceptable window
        %(which we determine as < absolute value of the error in the ssa + some
        %percentage, x, of the error in the ssa) then we pass, otherwise we
        %fail.  The constraint must hold for each variable
        thisSsaDiff = sum(abs(ode(:,2:end)-ssa(:,2:end)));
        thisNfDiff = sum(abs(ode(:,2:end)-nf(:,2:end)));
        
        if runs==1
            ssaDiff= thisSsaDiff;
            nfDiff= thisNfDiff;
        else
            ssaDiff= ssaDiff+thisSsaDiff;
            nfDiff= nfDiff+thisNfDiff;
        end;
        
    end;
    ssaDiff=ssaDiff./maxRuns;
    nfDiff=nfDiff./maxRuns;
    
    
    passed = true;
    badObservables = [];
    for i=1:length(ssaDiff)
        if nfDiff(i)>(ssaDiff(i)+tol*ssaDiff(i))
            passed = false;
            badObservables = [badObservables;i]; %#ok<AGROW>
        end;
    end;
    
    
    
    % some final output messaging
    if passed,
        fprintf(['[pass]   ',files(f).name]);
        x = length(files(f).name);
        while x<13,fprintf(' '); x=x+1; end;
        fprintf([' ',description,'\n']);
    else
        fprintf(1,'FAIL!!'); fprintf(['   ',files(f).name]);
        x = length(files(f).name);
        while x<13,fprintf(' '); x=x+1; end;
        fprintf([' ',description,'\n']);
        
        fprintf('   ---> failing observables are: ');
        for i=1:length(badObservables)
            fprintf([' ',num2str(badObservables(i))]);
        end;
        fprintf('\n');
        
        
        for i=1:length(badObservables)
            fprintf(['   ---> [',num2str(badObservables(i)),']: ssaDiff = ']);
            fprintf([num2str(ssaDiff(badObservables(i))),', nfDiff = ']);
            fprintf([num2str(nfDiff(badObservables(i))),', accepted tolerance = ']);
            fprintf([num2str(ssaDiff(badObservables(i))+tol*ssaDiff(badObservables(i))),'\n']);
        end;
        fprintf('\n');
        
    end;
    
    
    if ~passed
        % display the failed 
        figure; hold on; set(gcf,'color','white');
        t=title(['Error detected in model: ',files(f).name]);
        plot(ssa(:,1),ssa(:,2:end),'-');
        plot(ode(:,1),ode(:,2:end),'k-');
        plot(nf(:,1),nf(:,2:end),'-');
        plot(nf(:,1),nf(:,2:end),'.');
        l=legend(headers(2:end,:),'Location','NorthEastOutside');
        set(l,'Interpreter','none');
        set(t,'Interpreter','none');
        xlabel('Time');
        ylabel('Molecule Number');
        drawnow;
        
    elseif viewAllPlots
         % display if requested
        figure; hold on; set(gcf,'color','white');
        t=title(['Passed model: ',files(f).name]);
        plot(ssa(:,1),ssa(:,2:end),'-');
        plot(ode(:,1),ode(:,2:end),'k-');
        plot(nf(:,1),nf(:,2:end),'-');
        plot(nf(:,1),nf(:,2:end),'.');
        l=legend(headers(2:end,:),'Location','NorthEastOutside');
        set(l,'Interpreter','none');
        set(t,'Interpreter','none');
        xlabel('Time');
        ylabel('Molecule Number');
        drawnow;
    end;
    
end;




%remove temporary files
if deleteTempFiles,
    system(['rm -f ',modelDirName,'/*.cdat']);
    system(['rm -f ',modelDirName,'/*.gdat']);
    system(['rm -f ',modelDirName,'/*.net']);
    system(['rm -f ',modelDirName,'/*.xml']);
end;


