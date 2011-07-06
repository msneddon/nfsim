function [value] = evaluateTLBRparams(parameters,xdata)
%
%  Before you modify this file, be sure to look at runTLBRfit.m to see how
%  this function can be called to fit parameter values!  The purpose of
%  this function is to run the NFsim simulations that you need for a given
%  set of parameter values, then return the result so that they can be
%  compared (by Matlab) to the experimental data to constrain your results
%  and estimate your parameter values.
%
%
%



% first things first - get the experimental results
global ydata;
xdata = exp(xdata);



% set to 1 to display the fit as we run, 0 otherwise.  It is often useful
% to display the fit as you go so you can see if your fitting options are
% working or not
displayFit = 1;
fprintf([' params: [K1=',num2str(parameters(1),10),',K2=',num2str(parameters(2),10),']:\n']);


% set the default parameters that tell us the model and the paths to run
% NFsim.  This is standard for calling the runNFsimOnce routine.
pathToModel = 'tlbrExample/';
bnglFileName = 'tlbr.bngl';
pathToNFsim = '../../';
pathToOutput = 'tlbrExample/output/';



% Run the NFsim model with the parameters given to us by Matlab.  These
% parameters are stored in the parameters array.  Here, for the TLBR model,
% we actually have to perform a parameter scan as well across the ligand
% concentrations.  You will have to modify this code yourself, depending on
% your specific model and what results from the runs you need to extract.
tic;
for k=1:length(xdata)
    fprintf('.');
    runNFsimOnce(pathToModel,bnglFileName,pathToNFsim,pathToOutput, ... 
         k,{'Lig_conc','K1','K2'},[xdata(k);parameters(1);parameters(2)]);
end;
toc;



% Parse the results of this fit, and save them to these vectors that can be
% compared to the experimental data by Matlab.  Here we are looking at the
% fitted data and the crosslinked number of receptors.  We want the value
% at steady-state, so we can average over the last time points and 
fittedData = zeros(size(xdata(:,1)));
crossLinked = zeros(size(xdata(:,1)));
for k=1:length(xdata)
    [x,header] = tblread([pathToOutput,'PS_tlbr_',num2str(k),'.gdat']);
    fittedData(k) = mean((x(75:end,3)-(x(75:end,2)./6)) ./ (2*300));
    crossLinked(k) = mean( x(75:end,4) ./ (4*300) );
end;

%, finally, set the output value the minimizer will use
value = fittedData;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if you selected earlier to display the fit, this is where it is done.
if displayFit,
    
    fontsize = 14;
    figure(1); hold on; 
    cla; box on;
    set(gcf,'color','white');

    plot(xdata,ydata,'ro');
    plot(xdata,fittedData,'k-');
    %plot(xdata,crossLinked,'b-');


    set(gca,'XScale','log');
    axis([10^-3.5,10^2.8,-0.05,1]);

    ylabel('Fraction of bound ligand','FontName','Arial','fontSize',fontsize);
    xlabel('Ligand concentration (nM)','FontName','Arial','fontSize',fontsize);
    legend('Experimental data','Fitted simulation curve','Location','NorthWest');
    legend boxoff;
    
    set(gca,'FontName','Arial');
    set(gca,'fontSize',fontsize);
    
    % be sure to add the drawnow command, so the figure is updated even
    % though matlab is still running the fitting routine elsewhere.
    drawnow;
    
    
    
end;

    