% RUNTLBRFIT
%
%   This script provides an example that can be modified for other models
%   of how the set of Matlab routines can be used to run an NFsim model,
%   scan parameters of an NFsim model, and fit the parameters to some
%   experimental data.  Because data fitting is a problem that is problem
%   specific (meaning that the type and format of the experimental data is
%   specific to a problem), this is not a generalized script.  Rather, it
%   shows what steps you have to take to modify this code for your model.
%
%
%   The basic steps of defining a fitting routine in Matlab to run an NFsim
%   model is to 1) write an evaluate method that runs the NFsim simulation
%   you need, reads in the results, and compares the results to data (see
%   evaluateTLBRparams.m file in this directory).  Then, determine the
%   algorithm and options you want to run your fit.  Finally, use the
%   scripts below to actually execute the fit.
%
%   Follow and modify the outline below to run your fit.
%
%
%   created by Michael Sneddon, 8/19/2010






% first we clear the console and memory, sot that we can start fresh
clear; clc;

% for comparing the NFsim simulations to the TLBR data, we are going to
% create two global variables, defined as follows, to store the
% information.  We do this so that we do not have to read in the data fresh
% everytime the evaluateTLBRparams.m function is called - instead we can
% read it once here and refer to it in any of the other parameter scanning
% and fitting scripts
global ydata;
global xdata;



% In the TLBR example, we fit 2 free parameters at the same time.  Those
% parameters are K1 and K2 that control initial binding of a free ligand to
% a receptor and the subsequent crosslinking reaction.  As for any fit, you
% must choose the initial starting values for constraing your parameters.
% This is something you will have to play around with.  Remember that for
% some models, not all parameters may be properly constrained!  So you may
% converge to a different value depending on your initial conditions.  This
% is a general problem of fitting these types of models.  For some ideas in
% handling this, check out this paper: Gutenkunst, et al. Universally
% Sloppy Parameter Sensitivities in Systems Biology Models. PLoS Comp.
% Biol., 2007.
    
    
% Here are the default parameters for the TLBR model - note that it is
% sometimes useful to start with random initial values within some range,
% then run many fitting trials.  For this model, these values are in units
% of per nanomolar, because they are binding constants.  The parameters are
% then saved to an array named parameters.
K1 = 2; %rand(1)*10;
K2 = 50; %rand(1)*100;
parameters = [K1,K2];

% read in the experimental data to fit.  This data is used here courtesy of
% Michael Monine, originally published in Monine, et al. Modeling
% Multivalent Ligand-Receptor Interactions with Steric Constraints on
% Configurations of Cell-Surface Receptor Aggregates. Biophys. J.
% 98(1):48-56, 2010.
rawData=dlmread('tlbrExample/data/exp_lambda_TLBR.txt');

% in this example, we use the raw data to generate which runs of NFsim we
% would like.  Because there are more data points than we want to use to
% compare to, we can put in a step that allows us to fit to every other
% point, which will give us quicker results, but will be not quite as good
% a fit in the end.
stepOfRawPoints = 1;
dataToFit = rawData(1:stepOfRawPoints:end,:);

% save the raw data to our global functions.  Note that in the TLBR
% example, we have to scale the data by 0.85 as was done in the original
% data in order to directly compare the binding data.
xdata = dataToFit(:,1);
ydata = dataToFit(:,2).*0.85;


% here is where we set the options for the fit.  These were adjusted to
% create good fits for the TLBR data, but in general should be tuned for
% your specific model.  The choice of fitting method depends on the data to
% fit, the type of model you are running, the parameter values you are
% fitting, and ultimately what precision you need to fit your data.  See
% matlab documentation for the meaning of these options, and other options
% that are available to you.
options = optimset();
options = optimset(options,'Display','iter');  % [ off | iter | notify | final ]
options = optimset(options,'algorithm','trust-region-reflective');
options = optimset(options,'MaxIter',50);
options = optimset(options,'TolX',1e-6);
options = optimset(options,'TolFun',1e-6);
%options = optimset(options,'LevenbergMarquardt','on');
%options = optimset(options,'PrecondBandWidth',1);

options = optimset(options,'DiffMinChange',0.4);
options = optimset(options,'DiffMaxChange',5);
options = optimset(options,'FinDiffType','central'); % or forward


% here we can set the upper and lower bounds for the parameter values that
% we are trying to fit.  These two arrays must be the same length as the
% parameters array.  For the TLBR model, these are binding constants so
% they should be positive values (lower bound set to zero), but could be as
% high as we need, so the upperBound is set to infinity (inf).
lowerBound = zeros(length(parameters),1);
upperBound = zeros(length(parameters),1)+inf;



%  finally, we can actually run the fit using the lsqcurvefit routine.
%  LSQCURVEFIT is a routine in the optimization toolbox, which most users
%  who have access to matlab should have.  If you do not have this
%  function, there are other fitting routines you can use that may be
%  suffiecient for some cases, but we highly recommend you get the
%  optimization toolbox because it is designed to fit data!

% notice here also that we use a reference the evaluateTLBRparams file that
% is located in this directory.  You will have to modify  that file too in
% order to perform your fit.
[paramhat,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(  @evaluateTLBRparams, ... 
                      parameters, ...
                      log(xdata), ...
                      ydata, ...
                      lowerBound, ...
                      [], ...
                      options);
                  
% The fitted parameter values are saved to the output array called paramHat
% and can then be accessed by you
fprintf('\n\nparamHat\n');
for i=1:length(paramhat)
    fprintf([num2str(paramhat(i)),';\n']);
end;

% lsqcurvefit gives a number of other output options so that you can plot the
% confidence intervals of your measurements.  These methods were developed
% for continuous models, however, so may not be correctly estimating the
% error!  (This is particularly a problem when it is the ratio of two
% parameters that control most of the variability, not just a single
% parameter alone).  For instance, these parameter intervals are WRONG for
% the TLBR fitting example.  But, for some cases, they may be ok, so here
% is how you access them.
fprintf('\n\nConfidence Intervals\n');
confInt = nlparci(paramhat,residual,'jacobian',jacobian);
for i=1:length(paramhat);
   fprintf([num2str(paramhat(i)),' interval: ',num2str(confInt(i,1)),' to ',num2str(confInt(i,2)),'\n']);
end;

fprintf(['\n\nNorm of Residuals:   ', num2str(sum(resnorm)),'\n']);



return;



% Finally, here is some code that allows you to plot the final fit for the
% TLBR example.


% get the data to fit
rawData=dlmread('tlbrExample/data/exp_lambda_TLBR.txt');
dataToFit = rawData(:,:);
xdata = dataToFit(:,1);
ydata = dataToFit(:,2).*0.85;    

% Here are the estimated parameter values.  You can set this to whatever
% values you need
K1 = 0.327;
K2 = 52.23;

% set the default parameters that tell us the model and the paths
pathToModel = 'tlbrExample/';
bnglFileName = 'tlbr.bngl';
pathToNFsim = '../../../NFsim_v1.06';
pathToOutput = 'tlbrExample/output/';


% Create a small loop that runs NFsim on the model for multiple
% iteratations.  This is actually an example of a parameter scan (we are
% scanning the Ligand concentration, to create a dose response curve).
tic;
for k=1:length(xdata)
    fprintf('.');
    runNFsimOnce(pathToModel,bnglFileName,pathToNFsim,pathToOutput, ... 
         k,{'Lig_conc','K1','K2'},[xdata(k);K1;K2]);
end;
toc;

% Parse the results of this fit, and save them to these vectors.  When we
% look at the output, we want to plot the steady-state value of bound
% receptors to ligand, normalized by the total number of binding sites
% available to the ligand.  For that, we average from output step 75 to the
% end, and divide by the number of available binding site.  This is what is
% meant by problem specific comparison to data.
fittedData = zeros(size(xdata(:,1)));
crossLinked = zeros(size(xdata(:,1)));
for k=1:length(xdata)
    [x,header] = tblread([pathToOutput,'PS_tlbr_',num2str(k),'.gdat']);
    fittedData(k) = mean((x(75:end,3)-(x(75:end,2)./6)) ./ (2*300));
    crossLinked(k) = mean( x(75:end,4) ./ (4*300) );
end;



% finally, here is the script to plot the results
fontsize = 14;
figure; hold on; box on;
set(gcf,'color','white');

plot(xdata,ydata,'ro');
plot(xdata,fittedData,'b-');

set(gca,'xtick',[10^-3,10^-2,10^-1,1,10,100]);
set(gca,'XScale','log');
axis([10^-3.5,10^2.8,-0.05,1]);
set(gca,'TickLength',[0.02 0.025]);

ylabel('Fraction of bound ligand','FontName','Arial','fontSize',fontsize);
xlabel('Ligand concentration (nM)','FontName','Arial','fontSize',fontsize);
legend('Experimental data','Fitted simulation curve','Location','NorthWest');
legend boxoff;

set(gca,'FontName','Arial');
set(gca,'fontSize',fontsize);
   

    
    
    
    











