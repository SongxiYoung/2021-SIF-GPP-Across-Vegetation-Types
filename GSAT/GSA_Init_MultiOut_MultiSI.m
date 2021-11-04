%% GSA_Init_MultiOut_MultiSI: initialize the variables used in the GSA computation
%
% Usage:
%   pro = GSA_Init_MultiOut_MultiSI(pro)
%
% Inputs:
%    pro                project structure
%
% Output:
%    pro                updated project structure
%
% ------------------------------------------------------------------------
% Citation: Cannavo' F., Sensitivity analysis for volcanic source modeling quality assessment and model selection, Computers & Geosciences, Vol. 44, July 2012, Pages 52-59, ISSN 0098-3004, http://dx.doi.org/10.1016/j.cageo.2012.03.008.
% See also
%
% Author : Flavio Cannavo'
% e-mail: flavio(dot)cannavo(at)gmail(dot)com
% Release: 1.0
% Date   : 15-02-2011
%
% History:
% 1.0  15-01-2011  First release.
%      06-01-2014  Added comments.
%%

function pro = GSA_Init_MultiOut_MultiSI(pro)

% get two sets of samples of the input variables
[E, T] = fnc_SampleInputs(pro);

% Store the sample sets in the project structure array
pro.SampleSets.E = E;
pro.SampleSets.T = T;

% get the number of input variables
n = length(pro.Inputs.pdfs);

% set the number of possible combinations of the input variables
L = 2^n;
N = pro.N;

% Parfor loop requirements
model_function_handle = pro.Model.handle;
fE_cell = cell(N,1);

% Evaluate the model at the sample points in set E
parfor j=1:N
    fE_cell{j} = feval(model_function_handle,E(j,:)); % the function output must be a single variable that takes the form of a row vector
end

% Convert the cell array to a numeric array, and store the simulation 
% results in the project structure array
pro.GSA.fE = cell2mat(fE_cell); % this will only work if all output generated by the function is numeric

% Number of output variables
output = size(pro.GSA.fE,2);

% calculate the mean value for all the outcome variables
pro.GSA.mfE = nanmean(pro.GSA.fE);

pro.GSA.mfE(isnan(pro.GSA.mfE)) = 0;

% calculate the mean value of the differences from the mean (it will be 0)
pro.GSA.f0 = nanmean(pro.GSA.fE - repmat(pro.GSA.mfE,size(pro.GSA.fE,1),1));

% calculate the total variance of the model outcomes
pro.GSA.D  = nanmean((pro.GSA.fE - repmat(pro.GSA.mfE,size(pro.GSA.fE,1),1)).^2) - pro.GSA.f0.^2;

% approximate the error of the mean value
pro.GSA.ef0 = 1.96*sqrt(pro.GSA.D./sum(~isnan(pro.GSA.fE - repmat(pro.GSA.mfE,size(pro.GSA.fE,1),1))));

% prepare the structures for the temporary calculations of sensitivity
% coefficients
pro.SampleSets.H = cell(1,L-1);

pro.GSA.fH = cell(1,L-1);

pro.GSA.Dmi = nan(output,L-1);
pro.GSA.eDmi = nan(output,L-1);

pro.GSA.Di = nan(output,L-1);
pro.GSA.eDi = nan(output,L-1);

pro.GSA.GSI = nan(output,L-1);
pro.GSA.eGSI = nan(output,L-1);