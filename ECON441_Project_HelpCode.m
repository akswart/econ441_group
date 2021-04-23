% Parameters: alpha, beta, sigma, q
% Platform: fixed fee T, ad valorem tax t
% alpha: price-quality tradeoff in utility
% beta: price-quality tradeoff in platform choice
% sigma: technology parameter about how tightly the platform can make the
% choice

% Parameters

qstart      =           % enter fixed value of q1, ignore the normalized q2
alphastart  = 
sigmastart  =
fees        =           % ordered as a 2x1 vector of T first, then t.
cstart      =           % ordered as a 2x1 vector of c1 first, then the normalized c2.


%% Part 1 - Recreate the six figures of how various values differ across
% values of beta


% Beta evaluation points for six figures: YOU NEED TO CREATE
% Evaluate 1,000 points evenly spaced from 1/100 to 10. 
% Below I refer to the number of beta points as 'numbeta' and the
% evaluation points themselves 'numvec'


% create placeholder vectors of zeros where you will fill in the resulting
% value for each beta you're evaluating.
p1holder = zeros(numbeta,1);    % P1 placeholder
p2holder = zeros(numbeta,1);    % P2 placeholder
objholder = zeros(numbeta,1); 	% Objective function (satisfy optimality) placeholder
csholder = zeros(numbeta,1);    % Consumer Surplus placeholder
piholder = zeros(numbeta,1);    % Profit placeholder

% Loop over each value of beta
for ii=1:numbeta
    
vectry = [qstart;alphastart;numvec(ii);sigmastart;fees]; % arrange the vector of parameters you're evaluating at: q/alpha/beta/sigma/T/t <-- MUST be in th is order because subsequent functions assume this order.

h = @(x) func_foc_costs(x,vectry,cstart); % define objective and performance to solve FOC=0 with given parameters and compute corresponding Consumer Surplus and profit associated with it. x is pair of prices, vectry is the params considering, cstart is costs
tempp = func_find_prices(1000,vectry,cstart); % function solves for equilibrium prices [p1,p2] given parameters

[tempobj,tempcs,temppi] = h(tempp); % given the equilibrium prices (tempp), compute the outputs to store and compare across betas.

p1holder(ii)    = tempp(1);
p2holder(ii)    = tempp(2);
objholder(ii)   = tempobj;
csholder(ii)    = tempcs;
piholder(ii)    = temppi;
end


% Example of how to plot to look like in assignment (for p_1)
figure
plot(numvec,p1holder,'LineWidth',4)
title(['\fontsize{20}P1 by \beta (q=1, \alpha=0.5, \sigma=1, c=0.5)'])
xlabel(['\fontsize{20}\beta'])
ylabel(['\fontsize{20}P1'])
set(gca,'fontsize',16)
saveas(gcf,'p1_by_beta.png')




%% Part 2: Recreate four subfigures of optimal beta (maximizes CS).
% Vary q, alpha, sigma, c and find the optimal beta for each. Hold 
% the other paramters fixed at their standard values

% This time, evaluate 800 points of beta, evenly spaced from 1/100 to 8.



% alpha - evaluate from .1 to .6 with step size of .1.

alphastore  = zeros(length(alphavec),numbeta2); % for each pair of evaluation points (alpha,beta), along with all other params standard, save the CS here.
% You'll want to create a double for loop to do this.

% Then, for each value of alpha, find the optimal beta (in terms of
% maximizing CS) and save that to this next object.
alphaoptbeta = zeros(length(alphavec),1);


% c - vary c from 0 to 1 with step size of .1. Then follow similar approach
% as alpha
cstore      = zeros(length(cvec),numbeta2);
coptbeta = zeros(length(cvec),1);

% q - vary from 0 to 5 with step size of .1. Then follow similar approach
% as above. Note this is a little slower as it is more evaluation points.
qcstore     = zeros(length(qcvec),numbeta2);
qcoptbeta   = zeros(length(qcvec),1);


% sigma - vary from .1 to 1.5 with step size of .1. Then follow similar
% approach as above.
sigmastore  = zeros(length(sigmavec),numbeta2);
sigmaoptbeta = zeros(length(sigmavec),1);



