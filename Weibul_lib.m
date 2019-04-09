function wbl=Weibul_lib()
wbl.weibul_fit=@weibul_fit;                                 %[t,sd]=weibul_fit(q)
wbl.const_weibul_fit=@const_weibul_fit;                     %t=const_weibul_fit(q)
wbl.Weibul_logP=@Weibul_logP;                               %LogP=Weibul_logP(q,t,sd)
wbl.init_Weibul_link_functions=@init_Weibul_link_functions; %init_Weibul_link_functions(q);

%*** The following functions assume Weib_glb beeing initialized, either by calling init_Weibul_link_functions, weibul_fit, 
%*** const_weibul_fit, or Weibul_logP

wbl.Weibul_cdf=@Weibul_cdf;    %p=Weibul_cdf(x)
wbl.Weibul_inv=@Weibul_inv;    %x=Weibul_inv(p)
wbl.Weibul_pdf=@Weibul_pdf;        %dp=Weibul_pdf(x)

end



%************************ Exported :*********************************************************************
function [t,sd,sd_t,sd_sd]=weibul_fit(q,plot_logP)
global Weib_glb
if isfield(q,'trialCount') & length(q.intensity)>q.trialCount,
   q.intensity=q.intensity(1:q.trialCount);
   q.response=q.response(1:q.trialCount);
end;
if nargin<2,
   plot_logP=[];
end;
if isempty(plot_logP),
   plot_logP=true;
end;
init_Weibul_link_functions(q);
Weib_glb.plot_logP=plot_logP;

wstate = warning('off', 'stats:glmfit:BadScaling');
[b,dev,stats]=glmfit_mystartvalue(q.intensity',q.response','binomial','link',{@Weibul_inv, @dlink_f, @Weibul_cdf},'init_fun',@my_startingVals,'iterLim',1000);
warning(wstate);

t=-b(1)/b(2);
sd=1/b(2);

init_val=Weib_glb.init_val;  % save init_val before Weib_glb is cleared in Weibul_logP!!
LogP=Weibul_logP(q,t,sd);
if init_val.logP>LogP || sd<0,
   t=init_val.t;
   sd=init_val.sd;
   sd_t=NaN;
   sd_sd=NaN;
else
   df1=-1/b(2);
   df2=b(1)/b(2)^2;
   sd_t=(stats.covb(1,1)*df1^2+stats.covb(2,2)*df2^2+2*stats.covb(1,2)*df1*df2)^0.5;
   df2=1/b(2)^2;
   sd_sd=stats.se(2)*df2;
end;

end


function [t,sd_t]=const_weibul_fit(q,plot_logP)
global Weib_glb
if isfield(q,'trialCount') & length(q.intensity)>q.trialCount,
   q.intensity=q.intensity(1:q.trialCount);
   q.response=q.response(1:q.trialCount);
end;
if nargin<2,
   plot_logP=[];
end;
if isempty(plot_logP),
   plot_logP=true;
end;
init_Weibul_link_functions(q);
Weib_glb.plot_logP=plot_logP;

wstate = warning('off', 'stats:glmfit:BadScaling');
[b,dev,stats]=glmfit_mystartvalue([],q.response','binomial','link',{@Weibul_inv, @dlink_f, @Weibul_cdf} ...
   ,'init_fun',@my_const_startingVals,'iterLim',1000,'offset',q.intensity'*q.beta);
warning(wstate);

t=-b(1)/q.beta;
sd=1/q.beta;

init_val=Weib_glb.init_val;  % save init_val before Weib_glb is cleared in Weibul_logP!!
LogP=Weibul_logP(q,t,sd);
if init_val.logP>LogP || sd<0,
   t=init_val.t;
   sd_t=NaN;
else
   df1=-1/q.beta;
   sd_t=(stats.covb(1,1)*df1^2)^0.5;
end;

end

function LogP=Weibul_logP(q,t,sd)
if isfield(q,'trialCount') & length(q.intensity)>q.trialCount,
   q.intensity=q.intensity(1:q.trialCount);
   q.response=q.response(1:q.trialCount);
end;
init_Weibul_link_functions(q);

y=logical(q.response');
x=q.intensity';
LogP=sum(log(Weibul_cdf((x(y)-t)/sd)));
LogP=LogP+sum(log(1-Weibul_cdf((x(~y)-t)/sd)));
end


%*** The following functions assume Weib_glb beeing initialized, either by calling weibul_fit, or Weibul_logP
function p=Weibul_cdf(x)
global Weib_glb
x(x<Weib_glb.lowerBnd) = Weib_glb.lowerBnd;
x(x>Weib_glb.upperBnd) = Weib_glb.upperBnd;
p=Weib_glb.q.delta*Weib_glb.q.gamma+(1-Weib_glb.q.delta)*(1-(1-Weib_glb.q.gamma)*exp(-10.^(x)));
end
function x=Weibul_inv(p)
global Weib_glb
   p(p<Weib_glb.q.gamma+eps(class(p)))=Weib_glb.q.gamma+eps(class(p));
   p(p>1-Weib_glb.q.delta*(1-Weib_glb.q.gamma)-eps(class(p)))=1-Weib_glb.q.delta*(1-Weib_glb.q.gamma)-eps(class(p));
   x=log(-log((1-(p-Weib_glb.q.delta*Weib_glb.q.gamma)/(1-Weib_glb.q.delta))/(1-Weib_glb.q.gamma)))/log(10);
end
function dp=Weibul_pdf(x)
global Weib_glb
x(x<Weib_glb.lowerBnd) = Weib_glb.lowerBnd;
x(x>Weib_glb.upperBnd) = Weib_glb.upperBnd;
%**** compute the derivative of the shifted Weibul-CDF. The function is shifted such that Weibul-CDF(0)=q.pThreshold ************
dp=(1-Weib_glb.q.delta)*(1-Weib_glb.q.gamma)*exp(-10.^(x)).*10.^(x)*log(10);
end



%************************ Private :*********************************************************************
function init_Weibul_link_functions(q)
global Weib_glb
if isfield(q,'trialCount') & length(q.intensity)>q.trialCount,
   q.intensity=q.intensity(1:q.trialCount);
   q.response=q.response(1:q.trialCount);
end;
x=1;
Weib_glb=[];
Weib_glb.q=q;
Weib_glb.upperBnd = Weibul_inv(1-Weib_glb.q.delta*(1-Weib_glb.q.gamma)-eps(class(x))); 
Weib_glb.lowerBnd = Weibul_inv(Weib_glb.q.gamma+eps(class(x)));
end


%*** The following functions assume Weib_glb beeing initialized, either by calling weibul_fit, or Weibul_logP

function mu=my_startingVals(distr,y,N)
global Weib_glb
t=QuestMean(Weib_glb.q);
sd=QuestSd(Weib_glb.q);


plot_logP=Weib_glb.plot_logP;  % save plot_logP before Weib_glb is cleared in Weibul_logP!!

x=t+(-10:0.5:10)';
yrange=[min(1/Weib_glb.q.beta/2,0.01),max(2/Weib_glb.q.beta,0.04)];
y=yrange(1)+(0:40)'*diff(yrange)/40;
X=repmat(x',length(y),1);
Y=repmat(y,1,length(x));
P=zeros(length(y),length(x));
for j=1:length(y),
for i=1:length(x),
P(j,i)=Weibul_logP(Weib_glb.q,X(j,i),Y(j,i));
end;
end;

if plot_logP,
   FS=12;
   
   fignr=gcf;
   figure(100);
   clf
   set(gca,'Linewidth',1.5,'Fontsize',FS);
   surf(X,Y,P);
   set(get(gca,'Xlabel'),'string','t','Fontsize',FS);
   set(get(gca,'Ylabel'),'string','sd','Fontsize',FS);
   set(get(gca,'Zlabel'),'string','logP','Fontsize',FS);
   
   
   figure(fignr);
end;
i=find(P(:)==max(P(:)),1,'first');
t0=X(i);
sd0=Y(i);
fprintf('      StartingValue:\n');
fprintf('t=%10.4f%11s; sd=%10.4f%11s; LogP=%10.4f\n',t0,' ',sd0,' ',P(i));

mu=Weibul_cdf((Weib_glb.q.intensity'-t0)/sd0);
Weib_glb.init_val.t=t0;
Weib_glb.init_val.sd=sd0;
Weib_glb.init_val.logP=P(i);

end


function mu=my_const_startingVals(distr,y,N)
global Weib_glb
t=QuestMean(Weib_glb.q);
sd=1/Weib_glb.q.beta;


plot_logP=Weib_glb.plot_logP;  % save plot_logP before Weib_glb is cleared in Weibul_logP!!

x=t+(-10:0.5:10)';
P=zeros(1,length(x));
for i=1:length(x),
P(i)=Weibul_logP(Weib_glb.q,x(i),sd);
end;

if plot_logP,
   FS=12;
   
   fignr=gcf;
   figure(100);
   clf
   set(gca,'Linewidth',1.5,'Fontsize',FS);
   plot(x,P);
   set(get(gca,'Xlabel'),'string','t','Fontsize',FS);
   set(get(gca,'Ylabel'),'string','logP','Fontsize',FS);
   
   
   figure(fignr);
end;
i=find(P(:)==max(P(:)),1,'first');
t0=x(i);
fprintf('      StartingValue:\n');
fprintf('t=%10.4f%11s; LogP=%10.4f\n',t0,' ',P(i));

mu=Weibul_cdf((Weib_glb.q.intensity'-t0)/sd);
Weib_glb.init_val.t=t0;
Weib_glb.init_val.sd=sd;
Weib_glb.init_val.logP=P(i);

end

function y=dlink_f(mu)
y=1 ./ Weibul_pdf(Weibul_inv(mu));
end





%***************************************************************************************************
%***** all subsequent functions belong to a modified version of the matlab tool glmfit *************
%***************************************************************************************************
function [bb,dev,stats] = glmfit_mystartvalue(x,y,distr,varargin)
%GLMFIT Fit a generalized linear model.
%   B = GLMFIT(X,Y,DISTR) fits a generalized linear model using the
%   predictor matrix X, response Y, and distribution DISTR.  The result B
%   is a vector of coefficient estimates.  Acceptable values for DISTR are
%   'normal', 'binomial', 'poisson', 'gamma', and 'inverse gaussian'.  The
%   distribution is fit using the canonical link corresponding to DISTR.
%
%   X is a matrix with rows corresponding to observations, and columns to
%   predictor variables.  GLMFIT automatically includes a constant term in the
%   model (do not enter a column of ones directly into X).  Y is a vector of
%   response values.  If DISTR is 'binomial' Y may a binary vector indicating
%   success/failure, and the total number of trials is taken to be 1 for all
%   observations.  If DISTR is 'binomial', Y may also be a two column matrix,
%   the first column containing the number of successes for each observation,
%   and the second containing the total number of trials.
%
%   GLMFIT treats NaNs in X and Y as missing data, and removes the
%   corresponding observations.
%
%   B = GLMFIT(X,Y,DISTR,'PARAM1',val1,'PARAM2',val2,...) allows you to
%   specify optional parameter name/value pairs to control the model fit.
%   Parameters are:
%
%      'link' - the link function to use in place of the canonical link.
%         The link function defines the relationship f(mu) = x*b
%         between the mean response mu and the linear combination of
%         predictors x*b.  Specify the link parameter value as one of
%            - the text strings 'identity', 'log', 'logit', 'probit',
%              'comploglog', 'reciprocal', 'loglog', or
%            - an exponent P defining the power link, mu = (x*b)^P for
%              x*b > 0, or
%            - a cell array of the form {FL FD FI}, containing three
%              function handles, created using @, that define the link (FL),
%              the derivative of the link (FD), and the inverse link (FI).
%
%      'estdisp' - specify as 'on' to estimate a dispersion parameter for
%         the binomial or Poisson distribution in computing standard
%         errors, or 'off' (the default) to use the theoretical dispersion
%         value.  GLMFIT always estimates the disperson for other
%         distributions.
%
%      'offset' - a vector to use as an additional predictor variable, but
%         with a coefficient value fixed at 1.0.
%
%      'weights' - a vector of prior weights, such as the inverses of the
%         relative variance of each observation.
%
%      'constant' - specify as 'on' (the default) to include a constant
%         term in the model, or 'off' to omit it.  The coefficient of the
%         constant term is the first element of B.
%
%   [B,DEV] = GLMFIT(...) returns the deviance of the fit.
%
%   [B,DEV,STATS] = GLMFIT(...) returns a structure that contains the
%   following fields:
%       'dfe'       degrees of freedom for error
%       's'         theoretical or estimated dispersion parameter
%       'sfit'      estimated dispersion parameter
%       'se'        standard errors of coefficient estimates B
%       'coeffcorr' correlation matrix for B
%       'covb'      estimated covariance matrix for B
%       't'         t statistics for B
%       'p'         p-values for B
%       'resid'     residuals
%       'residp'    Pearson residuals
%       'residd'    deviance residuals
%       'resida'    Anscombe residuals
%
%   Example:  Fit a probit regression model for y on x.  Each y(i) is the
%   number of successes in n(i) trials.
%
%       x = [2100 2300 2500 2700 2900 3100 3300 3500 3700 3900 4100 4300]';
%       n = [48 42 31 34 31 21 23 23 21 16 17 21]';
%       y = [1 2 0 3 8 8 14 17 19 15 17 21]';
%       b = glmfit(x, [y n], 'binomial', 'link', 'probit');
%       yfit = glmval(b, x, 'probit', 'size', n);
%       plot(x, y./n, 'o', x, yfit./n, '-')
%
%   See also GLMVAL, REGSTATS, REGRESS.

%   References:
%      [1] Dobson, A.J. (2001) An Introduction to Generalized Linear
%          Models, 2nd edition, Chapman&Hall/CRC Press.
%      [2] McCullagh, P., and J.A. Nelder (1990) Generalized Linear
%          Models, 2nd edition, Chapman&Hall/CRC Press.
%      [3] Collett, D. (2002) Modelling Binary Data, 2nd edition,
%          Chapman&Hall/CRC Press.

%   The built-in power link is intended for cases where the response data y,
%   and the linear predictor x*b, are positive, and so the link function
%   is calculated as eta = max(mu^p,delta1), where delta1 is a small positive
%   value.  Similarly, the inverse link is mu = max(eta^(1/p),delta2).  It is
%   also possible to define a custom link as
%
%      FL = @(x) sign(x).*abs(x).^p;
%      FD = @(x) p.*abs(x).^(1-p);
%      FI = @(x) sign(x).*abs(x).^(1/p);
%      link = {FL FD FI};
%
%   which may be useful in cases where the data are not positive.

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.3.2.10 $  $Date: 2009/05/07 18:31:08 $

if nargin < 2
    error('stats:glmfit:TooFewInputs','At least two arguments are required');
end

if nargin < 3 || isempty(distr), distr = 'normal'; end

% Determine the syntax.
if nargin < 4
    newSyntax = false;
else
    arg = varargin{1};
    if ischar(arg) % either a link name (old syntax), or a parameter name
        newSyntax = isempty(strmatch(lower(arg), ...
            {'identity' 'log' 'logit' 'probit' 'comploglog' 'reciprocal' 'logloglink'},'exact'));
    else % power link exponent, or custom link, but not a parameter name
        newSyntax = false;
    end
end

% Process optional name/value pairs.
if newSyntax
    paramNames = {     'link' 'estdisp' 'offset' 'weights' 'constant','init_fun','iterlim'};
    paramDflts = {'canonical'     'off'      []        []        'on',@startingVals,100};
    [errid,errmsg,link,estdisp,offset,pwts,const,startingVals_p,iterLim] = ...
                           ... internal.stats.getargs(paramNames, paramDflts, varargin{:});
                           internal_stats_getargs(paramNames, paramDflts, varargin{:});
    if ~isempty(errid)
        error(sprintf('stats:glmfit:%s',errid),errmsg);
    end

else % the old syntax glmfit(x,y,distr,link,estdisp,offset,pwts,const)
    link = 'canonical';
    estdisp = 'off';
    offset = [];
    pwts = [];
    const = 'on';
    if nargin > 3 && ~isempty(varargin{1}), link = varargin{1}; end
    if nargin > 4 && ~isempty(varargin{2}), estdisp = varargin{2}; end
    if nargin > 5 && ~isempty(varargin{3}), offset = varargin{3}; end
    if nargin > 6 && ~isempty(varargin{4}), pwts = varargin{4}; end
    if nargin > 7 && ~isempty(varargin{5}), const = varargin{5}; end
end

% Set distribution-specific defaults.
N = []; % needed only for binomial
switch distr
case 'normal'
    sqrtvarFun = @(mu) ones(size(mu));
    devFun = @(mu,y) (y - mu).^2;
    if isequal(link, 'canonical'), link  = 'identity'; end
    estdisp = 'on';
case 'binomial'
    if size(y,2) == 1
        % N will get set to 1 below
        if any(y < 0 | y > 1)
            error('stats:glmfit:BadData', ...
                  'For the binomial distribution, Y must be a binary vector or\na matrix with two columns with the number of trials in the second column.');
        end
    elseif size(y,2) == 2
        y(y(:,2)==0,2) = NaN;
        N = y(:,2);
        y = y(:,1) ./ N;
        if any(y < 0 | y > 1)
            error('stats:glmfit:BadData', 'Y must contain values in the interval [0,N] for the binomial distribution.');
        end
    else
        error('stats:glmfit:MatrixOrBernoulliRequired',...
              'Y must be a two column matrix or a vector for the binomial distribution.');
    end
    % Wait until N has NaNs removed to define variance function p*(1-p)/N and
    % the deviance function 2*(y*log(y/mu) + (N-y)*log((N-y)/(N-mu))).
    if isequal(link, 'canonical'), link = 'logit'; end
case 'poisson'
    if any(y < 0)
        error('stats:glmfit:BadData', 'Y must contain non-negative values for the Poisson distribution.');
    end
    sqrtvarFun = @(mu) sqrt(mu);
    devFun = @(mu,y) 2*(y .* (log((y+(y==0)) ./ mu)) - (y - mu));
    if isequal(link, 'canonical'), link = 'log'; end
case 'gamma'
    if any(y <= 0)
        error('stats:glmfit:BadData', 'Y must contain positive values for the gamma distribution.');
    end
    sqrtvarFun = @(mu) mu;
    devFun = @(mu,y) 2*(-log(y ./ mu) + (y - mu) ./ mu);
    if isequal(link, 'canonical'), link = 'reciprocal'; end
    estdisp = 'on';
case 'inverse gaussian'
    if any(y <= 0)
        error('stats:glmfit:BadData', 'Y must contain positive values for the inverse Gaussian distribution.');
    end
    sqrtvarFun = @(mu) mu.^(3/2);
    devFun = @(mu,y) ((y - mu)./mu).^2 ./  y;
    if isequal(link, 'canonical'), link = -2; end
    estdisp = 'on';
otherwise
    error('stats:glmfit:BadDistribution', 'Distribution name is invalid.');
end

% Remove missing values from the data.  Also turns row vectors into columns.
[anybad,wasnan,y,x,offset,pwts,N] = statremovenan(y,x,offset,pwts,N);
if anybad > 0
    switch anybad
    case 2
        error('stats:glmfit:InputSizeMismatch', 'Number of observations in X and Y must match.')
    case 3
        error('stats:glmfit:InputSizeMismatch', 'Lengths of OFFSET and Y must match.')
    case 4
        error('stats:glmfit:InputSizeMismatch', 'Lengths of PWTS and Y must match.')
%   case 5
        % N is empty, or was created from y (so its length must match)
    end
end

if isequal(const,'on')
    x = [ones(size(y,1),1) x];
end
dataClass = superiorfloat(x,y);
x = cast(x,dataClass);
y = cast(y,dataClass);

% If x is rank deficient (perhaps because it is overparameterized), we will
% warn and remove columns, and the corresponding coefficients and std. errs.
% will be forced to zero.
[n,ncolx] = size(x);
if isempty(pwts)
    [Q,R,perm] = qr(x,0);
else
    [Q,R,perm] = qr(x .* pwts(:,ones(1,ncolx)),0);
end
rankx = sum(abs(diag(R)) > abs(R(1))*max(n,ncolx)*eps(class(R)));
if rankx < ncolx
    warning('stats:glmfit:IllConditioned', ...
            ['X is ill conditioned, or the model is overparameterized, and\n' ...
             'some coefficients are not identifiable.  You should use caution\n' ...
             'in making predictions.']);
    perm = perm(1:rankx);
    x = x(:,perm);
else
    perm = 1:ncolx;
end

% Number of observations after removing missing data, number of coeffs after
% removing dependent cols and (possibly) adding a constant term.
[n,p] = size(x);

if isempty(pwts)
    pwts = 1;
elseif any(pwts == 0)
    % A zero weight means ignore the observation, so n is reduced by one.
    % Residuals will be computed, however.
    n = n - sum(pwts == 0);
end
if isempty(offset), offset = 0; end
if isempty(N), N = 1; end

% Define variance and deviance for binomial, now that N has NaNs removed.
if isequal(distr, 'binomial')
    sqrtN = sqrt(N);
    sqrtvarFun = @(mu) sqrt(mu).*sqrt(1-mu) ./ sqrtN;
    devFun = @(mu,y) 2*N.*(y.*log((y+(y==0))./mu) + (1-y).*log((1-y+(y==1))./(1-mu)));
end

% Instantiate functions for one of the canned links, or validate a
% user-defined link specification.
[emsg,linkFun,dlinkFun,ilinkFun] = stattestlink(link,dataClass);
if ~isempty(emsg)
    error('stats:glmfit:BadLink',emsg);
end

% Initialize mu and eta from y.
mu = startingVals_p(distr,y,N);
eta = linkFun(mu);

% Set up for iterations
iter = 0;
%iterLim = 1000;  % replaced by optional argument iterLim
warned = false;
seps = sqrt(eps);
convcrit = 1e-6;
b = zeros(p,1,dataClass);

% Enforce limits on mu to guard against an inverse link that doesn't map into
% the support of the distribution.
switch distr
case 'binomial'
    % mu is a probability, so order one is the natural scale, and eps is a
    % reasonable lower limit on that scale (plus it's symmetric).
    muLims = [eps(dataClass) 1-eps(dataClass)];
case {'poisson' 'gamma' 'inverse gaussian'}
    % Here we don't know the natural scale for mu, so make the lower limit
    % small.  This choice keeps mu^4 from underflowing.  No upper limit.
    muLims = realmin(dataClass).^.25;
end

while iter <= iterLim
    iter = iter+1;

    % Compute adjusted dependent variable for least squares fit
    deta = dlinkFun(mu);
    z = eta + (y - mu) .* deta;

    % Compute IRLS weights the inverse of the variance function
    sqrtw = sqrt(pwts) ./ (abs(deta) .* sqrtvarFun(mu));

    % If the weights have an enormous range, we won't be able to do IRLS very
    % well.  The prior weights may be bad, or the fitted mu's may have too
    % wide a range, which is probably because the data do as well, or because
    % the link function is trying to go outside the distribution's support.
    wtol = max(sqrtw)*eps(dataClass)^(2/3);
    t = (sqrtw < wtol);
    if any(t)
        t = t & (sqrtw ~= 0);
        if any(t)
            sqrtw(t) = wtol;
            if ~warned
                warning('stats:glmfit:BadScaling', ...
                       ['Weights are ill-conditioned.   Data may be badly ' ...
                        'scaled, or\nthe link function may be inappropriate.']);
            end
            warned = true;
        end
    end

    % Compute coefficient estimates for this iteration - the IRLS step
    b_old = b;
    [b,R] = wfit(z - offset, x, sqrtw);

    % Form current linear predictor, including offset
    eta = offset + x * b;

    % Compute predicted mean using inverse link function
    mu = ilinkFun(eta);

    % Force mean in bounds, in case the link function is a wacky choice
    switch distr
    case 'binomial'
        if any(mu < muLims(1) | muLims(2) < mu)
        mu = max(min(mu,muLims(2)),muLims(1));
        end
    case {'poisson' 'gamma' 'inverse gaussian'}
        if any(mu < muLims(1))
        mu = max(mu,muLims(1));
        end
    end

    % Check stopping conditions
    if (~any(abs(b-b_old) > convcrit * max(seps, abs(b_old)))), break; end
end
if iter > iterLim
    warning('stats:glmfit:IterationLimit','Iteration limit reached.');
end

bb = zeros(ncolx,1,dataClass); bb(perm) = b;

if nargout > 1
    % Sum components of deviance to get the total deviance.
    di = devFun(mu,y);
    dev = sum(pwts .* di);
end

% Return additional statistics if requested
if nargout > 2
    % Compute the sum of squares used to estimate dispersion, and the
    % Anscombe residuals.
    switch(distr)
    case 'normal'
        ssr = sum(pwts .* (y - mu).^2);
        anscresid = y - mu;
    case 'binomial'
        ssr = sum(pwts .* (y - mu).^2 ./ (mu .* (1 - mu) ./ N));
        t = 2/3;
        anscresid = beta(t,t) * ...
            (betainc(y,t,t)-betainc(mu,t,t)) ./ ((mu.*(1-mu)).^(1/6) ./ sqrt(N));
    case 'poisson'
        ssr = sum(pwts .* (y - mu).^2 ./ mu);
        anscresid = 1.5 * ((y.^(2/3) - mu.^(2/3)) ./ mu.^(1/6));
    case 'gamma'
        ssr = sum(pwts .* ((y - mu) ./ mu).^2);
        anscresid = 3 * (y.^(1/3) - mu.^(1/3)) ./ mu.^(1/3);
    case 'inverse gaussian'
        ssr = sum(pwts .* ((y - mu) ./ mu.^(3/2)).^2);
        anscresid = (log(y) - log(mu)) ./ mu;
    end

    % Compute residuals, using original count scale for binomial
    if (isequal(distr, 'binomial'))
        resid = (y - mu) .* N;
    else
        resid  = y - mu;
    end

    dfe = max(n - p, 0);
    stats.beta = bb;
    stats.dfe = dfe;
    if dfe > 0
        stats.sfit = sqrt(ssr / dfe);
    else
        stats.sfit = NaN;
    end
    if isequal(estdisp,'off')
        stats.s = 1;
        stats.estdisp = false;
    else
        stats.s = stats.sfit;
        stats.estdisp = true;
    end

    % Find coefficient standard errors and correlations
    if ~isnan(stats.s) % dfe > 0 or estdisp == 'off'
        RI = R\eye(p);
        C = RI * RI';
        if isequal(estdisp,'on'), C = C * stats.s^2; end
        se = sqrt(diag(C)); se = se(:);   % insure vector even if empty
        stats.covb = zeros(ncolx,ncolx,dataClass);
        stats.covb(perm,perm) = C;
        C = C ./ (se * se');
        stats.se = zeros(ncolx,1,dataClass); stats.se(perm) = se;
        stats.coeffcorr = zeros(ncolx,ncolx,dataClass);
        stats.coeffcorr(perm,perm) = C;
        stats.t = NaN(ncolx,1,dataClass); stats.t(perm) = b ./ se;
        if isequal(estdisp,'on')
            stats.p = 2 * tcdf(-abs(stats.t), dfe);
        else
            stats.p = 2 * normcdf_TE(-abs(stats.t));
        end
    else
        stats.se = NaN(size(bb),class(bb));
        stats.coeffcorr = NaN(length(bb),class(bb));
        stats.t = NaN(size(bb),class(bb));
        stats.p = NaN(size(bb),class(bb));
    end

    stats.resid  = statinsertnan(wasnan, resid);
    stats.residp = statinsertnan(wasnan, (y - mu) ./ (sqrtvarFun(mu) + (y==mu)));
    stats.residd = statinsertnan(wasnan, sign(y - mu) .* sqrt(max(0,di)));
    stats.resida = statinsertnan(wasnan, anscresid);
end
end

function [b,R] = wfit(y,x,sw)
% Perform a weighted least squares fit
[n,p] = size(x);
yw = y .* sw;
xw = x .* sw(:,ones(1,p));
% No pivoting, no basic solution.  We've removed dependent cols from x, and
% checked the weights, so xw should be full rank.
[Q,R] = qr(xw,0);
b = R \ (Q'*yw);
end


function mu = startingVals(distr,y,N)
% Find a starting value for the mean, avoiding boundary values
switch distr
case 'poisson'
    mu = y + 0.25;
case 'binomial'
    mu = (N .* y + 0.5) ./ (N + 1);
    mu= (mu-0.5)*0.8+0.5;
case {'gamma' 'inverse gaussian'}
    mu = max(y, eps(class(y))); % somewhat arbitrary
otherwise
    mu = y;
end
end


function [badin,wasnan,varargout]=statremovenan(varargin)
%STATREMOVENAN Remove NaN values from inputs

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 1.3.2.2 $  $Date: 2005/05/31 16:45:15 $

badin = 0;
wasnan = 0;
n = -1;

% Find NaN, check length, and store outputs temporarily
varargout = cell(nargout,1);
for j=1:nargin
   y = varargin{j};
   if (size(y,1)==1) && (n~=1) 
       y =  y';
   end

   ny = size(y,1);
   if (n==-1)
      n = ny;
   elseif (n~=ny && ny~=0)
      if (badin==0), badin = j; end
   end
   
   varargout{j} = y;

   if (badin==0 && ny>0)
       wasnan = wasnan | any(isnan(y),2);
   end
end

if (badin>0), return; end

% Fix outputs
if (any(wasnan))
   t = ~wasnan;
   for j=1:nargin
      y = varargout{j};
      if (length(y)>0), varargout{j} = y(t,:); end
   end
end
end


function [emsg,link,dlink,ilink] = stattestlink(linkArg,dataClass)
%STATTESTLINK Test link function for GLMFIT and GLMVAL

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.5.4.3 $  $Date: 2009/05/07 18:34:04 $

emsg = '';
link = []; dlink = []; ilink = [];

% Recognize the power link.
if isfloat(linkArg) && isscalar(linkArg)
    if linkArg == 0 % equivalent to the log link
        linkArg = 'log';
    else
        linkExponent = linkArg;
        linkArg = 'power';
    end
end

% Some inverse links are defined with limits on eta so that we always get a
% [finite|positive|(0,1)] value for mu, and so we can then evaluate the link
% and its derivative with such a mu and get a finite eta.  The link functions
% themselves all have the corresponding property intrinsically, except power
% with abs(exponent) > 1.
%
% The links that map from [0,1] have order one as the natural scale, and eps is
% a reasonable lower limit on that scale (plus it's symmetric).  We don't know
% the natural scale for the other cases, so make the limits wide.
tiny = realmin(dataClass)^.25; % keep fourth powers from under/overflowing
if ischar(linkArg)
    switch linkArg
    case 'identity'
        link = @(mu) mu;
        dlink = @(mu) ones(size(mu));
        ilink = @(eta) eta;
    case 'log'
        link = @(mu) log(mu);
        dlink = @(mu) 1 ./ mu;
        % keep mu = ilink(eta) in [tiny, 1/tiny];
        lowerBnd = log(tiny); upperBnd = -lowerBnd;
        ilink = @(eta) exp(constrain(eta,lowerBnd,upperBnd));
    case 'logit'
        link = @(mu) log(mu ./ (1-mu));
        dlink = @(mu) 1 ./ (mu .* (1-mu));
        % keep mu = ilink(eta) in (approx) [eps, (1-eps)];
        lowerBnd = log(eps(dataClass)); upperBnd = -lowerBnd;
        ilink = @(eta) 1 ./ (1 + exp(-constrain(eta,lowerBnd,upperBnd)));
    case 'probit'
        link = @(mu) norminv(mu);
        dlink = @(mu) 1 ./ normpdf(norminv(mu));
        % keep mu = ilink(eta) in [eps, (1-eps)];
        lowerBnd = norminv(eps(dataClass)); upperBnd = -lowerBnd;
        ilink = @(eta) normcdf(constrain(eta,lowerBnd,upperBnd));
    case 'comploglog'
        link = @(mu) log(-log1p(-mu));
        dlink = @(mu) 1 ./ -((1-mu) .* log1p(-mu));
        % keep mu = ilink(eta) in [eps, (1-eps)];
        lowerBnd = log(-log1p(-eps(dataClass))); upperBnd = log(-log(eps(dataClass)));
        ilink = @(eta) -expm1(-exp(constrain(eta,lowerBnd,upperBnd)));
    case {'loglog', 'logloglink'}
        link = @(mu) log(-log(mu));
        dlink = @(mu)  1 ./ (mu .* log(mu));
        % keep mu = ilink(eta) in [eps, (1-eps)];
        lowerBnd = log(-log1p(-eps(dataClass))); upperBnd = log(-log(eps(dataClass)));
        ilink = @(eta) exp(-exp(constrain(eta,lowerBnd,upperBnd)));
    case 'reciprocal'
        link = @(mu) 1 ./ mu;
        dlink = @(mu) -1 ./ mu.^2;
        % keep mu = ilink(eta) in [tiny, 1/tiny];
        lowerBnd = tiny; upperBnd = 1/lowerBnd;
        ilink = @(eta) 1 ./ constrain(eta,lowerBnd,upperBnd);
    case 'power' % linkExponent==0 (equivalent to 'log') has been weeded out already
        % keep eta = link(mu) in [tiny, 1/tiny];
        lowerBnd = tiny.^min(abs(1/linkExponent),1); upperBnd = 1/lowerBnd;
        link = @(mu) constrain(mu,lowerBnd,upperBnd).^linkExponent;
        dlink = @(mu) linkExponent * mu.^(linkExponent-1);
        % keep mu = ilink(eta) in [tiny, 1/tiny];
        lowerBnd = tiny.^min(abs(linkExponent),1); upperBnd = 1/lowerBnd;
        ilink = @(eta) constrain(eta,lowerBnd,upperBnd) .^ (1/linkExponent);
    otherwise
        emsg = 'Unrecognized link function name.';
    end

elseif iscell(linkArg)
    % A cell array of three functions is okay
    if numel(linkArg) ~= 3
        emsg = 'LINK cell array must have three components';
        return
    end

    link = linkArg{1};
    if ischar(link) && ~isempty(which(link))
        name = link; link = @(mu) feval(name,mu);
    elseif ~isa(link,'function_handle') && ~isa(link,'inline')
        emsg = 'LINK function is not valid';
        return
    end

    dlink = linkArg{2};
    if ischar(dlink) && ~isempty(which(dlink))
        name = dlink; dlink = @(mu) feval(name,mu);
    elseif ~isa(dlink,'function_handle') && ~isa(dlink,'inline')
        emsg = 'LINK function derivative is not valid';
        return
    end

    ilink = linkArg{3};
    if ischar(ilink) && ~isempty(which(ilink))
        name = ilink; ilink = @(eta) feval(name,eta);
    elseif ~isa(ilink,'function_handle') && ~isa(ilink,'inline')
        emsg = 'LINK function inverse is not valid';
        return
    end

else
    emsg = 'LINK argument is not valid.';
end
end

%-----------------------------------------------------------------------
function x = constrain(x,lower,upper)
% Constrain between upper and lower limits, and do not ignore NaN
x(x<lower) = lower;
x(x>upper) = upper;
end

function [varargout]=statinsertnan(wasnan,varargin)
%STATINSERTNAN Insert NaN, space, '' or undefined value into inputs.
%   X1 = STATINSERTNAN(WASNAN, Y1) inserts missing values in Y1 and returns
%   it as X1. WASNAN is a logical column vector and the output of
%   STATREMOVENAN. Its TRUE values indicate the rows in X1 that will
%   contain missing values. Y1 is a column vector or a matrix. The type of
%   Y1 can be:
%   Categorical       - X1 is categorical, undefined values represents
%                       missing values.
%   Double            - X1 is double. NaN values represents missing values.
%   Single            - X1 is single. NaN values represents missing values.
%   Character matrix  - X1 is a character matrix. Space represents missing
%                       values.
%   Cell              - X1 is a cell array. empty string '' represents
%                       missing values.
%
%  [X1,X2,...] = STATINSERTNAN(WASNAN,Y1,Y2,...) accepts any number of
%  input variables Y1,Y2,Y3,.... STATINSERTNAN inserts missing values in
%  Y1, Y2,...  and returns them as X1, X2,... respectively.
%
%  This utility is used by some Statistics Toolbox functions to handle
%  missing values.
%
%  See also STATREMOVENAN.


%   Copyright 1993-2008 The MathWorks, Inc.
%   $Revision: 1.4.2.3 $  $Date: 2008/10/31 07:41:07 $

if ~any(wasnan)
     varargout = varargin;
     return;
end

ok = ~wasnan;
len = length(wasnan);
for j=1:nargin-1
    y = varargin{j};
    if (size(y,1)==1) && sum(ok) > 1
        y =  y';
    end
    
    [n,p] = size(y);
    
    if ischar(y)
        x = repmat(' ', [len,p]);
    elseif isa(y, 'nominal')
        x = nominal(NaN([len,p]));
        x = addlevels(x,getlabels(y));
    elseif isa(y, 'ordinal')
        x = ordinal(NaN([len,p]));
        x = addlevels(x,getlabels(y));
    elseif iscell(y)
        x = repmat({''},[len,p]);
    elseif isfloat(y)
            x = nan([len,p],class(y));
    elseif islogical(y)
        error('stats:statinsertnan:InputTypeIncorrect',...
            ['Logical input is not allowed because it can''t '...
            'present missing values. Use CATEGORICAL variable instead.']);
    else
        error('stats:statinsertnan:InputTypeIncorrect',...
            ['Y must be categorical, double, single '...
            ' cell array; or a 2D character array.']);
    end
    
    x(ok,:) = y;
    
    varargout{j} = x;
end
end

function [eid,emsg,varargout]=internal_stats_getargs(pnames,dflts,varargin)
%GETARGS Process parameter name/value pairs for statistics functions
%   [EID,EMSG,A,B,...]=GETARGS(PNAMES,DFLTS,'NAME1',VAL1,'NAME2',VAL2,...)
%   accepts a cell array PNAMES of valid parameter names, a cell array
%   DFLTS of default values for the parameters named in PNAMES, and
%   additional parameter name/value pairs.  Returns parameter values A,B,...
%   in the same order as the names in PNAMES.  Outputs corresponding to
%   entries in PNAMES that are not specified in the name/value pairs are
%   set to the corresponding value from DFLTS.  If nargout is equal to
%   length(PNAMES)+1, then unrecognized name/value pairs are an error.  If
%   nargout is equal to length(PNAMES)+2, then all unrecognized name/value
%   pairs are returned in a single cell array following any other outputs.
%
%   EID and EMSG are empty if the arguments are valid.  If an error occurs,
%   EMSG is the text of an error message and EID is the final component
%   of an error message id.  GETARGS does not actually throw any errors,
%   but rather returns EID and EMSG so that the caller may throw the error.
%   Outputs will be partially processed after an error occurs.
%
%   This utility is used by some Statistics Toolbox functions to process
%   name/value pair arguments.
%
%   Example:
%       pnames = {'color' 'linestyle', 'linewidth'}
%       dflts  = {    'r'         '_'          '1'}
%       varargin = {{'linew' 2 'nonesuch' [1 2 3] 'linestyle' ':'}
%       [eid,emsg,c,ls,lw] = internal.stats.getargs(pnames,dflts,varargin{:})    % error
%       [eid,emsg,c,ls,lw,ur] = internal.stats.getargs(pnames,dflts,varargin{:}) % ok

%   Copyright 1993-2009 The MathWorks, Inc. 
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:18:35 $ 

% We always create (nparams+2) outputs:
%    one each for emsg and eid
%    nparams varargs for values corresponding to names in pnames
% If they ask for one more (nargout == nparams+3), it's for unrecognized
% names/values

% Initialize some variables
emsg = '';
eid = '';
nparams = length(pnames);
varargout = dflts;
unrecog = {};
nargs = length(varargin);

% Must have name/value pairs
if mod(nargs,2)~=0
    eid = 'WrongNumberArgs';
    emsg = 'Wrong number of arguments.';
else
    % Process name/value pairs
    for j=1:2:nargs
        pname = varargin{j};
        if ~ischar(pname)
            eid = 'BadParamName';
            emsg = 'Parameter name must be text.';
            break;
        end
        i = strmatch(lower(pname),pnames);
        if isempty(i)
            % if they've asked to get back unrecognized names/values, add this
            % one to the list
            if nargout > nparams+2
                unrecog((end+1):(end+2)) = {varargin{j} varargin{j+1}};
                
                % otherwise, it's an error
            else
                eid = 'BadParamName';
                emsg = sprintf('Invalid parameter name:  %s.',pname);
                break;
            end
        elseif length(i)>1
            i = strmatch(lower(pname),pnames,'exact');
            if length(i)~=1
                eid = 'BadParamName';
                emsg = sprintf('Ambiguous parameter name:  %s.',pname);
                break;
            else
                varargout{i} = varargin{j+1};
            end
        else
            varargout{i} = varargin{j+1};
        end
    end
end

varargout{nparams+1} = unrecog;
end


function [p,plo,pup] = normcdf_TE(x,mu,sigma,pcov,alpha)
%NORMCDF Normal cumulative distribution function (cdf).
%   P = NORMCDF(X,MU,SIGMA) returns the cdf of the normal distribution with
%   mean MU and standard deviation SIGMA, evaluated at the values in X.
%   The size of P is the common size of X, MU and SIGMA.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [P,PLO,PUP] = NORMCDF(X,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%   for P when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  PLO and PUP are arrays of the same size as P containing the lower
%   and upper confidence bounds.
%
%   See also ERF, ERFC, NORMFIT, NORMINV, NORMLIKE, NORMPDF, NORMRND, NORMSTAT.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:26:18 $

if nargin<1
    error(message('stats:normcdf:TooFewInputsX'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>1
   if nargin<4
      error(message('stats:normcdf:TooFewInputsCovariance'));
   end
   if ~isequal(size(pcov),[2 2])
      error(message('stats:normcdf:BadCovarianceSize'));
   end
   if nargin<5
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:normcdf:BadAlpha'));
   end
end

try
    z = (x-mu) ./ sigma;
catch
    error(message('stats:normcdf:InputSizeMismatch'));
end

% Prepare output
p = NaN(size(z),class(z));
if nargout>=2
    plo = NaN(size(z),class(z));
    pup = NaN(size(z),class(z));
end

% Set edge case sigma=0
p(sigma==0 & x<mu) = 0;
p(sigma==0 & x>=mu) = 1;
if nargout>=2
    plo(sigma==0 & x<mu) = 0;
    plo(sigma==0 & x>=mu) = 1;
    pup(sigma==0 & x<mu) = 0;
    pup(sigma==0 & x>=mu) = 1;
end

% Normal cases
if isscalar(sigma)
    if sigma>0
        todo = true(size(z));
    else
        return;
    end
else
    todo = sigma>0;
end
z = z(todo);

% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
p(todo) = 0.5 * erfc(-z ./ sqrt(2));

% Compute confidence bounds if requested.
if nargout>=2
   zvar = (pcov(1,1) + 2*pcov(1,2)*z + pcov(2,2)*z.^2) ./ (sigma.^2);
   if any(zvar<0)
      error(message('stats:normcdf:BadCovarianceSymPos'));
   end
   normz = -norminv(alpha/2);
   halfwidth = normz * sqrt(zvar);
   zlo = z - halfwidth;
   zup = z + halfwidth;

   plo(todo) = 0.5 * erfc(-zlo./sqrt(2));
   pup(todo) = 0.5 * erfc(-zup./sqrt(2));
end
end
