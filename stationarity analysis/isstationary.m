% +------------------------------------------------------+
% |      Signal wide-sense stationarity estimation       |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        08/22/22 | 
% +------------------------------------------------------+
% 
% function: [wss, mean_stat, var_stat, cov_stat] = isstationary(x, gamma)
%
% Input:
% x - signal in the time domain;
% gamma - confidence level (e.g., 0.9, 0.95) for the hypothesis that
%         the signal under consideration is stationary .
% 
% Output:
% wss - a Boolean flag showing whether the signal is wide-sense   
%       stationary (WSS), that is, simultaneously stationary about   
%       its mean, variance and autocovariance;
% mean_stat - a Boolean flag showing whether the signal is stationary
%             about its mean;
% var_stat - a Boolean flag showing whether the signal is stationary
%            about its variance; 
% cov_stat - a Boolean flag showing whether the signal is stationary
%            about its autocovariance. 
%
% Note: a signal (e.g., time series) is said to be weakly stationary or
% wide sense stationary (WSS) if its (time-localized) mean and variance 
% are constant over time and if its autocovariance function Cxx(t1, t2) 
% depends only on the difference t2-t1, but not on their particular values. 
% Be aware - one does not consider the underling process or the entire
% population here! The function estimates the WSS of the signal "as it is"!
function [wss, mean_stat, var_stat, cov_stat] = isstationary(x, gamma)
% input validation
validateattributes(x, {'single', 'double'}, ...
                      {'vector', 'real', 'nonnan', 'nonempty', 'finite'}, ...
                      '', 'x', 1)
validateattributes(gamma, {'single', 'double'}, ...
                          {'scalar', 'real', 'nonnan', 'nonempty', 'finite', ...
                           '>', 0, '<', 1}, ...
                           '', 'gamma', 2)
% split the signal into three equally length parts with 50% overlapping
N = length(x);
[X, ~] = buffer(x, floor(N/2), ceil(N/4), 'nodelay');   
% check for the statistical identity of the three partial signals
x1 = X(:, 1); x2 = X(:, 2); x3 = X(:, 3);
[mean_ident_12, var_ident_12, cov_ident_12] = isstatid(x1, x2, gamma);   
[mean_ident_23, var_ident_23, cov_ident_23] = isstatid(x2, x3, gamma);   
[mean_ident_13, var_ident_13, cov_ident_13] = isstatid(x1, x3, gamma);
% check for stationarity of the signal under consideration
mean_stat = mean_ident_12  && mean_ident_23  && mean_ident_13;
var_stat  =  var_ident_12  &&  var_ident_23  &&  var_ident_13;
cov_stat  =  cov_ident_12  &&  cov_ident_23  &&  cov_ident_13;
wss       =  mean_stat     &&  var_stat      &&  cov_stat;  
                      
end
function [mean_ident, var_ident, cov_ident] = isstatid(x1, x2, gamma)                                       
% set the significance level (alpha)
% Note: actually, alpha is the predefined threshold level under which the
% corresponding condition should be considered logical truth. It is a
% direct equivalent to the significance level � the probability of making a
% wrong decision to reject the null hypothesis when it is actually true (in
% statistics, it is called a false positive error or type I error).
% Particularly, the null hypothesis is that the signal under test is
% stationary. The opposite quantity is the confidence level (gamma) � the
% probability of not rejecting the null hypothesis when it is true. The
% significance level (alpha) complements the confidence level (gamma) to 1.
alpha = 1-gamma;
% check for the statistical significance of the first moment (i.e., the
% mean) using the inverse form factor (i.e., the MEAN to STD ratio) and
% then empirically test for first moment identity
m1 = mean(x1); m2 = mean(x2);
s1 = std(x1); s2 = std(x2);
if (abs(m1)/s1 <= alpha) && (abs(m2)/s2 <= alpha)                             
    % when both mean values are not statistically significant...
    mean_ident = true;   
else
    % when at least one of the mean values is statistically significant...
    mean_ident = abs((m1 - m2)/min(m1, m2)) <= alpha;                    
end
          
% empirically test for second moment (i.e., variance) identity   
v1 = var(x1); v2 = var(x2);
var_ident = abs(v1 - v2)/min(v1, v2) <= alpha;                           
% empirically test for autocovariance identity via: (i) computation of
% the autocovariances of the two signal halves; (ii) estimation of the
% similarity between them using the sample Pearson correlation coefficient
% and (iii) additional comparison of the variances of the autocovariance
% sequences of the two signal halves, since (ii) does not compare the
% scales of the two sequences but only their structure/pattern.
% Note: the p-value obtained by the corrcof function ranges from 0 to 1,
% where values close to 0 (less than alpha) correspond to a significant
% correlation and a low probability of observing the null hypothesis that
% there is no relationship between the sequences. For more information
% about the computation of the p-value see the script of the built-in
% Matlab corrcoef function (e.g., type "edit corrcoef" in the Command
% Window).
c1 = xcov(x1, 'coeff'); c2 = xcov(x2, 'coeff'); 
vc1 = var(c1); vc2 = var(c2);
[~, p] = corrcoef(c1, c2);
cov_ident = abs(vc1 - vc2)/min(vc1, vc2) < alpha && p(1, 2) <= alpha;    
                      
end