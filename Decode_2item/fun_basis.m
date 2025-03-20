function [c, chan_th] = fun_basis(s,n_chans)
% update to fun_basis.m used by van Bergen et al to match basis set used in
% Serences/Sprague/Brouwer papers, where exponent scales w/ # of channels
% to ensure filter 'steerability', which can avoid biases in decoding
% performance
%
% TCS 8/1/2018

if nargin < 2
    n_chans = 8;
end

chan_th = linspace(0,2*pi-(2*pi/n_chans),n_chans);

%TuningCentres = 0:2*pi/8:2*pi-0.001; 
assert(size(s,2)==1, 's must be a scalar or column vector');        

% utility function to compute distance between two angles
ang_dist = @(a,b) min(mod(a-b,2*pi),mod(b-a,2*pi));
chan_size = pi;

chan_pow = n_chans-mod(n_chans,2);

c = (cos( pi * ang_dist(s,chan_th) ./ (2*chan_size) ).^chan_pow) .* (ang_dist(s,chan_th)<=chan_size) ;


 %   c = max(0, cos(bsxfun(@minus, s, TuningCentres))).^5;        
end