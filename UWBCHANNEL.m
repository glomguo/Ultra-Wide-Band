function [t,h,Tbar,Tba,RMS] = UWBCHANNEL(Lam,lambda,Gam,gamma,std_ln_1,std_ln_2,nlos,std_shdw,num_channels)

% IEEE 802.15.3a UWB channel model for PHY proposal evaluation
%      continuous-time realization of modified S-V channel model
% Input parameters:
%   Lam      Cluster arrival rate in GHz (avg # of clusters per nsec)
%   lambda   Ray arrival rate in GHz (avg # of rays per nsec)
%   Gam      Cluster decay factor (time constant, nsec)
%   gamma    Ray decay factor (time constant, nsec)
%   std_ln_1 Standard deviation of log-normal variable for cluster fading
%   std_ln_2 Standard deviation of log-normal variable for ray fading
%   nlos     Flag to specify generation of Non Line Of Sight channels
%   std_shdw  Standard deviation of log-normal shadowing of entire impulse response
%   num_channels  number of random realizations to generate
% Outputs

%   h is returned as a matrix with num_channels columns, each column 
%     holding a random realization of the channel model (an impulse response)
%   t is organized as h, but holds the time instances (in nsec) of the paths whose
%     signed amplitudes are stored in h
%   t0 is the arrival time of the first cluster for each realization
%   np is the number of paths for each realization.
% Thus, the k'th realization of the channel impulse response is the sequence 
%   of (time,value) pairs given by (t(1:np(k),k), h(1:np(k),k))

% initialize and precompute some things
std_L = 1/sqrt(2*Lam);      % std dev (nsec) of cluster arrival spacing
std_lam = 1/sqrt(2*lambda); % std dev (nsec) of ray arrival spacing
mu_const = (std_ln_1^2+std_ln_2^2)*log(10)/20;  % pre-compute for later
h_len = 100;  % there must be a better estimate of # of paths than this
ngrow = 100; % amount to grow data structure if more paths are needed
h = zeros(h_len,num_channels);
t = zeros(h_len,num_channels);
t0 = zeros(1,num_channels);
np = zeros(1,num_channels);

for k = 1:num_channels       % loop over number of channels
  tmp_h = zeros(size(h,1),1);
  tmp_t = zeros(size(h,1),1);
  if nlos,
    Tc = (std_L*randn)^2 + (std_L*randn)^2;  % First cluster random arrival
  else
    Tc = 0;                   % First cluster arrival occurs at time 0
  end
  t0(k) = Tc;
  
  %%%%%%%%%%%%%%%%%%%  The following line in rev. 1 was in error
  % ln_xi = std_ln_1*randn;   % set cluster fading
  %%%%%%%%%%%%%%%%%%%
  
  path_ix = 0;
  while (Tc < 10*Gam)
    % Determine Ray arrivals for each cluster
    Tr = 0;  % first ray arrival defined to be time 0 relative to cluster
    ln_xi = std_ln_1*randn;   % set cluster fading (new line added in rev. 1)
    while (Tr < 10*gamma)
      t_val = (Tc+Tr);  % time of arrival of this ray
      mu = (-10*Tc/Gam-10*Tr/gamma)/log(10) - mu_const;
      ln_beta = mu + std_ln_2*randn;
      pk = 2*round(rand)-1;
      h_val = pk * 10^((ln_xi+ln_beta)/20);  % signed amplitude of this ray
      path_ix = path_ix + 1;  % row index of this ray
      if path_ix > h_len,
        % grow the output structures to handle more paths as needed
%        fprintf(1,'Growing CIR length from %d paths to %d\n', length(tmp_h)+[0 ngrow]);
        tmp_h = [tmp_h; zeros(ngrow,1)];
        tmp_t = [tmp_t; zeros(ngrow,1)];
        h = [h; zeros(ngrow,num_channels)];
        t = [t; zeros(ngrow,num_channels)];
        h_len = h_len + ngrow;
      end
      tmp_h(path_ix) = h_val;
      tmp_t(path_ix) = t_val;
      Tr = Tr + (std_lam*randn)^2 + (std_lam*randn)^2;
    end
    Tc = Tc + (std_L*randn)^2 + (std_L*randn)^2;
  end
  np(k) = path_ix;  % number of rays (or paths) for this realization
  [sort_tmp_t,sort_ix] = sort(tmp_t(1:np(k)));    % sort in ascending time order
  t(1:np(k),k) = sort_tmp_t;
  h(1:np(k),k) = tmp_h(sort_ix(1:np(k)));
  % now impose a log-normal shadowing on this realization
  fac = 10^(std_shdw*randn/20) / sqrt( h(1:np(k),k)' * h(1:np(k),k) );
  h(1:np(k),k) = h(1:np(k),k) * fac;
end

% plot(t,h);
%plot(h);
% MEAN EXCESS DELAY CALCULATION.

[a b] = size(h);
for ii = 1:1:b

Tbar(ii) = sum(h((1:100),ii).^2.*t((1:100),ii))/(sum(h((1:100),ii).^2));

Tba(ii) = sum(h((1:100),ii).^2.*t((1:100),ii).^2)/(sum(h((1:100),ii).^2));

RMS(ii) = sqrt(Tba(ii) - (Tbar(ii)^2));

end


%{
Points_Gain = 20;
h_channel = h.';
h_channel_normalized = zeros(1,Points_Gain*h_len);
for i = 1:h_len
    for j = 1:Points_Gain
        h_channel_normalized(1,Points_Gain*(i-1)+1:1:Points_Gain*(i-1)+j) = h_channel(1,i);
    end
end

h= h_channel_normalized;
%}

end

