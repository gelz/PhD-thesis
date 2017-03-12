function [x_p,hyper_par] = initialize_particles_sin(data_input,N_p,N_number,T0)


% define x_p structure;

% model discription:  z = ARmodel (conv) source_input;

% source is modelled as a time-varying fourier series; 
% source = sum(Amplitude*cos(k*phase);  phase = phase + period;


% x_p =  [AR_par, source_input_par];

% AR_par = [AR1, AR2, ... AR_p];
% source_input_par = [harmonics_number, Amplitude_1, Amp_2, ..., Amp_n,
% period];

% hyper_par includes the variance of period and amplitude changes;

% define the AR parameter:
fs = 16000;
p_order = 11;
K_order = 40;   % define the number of harmonic components used;
if nargin <4
    T0 = 188;
end
% define the first guess of AR coefficients and source parameters, using
% algorithm in 'Joint source-filter optimization for robust glottal source estimation
% in the presence of shimmer and jitter',Prasanta Kumar Ghosh *, Shrikanth
% S. Narayanan, Speech Communication, 2011
% [data_input,fs] = data_read();

x_p_est = joint_est_para_sin(p_order,K_order,data_input,fs,T0);   % use the first 3 periods to esimate the initial parameters;

% x_p_est_plus = 0.05*randn(N_p,length(x_p_est));
% x_p = (x_p_est_plus+1)* x_p_est;
% x_p_est = 0.1*randn(p_order+K_order,1);
x_p = zeros(N_p,N_number,length(x_p_est)+1);

for i = 1 : N_p
    x_p(i,1,1:(end-1)) = (0.005*randn(1,1,length(x_p_est))+1).*reshape(x_p_est,1,1,length(x_p_est));    %%% this time complete here ( 08/04/2013  16:49)
    x_p(i,1,end) = T0;
    x_p(i,2,1:(end-1)) = (0.005*randn(1,1,length(x_p_est))+1).*reshape(x_p_est,1,1,length(x_p_est));    %%% this time complete here ( 08/04/2013  16:49)
    x_p(i,2,end) = T0;
    x_p(i,3,1:(end-1)) = (0.005*randn(1,1,length(x_p_est))+1).*reshape(x_p_est,1,1,length(x_p_est));    %%% this time complete here ( 08/04/2013  16:49)
    x_p(i,3,end) = T0;
    x_p(i,4,1:(end-1)) = (0.005*randn(1,1,length(x_p_est))+1).*reshape(x_p_est,1,1,length(x_p_est));    %%% this time complete here ( 08/04/2013  16:49)
    x_p(i,4,end) = T0;        
end
var_period =10;
var_amplitude = 0.02;
var_amplitude_ar = 0.02;

hyper_par = [var_period, var_amplitude,var_amplitude_ar];













