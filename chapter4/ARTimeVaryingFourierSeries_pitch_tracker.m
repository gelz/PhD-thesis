 function [T_out2,std_T_out2] = ARTimeVaryingFourierSeries_pitch_tracker(data_input,T_0,B_size,SNR_in)

if nargin<4
    SNR_in = 80;
end

if nargin<3
    B_size = 128;
end

data_sample = data_input;
SNR_desired_input = SNR_in;

noise = randn(1,length(data_sample));
prior_SNR = 10*log10((sum(data_sample.^2)./sum((noise).^2)));
k = 10^((prior_SNR-SNR_desired_input)/20);

k = 0;

y_ob = data_sample + k*noise;
y = y_ob;


% set particle filter parameters;

N_p = 1000;  %  number of particles used;
N_number = 100;  % define the maximum number of periods to be estimated;

[x_p,hyper_par] = initialize_particles_sin(y,N_p,N_number,T_0);


w0 = ones(N_p,1)./N_p;  % assign weights to each particles;

% the standard paticle filter process:
%   assume  z_k = f(x_k,v_k);  x_k denotes the hidden states, while v_k
%   denotes the watching noise;
%   and x_k  = g(x_(k-1), n_k);  this is the evolving process of states;

% model description:
% for the speech model, here we use the AR model driven by a time-varying
% fourier series; 
% x_p =  [AR_par, source_input_par];
% AR_par = [p_order, AR1, AR2, ... AR_p];
% source_input_par = [harmonics_number, Amplitude_1, Amp_2, ..., Amp_n,
% period];
% hyper_par includes the variance of period and amplitude changes;


% to use blocks, set block size = B_size
% B_size = 128;
Win = window(@hamming,B_size);
INC = 0.5;

T_out_size = floor(length(data_sample)./B_size/INC)-2;
T_out2 = zeros(T_out_size,1);
std_T_out2 = zeros(T_out_size,1);
w_old = w0;
w_new = w0;
noise_power = 0.1;

for t = 1 : (floor(length(data_sample)./B_size/INC)-2)
  
    prior = w_old;
    data_pointer = ceil(((t-1)*B_size*INC+1)):ceil(((t-1)*B_size*INC+B_size));
    new_data = Win.*data_sample(data_pointer)';
    
    for k = 1:N_p
        
        model_data = model_data_ARTimeVaryingFourierSeries_sin(x_p(k,:,:),data_pointer,data_sample,B_size);
        
%         log_likelihood(k) = sum(log10(pdf('norm',model_data-new_data',0,noise_power)));
        log_likelihood_temp = -(model_data - new_data').^2/2/(noise_power.^2)/log(10);
        log_likelihood(k) = sum(log_likelihood_temp);

        log_prior(k) = log10(prior(k));
        log_posterior(k) = log_prior(k) + log_likelihood(k);
    end
    
    w_new = 10.^(log_posterior - max(log_posterior));
    
    w_new = w_new./(sum(w_new));
    
    [T_out2(t),std_T_out2(t)] = mean_out_TimeVary_harmonics(w_new,x_p,data_pointer);   
    
    
    ess = 1./(sum(w_new.^2));
    
    if ess > N_p/2
        resample_decision = 0;
    else resample_decision = 1;
    end
    
    if resample_decision == 0
        w_old = w_new;
    else
        flags = resample([1:N_p]',w_new);
        for k = 1:N_p
            x_p(k,:,:) = x_p(flags(k),:,:);
        end
        w_old = 1./N_p.*ones(N_p,1);
    end
    
    for k = 1:N_p
        x_p(k,:,:) = update_rule1_harmonics_sin(x_p(k,:,:),data_pointer,hyper_par);
    end
%     t,
  

    
    
  
end