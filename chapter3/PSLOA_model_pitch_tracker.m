function [T_out2,std_T_out2] = PSLOA_model_pitch_tracker(data_input,T_0,clean_data,B_size,SNR_in)



if nargin<5
    SNR_in = 80;
end
if nargin<4
    B_size = 64;
end

data_sample = data_input;
SNR_desired_input = SNR_in;

noise = randn(1,length(data_sample));
prior_SNR = 10*log10((sum(data_sample.^2)./sum((noise).^2)));
k = 10^((prior_SNR-SNR_desired_input)/20);

k = 0;

y_ob = data_sample + k*noise;
y = y_ob;



N_number = 100;  % define the maximum number of states;
N_p = 1000;     % define the number of particles;

% seita = [seita1, seita2];  seita1 = Period move range; seita2 = amplitude
% move range;

% B_size = 64;

% to initialize the parameters;

T_out_size = floor(length(data_sample)./B_size);
T_out = zeros(N_number,1);
T_out2 = zeros(T_out_size,1);
std_T_out2 = zeros(T_out_size,1);

% T_0 = 80;
amplitude_0 = 0.99;
T_variance = 10;
amplitude_variance = 0.15;
noise_power_2 = 0.001;

x_0 = [T_0,amplitude_0];
seita_0 = [T_variance, amplitude_variance,noise_power_2];
s0 = clean_data(1:3*T_0);
noise_power = 0.1;
D_store = zeros(N_p,length(data_input)+300);
% set the proposal function: q(x_k1|x_k1-1) = x_k1-1 +
% uniform(-seita_01,seita_01); q(x_k2|x_k2-1) = x_k2-1 +
% uniform(-seita_02,seita_02);

% initialize particles;
x_p = zeros(N_p,N_number,length(x_0));
w_old = 1./N_p.*ones(N_p,1);
w_new = w_old;
log_likelihood = zeros(N_p,1);
log_prior = zeros(N_p,1);
log_posterior = zeros(N_p,1);

% below is initiliazation for x_p; to make it robust to data, initialize 5
% periods;
for i = 1:5
    for k = 1:N_p
        if i==1
            x_p(k,i,1) = round(x_0(1) - seita_0(1) + 2*seita_0(1)*rand(1));
            x_p(k,i,2) = x_0(2) - seita_0(2) + 2*seita_0(2)*rand(1);
        else
            x_p(k,i,1) = round(x_p(k,i-1,1) - seita_0(1) + 2*seita_0(1)*rand(1));
            if (x_p(k,i,1)<40)||(x_p(k,i,1)>200)
                x_p(k,i,1) = x_p(k,i-1,1);
            end
            x_p(k,i,2) = x_p(k,i-1,2) - seita_0(2) + 2*seita_0(2)*rand(1);
            if (x_p(k,i,2)<0.4)||(x_p(k,i,2)>1.5)
                x_p(k,i,2) = x_p(k,i-1,2);
            end            
        end
    end
end


for t = 1 : floor(length(data_sample)./B_size)
  
    prior = w_old;
    data_pointer = ((t-1)*B_size+1):t*B_size;
    new_data = data_sample(data_pointer);
    
    for k = 1:N_p
        
        [model_data,D_store] = model_PSLOA(x_p(k,:,:),data_pointer,data_sample,B_size,s0,D_store,k,seita_0);    % calculate the model prediction data x_block(t)use parameters particles k; 
        
        log_likelihood(k) = sum(log10(pdf('norm',model_data-new_data,0,noise_power)));  % calculate the likelihood p(y_block(t) | x_block(t) )
        log_prior(k) = log10(prior(k));
        log_posterior(k) = log_prior(k) + log_likelihood(k);
    end
    
    w_new = 10.^(log_posterior - max(log_posterior));
    
    w_new = w_new./(sum(w_new));
    
    [T_out2(t),std_T_out2(t)] = mean_out_PSLOA(w_new,x_p,data_pointer);
    
    
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
            D_store(k,:,:) = D_store(flags(k),:,:);
        end
        w_old = 1./N_p.*ones(N_p,1);
    end
    
    for k = 1:N_p
        x_p(k,:,:) = update_rule1(x_p(k,:,:),data_pointer,seita_0);
    end
    t,
      
  
end
        
x_p_out = mean(x_p,1);        
T_out = x_p_out(1,:,1);

end
