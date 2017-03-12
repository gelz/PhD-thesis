function [T_out2,std_T_out2] = RB_pitch_tracker_modified_fast(data_input,k_level,T0,SNR_desired_input)

% data_sample = clean_data;
fs = 16000;
% 
 noise = randn(length(data_input),1);
% prior_SNR = 10*log10((sum(data_sample.^2)./sum((noise).^2)));
% k = 10^((prior_SNR-SNR_input)/20);
% 
% y_ob = data_sample + k*noise;
% y = y_ob;
% data_sample2 = data_sample;
% data_sample = y;

y = data_input';
data_sample = data_input';

noise_level = mean(k_level*abs(noise));
c_test = randn(100,1);
noise_compare = mean(abs(c_test));

noise_level = noise_level./noise_compare;


% N_p =500;  %  number of particles used;
N_p = 1000;
N_number = 100;  % define the maximum number of periods to be estimated;

p_order = 8;
K_order = 10;   % define the number of harmonic components used;

data_init = data_sample;
%[x_p,hyper_par] = initialize_particles_sin_RB2_temp(N_p,N_number,T0,p_order,K_order,data_init,fs);
[x_p,hyper_par] = initialize_particles_sin_RB2_temp_randomP(N_p,N_number,T0,p_order,K_order,data_init,fs);


hyper_par2 = 0.0001; % the second part of ar coefficent variance;
hyper_par = [hyper_par,hyper_par2];

w0 = ones(N_p,1)./N_p;  % assign weights to each particles;
w_old = w0;
w_new = w0;
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
% the updating scheme:
%     update Kalman filter parameters every point; 
%     then resampling only after a block has been processed; 
block_size = 128;


% model_process_noise_power = 0.001;
% model_process_noise_power = 0.002;
model_process_noise_power = 0.000009;

sigma_noise = model_process_noise_power;
ob_noise_power = sigma_noise;
% if k_level <0.15
%     k_level = 0.15;
% end
% 
%  ob_noise_power =  1*k_level;
%  ob_noise_power2 = 1*k_level;
 
 
%  if SNR_desired_input<=0
% %      if SNR_desired_input > -5
% %       ob_noise_power = 2*k_level;
% %       ob_noise_power2 = 2*k_level;
% %      else
%       ob_noise_power = 2*k_level;
%       ob_noise_power2 = 2*k_level;
% %      end
%  end
%  ob_noise_power = 2*noise_level;
%  ob_noise_power2 = 2*noise_level;

noise_level = k_level;
%  if SNR_desired_input<5
% %       ob_noise_power = 2*noise_level;
%       ob_noise_power2 = 2*noise_level;
%  else
% %      ob_noise_power = 0.05;
%      ob_noise_power2 = 0.05;
%  end
% ob_noise_power2 = k_level.^2;
% ob_noise_power2 = 1*k_level.^2;  % 2*kevel works for -5 dB;  4*k_level works for 0 dB;
ob_noise_power2 =  1*k_level.^2;
if ob_noise_power2 < 0.01
    ob_noise_power2 = 0.01;  
end

% if ob_noise_power2 < 0.0001
%     ob_noise_power2 = 0.0001;
% end


speech_est = zeros(N_p,length(data_sample));
speech_est(:,1:12) = ones(N_p,1)*data_sample(1:12)';
drive_est = zeros(N_p,length(data_sample));
ar_est = zeros(N_p,length(data_sample),p_order);
K1 = cell(N_p,1);
K2 = cell(N_p,1);


var_coef = zeros(p_order+K_order,1);
var_coef2 = zeros(p_order+K_order,1);

for i = 1:p_order
% var_coef(i) = 1.5*abs(squeeze(x_p(1,1,i)));    % optimal in 5 dB and
% below
% var_coef(i) = 0.5*abs(squeeze(x_p(1,1,1)));    % between period;
% var_coef(i) = (0.145*abs(squeeze(x_p(1,1,1)))).^2;    % between period;
% var_coef(i) = 0.001;
var_coef(i) = 0.001;

var_coef2(i) = (0.000*abs(squeeze(x_p(1,1,1)))).^2;   % within period;
end
for j = (p_order+1 : p_order+K_order)
%         var_coef(j) = 0.1*abs(squeeze(x_p(1,1,j)));  % better in high snr dB;
% var_coef(i) = 0.01;
%         var_coef(j) = (0.05*abs(squeeze(x_p(1,1,j)))).^2;  % between period;
                var_coef(j) = 0.0005;  % between period;

var_coef2(i) = (0.0000*abs(squeeze(x_p(1,1,1)))).^2;   % within period;

end
% for i = 1:N_p
%     K1{i,1} = diag(var_coef2);
% end
for i = 1:N_p
    K1{i,1} = diag(var_coef2);
end
posterior_prob = ones(N_p,length(data_sample));

% Pointer_num = zeros(N_p,1);
% Pointer_num2 = zeros(N_p,1);
c = zeros(N_p, K_order);
T_out2 = nan(length(data_sample)-1,1);
std_T_out2 = T_out2;
Pointer_num = nan(N_p,1);
Pointer_previous = ones(N_p,1);
Pointer_current = nan(N_p,1);


for k = 1:N_p
    Pointer_num(k) = length(find(x_p(k,:,end)~=0));
    Pointer_current(k) = sum(x_p(k,:,end));
end
C = nan(1,K_order+p_order);

for i = 12 : (length(data_sample)-1)
    

    w_new = w_old;
    
     A =  eye(p_order+K_order);                   % A_i,k   
     G = diag(var_coef2.*ones(1,p_order+K_order)');

%      w_k_noise = randn(p_order+K_order,1);    
    
    
    
    for k = 1: N_p
        
%         Pointer_num(k) = length(find(x_p(k,:,end)~=0));  % determine the current number of periods for particle k;                    
%         Pointer_num2(k) = Pointer_num(k);
     if     Pointer_current(k)<i                                        % determine if N_t is complete;
            x_p(k,:,:) = update_rule2_harmonics_sin_RB(x_p(k,:,:),i,hyper_par);  % add rao-blackwellization method; update nonlinear part;
            Pointer_num(k) = Pointer_num(k)+1;
           
            T0 = round(x_p(k,Pointer_num(k),end));
            Pointer_previous(k) = Pointer_current(k);
            Pointer_current(k) = Pointer_current(k)+T0 ;

            
            
            unit_freq = 2*pi/T0;
            unit_phase = unit_freq *(i- Pointer_previous(k));
            for m = 1 : (K_order/2)
                c(k,m) = cos((m-1)*unit_phase);
            end
            for m = ((K_order/2)+1) : K_order
                c(k,m) = sin((m-K_order/2)*unit_phase);
            end

        if i < (p_order+1)    
            C =  [zeros(1,p_order+1-i),-speech_est(k,i-1:-1:1),c(k,:)];
        else
            C =  [-speech_est(k,i-1:-1:i-p_order),c(k,:)];
        end

     sigma_v_k = ob_noise_power;
%      v_k = randn(1);
     G = diag(var_coef.*ones(1,p_order+K_order)');
     
     Temp = reshape((x_p(k,Pointer_num(k)-1,1:(p_order+K_order))),(p_order+K_order),1);
     x_temp = A*Temp;     
     x_p(k,Pointer_num(k),1:(p_order+K_order)) = reshape(x_temp,1,1,(p_order+K_order));
     
     K2{k,1} = G*G' + A*K1{k,1}*A';      
      
     C_st = C*K2{k,1}*C' + sigma_v_k;
%      u_st = C*x_temp;
     u_st = C*x_temp;
%      speech_est(k,i) = u_st + sqrt(C_st)*randn(1);
     
     
     phi_t = 1/ob_noise_power2 + 1/C_st;
     theta_t = data_sample(i)/ob_noise_power2 + u_st/C_st;
     s_est_t = theta_t / phi_t;
     speech_est(k,i) = s_est_t + sqrt(1/phi_t)*randn(1);    
%      speech_est(k,i) = s_est_t;

 
     
     TT = sigma_v_k + C*K2{k,1}*C';
     
     J = K2{k,1}*C'*TT^(-1);
     
     K1{k,1} = (eye(p_order+K_order)-J*C)*K2{k,1};
     
     x_temp2 = x_temp + J*(speech_est(k,i)-C*x_temp);
     x_p(k,Pointer_num(k),1:(p_order+K_order)) = reshape(x_temp2,1,1,(p_order+K_order));
%      posterior_prob(k,i) = 1/realsqrt(2*pi)/realsqrt(phi_t)*exp(-(data_sample(i)-s_est_t).^2/2/(ob_noise_power2)-(s_est_t-u_st).^2/2/C_st)/realsqrt(ob_noise_power2*C_st);

    posterior_prob(k,i) = -0.5*log10(phi_t)-0.5*log10(ob_noise_power2)-0.5*log10(2*pi)+1/log(10)*(-(data_sample(i)-s_est_t).^2/2/(ob_noise_power2)-(s_est_t-u_st).^2/2/C_st)- 0.5*log10(C_st);

%      posterior_prob(k,i) = sqrt(2*pi)/sqrt(phi_t)*normpdf(data_sample(i),s_est_t,sqrt(ob_noise_power2))*normpdf(s_est_t,u_st,sqrt(C_st));
 %      speech_est(k,i) = C*x_temp+sigma_noise*randn(1);
%      drive_est(k,i) = C(p_order+1:end)*x_temp(p_order+1:end);
%      ar_est(k,i,:) = reshape(x_temp(1:p_order),1,1,p_order);
%      posterior_prob(k,i) = normpdf(data_sample(i),speech_est(k,i),ob_noise_power2);
%      
   
        
        else

%      w_k_noise = randn(p_order+K_order,1);
            T0 = round(x_p(k,Pointer_num(k),end));
            unit_freq = 2*pi/T0;
            unit_phase = unit_freq *(i- Pointer_previous(k));
            for m = 1 : (K_order/2)
                c(k,m) = cos((m-1)*unit_phase);
            end
            for m = ((K_order/2)+1) : K_order
                c(k,m) = sin((m-K_order/2)*unit_phase);
            end
        if i < (p_order+1)    
            C =  [zeros(1,p_order+1-i),-speech_est(k,i-1:-1:1),c(k,:)];
        else
            C =  [-speech_est(k,i-1:-1:i-p_order),c(k,:)];
        end

     sigma_v_k = ob_noise_power;

%      G = diag(var_coef2.*ones(1,p_order+K_order)');  % G =zeros;

     x_temp = reshape((x_p(k,Pointer_num(k),1:(p_order+K_order))),(p_order+K_order),1);
%      x_temp = Temp;     
%      x_p(k,Pointer_num(k),1:(p_order+K_order)) = reshape(x_temp,1,1,(p_order+K_order));

% 
%      sigma_v_k = ob_noise_power;
%      v_k = randn(1);
     G = diag(var_coef.*ones(1,p_order+K_order)');
%        G = diag(var_coef2.*ones(1,p_order+K_order)');
%      Temp = reshape((x_p(k,Pointer_num2(k),1:(p_order+K_order))),(p_order+K_order),1);
%      x_temp = A*Temp;     
%      x_p(k,Pointer_num(k),1:(p_order+K_order)) = reshape(x_temp,1,1,(p_order+K_order));
%      
%      K2{k,1} = K1{k,1};   
     k1 = K1{k,1};
     k2 = k1+ G*G';
%      k2 = k1;
     C_st = C*k2*C' + sigma_v_k;
%      u_st = C*x_temp;
     u_st = C*x_temp;
%      speech_est(k,i) = u_st + sqrt(C_st)*randn(1);
     
     
     phi_t = 1/ob_noise_power2 + 1/C_st;
     theta_t = data_sample(i)/ob_noise_power2 + u_st/C_st;
     s_est_t = theta_t / phi_t;
     speech_est(k,i) = s_est_t + sqrt(1/phi_t)*randn(1);    
%      speech_est(k,i) = s_est_t;
 
     
     TT = sigma_v_k + C*k2*C';
     
     J = k2*C'*TT^(-1);
     
     k1 = (eye(p_order+K_order)-J*C)*k2;
     
     x_temp2 = x_temp + J*(speech_est(k,i)-C*x_temp);
     x_p(k,Pointer_num(k),1:(p_order+K_order)) = reshape(x_temp2,1,1,(p_order+K_order));
     
     K1{k,1} = k1;
     
%      posterior_prob(k,i) = realsqrt(2*pi)/realsqrt(phi_t)*normpdf(data_sample(i),s_est_t,sqrt(ob_noise_power2))*normpdf(s_est_t,u_st,sqrt(C_st));
%     posterior_prob(k,i) = 1/realsqrt(2*pi)/realsqrt(phi_t)*exp(-(data_sample(i)-s_est_t).^2/2/(ob_noise_power2)-(s_est_t-u_st).^2/2/C_st)/realsqrt(ob_noise_power2*C_st);
    posterior_prob(k,i) = -0.5*log10(phi_t)-0.5*log10(ob_noise_power2)-0.5*log10(2*pi)+1/log(10)*(-(data_sample(i)-s_est_t).^2/2/(ob_noise_power2)-(s_est_t-u_st).^2/2/C_st)- 0.5*log10(C_st);
%      posterior_prob(k,i) = sqrt(2*pi)/sqrt(phi_t)*normpdf(data_sample(i),s_est_t,sqrt(ob_noise_power2))*normpdf(s_est_t,u_st,sqrt(C_st));

%      K1{k,1} =  matrix_temp;     
     end

    end
    

%     if ((i<180)&&(mod(i,block_size) == 0))||((i>180)&&(mod(i,block_size/4) == 0))
    if ((i<180)&&(mod(i,block_size) == 0))||((i>180)&&(mod(i,block_size/2) == 0))

        log_posterior = sum((posterior_prob(:,i-block_size+12:i)),2);
%         log_posterior = sum(log10(posterior_prob(:,i-block_size+12:i)),2);
        
        w_new = 10.^(log_posterior - max(log_posterior));
        w_new = w_new.*w_old;
        w_new = w_new./(sum(w_new));       
        
        
        ess = 1./(sum(w_new.^2));
        if (ess > N_p/2)||(i<15)
            resample_decision = 0;
        else resample_decision = 1;
        end
    
        if resample_decision == 0
            w_old = w_new;
        else
        w_old = w_new;
        flags = resample([1:N_p]',w_old);
            for k = 1:N_p
                x_p(k,:,:) = x_p(flags(k),:,:);
                K1{k,1} = K1{flags(k),1};
                speech_est(k,:)=speech_est(flags(k),:);
                drive_est(k,:)=drive_est(flags(k),:);
                Pointer_num(k) = Pointer_num(flags(k)) ;
                Pointer_previous(k) = Pointer_previous(flags(k));
                Pointer_current(k) = Pointer_current(flags(k)) ;
            end
        w_old = 1./N_p.*ones(N_p,1);
        end
% %     i,
    i, 
%     sum(w_new.^2),   
    
    end

   
    
   [T_out2(i),std_T_out2(i)] = mean_out_TimeVary_harmonics_RB(w_new,x_p,i);   

%%      comment out the belowing to avoid output;
%     i, 
%     sum(w_new.^2),    
end
%     figure();
%     plot(mean(speech_est(:,:),1));
%     figure();
%     plot(squeeze(mean(x_p(:,:,1:p_order))));
%     figure();
%     plot(squeeze(mean(x_p(:,:,(1+p_order):(K_order+p_order)))));
% i = 1;   
% figure();
% plot(squeeze(x_p(i,:,1)));
% hold on;
% plot(mean(squeeze(x_p(:,:,1)),1),'g');
% hold on;
% plot(squeeze(x_p(i,:,2)),'r');
% hold on;
% plot(squeeze(x_p(i,:,3)),'k');
% figure();
% plot(squeeze(x_p(i,:,19)));
% hold on;
% plot(mean(squeeze(x_p(:,:,19)),1),'g');
% hold on;
% plot(squeeze(x_p(i,:,20)),'r');
% hold on;
% plot(squeeze(x_p(i,:,23)),'k');
% figure();
% plot(speech_est(1,:));
% hold on;
% plot(data_input,'r');

end


function [x_p,hyper_par] = initialize_particles_sin_RB2_temp_randomP(N_p,N_number,T0,p_order,K_order,data_input,fs)


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
% 
% p_order = 11;
% K_order = 40;   % define the number of harmonic components used;
if nargin <3
    T0 = 188;
end
% define the first guess of AR coefficients and source parameters, using
% algorithm in 'Joint source-filter optimization for robust glottal source estimation
% in the presence of shimmer and jitter',Prasanta Kumar Ghosh *, Shrikanth
% S. Narayanan, Speech Communication, 2011
% [data_input,fs] = data_read();

T02 = T0;
% T0 = 40:5:100;
% if T02 > 70
%     T0 = [(T02-6*5):5:(T02+6*5)];
% else
%     T0 = [(T02-4*5):5:(T02+4*5)];
% end
if T02 > 70
    T0 = [(T02-2*5):5:(T02+2*5)];
else
    T0 = [(T02-2*5):5:(T02+2*5)];
end
N_T = length(T0);
% N_T,
kk = 1;

for i = 1:N_T
x_p_est(:,i) = joint_est_para_sin_v2(p_order,K_order,data_input,fs,T0(i));   % use the first 3 periods to esimate the initial parameters;
end

%     for j = 1:p_order
%         if abs(x_p_est(j))>0.8
%            x_p_est(j) = sign(x_p_est(j))*0.7;
%         end
%     end
% x_p_est_plus = 0.05*randn(N_p,length(x_p_est));
% x_p = (x_p_est_plus+1)* x_p_est;
% x_p_est = 0.1*randn(p_order+K_order,1);
x_p = zeros(N_p,N_number,size(x_p_est,1)+1);
var_period =8;
for i = 1 : N_p
    ccc = floor(N_p./N_T);
    dd = floor(i./ccc)+1;
    if dd>N_T
        dd= N_T;
    end
    x_p(i,1,1:(end-1)) = reshape(x_p_est(:,dd),1,1,size(x_p_est,1));    %%% this time complete here ( 08/04/2013  16:49)
    x_p(i,1,end) = T0(dd)+round(var_period/3*randn(1));
%     x_p(i,2,1:(end-1)) = (0.05*randn(1,1,length(x_p_est))+1).*reshape(x_p_est,1,1,length(x_p_est));    %%% this time complete here ( 08/04/2013  16:49)
%     x_p(i,2,end) = T0;
%     x_p(i,3,1:(end-1)) = (0.05*randn(1,1,length(x_p_est))+1).*reshape(x_p_est,1,1,length(x_p_est));    %%% this time complete here ( 08/04/2013  16:49)
%     x_p(i,3,end) = T0;
%     x_p(i,4,1:(end-1)) = (0.05*randn(1,1,length(x_p_est))+1).*reshape(x_p_est,1,1,length(x_p_est));    %%% this time complete here ( 08/04/2013  16:49)
%     x_p(i,4,end) = T0;        
end

var_amplitude = 0.001;
var_amplitude_ar = 0.01;

hyper_par = [var_period, var_amplitude,var_amplitude_ar];

end











