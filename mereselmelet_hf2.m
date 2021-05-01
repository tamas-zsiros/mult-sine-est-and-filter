%% base parameters
M=100;
N=4*M+1;
r=0.87;
q=0.25;
P=21;

%% 1,
% generate multi-sine
rand_phase = 2*pi*randn(M,1);
multi_sine = zeros(N,1);
multi_sine_zero_phase = zeros(N,1); 
for j = 1:N
   multi_sine(j) = 1;
   multi_sine_zero_phase(j) = 1;
   for k = 1:M
       c = exp(2*pi*1i/N*j*(2*k+1));
       x = exp(1i*rand_phase(k));
       x_0 = exp(1i*0); %% = 1
       multi_sine(j) = multi_sine(j) + 0.5 * (c*x + conj(c)*conj(x));
       multi_sine_zero_phase(j) = multi_sine_zero_phase(j) + 0.5 * (c*x_0 + conj(c)*conj(x_0));
   end
end

figure(1);
plot(1:N, multi_sine);
hold on;
plot(1:N, multi_sine_zero_phase);
hold off

% analyze the signal
multi_sine = [multi_sine;multi_sine;multi_sine];   %% measure more periods
[y,x_pred] = signal_analyzer(multi_sine, N, M);

error_signal = multi_sine - y;
figure(2)
subplot(3,1,1); 
plot(multi_sine); title('multi sine');
subplot(3,1,2); 
plot(y); title('multi sine analyzer');
subplot(3,1,3); 
plot(error_signal); title('error signal');

%% use the given system with the signal
num_of_sys = [0, (1-r), (1-r), 0, 0];
den_of_sys = [2, 0, 0, 0,-2*r];

sys = tf(num_of_sys, den_of_sys, 1/N) % safety check system
sys_output = filter(num_of_sys,den_of_sys, multi_sine);   % inverz z transform and convolution in one step
figure(3)
subplot(2,1,1); 
plot(multi_sine); title('multi sine');
subplot(2,1,2); 
plot(sys_output); title('sys output');

[sys_output_est, sys_x_0] = signal_analyzer_with_all_freqs(sys_output, N, M);
% Get the amplitudes and phases of components 
sys_magnitudes = abs(sys_x_0(:, end));
sys_phases = rad2deg(angle(sys_x_0(:, end)));

variance = max(multi_sine) * 0.01;
additive_noise = randn(length(multi_sine), 1) * variance^2;

noisy_sys_output = filter(num_of_sys,den_of_sys, multi_sine) + additive_noise;
[noisy_sys_output_est, noisy_sys_x_0] = signal_analyzer_with_all_freqs(noisy_sys_output, N, M);
noisy_sys_magnitudes = abs(noisy_sys_x_0(:, end));
noisy_sys_phases = rad2deg(angle(noisy_sys_x_0(:, end)));

amplitude_diff_summ = 0;
for i = 1:M*2
    if(mod(i,2) == 0)
        amplitude_diff_summ = amplitude_diff_summ + abs(sys_magnitudes(i) - noisy_sys_magnitudes(i));
    end
end

amplitude_diff = zeros(M,1);
for i = 1:M*2
    if(mod(i,1) == 0)
        amplitude_diff(i) = abs(sys_magnitudes(i) - noisy_sys_magnitudes(i));
    end
end

figure(4)
subplot(3,1,1); 
plot(1:M*2, sys_magnitudes); title('sys output magnitudes');
hold on;
plot(1:M*2, noisy_sys_magnitudes);
hold off;
subplot(3,1,2);
plot(1:M*2, sys_phases); title('phases');
hold on;
plot(1:M*2, noisy_sys_phases);
hold off;
subplot(3,1,3);
plot(1:M*2, amplitude_diff); title('magnitude diffs');
%% 1.4
[W, prev_W] = recursive_LS(multi_sine,sys_output, P, 300);
plot_est_sys_output(multi_sine,sys_output, W, P, 5)
figure(6)
plot(1:length(prev_W), prev_W(:, 1:end)'); title('konvergencia LS diagram');
%% 1.5

%% 2.0
%LMS

[W_LMS, prev_W_LMS] = LMS(multi_sine, sys_output, P, 300,zeros(P,1), 0.00001);
prev_W_LMS = prev_W_LMS(:, 2:end);  % cut the first - invalid - value
prev_W_LMS = [prev_W_LMS(21, 1:end); prev_W_LMS(20, 1:end);prev_W_LMS(17, 1:end);prev_W_LMS(16, 1:end);prev_W_LMS(13, 1:end)];  %cut abs. largest values
figure(7);
plot(1:length(prev_W_LMS), prev_W_LMS(:, 1:end)'); title('konvergencia LMS diagram');

plot_est_sys_output(multi_sine, sys_output, W_LMS, P, 8); % plot comparison between estimate and sys output

% modify system, repeat
num_of_sys_mod = [0, (1-r-q), (1-r-q), 0, 0];
den_of_sys_mod = [2, 0, 0, 0,-2*(r-q)];
sys_output_mod = filter(num_of_sys_mod,den_of_sys_mod, multi_sine);
[W_LMS_mod, prev_W_LMS_mod] = LMS(multi_sine, sys_output_mod, P, 300,W_LMS, 0.00001);
prev_W_LMS_mod = prev_W_LMS_mod(:, 2:end);
prev_W_LMS_mod = [prev_W_LMS_mod(21, 1:end); prev_W_LMS_mod(20, 1:end);prev_W_LMS_mod(17, 1:end);prev_W_LMS_mod(16, 1:end);prev_W_LMS_mod(13, 1:end)];
figure(9);
plot(1:length(prev_W_LMS_mod), prev_W_LMS_mod(:, 1:end)'); title('konvergencia LMS diagram módosítás után');
plot_est_sys_output(multi_sine, sys_output_mod, W_LMS_mod, P, 10) % plot comparison between estimate and sys output
%% 3.0
[W_EE, prev_W_EE] = EE(multi_sine, sys_output, 300, zeros(7,1),0.0003); % constants: 300: starting point, 0.0003 convergence rate
figure(11);
plot(1:length(prev_W_EE), prev_W_EE(:, 1:end)'); title('konvergencia EE diagram');
plot_EE(multi_sine, sys_output, W_EE, 12);

[W_EE_mod, prev_W_EE_mod] = EE(multi_sine, sys_output_mod, 300, W_EE, 0.0001);
figure(13);
plot(1:length(prev_W_EE_mod), prev_W_EE_mod(:, 1:end)'); title('konvergencia EE diagram mod után');
plot_EE(multi_sine, sys_output_mod, W_EE_mod, 14);

%% 4.0
[A,B,C,D] = tf2ss(num_of_sys, den_of_sys);  %create state space variables
[A_mod,B_mod,C_mod,D_mod] = tf2ss(num_of_sys_mod, den_of_sys_mod);
system_noise = 0.01 * max(multi_sine);
R = variance ^2 * eye(1); % scalar output, variance previously calculated in 1.3
Q = system_noise ^2 * eye(4); % number of state variables
state_noise = mvnrnd([0 0 0 0], Q, length(noisy_sys_output));
[X, P, traces, prev_X] = kalman_predict(noisy_sys_output,multi_sine, A, B, C, additive_noise, state_noise); %run kalman prediction
%plot
figure(15);
plot(1:200, traces(1:200)); title("traces");
figure(16);
subplot(3,1,1);
plot(1:200, prev_X(:, 1:200)'); title('prev X');
predicted_signal = C * prev_X;
subplot(3,1,2);
plot(1:200, predicted_signal(:, 1:200)'); title('predicted signal');
subplot(3,1,3);
plot(1:200, noisy_sys_output(1:200)); title('sys output');
%subplot(4,1,4);
%error = noisy_sys_output(1:length(predicted_signal)) - predicted_signal(1:end)';
%plot(1:length(predicted_signal), error); title('error');

noisy_sys_output_mod = sys_output_mod + additive_noise;
[X_mod, P_mod, traces_mod, prev_X_mod] = kalman_predict(noisy_sys_output_mod, multi_sine, A_mod, B_mod, C_mod, additive_noise, state_noise);
traces_mod = rmoutliers(traces_mod);

figure(17);
plot(1:500, traces_mod(1:500)); title("traces mod");
figure(18);
subplot(3,1,1);
plot(1:500, prev_X_mod(:, 1:500)'); title('prev X mod');
predicted_signal_mod = C * prev_X_mod;
subplot(3,1,2);
plot(1:500, predicted_signal_mod(:, 1:500)'); title('predicted mod signal');
subplot(3,1,3);
plot(1:500, noisy_sys_output_mod(1:500)); title('sys output mod');

%% functions

function [signal_estimate, x_0] = signal_analyzer(signal, N, M)
    signal_estimate = zeros(3*N,1);
    x_0 = zeros(M, 3*N);
    for j=1:1:length(signal)
        for k=1:M
            c = exp(2*pi*1i/N*j*(2*k+1));
            signal_estimate(j) = signal_estimate(j) + (x_0(k,j) * c + conj(x_0(k,j)) * conj(c));
        end
    
        for l=1:M
        g = 1/N*exp(-2*pi*1i/N*j*(2*l+1));
        x_0(l,j+1) = (signal(j) - signal_estimate(j)) * g + x_0(l,j);
        end
    end
end

function [signal_estimate, x_0] = signal_analyzer_with_all_freqs(signal, N, M)
    signal_estimate = zeros(3*N,1);
    M = M*2;
    x_0 = zeros(M, 3*N);
    for j=1:1:length(signal)
        for k=1:M
            c = exp(2*pi*1i/N*j*(k));
            signal_estimate(j) = signal_estimate(j) + (x_0(k,j) * c + conj(x_0(k,j)) * conj(c));
        end
    
        for l=1:M
        g = 1/N*exp(-2*pi*1i/N*j*(l));
        x_0(l,j+1) = (signal(j) - signal_estimate(j)) * g + x_0(l,j);
        end
    end
end

function [W, prev_W] = recursive_LS(signal,sys_output, window_size, start_point)
    W = zeros(window_size,1);
    P = eye(window_size);
    prev_W = ones(window_size,1);
    counter = 0;
    while counter < 1200-start_point    %1200: the current length of the signal
        X = signal(start_point + counter - window_size:start_point + counter -1)';
        z = sys_output(start_point + counter);
        G = (P*X')/(1 + X * P * X');
        prev_W(:, counter+1) = W;
        W = W + G*(z - X*W);
        P = (eye(window_size) -G*X)*P;
        counter = counter+1;
        
    end
end

function [W_LMS, prev_W_LMS] = LMS(signal, sys_output, window_size, start_point, start_weights, learning_rate)
    W_LMS = start_weights;
    prev_W_LMS = ones(window_size,1); % first column to be deleted
    counter = 0;
    prev_ee_counter = 1;
    while counter < 6*1200-start_point      %iterate over more periods 
        offset = mod(start_point + counter, length(sys_output)+1);  % convert index to array dimension
        if(offset == 0)
            counter = counter+window_size+1;
            continue;
        end
        X = signal(offset - window_size:offset -1);     %previous input signal values
        z = sys_output(offset);                         %current system output
        counter = counter+1;
        prev_W_LMS(:,prev_ee_counter) = W_LMS;          %store previous weights
        prev_ee_counter = prev_ee_counter +1;
        e = z - X' * W_LMS;
        W_LMS = W_LMS + (2 * learning_rate * X' * e)';  %calcuate new weights
        
    end
end

function [W_EE, prev_EE] = EE(signal, sys_output,start_point, start_weights, learning_rate)
    W_EE = start_weights;
    prev_EE = ones(7,1);
    counter = 0;
    prev_ee_counter = 1;
    while counter < 3*1200-start_point      %iterate over more periods 
        offset = mod(start_point + counter, length(sys_output)+1);      % convert index to array dimension
        if(offset == 0)
            counter = counter+7;
            continue;
        end
        X = signal(offset - 3:offset -1);           %previous input signal values
        X_2 = sys_output(offset - 4:offset -1);     %previous system output values
        X = [X;X_2];                                %concat previous values
        z = sys_output(offset);                     %current system output
        prev_EE(:, prev_ee_counter) = W_EE;         %store previous weights
        counter = counter + 1;
        prev_ee_counter = prev_ee_counter + 1;
        e = z - X'*W_EE;
        W_EE = W_EE + (2 * learning_rate * X' * e)';    %calcuate new weights
    end
end

function [X, P, traces, prev_X] = kalman_predict (sys_output, multi_sine, A, B, C, R, Q)
    X = zeros(4,1);
    P = diag(ones(1,4));
    traces = trace(P);
    counter = 1;
    while counter < length(sys_output)
       G = A * P * C' / (C*P*C' + R(counter)) ;
       X = A * X + G * (sys_output(counter) - C * X) + B*multi_sine(counter);
       P = (A - G*C)*P*A' + diag(Q(counter,:)');
       traces(counter) = trace(P);
       prev_X(:, counter) = X;
       counter = counter + 1;
       
    end
end

function [] = plot_est_sys_output (signal, sys_output, weights, P, figure_counter)
error_LS = zeros(length(sys_output),1);
LMS_signal = zeros(length(sys_output),1);
for j=1+P:length(sys_output)
    LMS_signal(j) = dot(weights, signal(j - P:j-1)'); 
    error_LS(j) = dot(weights, signal(j - P:j-1)') - sys_output(j);
end
figure(figure_counter);
subplot(3,1,1);
plot(1:length(sys_output), sys_output); title('sys output');
subplot(3,1,2);
plot(1:length(sys_output), LMS_signal); title('estimation');
subplot(3,1,3);
plot(1:length(sys_output), error_LS); title('error ');
end

function [] = plot_EE (signal, sys_output, weights, figure_counter)
error_LS = zeros(length(sys_output),1);
EE_signal = zeros(length(sys_output),1);
for j=17:length(sys_output)
    EE_signal(j) = dot(weights(1:3), signal(j - 3:j-1)') + dot(weights(4:7), sys_output(j-4:j-1)'); 
    error_LS(j) = EE_signal(j) - sys_output(j);
end
figure(figure_counter);
subplot(3,1,1);
plot(1:length(sys_output), sys_output); title('sys output');
subplot(3,1,2);
plot(1:length(sys_output), EE_signal); title('estimation');
subplot(3,1,3);
plot(1:length(sys_output), error_LS); title('error ');
end
    