clear all;
close all;
clc;

% Simulation parameters
N = 10^5;          % Number of symbols
M = 16;            % Constellation size (M-PSK)
k = log2(M);       % Bits per symbol
m_fading = [1 4 10];             % Nakagami fading parameter
Eb_N0_dB = 0:25;   % Eb/N0 range in dB
Es_N0_dB = Eb_N0_dB + 10*log10(k);  % Converting to Es/N0

% Gray coding mapping
ref = 0:M-1;
map = bitxor(ref, floor(ref/2));
[~, ind] = sort(map);

% Preallocate arrays
nSymErr = zeros(1, length(Eb_N0_dB));
ser_theoretical = zeros(1, length(Eb_N0_dB));

colors = ['r', 'g', 'b', 'm'];
markers = ['o', 's', '^'];
lineStyles = {'-', '--'};
 
idx_m = 1;

figure;

for m = m_fading
    % Main simulation loop
    for ii = 1:length(Eb_N0_dB)
        % Symbol generation
        ipDec = randi([0 M-1], 1, N);  % Random symbols
        ipPhase = ipDec * 2*pi/M;  % Phase mapping
        s = exp(1j * ipPhase);  % Complex symbols
        
        % Nakagami-m fading channel
        omega = 1;  % Average power
        gamma = gamrnd(m, omega/m, 1, N);  % Generate gamma-distributed random variables
        h_amp = sqrt(gamma);  % Nakagami-m amplitude
        h_phase = 2*pi*rand(1,N);  % Random phase
        h = h_amp .* exp(1j*h_phase);  % Complex fading channel
        
        % Normalize fading to ensure unit average power
        h = h / sqrt(mean(abs(h).^2));
        
        % Add AWGN noise
        snr = 10^(Es_N0_dB(ii)/10);
        noise_var = 1/(2*snr);  % Noise variance per dimension
        n = sqrt(noise_var) * (randn(1,N) + 1j*randn(1,N));
        
        % Received signal
        y = h .* s + n;
        
        % Perfect channel estimation and equalization
        y_eq = y ./ h;
        
        % Phase detection
        rx_phase = angle(y_eq);
        rx_phase(rx_phase < 0) = rx_phase(rx_phase < 0) + 2*pi;
        
        % Symbol decision
        est_phase = 2*pi/M * round(rx_phase/(2*pi/M));
        est_phase(est_phase == 2*pi) = 0;
        est_dec = round(est_phase * M/(2*pi));
        
        % Count symbol errors
        nSymErr(ii) = sum(ipDec ~= est_dec);
        
        % Theoretical SER calculation for M-PSK over Nakagami-m fading
        % Using MGF-based approach
        gamma_s = 10^(Es_N0_dB(ii)/10);  % Average SNR per symbol
        
        theta = linspace(0, pi*(M-1)/M, 1000);
        integrand = zeros(size(theta));
            
        for j = 1:length(theta)
            % MGF-based approach for M-PSK
            gPSK = (sin(pi/M))^2;
            temp = (m/(m + gamma_s * gPSK/(sin(theta(j))^2)))^m;
            integrand(j) = temp;
        end
        
        % Calculate SER using numerical integration
        ser_theoretical(ii) = 1/(pi) * trapz(theta, integrand);
    end

    simSer = nSymErr/N;
    
    semilogy(Eb_N0_dB, ser_theoretical, [colors(idx_m) lineStyles{1} markers(1)], 'LineWidth', 2);
    hold on;
    semilogy(Eb_N0_dB, simSer, [colors(idx_m) lineStyles{2} markers(3)], 'LineWidth', 2);
    hold on

    idx_m = idx_m + 1;
end

% AWGN channel (no fading)
for ii = 1:length(Eb_N0_dB)
    % Symbol generation
    ipDec = randi([0 M-1], 1, N);  % Random symbols
    ipPhase = ipDec * 2*pi/M;  % Phase mapping
    s = exp(1j * ipPhase);  % Complex symbols
    
    % Add AWGN noise
    snr = 10^(Es_N0_dB(ii)/10);
    noise_var = 1/(2*snr);  % Noise variance per dimension
    n = sqrt(noise_var) * (randn(1,N) + 1j*randn(1,N));
    
    % Received signal
    y = s + n;
    
    % Phase detection
    rx_phase = angle(y);
    rx_phase(rx_phase < 0) = rx_phase(rx_phase < 0) + 2*pi;
    
    % Symbol decision
    est_phase = 2*pi/M * round(rx_phase/(2*pi/M));
    est_phase(est_phase == 2*pi) = 0;
    est_dec = round(est_phase * M/(2*pi));
    
    % Count symbol errors
    nSymErr(ii) = sum(ipDec ~= est_dec);
end

% Theoretical SER for AWGN channel
gamma_s = 10.^(Es_N0_dB/10);  % SNR per symbol
ser_theoretical = 2*qfunc(sqrt(2*gamma_s)*sin(pi/M));

simSer = nSymErr/N;

semilogy(Eb_N0_dB, ser_theoretical, [colors(idx_m) lineStyles{1} markers(1)], 'LineWidth', 2);
hold on;
semilogy(Eb_N0_dB, simSer, [colors(idx_m) lineStyles{2} markers(3)], 'LineWidth', 2);

xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate');
title(['SER vs Eb/N0 for ' num2str(M) '-PSK']);
legend('Theoretical Nakagami m=1', 'Simulation Nakagami m=1', ...
       'Theoretical Nakagami m=4', 'Simulation Nakagami m=4', ...
       'Theoretical Nakagami m=10', 'Simulation Nakagami m=10', ...
       'Theoretical AWGN', 'Simulation AWGN');
axis([0 25 1e-5 1]);
grid on;
hold off;