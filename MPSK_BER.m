clear all;
close all;
clc;

% Simulation parameters
N = 10^6;          % Number of symbols
M = 4;            % Constellation size (16-PSK)
k = log2(M);       % Bits per symbol
m_fading = [1 4 10];             % Nakagami fading parameter
Eb_N0_dB = 0:25;   % Eb/N0 range in dB
Es_N0_dB = Eb_N0_dB + 10*log10(k);  % Converting to Es/N0

% Gray coding mapping
ref = 0:M-1;
map = bitxor(ref, floor(ref/2));
[~, ind] = sort(map);

% Preallocate arrays
nBitErr = zeros(1, length(Eb_N0_dB));
ber_theoretical = zeros(1, length(Eb_N0_dB));

colors = ['r', 'g', 'b', 'm'];
markers = ['o', 's', '^'];
lineStyles = {'-', '--'};
 
idx_m = 1;

figure;

for m = m_fading


    % Main simulation loop
    for ii = 1:length(Eb_N0_dB)
        % Symbol generation
        ipBit = randi([0 1], 1, N*k);  % Random binary data
        ipBitReshape = reshape(ipBit, k, N)';  % Reshape to k bits per symbol
        
        % Binary to decimal conversion
        bin2DecMatrix = ones(N,1) * (2.^((k-1):-1:0));  
        grayRef = sum(ipBitReshape .* bin2DecMatrix, 2)';
        
        % Gray coded constellation mapping
        ipDec = ind(grayRef + 1) - 1;  
        ipPhase = ipDec * 2*pi/M;  % Phase mapping
        s = exp(1j * ipPhase);  % Complex symbols
        
        % Nakagami-m fading channel
        % Generate Nakagami-m envelope with unit mean power (Ω = 1)
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
        
        % Gray decoding
        rx_gray = map(est_dec + 1);
        
        % Convert decisions to binary
        rx_bits = zeros(1, N*k);
        for idx = 1:N
            temp = dec2bin(rx_gray(idx), k) - '0';
            rx_bits((idx-1)*k + 1 : idx*k) = temp;
        end
        
        % Count bit errors
        nBitErr(ii) = sum(ipBit ~= rx_bits);
        
        % Theoretical BER calculation for M-PSK over Nakagami-m fading
        % Using MGF-based approach
        gamma_b = 10^(Eb_N0_dB(ii)/10);  % Average SNR per bit
        gamma_s = gamma_b * k;            % Average SNR per symbol
        
        theta = linspace(0, pi*(M-1)/M, 1000);
        integrand = zeros(size(theta));
            
        for j = 1:length(theta)
            % MGF-based approach for M-PSK
            gPSK = (sin(pi/M))^2;
            temp = (m/(m + gamma_s * gPSK/(sin(theta(j))^2)))^m;
            integrand(j) = temp;
        end
            
            % Calculate BER using numerical integration
            Ps = 1/(pi) * trapz(theta, integrand);
            % Convert SER to BER (Gray coding approximation)
            ber_theoretical(ii) = Ps/k;
    end

    simBer = nBitErr/(N*k);
    
    semilogy(Eb_N0_dB, ber_theoretical, [colors(idx_m) lineStyles{1} markers(1)], 'LineWidth', 2);
    hold on;
    semilogy(Eb_N0_dB, simBer, [colors(idx_m) lineStyles{2} markers(3)], 'LineWidth', 2);
    hold on

    idx_m = idx_m + 1;

end





    % AWGN channel (no fading)
    for ii = 1:length(Eb_N0_dB)
        % Symbol generation
        ipBit = randi([0 1], 1, N*k);  % Random binary data
        ipBitReshape = reshape(ipBit, k, N)';  % Reshape to k bits per symbol
        
        % Binary to decimal conversion
        bin2DecMatrix = ones(N,1) * (2.^((k-1):-1:0));  
        grayRef = sum(ipBitReshape .* bin2DecMatrix, 2)';
        
        % Gray coded constellation mapping
        ipDec = ind(grayRef + 1) - 1;  
        ipPhase = ipDec * 2*pi/M;  % Phase mapping
        s = exp(1j * ipPhase);  % Complex symbols
        
        % Nakagami-m fading channel
        % Generate Nakagami-m envelope with unit mean power (Ω = 1)
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
        y = s + n;
        
        % Perfect channel estimation and equalization
        y_eq = y;
        
        % Phase detection
        rx_phase = angle(y_eq);
        rx_phase(rx_phase < 0) = rx_phase(rx_phase < 0) + 2*pi;
        
        % Symbol decision
        est_phase = 2*pi/M * round(rx_phase/(2*pi/M));
        est_phase(est_phase == 2*pi) = 0;
        est_dec = round(est_phase * M/(2*pi));
        
        % Gray decoding
        rx_gray = map(est_dec + 1);
        
        % Convert decisions to binary
        rx_bits = zeros(1, N*k);
        for idx = 1:N
            temp = dec2bin(rx_gray(idx), k) - '0';
            rx_bits((idx-1)*k + 1 : idx*k) = temp;
        end
        
        % Count bit errors
        nBitErr(ii) = sum(ipBit ~= rx_bits);
        
        % Theoretical BER calculation for M-PSK over Nakagami-m fading
        % Using MGF-based approach
        gamma_b = 10^(Eb_N0_dB(ii)/10);  % Average SNR per bit
        gamma_s = gamma_b * k;            % Average SNR per symbol
        
        % Integration for BER calculation
        theta = linspace(0, pi/2, 1000);
        gamma_psk = (sin(pi/M))^2;
        mgf_integrand = (1 + (gamma_s * gamma_psk)./(m * sin(theta).^2)).^(-m);
        
        if(M<16)
        ber_theoretical(ii) = qfunc(sqrt(2*gamma_b));
        end
    
    end

if (M==16)
ber_theoretical = (1/k)*erfc(sqrt(k*10.^(Eb_N0_dB/10))*sin(pi/M));
end

    simBer = nBitErr/(N*k);
    
    semilogy(Eb_N0_dB, ber_theoretical, [colors(idx_m) lineStyles{1} markers(1)], 'LineWidth', 2);
    hold on;
    semilogy(Eb_N0_dB, simBer, [colors(idx_m) lineStyles{2} markers(3)], 'LineWidth', 2);
    hold on

xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
title(['BER vs Eb/N0 for ' num2str(M) 'PSK']);
legend('Theoretical Nakagami m=1', 'Simulation Nakagami m=1', 'Theoretical Nakagami m=4', 'Simulation Nakagami m=4', 'Theoretical Nakagami m=10', 'Simulation Nakagami m=10', 'Theoretical Nakagami AWGN', 'Simulation Nakagami AWGN');
axis([0 25 1e-5 1]);
hold off;
