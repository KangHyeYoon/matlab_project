clc; clear; close all;

%% BPSK BER simulation
%% Bit generation
Bp_Nbits = 10^5; % BPSK의 Bit의 수를 10^5개로 설정 
Bp_Tx_bits = randi([0 1],1,Bp_Nbits);  % Bp_Nbits개의 0 또는 1의 랜덤한 값을 갖는 BPSK 송신Bits 생성
snr_idx = 0; 
for snr_dB = -10:1:10 % BPSK의 SNR의 범위를 –10부터 10까지 설정
    snr_idx = snr_idx + 1; % BPSK의 SNR Index 생성
    snr_set(snr_idx) = snr_dB; % SNR Value Storage Variables % BPSK의 SNR 값을 벡터로 생성

    %% BPSK modulation 
    Bp_M = 2; % Modulation order
    Bp_k = log2(Bp_M); % BPSK의 k = 1
    Bp_Nsyms = Bp_Nbits / Bp_k; % Number of modulated symbols % Bp_Nsyms = 10^5 / 1 = 10^5
    Bp_Tx_mod = 1 – 2 * (Bp_Tx_bits==0); % BPSK modulation  % 연산을 통해 송신 Bits를 –1 또는 1로 변조

    %% AWGN channel passing    
    Bp_EbN0 = 10^(snr_dB/10); % BPSK EbN0
    Bp_EsN0 = Bp_EbN0 * Bp_k; % BPSK EsN0
    

    Bp_EbN0_set(snr_idx) = Bp_EbN0; % BPSK의 Eb/N0 값을 벡터로 생성

    % Noise power
    Bp_Eb_N0 = 1/Bp_EbN0; % Eb/N0 BPSK의 Noise Power 생성
    Bp_Es_N0 = 1/Bp_EsN0;  % Es/N0 BPSK의 Noise Power 생성

    % Noise generation
    Bp_Eb_n = sqrt(Bp_Eb_N0/2) * randn(1,Bp_Nsyms); % Eb/N0 BPSK의 Noise 생성
    Bp_Es_n = sqrt(Bp_Es_N0/2) * randn(1,Bp_Nsyms); % Es/N0 BPSK의 Noise 생성

    % Es/N0 BPSK Received signal
    Bp_Eb_Rx = Bp_Tx_mod + Bp_Eb_n; % Eb/N0 QPSK 변조신호에 Noise를 더한 수신신호 생성
    Bp_Es_Rx = Bp_Tx_mod + Bp_Es_n; % Es/N0 QPSK 변조신호에 Noise를 더한 수신신호 생성

    %% BPSK demodulation
    DecisionBoundary = 0; % Bits를 나누는 기준 = 0
    Bp_Eb_Tx_bits_hat = Bp_Eb_Rx > DecisionBoundary; % Eb/N0 QPSK 수신신호가 0보다 크면 1 작으면 0을 저장
    Bp_Es_Tx_bits_hat = Bp_Es_Rx > DecisionBoundary; % Es/N0 QPSK 수신신호가 0보다 크면 1 작으면 0을 저장

    %% Bit error rate
    Bp_Eb_err = sum(Bp_Eb_Tx_bits_hat ~= Bp_Tx_bits); % Eb/N0 QPSK 수신Bits와 송신Bits를 비교 후 합
    Bp_Eb_BER(snr_idx) = Bp_Eb_err / Bp_Nbits; % 더한 값을 bit 수 (Bp_Nbits) 로 나누어 BER 계산
    Bp_Es_err = sum(Bp_Es_Tx_bits_hat ~= Bp_Tx_bits); % Es/N0 QPSK 수신Bits와 송신Bits를 비교 후 합
    Bp_Es_BER(snr_idx) = Bp_Es_err / Bp_Nbits; % 더한 값을 bit 수 (Bp_Nbits) 로 나누어 BER 계산
end

%% QPSK BER simulation
% Bit generation
Qp_Nbits = 10^5; % QPSK의 Bit의 수를 10^5개로 설정 
Qp_Tx_bits = randi([0 1],1,Qp_Nbits); % Bp_Nbits개의 0 또는 1의 랜덤한 값을 갖는 QPSK 송신Bits 생성

snr_idx = 0;
for snr_dB = -10:1:10 % QPSK의 SNR의 범위를 –10부터 10까지 설정
    snr_idx = snr_idx + 1; % QPSK의 SNR Index 생성
    snr_set(snr_idx) = snr_dB; % QPSK의 SNR 값을 벡터로 생성

    %% QPSK modulation 
    Qp_M = 4; % Modulation order
    Qp_k = log2(Qp_M); % QPSK의 k = 2
    Qp_Nsyms = Qp_Nbits/ Qp_k; % Number of modulated Symbols % Bp_Nsyms = 10^5 / 2 = 50000

    % Eb/N0 QPSK Gray code & Equally interval modulation % 연산을 통해 Eb/N0 QPSK 변조
    Qp_Eb_Tx_mod = ((2*Qp_Tx_bits(1:2:end)-1) + 1j*(2*Qp_Tx_bits(2:2:end)-1)); 
    % Es/N0 QPSK Gray code & Equally interval modulation % 연산을 통해 Es/N0 QPSK 변조
    Qp_Es_Tx_mod = ((2*Qp_Tx_bits(1:2:end)-1) + 1j*(2*Qp_Tx_bits(2:2:end)-1)) / sqrt(2); 

    % Eb/N0 QPSK None eq.interval code modulation % 연산을 통해 None eq.interval Eb/N0 QPSK 변조
    Qp_Eb_Tx_ne_mod = (1/sqrt(2)*(2*Qp_Tx_bits(1:2:end)-1) + 1j*(sqrt(3/2)*(2*Qp_Tx_bits(2:2:end)-1)));
    % Es/N0 QPSK None eq.interval code modulation % 연산을 통해 None eq.interval Es/N0 QPSK 변조
    Qp_Es_Tx_ne_mod = (1/2*(2*Qp_Tx_bits(1:2:end)-1) + 1j*(sqrt(3)/2*(2*Qp_Tx_bits(2:2:end)-1))); 

    real_Qp_Tx_ng_bits = Qp_Tx_bits(1:2:end); % Real part of the Transmit Bits % 실수 부분만 변수에 저장
    imag_Qp_Tx_ng_bits = Qp_Tx_bits(2:2:end); % Imag part of the Transmit Bits % 복소수 부분만 변수에 저장

    for i = 1 : Qp_Nsyms
        if real_Qp_Tx_ng_bits(i) == imag_Qp_Tx_ng_bits(i)
      % Eb/N0 None Gray Code QPSK modulation
            Qp_Eb_Tx_ng_mod(i) = ((-1) + 1j*(2*imag_Qp_Tx_ng_bits(i)-1)); 
      % Es/N0 None Gray Code QPSK modulation
            Qp_Es_Tx_ng_mod(i) = ((-1) + 1j*(2*imag_Qp_Tx_ng_bits(i)-1)) / sqrt(2); 
        elseif real_Qp_Tx_ng_bits(i) ~= imag_Qp_Tx_ng_bits(i)
      % Eb/N0 None Gray Code QPSK modulation
            Qp_Eb_Tx_ng_mod(i) = ((1) + 1j*(2*imag_Qp_Tx_ng_bits(i)-1)); 
      % Es/N0 None Gray Code QPSK modulation
            Qp_Es_Tx_ng_mod(i) = ((1) + 1j*(2*imag_Qp_Tx_ng_bits(i)-1)) / sqrt(2); 
        end
    end % 반복문(for)과 조건문(if)을 이용하여 None gray code Eb/N0, Es/N0 QPSK 변조

    %% Eb/N0 AWGN channel passing
    Qp_EbN0 = 10^(snr_dB/10); % QPSK EbN0
    Qp_EbN0_set(snr_idx) = Qp_EbN0; % QPSK의 Eb/N0 값을 벡터로 생성

    Qp_Eb_N0 = 1/Qp_EbN0; % Noise power  % Eb/N0 QPSK의 Noise Power 생성

    % noise generation % Eb/N0 QPSK의 Noise 생성
    Qp_Eb_n = sqrt(Qp_Eb_N0/2) * randn(1,Qp_Nsyms) + 1j * sqrt(Qp_Eb_N0/2) * randn(1,Qp_Nsyms); 

    % Gray code & Equally interval modulation Received signal Eb/N0
    Qp_Eb_Rx = Qp_Eb_Tx_mod + Qp_Eb_n; 
    % Eb/N0 QPSK 변조신호에 Noise를 더한 수신신호 생성
    % None Equally interval code modulation Received signal Eb/N0
    Qp_Eb_Rx_ne = Qp_Eb_Tx_ne_mod + Qp_Eb_n; 
    %  None Equally interval Eb/N0 QPSK 변조신호에 Noise를 더한 수신신호 생성
    % None Gray code Received signal Eb/N0
    Qp_Eb_Rx_ng = Qp_Eb_Tx_ng_mod + Qp_Eb_n;
    %  None Gray code Eb/N0 QPSK 변조신호에 Noise를 더한 수신신호 생성

    %% Es/N0 AWGN channel passing
    Qp_EsN0 = Qp_EbN0 * Qp_k; % QPSK EsN0

    % Noise power
    Qp_Es_N0 = 1/Qp_EsN0; % Es/N0 QPSK의 Noise Power 생성

    % noise generation % Es/N0 QPSK의 Noise 생성
    Qp_Es_n = sqrt(Qp_Es_N0/2) * randn(1,Qp_Nsyms) + 1j * sqrt(Qp_Es_N0/2)*randn(1,Qp_Nsyms); 

    % Gray code & Equally interval modulation Received signal Es/N0
    Qp_Es_Rx = Qp_Es_Tx_mod + Qp_Es_n; % Es/N0 QPSK 변조신호에 Noise를 더한 수신신호 생성
    % None Equally interval modulation Received signal Es/N0
    Qp_Es_Rx_ne = Qp_Es_Tx_ne_mod + Qp_Es_n; 
    %  None Equally interval Es/N0 QPSK 변조신호에 Noise를 더한 수신신호 생성
    % None Gray code Received signal Es/N0
    Qp_Es_Rx_ng = Qp_Es_Tx_ng_mod + Qp_Es_n; 
    %  None Gray code Eb/N0 QPSK 변조신호에 Noise를 더한 수신신호 생성

    %% Gray code & Equally interval applications (Eb/N0) QPSK demodulation
    % Gray code & Equally interval Eb/N0 QPSK Symbol 
    Qp_Eb_Rx_Syms = (sign(real(Qp_Eb_Rx)) + 1j * sign(imag(Qp_Eb_Rx))); 
    % 생성된 수신 신호를 sign() 함수를 통해 Eb/N0 QPSK Symbol로 복조
    % 부호함수 sign(x) : x > 0 = 1, x = 0, x < 0 = - 1, x => 복소수  = x.abs(x)를 출력  

    % Demodulation Bit Storage Variable
    Qp_Eb_Rx_demod = zeros([1,Qp_Nbits]); % Eb/N0 복조 Bits를 저장하기 위한 영벡터 생성

    % Gray code & Equally interval Eb/N0 demodulation 
    Qp_Eb_Rx_demod(1:2:end) = (real(Qp_Eb_Rx_Syms) + 1) / 2; : % 홀수번째 열에 복조된 실수 Bits 저장
    Qp_Eb_Rx_demod(2:2:end) = (imag(Qp_Eb_Rx_Syms) + 1) / 2; % 짝수번째 열에 복조된 복소수 Bits 저장

    Qp_Eb_Rx_Bits = Qp_Eb_Rx_demod;

    %% Gray code & Equally interval applications (Es/N0) QPSK demodulation
    % Gray code & Equally interval Es/N0 QPSK Symbol
    Qp_Es_Rx_Syms = (sign(real(Qp_Es_Rx)) + 1j * sign(imag(Qp_Es_Rx))) / sqrt(2); 
    % 생성된 수신 신호를 sign() 함수를 통해 Es/N0 QPSK Symbol로 복조

    % Demodulation Bit Storage Variable
    Qp_Es_Rx_demod = zeros([1,Qp_Nbits]); % Es/N0 복조 Bits를 저장하기 위한 영벡터 생성

    % Gray code Es/N0 QPSK demodulation
    Qp_Es_Rx_demod(1:2:end) = (real(Qp_Es_Rx_Syms) * sqrt(2) + 1) / 2; % 홀수번째 열에 복조된 실수 Bits 저장
    Qp_Es_Rx_demod(2:2:end) = (imag(Qp_Es_Rx_Syms) * sqrt(2) + 1) / 2; % 짝수번째 열에 복조된 복소수 Bits 저장

    Qp_Es_Rx_Bits = Qp_Es_Rx_demod;
    
    %% None Equally interval applications (Eb/N0) QPSK demodulation
    % None Equally interval Eb/N0 QPSK Symbol
    % 생성된 수신 신호를 sign() 함수를 통해 None Equally interval Eb/N0 QPSK Symbol로 복조
    Qp_Eb_Rx_ne_Syms = 1/2*sign(real(Qp_Eb_Rx_ne)) + 1j * (sqrt(3)/2*sign(imag(Qp_Eb_Rx_ne))); 

    % Demodulation Bit Storage Variable     % None Equally interval Eb/N0 QPSK 복조된 Bits를 저장하기 위한 영벡터 생성
    Qp_Eb_Rx_ne_demod = zeros([1,Qp_Nbits]); 

    % None Equally interval Eb/N0 QPSK demodulation
    % 홀수번째 열에 None Equally interval Eb/N0 복조된 실수 Bits 저장
    Qp_Eb_Rx_ne_demod(1:2:end) = (real(Qp_Eb_Rx_ne_Syms) * 2 + 1) / 2; 
    % 짝수번째 열에 None Equally interval Eb/N0 복조된 실수 Bits 저장
    Qp_Eb_Rx_ne_demod(2:2:end) = (imag(Qp_Eb_Rx_ne_Syms) / (sqrt(3)/2) + 1) / 2; 

    Qp_Eb_Rx_ne_Bits = Qp_Eb_Rx_ne_demod;
    
    % None Equally interval applications (Eb/N0) QPSK demodulation
    % None Equally interval Es/N0 QPSK Symbol
    Qp_Es_Rx_ne_Syms = 1/sqrt(2)*sign(real(Qp_Es_Rx_ne)) + 1j * (sqrt(3/2)*sign(imag(Qp_Es_Rx_ne)));
    % 생성된 수신 신호를 sign() 함수를 통해 None Equally interval Eb/N0 QPSK Symbol로 복조

    % Demodulation Bit Storage Variable % None Equally interval Es/N0 복조 Bits를 저장하기 위한 영벡터 생성
    Qp_Es_Rx_ne_demod = zeros([1,Qp_Nbits]); 

    % None Equally interval Es/N0 QPSK demodulation
    % 홀수번째 열에 None Equally interval Es/N0 복조된 실수 Bits 저장
    Qp_Es_Rx_ne_demod(1:2:end) = (real(Qp_Es_Rx_ne_Syms) * sqrt(2) + 1) / 2;
    % 짝수번째 열에 None Equally interval Es/N0 복조된 실수 Bits 저장
    Qp_Es_Rx_ne_demod(2:2:end) = (imag(Qp_Es_Rx_ne_Syms) / sqrt(3/2) + 1) / 2; 

    Qp_Es_Rx_ne_Bits = Qp_Es_Rx_ne_demod;

    %% None Gray code applications (Eb/N0) QPSK demodulation
    % None Gray code Eb/N0 QPSK Symbol
    Qp_Eb_Rx_ng_Syms = (sign(real(Qp_Eb_Rx_ng)) + 1j * sign(imag(Qp_Eb_Rx_ng))); 
    % 생성된 수신 신호를 sign() 함수를 통해 None Gray code Eb/N0 QPSK Symbol로 복조
    
    % Real part of the demodulation signal % 복조된 None Gray code Eb/N0 QPSK Symbol의 실수 부분 저장
    real_Qp_Eb_Rx_ng_Syms = real(Qp_Eb_Rx_ng_Syms); 
    % Imag part of the demodulation signal % 복조된 None Gray code Eb/N0 QPSK Symbol의 복소수 부분 저장
    imag_Qp_Eb_Rx_ng_Syms = imag(Qp_Eb_Rx_ng_Syms); 

    % Demodulation Bit Real partial Storage Variable
    Qp_Eb_Rx_ng_Bit_1 = zeros(1,Qp_Nsyms);
    % None Gray code Eb/N0 QPSK 복조된 실수 Bits를 저장하기 위한 영벡터 생성
    % Demodulation Bit Imag partial Storage Variable
    Qp_Eb_Rx_ng_Bit_2 = zeros(1,Qp_Nsyms); 
    % None Gray code Eb/N0 QPSK 복조된 복소수 Bits를 저장하기 위한 영벡터 생성

    % None Gray Code Eb/N0 QPSK demodulation
    for i = 1 : Qp_Nsyms
        if (real_Qp_Eb_Rx_ng_Syms(i)> 0) && (imag_Qp_Eb_Rx_ng_Syms(i)> 0)
            Qp_Eb_Rx_ng_Bit_1(i) = 0;
            Qp_Eb_Rx_ng_Bit_2(i) = 1;
        elseif (real_Qp_Eb_Rx_ng_Syms(i)< 0) && (imag_Qp_Eb_Rx_ng_Syms(i)< 0)
            Qp_Eb_Rx_ng_Bit_1(i) = 0;
            Qp_Eb_Rx_ng_Bit_2(i) = 0;
        elseif (real_Qp_Eb_Rx_ng_Syms(i)> 0) && (imag_Qp_Eb_Rx_ng_Syms(i)< 0)
            Qp_Eb_Rx_ng_Bit_1(i) = 1;
            Qp_Eb_Rx_ng_Bit_2(i) = 0;
        elseif (real_Qp_Eb_Rx_ng_Syms(i)< 0) && (imag_Qp_Eb_Rx_ng_Syms(i) > 0)
            Qp_Eb_Rx_ng_Bit_1(i) = 1;
            Qp_Eb_Rx_ng_Bit_2(i) = 1;
        end
    end % 생성한 영벡터에 복조된 실수 및 복소수 Bits를 각각 저장

    % Demodulation Bit Storage Variable
    Qp_Eb_Rx_ng_Bit = zeros(1,Qp_Nbits); % 실수와 복소수를 하나의 벡터로 합치기 위한 영벡터 생성 

    % None Gary code Eb/N0 QPSK demodulation
    Qp_Eb_Rx_ng_Bit(1:2:end) = Qp_Eb_Rx_ng_Bit_1; % 홀수번째 열에 복조된 실수 Bits 저장
    Qp_Eb_Rx_ng_Bit(2:2:end) = Qp_Eb_Rx_ng_Bit_2; % 짝수번째 열에 복조된 복소수 Bits 저장

   
    %% None Gray code applications (Es/N0) QPSK demodulation
    % None Gray code Es/N0 QPSK Symbol
    Qp_Es_Rx_ng_Syms = (sign(real(Qp_Es_Rx_ng)) + 1j * sign(imag(Qp_Es_Rx_ng))) / sqrt(2); 
    % 생성된 수신 신호를 sign() 함수를 통해 None Gray code Es/N0 QPSK Symbol로 복조

    % Real part of the demodulation signal % 복조된 None Gray code Es/N0 QPSK Symbol의 실수 부분 저장
    real_Qp_Es_Rx_ng_Syms = real(Qp_Es_Rx_ng_Syms); 
    % Imag part of the demodulation signal % 복조된 None Gray code Es/N0 QPSK Symbol의 복소수 부분 저장
    imag_Qp_Es_Rx_ng_Syms = imag(Qp_Es_Rx_ng_Syms); 

    % Real part of the demodulation signal
    Qp_Es_Rx_ng_Bit_1 = zeros(1,Qp_Nsyms); 
    % None Gray code Es/N0 QPSK 복조된 실수 Bits를 저장하기 위한 영벡터 생성
    % Imag part of the demodulation signal
    Qp_Es_Rx_ng_Bit_2 = zeros(1,Qp_Nsyms); 
    % None Gray code Es/N0 QPSK 복조된 복소수 Bits를 저장하기 위한 영벡터 생성

    % None Gray Code Es/N0 QPSK demodulation
    for i = 1 : Qp_Nsyms
        if (real_Qp_Es_Rx_ng_Syms(i)> 0) && (imag_Qp_Es_Rx_ng_Syms(i)> 0)
            Qp_Es_Rx_ng_Bit_1(i) = 0;
            Qp_Es_Rx_ng_Bit_2(i) = 1;
        elseif (real_Qp_Es_Rx_ng_Syms(i)< 0) && (imag_Qp_Es_Rx_ng_Syms(i)< 0)
            Qp_Es_Rx_ng_Bit_1(i) = 0;
            Qp_Es_Rx_ng_Bit_2(i) = 0;
        elseif (real_Qp_Es_Rx_ng_Syms(i)> 0) && (imag_Qp_Es_Rx_ng_Syms(i)< 0)
            Qp_Es_Rx_ng_Bit_1(i) = 1;
            Qp_Es_Rx_ng_Bit_2(i) = 0;
        elseif (real_Qp_Es_Rx_ng_Syms(i)< 0) && (imag_Qp_Es_Rx_ng_Syms(i)> 0)
            Qp_Es_Rx_ng_Bit_1(i) = 1;
            Qp_Es_Rx_ng_Bit_2(i) = 1;
        end
    end % 생성한 영벡터에 복조된 실수 및 복소수 Bits를 각각 저장

    % Demodulation Bit Storage Variable
    Qp_Es_Rx_ng_Bit = zeros(1,Qp_Nbits); % 실수와 복소수를 하나의 벡터로 합치기 위한 영벡터 생성 

    % None Gary code Es/N0 QPSK demodulation
    Qp_Es_Rx_ng_Bit(1:2:end) = Qp_Es_Rx_ng_Bit_1; % 홀수번째 열에 복조된 실수 Bits 저장
    Qp_Es_Rx_ng_Bit(2:2:end) = Qp_Es_Rx_ng_Bit_2; % 짝수번째 열에 복조된 복소수 Bits 저장

    %% BER according to Eb/N0 
    % Gray code & Equally interval Eb/N0 QPSK BER % 최적의 Eb/N0 QPSK BER 계산
    Qp_Eb_Bits_err = sum(Qp_Tx_bits ~= Qp_Eb_Rx_Bits); % 송신 Bits와 수신 Bits를 비교하여 오류의 합 저장
    Qp_Eb_BER(snr_idx) = Qp_Eb_Bits_err/Qp_Nbits; % 오류의 합을 Bits 수로 나누어 BER값 저장
 
    % None Equally interval Eb/N0 QPSK BER % None Equally interval의 Eb/N0 QPSK BER 계산
    Qp_Eb_ne_Bits_err = sum(Qp_Tx_bits ~= Qp_Eb_Rx_ne_Bits); % 송신 Bits와 수신 Bits를 비교하여 오류의 합 저장
    Qp_Eb_ne_BER(snr_idx) = Qp_Eb_ne_Bits_err/Qp_Nbits; % 오류의 합을 Bits 수로 나누어 BER값 저장

    % None Gray code Eb/N0 QPSK BER % None Gray code의 Eb/N0 QPSK BER 계산
    Qp_Eb_ng_Bits_err = sum(Qp_Tx_bits ~= Qp_Eb_Rx_ng_Bit); % 송신 Bits와 수신 Bits를 비교하여 오류의 합 저장
    Qp_Eb_ng_BER(snr_idx) = Qp_Eb_ng_Bits_err/Qp_Nbits; % 오류의 합을 Bits 수로 나누어 BER값 저장

    %% BER according to Es/N0
    % Gray code & Equally interval Es/N0 QPSK BER % 최적의 Es/N0 QPSK BER 계산
    Qp_Es_Bits_err = sum(Qp_Tx_bits ~= Qp_Es_Rx_Bits); % 송신 Bits와 수신 Bits를 비교하여 오류의 합 저장
    Qp_Es_BER(snr_idx) = Qp_Es_Bits_err/Qp_Nbits; % 오류의 합을 Bits 수로 나누어 BER값 저장

    % None Equally interval Es/N0 QPSK BER % None Equally interval의 Es/N0 QPSK BER 계산
    Qp_Es_ne_Bits_err = sum(Qp_Tx_bits ~= Qp_Es_Rx_ne_Bits); % 송신 Bits와 수신 Bits를 비교하여 오류의 합 저장
    Qp_Es_ne_BER(snr_idx) = Qp_Es_ne_Bits_err/Qp_Nbits; % 오류의 합을 Bits 수로 나누어 BER값 저장


    % None Gray code Es/N0 QPSK BER % None Gray code의 Es/N0 QPSK BER 계산
    Qp_Es_ng_Bits_err = sum(Qp_Tx_bits ~= Qp_Es_Rx_ng_Bit); % 송신 Bits와 수신 Bits를 비교하여 오류의 합 저장
    Qp_Es_ng_BER(snr_idx) = Qp_Es_ng_Bits_err/Qp_Nbits; % 오류의 합을 Bits 수로 나누어 BER값 저장
end

%% BER Theory
% BPSK theory
BPSK_BER_theory = 1/2*erfc(sqrt(Bp_EbN0_set)); % BPSK BER 이론 수식에 따른 이론값
% QPSK theory
QPSK_BER_theory = 1/Qp_k * erfc(sqrt(Qp_k*Qp_EbN0_set)*sin(pi/Qp_M)); % QPSK BER 이론 수식에 따른 이론값

%% performance of optimal QPSK Figure 

% (1) performance of optimal constellation % 최적의 Eb/N0 QPSK 변조 방식에 따른 성상도
figure(1);
plot(real(Qp_Eb_Rx), imag(Qp_Eb_Rx),'o'); hold on; grid on; % performance of optimal + noise Eb/N0 Rx
plot(real(Qp_Eb_Tx_mod), imag(Qp_Eb_Tx_mod),'rx',LineWidth=2); % performance of optimal Eb/N0 Tx
xlabel('real'); ylabel('imag');
legend('performance of optimal + noise Eb/N0 Rx','performance of optimal Eb/N0 Tx');
title('performance of optimal Eb/N0 constellation');

% (2) performance of optimal Es/N0 constellation % 최적의 Es/N0 QPSK 변조 방식에 따른 성상도
figure(2);
plot(real(Qp_Es_Rx), imag(Qp_Es_Rx),'o'); hold on; grid on; % performance of optimal + noise Es/N0 Rx
plot(real(Qp_Es_Tx_mod), imag(Qp_Es_Tx_mod),'rx',LineWidth=2); % performance of optimal Es/N0 Tx
xlabel('real'); ylabel('imag');
legend('performance of optimal + noise Es/N0 Rx','performance of optimal Es/N0 Tx');
title('performance of optimal Es/N0 constellation');

% (3) BER-SNR of performance of optimal QPSK modulation scheme (AWGN channel) 
% 최적의 Eb/N0 & Es/N0 QPSK 변조 방식과 QPSK 이론값에 따른 성능 평가
figure(3);
semilogy(snr_set, Qp_Eb_BER, 'ro-', LineWidth=2); grid on; hold on; % performance of optimal Eb/N0 BER-SNR sim.
semilogy(snr_set, Qp_Es_BER, 'gs-', LineWidth=2); hold on; grid on; % performance of optimal Es/N0 BER-SNR sim.
semilogy(snr_set,QPSK_BER_theory,'k--',LineWidth=2); hold on; grid on; % QPSK BER-SNR Theory
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) P_{b}');
legend('performance of optimal Eb/N0 QPSK', 'performance of optimal Es/n0 QPSK','QPSK theory');
title('최적의 QPSK 변조 방식의 BER-SNR 성능(AWGN 채널)');

% (4) Eb/N0 BER-SNR vs. QPSK Eb/N0 Theoretical Values
% 최적의 Eb/N0 QPSK 변조 방식과 QPSK 이론값에 대한 성능 비교
figure(4);
semilogy(snr_set, Qp_Eb_BER, 'ro-', LineWidth=2); hold on; grid on; % performance of optimal Eb/N0 BER-SNR sim.
semilogy(snr_set,QPSK_BER_theory,'k--',LineWidth=2); % QPSK BER-SNR Theory
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('performance of optimal Eb/N0 QPSK', 'QPSK theory');
title('Eb/N0 QPSK BER-SNR vs. QPSK Eb/N0 이론값 성능 비교');

% (5) Es/N0 BER-SNR vs. QPSK Es/N0 Theoretical Values
% 최적의 Es/N0 QPSK 변조 방식과 QPSK 이론값에 대한 성능 비교
figure(5);
semilogy(snr_set, Qp_Es_BER, 'gs-', LineWidth=2); hold on; grid on; % performance of optimal Es/N0 BER-SNR sim.
semilogy(snr_set,QPSK_BER_theory,'k--',LineWidth=2); % QPSK BER-SNR Theory
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('gray code Es/N0 QPSK', 'QPSK theory');
title('Es/N0 QPSK BER-SNR vs. QPSK Eb/N0 이론값 성능 비교');

% (6) Compare performance with BPSK for each result
% 최적의 Eb/N0 & Es/N0 QPSK 변조 방식과 QPSK 이론값, BPSK 이론값에 대한 성능 비교
figure(6);
semilogy(snr_set, Qp_Eb_BER, 'ro-', LineWidth=2); grid on; hold on; % performance of optimal Eb/N0 BER-SNR sim.
semilogy(snr_set, Qp_Es_BER, 'ms-', LineWidth=2); % performance of optimal  Es/N0 BER-SNR sim.
semilogy(snr_set, Bp_Es_BER, 'gs-', LineWidth=2); % BPSK sim.
semilogy(snr_set, BPSK_BER_theory,'k--'); hold on; grid on; % BPSK BER-SNR Theory
semilogy(snr_set,QPSK_BER_theory,'bx', LineWidth=2);  % QPSK BER-SNR Theory
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('performance of optimal Eb/N0 QPSK', 'performance of optimal Es/n0 QPSK', 'BPSK simulation', 'QPSK theory', 'BPSK theory');
title('각 결과에 대한 BPSK와 성능 비교');

%% None eq.interval QPSK Figure

% (7) N-eq.interval code Eb/N0 constellation % None Equally interval Eb/N0 QPSK 변조 방식에 따른 성상도
figure(7);
plot(real(Qp_Eb_Rx_ne), imag(Qp_Eb_Rx_ne),'o'); hold on; grid on; % N-eq.interval code + noise Eb/N0 Rx
plot(real(Qp_Eb_Tx_ne_mod), imag(Qp_Eb_Tx_ne_mod),'rx',LineWidth=2); % N-eq.interval code Eb/N0 Tx
xlabel('real'); ylabel('imag');
legend('N-eq.interval code + noise Eb/N0 Rx', 'N-eq.interval code Eb/N0 Tx');
title('N-eq.interval code Eb/N0 constellation');

% (8) N-eq.interval code Es/N0 constellation % None Equally interval Es/N0 QPSK 변조 방식에 따른 성상도
figure(8);
plot(real(Qp_Es_Rx_ne), imag(Qp_Es_Rx_ne),'o'); hold on; grid on; % N-eq.interval code + noise Es/N0 Rx
plot(real(Qp_Es_Tx_ne_mod), imag(Qp_Es_Tx_ne_mod),'rx',LineWidth=2); % N-eq.interval code Es/N0 Tx
xlabel('real'); ylabel('imag');
legend('N-eq.interval code + noise Es/N0 Rx', 'N-eq.interval code Es/N0 Tx');
title('N-eq.interval code Es/N0 constellation');

% (9) BER-SNR of N-Equally interval QPSK modulation scheme (AWGN channel)
%  None Equally interval Eb/N0 & Es/N0 QPSK 변조 방식과 QPSK 이론값에 따른 성능 평가
figure(9);
semilogy(snr_set,QPSK_BER_theory,'k--',"lineWidth",2); hold on; grid on; % eq.interval Eb/N0 BER-SNR Theory
semilogy(snr_set, Qp_Eb_ne_BER, 'ro-', LineWidth=2); grid on; hold on; % None eq.interval Eb/N0 BER-SNR sim.
semilogy(snr_set, Qp_Es_ne_BER, 'gs--', LineWidth=2); % None eq.interval Es/N0 BER-SNR sim.
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('QPSK theory','N-eq.interval Eb/N0 QPSK','N-eq.interval Es/N0 QPSK');
title('N-Equally interval QPSK 변조 방식 BER-SNR 성능 (AWGN 채널)');

% (10) Eq interval QPSK Eb/N0 BER-SNR vs. N-Eq interval QPSK Eb/N0 BER-SNR
%  Equally interval Eb/N0 QPSK 변조 방식과 None Equally interval Eb/N0 QPSK 변조 방식에 대한 성능 비교
figure(10);
semilogy(snr_set, Qp_Eb_BER, 'ro-', LineWidth=2); hold on; grid on; % eq.interval Eb/N0 BER-SNR sim.
semilogy(snr_set, Qp_Eb_ne_BER, 'go-', LineWidth=2); % None eq.interval Eb/N0 BER-SNR sim.
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('eq.interval Eb/N0 QPSK','N-eq.interval Eb/N0 QPSK');
title('Eq interval QPSK Eb/N0 BER-SNR vs. N-Eq interval QPSK Eb/N0 BER-SNR 성능 비교');

% (11) Eq interval QPSK Es/N0 BER-SNR vs. N-Eq interval QPSK Es/N0 BER-SNR
%  Equally interval Es/N0 QPSK 변조 방식과 None Equally interval Es/N0 QPSK 변조 방식에 대한 성능 비교
figure(11);
semilogy(snr_set, Qp_Es_BER, 'rs-', LineWidth=2); hold on; grid on; % eq.interval Es/N0 BER-SNR sim.
semilogy(snr_set, Qp_Es_ne_BER, 'gs-', LineWidth=2);% None eq.interval Es/N0 BER-SNR sim.
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('eq.interval Es/N0 QPSK','N-eq.interval Es/N0 QPSK');
title('Eq interval QPSK Es/N0 BER-SNR vs. N-Eq interval QPSK Es/N0 BER-SNR 성능 비교');

%% None Gray Code QPSK Figure

% (12) N-gray code Eb/N0 constellation % None Gray code Eb/N0 QPSK 변조 방식에 따른 성상도
figure(12);
plot(real(Qp_Eb_Rx_ng), imag(Qp_Eb_Rx_ng),'o'); hold on; grid on; % N-gray code + noise Eb/N0 Rx
plot(real(Qp_Eb_Tx_ng_mod), imag(Qp_Eb_Tx_ng_mod),'rx',LineWidth=2); % N-gray code Eb/N0 Tx
xlabel('real'); ylabel('imag');
legend('N-gray code + noise Eb/N0 Rx', 'N-gray code Eb/N0 Tx');
title('N-gray code Eb/N0 constellation');

% (13) N-gray code Es/N0 constellation % None Gray code Es/N0 QPSK 변조 방식에 따른 성상도
figure(13);
plot(real(Qp_Es_Rx_ng), imag(Qp_Es_Rx_ng),'o'); hold on; grid on; % N-gray code + noise Es/N0 Rx
plot(real(Qp_Es_Tx_ng_mod), imag(Qp_Es_Tx_ng_mod),'rx',LineWidth=2); % N-gray code Es/N0 Tx
xlabel('real'); ylabel('imag');
legend('N-gray code + noise Es/N0 Rx', 'N-gray code Es/N0 Tx');
title('N-gray code Es/N0 constellation');

% (14) BER-SNR of N-Gay code QPSK modulation scheme (AWGN channel)
%  None Gray code Eb/N0 QPSK 변조 방식과 None Gray code Es/N0 QPSK 변조 방식에 따른 성능 평가
figure(14);
semilogy(snr_set, Qp_Eb_ng_BER, 'ro-', LineWidth=2); grid on; hold on; % None Gray code Eb/N0 BER-SNR sim.
semilogy(snr_set, Qp_Es_ng_BER, 'gs--', LineWidth=2); hold on; grid on; % None Gray code Es/N0 BER-SNR sim.
semilogy(snr_set,QPSK_BER_theory,'k--',"lineWidth",2); % QPSK BER-SNR Theory
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('N-gray code Eb/N0 QPSK', 'N-gray code Es/n0 QPSK', 'QPSK theory');
title('N-gray code QPSK 변조 방식의 BER-SNR 성능(AWGN 채널)');

% (15) gray code QPSK Eb/N0 BER-SNR vs. N-gray code QPSK Eb/N0 BER-SNR
%  Gray code Eb/N0 QPSK 변조 방식과 None Gray code Eb/N0 QPSK 변조 방식에 대한 성능 비교
figure(15);
semilogy(snr_set, Qp_Eb_BER, 'ro-', LineWidth=2); hold on; grid on; % None Gray code Eb/N0 BER-SNR sim.
semilogy(snr_set, Qp_Eb_ng_BER, 'go-', LineWidth=2); % Eb/N0 QPSK BER-SNR Theory
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('gray code Eb/N0 QPSK','N-gray code Eb/N0 QPSK');
title('gray code QPSK Eb/N0 BER-SNR vs. N-gray code QPSK Eb/N0 BER-SNR 성능 비교');

% (16) gray code QPSK Es/N0 BER-SNR vs. N-gray code QPSK Es/N0 BER-SNR
%  Gray code Es/N0 QPSK 변조 방식과 None Gray code Es/N0 QPSK 변조 방식에 대한 성능 비교
figure(16);
semilogy(snr_set, Qp_Es_BER, 'rs-', LineWidth=2); hold on; grid on; % None Gray code Es/N0 BER-SNR sim.
semilogy(snr_set, Qp_Es_ng_BER, 'gs-', LineWidth=2); % Es/N0 QPSK BER-SNR Theory
xlabel('SNR(dB)'); ylabel('Bits 에러 확률(BER) Pb');
legend('gray code Es/N0 QPSK','N-gray code Es/N0 QPSK');
title('gray code QPSK Es/N0 BER-SNR vs. N-gray code QPSK Es/N0 BER-SNR 성능 비교');