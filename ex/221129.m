%% OFDM 기초 실습

clc; clear; close all;

N_sym = 128;    % 디지털 변조 심볼 개수
X = 2 * randi([0 1], 1, N_sym) - 1;    % BPSK 심볼

%% OFDM 변조
% IDFT 수행
x = ifft(X) * sqrt(N_sym);    % IDFT (이산신호에 대하여 주파수 -> 신호 OFDM 변조)
% CP 생성
N_cp = 16;
% CP 삽입
x_cp = x(N_sym-16+1 : end );    % ICI 억제하기 (뒤에있는거 앞으로 갖다붙히기)

x_OFDM = [x_cp, x];

figure;
plot(real(x)); hold on;
plot(real(x_OFDM));
legend('OFDM w/o CP', 'OFDM w/ XCP');

%% OFDM 복조
% CP 제거
x_off = x_OFDM(N_cp+1 : end);
% DFT 수행
X_off = fft(x_off)/sqrt(N_sym);

figure;
stem(X_off);