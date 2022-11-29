%% OFDM 기초 실습

dt = 0.01;  % 시간 간격
t = 0:dt:20-dt;  % 시간축 생성

% 전송 신호 생성
T  = 10;    % 심볼 길이
x = zeros(1, numel(t));
x(1:numel(0:dt:T-dt)) = 1;

figure;
plot(t,x);

% 반송파 매핑
c_1  = cos(2*pi*(1/T)*t);
c_2  = cos(2*pi*(2/T)*t);
c_3  = cos(2*pi*(3/T)*t);
c_4  = cos(2*pi*(4/T)*t);

% 시간에 반송파 곱하기
x_1 = x.*c_1;
x_2 = x.*c_2;
x_3 = x.*c_3;
x_4 = x.*c_4;

% 각각 figure
figure;
subplot(4, 1, 1);
plot(t, x_1);
subplot(4, 1, 2);
plot(t, x_2);
subplot(4, 1, 3);
plot(t, x_3);
subplot(4, 1, 4);
plot(t, x_4);


% OFDM 심볼 생성
y = x_1 + x_2 + x_3 + x_4;     % 더함

figure;
plot(x,y);