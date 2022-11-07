%% 블록코드 실습
clc; clear; close all;

%% 블록코드 부호화
n = 7;
k = 4;
parity_matrix = [
    1 0 1
    1 1 0
    0 1 1
    1 1 1
];   % 패리티 비트
message = [1 1 0 0];    % 임의의 메세지 비트

Ik = [
    1 0 0 0
    0 1 0 0
    0 0 1 0
    0 0 0 1
];
general_matrix =  [
    Ik, parity_matrix
];  % 행렬 생성 (행을 기준으로 양옆에 각각 생성되어 합쳐짐)

% 부호화
c = message * general_matrix;   % [1 1 0 0 2 1 1]
% 2는 이진법에서는 0과 다를 게 없음 -> mod()를 통해 2->0으로 모들러 연산 해주자
c = mod(c,2);   %[1 1 0 0(메세지 비트) / 0 1 1(패리티 비트)]

%% 전송 오류 반영
error_bit = [0 0 0 0 1 0 0];    % 에러가 나면 1이 들어가도록 할 거임
e = c + error_bit;

%% 복호화

% 패리티 검사 행렬 생성
In_k = [
    1 0 0
    0 1 0
    0 0 1
];
H = [parity_matrix.', In_k]; % 행을 중심으로 왼쪽 4열은 parity_matrix를 행렬 뒤집은것(트랜스포즈), 오른쪽 3열은 In_k
H_T = H.';

% 신드롬 생성
syndrom = e * H_T;   % e 대신 c를 넣으면 신드롬도 0,0,0이 나온다
% 모듈러 연산도 꼭 해줘야 한다
syndrom = mod(s,2);    % [1,0,0]    

%% 오류 정정
for i = 1 : length(H_T)
    val(i) = sum(syndrom == H_T(i, :));    % H_T와 syndrom의 각 원소를 비교
                            % syndrom이 100이므로 [101, 110, ...]을 비트 단위로 돌면서 같은 것을 찾는다
                            % 비트가 같은 만큼 1이 나오고 그걸 더함. 결국 3이 나와야 완벽하게 같은 거임
end

err_idx = find(val==(n-k));     % 3이 나온 부분만 저장

% 수신 신호 r에서 에러가 난 위치만 고쳐줄거다
% 1이면 0으로, 0이면 1로 고쳐주면 되니까 +1 후 모듈러 연산 해주면 됨
c_hat = r;
c_hat(err_idx) = mod(c_hat(err_idx)+1,2);

%% 오류 정정 확인
syndrom_2 = mod(c_hat * H_T, 2);    % All 0가 뜨는지(오류가 없어졌는지) 확인