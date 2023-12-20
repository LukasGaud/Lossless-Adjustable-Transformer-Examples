% Quarter Car Example for the Thesis
% Clearing all the variables
clear all
% Defining symbols
syms s
ms = 380;
mu = 55;
kt = 350e3;
ks = 2e4;
cs = 10500;

% ms = 2;
% mu = 1;
% kt = 1;
% ks = 1;
% cs = 1;

% C1 = @(s)(-[
%     0, 1, 0;
%     1, 0, -1;
%      0, -1, 1]);
C1 = @(s)(-[
    0, -1, 0;
    1, 0, 1;
     0, 1, -1]);
% C1 = @(s)(-[
%     0, 1, 0;
%     1, 0, -1;
%      0, -1, 1]);
C2 = @(s)(-[
    0, 0;
    1, 1;
    -1, -1]);
Q1 = @(s)([
    cs, 0;
    0, ks]);
P1 = @(s)([
    1, 0;
    0, s]);
Q2 = @(s)([
    kt, 0 ,0;
    0, ms*s, 0;
    0, 0, mu*s]);
P2 = @(s)([
    s, 0, 0;
    0, 1, 0;
    0, 0, 1]);
I1 = [1, 0, 0; 0, 1, 0; 0, 0, 1];
fR = @(s)([ ...
    zeros(3,3), eye(3), zeros(3,2), C1(s)', zeros(3,2), zeros(3,3); ...
    zeros(2,3), zeros(2,3), P1(s), zeros(2,3), -Q1(s), zeros(2,3); ...
    zeros(3,3), zeros(3,3), zeros(3,2), Q2(s), zeros(3,2), -P2(s); ...
    zeros(2,3), zeros(2,3), zeros(2,2), C2(s)', eye(2), zeros(2,3); ...
    C1(s), zeros(3,3), -C2(s), zeros(3,3), zeros(3,2), eye(3)]);

R0 = fR(0);
Rs = fR(s);

% B = -R0(:, [3, 1, 5]);
% A = -R0(:, [2, 4, 6:end]);
B = -R0(:, [6, 1, 5]);
A = -R0(:, [2, 4, 3, 7:end]);
% A(:,1) = -A(:,1);
% E = double((Rs(:, [2, 4, 6:end]) + A)/s);
% Rsquare = Rs(:, [3, 2, 4, 7:end]);
% Rsquare(:,1) = -Rsquare(:,1);
% E = double((Rsquare + A)/s);
E = double((Rs(:, [2, 4, 3, 7:end]) + A)/s);
% C = [0, 0, -1, zeros(1, 10)];
C = [0, 0, 1, zeros(1, 10)];
D = 0;
[uE1i, uA1i, P, Q, rankTop, nullSize] = fKron(E, A);

% Form a transformation P1
P1i = eye(size(A));
for k = 1:length(P)
    P1i = P{k}*P1i;
end
% Form coordinate changes
Q1i = eye(size(A));
for k = 1:length(Q)
    Q1i = Q1i*Q{k};
end

% Constructing a block diagonal form
[blockA1, blockE1, blockU1, blockV1] = fBlockDiag(uA1i, uE1i, rankTop, nullSize);
fullQ = Q1i * blockV1;
fullP = blockU1 * P1i;
Cfull = C*fullQ;
uB1 = fullP*B;

% Eliminate infinite elementary divisors
E_tilde = blockE1(rankTop+1:end, rankTop+1:end);
A_tilde = blockA1(rankTop+1:end, rankTop+1:end);
B_tilde = uB1(rankTop+1:end, :);
C_tilde = Cfull(rankTop+1:end);

% Final State-Space
% D_bar = double(C_tilde*(s*E_tilde - A_tilde)^(-1)*B_tilde);
% Change of sign on the force coordinates
T = [
    1, 0, 0;
    0, -1 ,0;
    0, 0, -1];
D_bar = simplify(C_tilde*(s*E_tilde - A_tilde)^(-1)*B_tilde(:,1));
E_bar = blockE1(1:rankTop, 1:rankTop);
J_bar = blockA1(1:rankTop, 1:rankTop);
A_bar = T^(-1)*E_bar^(-1)*J_bar*T
B_bar = T^(-1)*E_bar^(-1)*uB1(1:rankTop, 1) 
F_bar = E_bar^(-1)*uB1(1:rankTop, [2,3]);
C_bar = Cfull(1:rankTop)*T
H_bar = simplify(C_tilde*(s*E_tilde - A_tilde)^(-1)*B_tilde(:,[2, 3]));
T = [
    1, 0, 0;
    0, -1 ,0;
    0, 0, -1];
% Checking coordinates
fullCoord = fullQ^(-1);