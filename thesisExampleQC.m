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

C1 = @(s)(-[
    0, 1, 0;
    1, 0, -1;
     0, -1, 1]);
C2 = @(s)(-[
    0, 0;
    -1, -1;
    1, 1]);
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

fR = @(s)([ ...
    zeros(3,3), eye(3), zeros(3,2), C1(s)', zeros(3,2), zeros(3,3); ...
    zeros(2,3), zeros(2,3), P1(s), zeros(2,3), -Q1(s), zeros(2,3); ...
    zeros(3,3), zeros(3,3), zeros(3,2), Q2(s), zeros(3,2), -P2(s); ...
    zeros(2,3), zeros(2,3), zeros(2,2), C2(s)', eye(2), zeros(2,3); ...
    -C1(s), zeros(3,3), -C2(s), zeros(3,3), zeros(3,2), eye(3)]);

R0 = fR(0);
Rs = fR(s);

B = -R0(:, [3, 1, 5]);
A = -R0(:, [2, 4, 6:end]);
E = double((Rs(:, [2, 4, 6:end]) + A)/s);
C = [0, 0, -1, zeros(1, 10)];
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
D_bar = double(C_tilde*(s*E_tilde - A_tilde)^(-1)*B_tilde);
E_bar = blockE1(1:rankTop, 1:rankTop);
J_bar = blockA1(1:rankTop, 1:rankTop);
A_bar = E_bar^(-1)*J_bar;
B_bar = E_bar^(-1)*uB1(1:rankTop, 1);
F_bar = E_bar^(-1)*uB1(1:rankTop, [2,3]);
C_bar = Cfull(1:rankTop);

% Checking coordinates
fullCoord = fullQ^(-1);
for k = 1:4
    index(k) = find(fullCoord(k, :) ~= 0);
end

% Transforming to standard coordinates and re-arranging the columsn with
% permutation matrix.
T1 = [
    -1, 0, 0, 0;
    0, -1, 0, 0;
    0, 0, 1/ks, 1/kt;
    0, 0, 0, 1/kt];
T2 = [
    0, 0, 1, 0;
    1, 0, 0, 0;
    0, 0, 0, 1;
    0, 1, 0, 0];
T = T2*T1;

A_can = T*A_bar*(T^(-1))
B_can = T*B_bar
C_can = C_bar*T^(-1)
D_can = D_bar
F_can = T*F_bar(:,1)
G_can = T*F_bar(:,2)

% Determinant Checks
D1 = det(eye(size(A_bar))*s - A_bar);
C1 = double(coeffs(D1, s, 'All'));
Roots1 = roots(C1);

% State-space model - derived by hand
ssA = [0, 1, 0, 0;...
        -ks/ms, -cs/ms, ks/ms, cs/ms;...
        0, 0, 0, 1;...
        ks/mu, cs/mu, -(kt+ks)/mu, -cs/mu];
Ds = det(eye(4)*s - ssA);
Cs = double(coeffs(Ds, s, 'All'));
Rs = roots(Cs);