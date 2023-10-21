syms ms mu ks kt cs p P positive real
syms s

% Representation 1
% x1 = zs-zr; x2 = zu-zr; x3 = zsDot; x4 = zuDot
M = [
    ms + p^2, -p^2;
    -p^2, mu + p^2];
N = [
    -ks, ks, -cs-p*(s*p), cs + p*(s*p);
    ks, -ks - kt, cs + p*(s*p), -cs - p*(s*p)];

A1 = [
    0, 0, 1, 0;
    0, 0, 0, 1;
    M^(-1)*N];
F1 = [1;0];
simplify(A1)

% Representation 2
A21 = [
    0, 0, -kt;
    0, 0, 0;
    1/(ms + mu), 0, 0];
C2 = [0, -ks, kt*ms/(ms+mu)];
B2 = [0; 1; -ms/(ms+mu)];
D21 = -cs;
D22 = -mu*ms/(ms + mu);
alpha = (1 - 1/p * D22 * 1/p)^(-1);

A2 = [
    A21, B2*1/p;
    alpha*1/p*C2, alpha*1/p*(D21*1/p - D22*(s*p/(p^2)))];

A3 = [
    A21, B2*1/p;
    alpha*1/p*C2, alpha*1/p*(D21*1/p - D22*(s*p/(P)))];
% simplify(A2)

% Transformation 
T = [
    0, 0, ms, mu;
    1, -1, 0, 0;
    0, 1, 0, 0;
    0, 0, p, -p];
Tdot = [
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, (s*p), (-s*p)];

simplify(T^(-1)*A2*T - T^(-1)*Tdot)
simplify(A2)
simplify(A3)