syms kt ks ms mu cs positive real
syms s

A = [0, 0, -kt; 0, 0, 0; 1/(ms + mu), 0, 0];
B = [0; 1; -ms/(ms+mu)];
% B = [0; ks; kt*ms/(ms + mu)];
C = [0, ks, -kt*ms/(ms + mu)];
% C = [0, 1, ms/(ms+mu)];
D1 = cs;
D2 = mu*ms/(ms + mu);
G1 = C*(s*eye(3)-A)^(-1)*B + D1 + s*D2;
simplify(1/G1);
A2 = [0, 1, 0, 0; -ks/ms, -cs/ms, ks/ms, cs/ms; 0, 0, 0, 1; ks/mu, cs/mu, -(kt + ks)/mu, -cs/mu];
C2 = [0, 1, 0, -1];
B2  =[0; 1/ms; 0; -1/mu];
G2 = C2*(s*eye(4) - A2)^(-1)*B2;
simplify(G2)
simplify(1/G1)