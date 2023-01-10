function [uA, uB, P, Q] = fBlockDiag(A, B, rankTop, nullSize)
    N = size(A, 1);
    P = eye(N);
    Q = eye(N);
    nullSum = cumsum(nullSize);
for i = 1:length(nullSize)-1
    n = rankTop + nullSum(i) + nullSize(i+1);
%     n = size(A, 1);
    B11 = B(1:rankTop, 1:rankTop);
    B21 = B(rankTop+nullSum(i)+1:n, 1:rankTop);
    U11 = -B21*B11^(-1);
    

    U = [
        [eye(rankTop), zeros(rankTop, N-rankTop);
        zeros(nullSize(i), rankTop), eye(nullSize(i)), zeros(nullSize(i), N - rankTop - nullSize(i));
        U11, zeros(nullSize(i+1), nullSize(i)), eye(nullSize(i+1)), zeros(nullSize(i+1), N - nullSize(i)-rankTop-nullSize(i+1))];
        zeros(N-nullSum(i+1)-rankTop, nullSize(i+1)+rankTop), eye(N-nullSum(i+1)-rankTop)];
    P = U*P;
%     U1 = [
%         eye(rankTop), zeros(rankTop, n-rankTop);
%         U11, eye(n-rankTop)];
%     U = [
%         U1, zeros(size(U1,1), N - size(U1,2)); ...
%         zeros(size(U1,2), N - size(U1,1)), eye(N - size(U1,2))];

    A = U*A;
    B = U*B;
    
    A21 = A(rankTop+nullSum(i)+1:n, 1:rankTop);
    A22 = A(rankTop+nullSum(i)+1:n, rankTop+nullSum(i)+1:n);
    V11 = -A22^(-1)*A21;

V = [
        [eye(rankTop), zeros(rankTop, N-rankTop);
        zeros(nullSize(i), rankTop), eye(nullSize(i)), zeros(nullSize(i), N - rankTop - nullSize(i));
        V11, zeros(nullSize(i+1), nullSize(i)), eye(nullSize(i+1)), zeros(nullSize(i+1), N - nullSize(i)-rankTop-nullSize(i+1))];
        zeros(N-nullSum(i+1)-rankTop, nullSize(i+1)+rankTop), eye(N-nullSum(i+1)-rankTop)];
Q = Q*V;

A = A*V;
B = B*V;
%     
% %     V1 = [eye(rankTop), zeros(rankTop, n-rankTop);
% %         V21, eye(n-rankTop)];
% %     V = [
% %         V1, zeros(size(V1,1), N - size(V1,2)); ...
% %         zeros(size(V1,2), N - size(V1,1)), eye(N - size(V1,2))];
%     V = [
%         eye(rankTop), zeros(rankTop, N-rankTop);
%         zeros(nullSize(i), rankTop+nullSum(i)), eye(nullSize(i)), zeros(nullSize(i), N - rankTop - nullSum(i));
%         V11, zeros(nullSize(i+1), nullSize(i)), eye(nullSize(i+1)), zeros(nullSize(i+1), N - nullSum(i+1)-rankTop);
%         zeros(N-nullSum(i+1)-rankTop, nullSize(i+1)+rankTop), eye(N-nullSum(i+1)-rankTop)];
% 
%     P = U*P;
%     Q = Q*V;
%     uA = hatA*V;
%     uB = hatB*V;
end
uA = A;
uB = B;

end