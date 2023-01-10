function [updateA, updateB, P, Q, rankTopBlock, nullSize] = fKron(A, B)
%  - i -> x-direction
%  - j -> y-direction
n(1) = size(A,1);
si(1) = 0;
sj(1) = 0;
d = size(A,1);
j = 1;
stop = true;
updateA = A;
updateB = B;
istartingBlock(1) = d;
jstartingBlock(1) = d;
while stop
    indexA = logical([]);
    indexB = logical([]);
    [aU, aSigma, aV] = svd(updateA(1:jstartingBlock(j),1:istartingBlock(j)));
    ri = rank(aSigma);
    sj(j+1) = n(j) - ri;
    istartingBlock(j+1) =  istartingBlock(j)-sj(j+1);
    if n(j) - ri == 0
        stop = false;
        l = j-1;
        break
    else
%     Ap = double(updateA(1:jstartingBlock(j),1:istartingBlock(j))*aV);
    aV = aV^(-1);
    for k = 1:size(aSigma,2)
%         indexA(k) = all(Ap(:,k) == 0);
        indexA(k) = all(aSigma(:,k) == 0);
    end
    aV = [aV(:,~indexA), aV(:,indexA)]';
    updateA(1:jstartingBlock(j), 1:istartingBlock(j)) = updateA(1:jstartingBlock(j), 1:istartingBlock(j))*aV;
    updateB(1:jstartingBlock(j), 1:istartingBlock(j)) = updateB(1:jstartingBlock(j), 1:istartingBlock(j))*aV;
    if j > 1
            updateA(jstartingBlock(j)+1:end, 1:istartingBlock(j)) = updateA(jstartingBlock(j)+1:end, 1:istartingBlock(j))*aV;
            updateB(jstartingBlock(j)+1:end, 1:istartingBlock(j)) = updateB(jstartingBlock(j)+1:end, 1:istartingBlock(j))*aV;
    end
    [bU, bSigma, bV] = svd(updateB(1:jstartingBlock(j),istartingBlock(j+1)+1:istartingBlock(j)));
%     Bp = double(bU*updateB(1:jstartingBlock(j),1:istartingBlock(j)));
    bU = bU^(-1);
    for k = 1:size(bSigma,1)
%         indexB(k) = all(Bp(k,:) == 0);
        indexB(k) = all(bSigma(k, :) == 0);
    end
    bU = [bU(indexB,:); bU(~indexB,:)];
    rj = rank(bSigma);  
    P{j} = [...
        bU, zeros(size(bU,2), size(A,2)-size(bU,2));
        zeros(size(A, 2)-size(bU, 2), size(bU,2)), eye(size(A,2)-size(bU,2))];
    Q{j} = [...
        aV, zeros(size(aV,2), size(A,2)-size(aV,2));
        zeros(size(A, 2)-size(aV, 2), size(aV,2)), eye(size(A,2)-size(aV,2))];
        jstartingBlock(j+1) =  jstartingBlock(j) - rj;
        updateB(1:jstartingBlock(j), 1:istartingBlock(j)) = bU* updateB(1:jstartingBlock(j), 1:istartingBlock(j));
        updateA(1:jstartingBlock(j), 1:istartingBlock(j)) = bU* updateA(1:jstartingBlock(j), 1:istartingBlock(j));
        n(j+1) = n(j) - sj(j+1);
        j = j+1;
    end
    U{j-1} = bU;
    V{j-1} = aV;
    rankTopBlock = ri;
    nullSize = [0, flip(sj(2:end))];
end







