function [C, U, R] = cur_deim(A, k)

%CUR_DEIM  DEIM incurred CUR decomposition
% function [icol, irow, M] = cur_deim(A, m)
%
% C = A(:,icol);  R = A(irow,:)
%
% Reference: Embree and Sorensen, 2016
%


if nargin < 2 || isempty(k), k = 2; end

[U, ~, V] = svds(full(A),k);

for i = 1:k
  [~, irow(i)] = max(abs(U(:,i)));
  [~, icol(i)] = max(abs(V(:,i)));
  U(:,i+1:k) = U(:,i+1:k) - U(:,1:i) * (U(irow,1:i) \ U(irow,i+1:k));
  V(:,i+1:k) = V(:,i+1:k) - V(:,1:i) * (V(icol,1:i) \ V(icol,i+1:k));
end
C=A(:,icol);
R=A(irow,:);
U = C \ (A / R);

