function [C, U, R] = cur_deim(A, k)

%CUR_DEIM  DEIM incurred CUR decomposition
% function [icol, irow, M] = cur_deim(A, k)
%
% C = A(:,icol);  R = A(irow,:)
%
% Reference: Embree and Sorensen, 2016
% 
% (C) Perfect Gidisu, Michiel Hochstenbach 2020


if nargin < 2 || isempty(k), k = 2; end

[U, ~, V] = svds(A,k);

irow=zeros(1,k);
icol=zeros(1,k);
for j = 1:k
  [~, irow(j)] = max(abs(U(:,j)));
  [~, icol(j)] = max(abs(V(:,j)));
  if j<k
   U(:,j+1) = U(:,j+1) - U(:,1:j) * (U(irow(1:j),1:j) \ U(irow(1:j),j+1));
   V(:,j+1) = V(:,j+1) - V(:,1:j) * (V(icol(1:j),1:j) \ V(icol(1:j),j+1));
  end
end
C=A(:,icol);
R=A(irow,:);
U = C \ (A / R);

