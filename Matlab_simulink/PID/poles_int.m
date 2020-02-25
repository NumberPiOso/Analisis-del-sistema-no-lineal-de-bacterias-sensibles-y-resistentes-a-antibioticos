function [K,L] = poles_int(S,P)

[p,m] = size(S.D); n=size(S.A); n=n(1);

if length(P) ~= (n+p)

    disp(['Debe asignar (n+p)= ' num2str(n+p) ' polos'])

    K = []; L=[];

    return

end

A = [S.A zeros(n,p); -S.C eye(p)];

B = [S.B; zeros(p,m)];

if p+m==2

    K1 = acker(A,B,P);

else

    K1 = place(A,B,P);

end

K = K1(:,1:n); L = K1(:,n+1:n+p);