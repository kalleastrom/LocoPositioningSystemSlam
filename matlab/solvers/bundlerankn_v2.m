function [U,V,res] = bundlerankn_v2(X,W,U0,V0)

lille = 1e-7;
[N,rk] = size(U0);
M = size(V0,2);

res0 = U0*V0-X;
res0 = res0(W~=0);

resgrad = calcresgrad_rankn_v2(W,U0,V0);

x_gn = -resgrad\res0;



U_gn = [zeros(rk,rk);reshape(x_gn(1:rk*(N-rk)),rk,(N-rk))'];
V_gn = reshape(x_gn(rk*(N-rk)+1:end),rk,M);

U = U0 + U_gn;
V = V0 + V_gn;

res = U*V-X;
res = res(W~=0);

while norm(res)>(lille + norm(res0)),
    U_gn = U_gn/10;
    V_gn = V_gn/10;
    U = U0 + U_gn;
    V = V0 + V_gn;

    res = U*V-X;
    res = res(W~=0);
end


