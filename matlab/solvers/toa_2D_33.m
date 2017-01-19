function [Rc,Sc]=toa_2D_33_v5(d)
% [Rc,Sc]=toa_2D_33_v4(d)
% Calculate possible placements of receivers R and senders S
% in the plane so that the distance between R(:,i) and S(:,j) is d(i,j)
% Input: 
%    d - 3x3 matrix
% Output:
%    Rc - cell array with possible receiver configurations R, i.e. R = Rc{k}.
%    Sc - cell array with possible sender configurations, i.e. S = S{k}
%    (Each R is of size 2x3 and S of size 2x3
%

%% Go from data, d, to, coefficients B1, B2, B3 of the three
% equations

Cm = [-1 1 0;-1 0 1];
dc = Cm*(d.^2)*Cm';
d11c = d(1,1).^2;
di1c = Cm*(d(:,1).^2);
d1jc = (d(1,:).^2)*Cm';

Rt = eye(2);
St = dc/(-2);

Rtt = [zeros(2,1) Rt];
Stt = [zeros(2,1) St];

K1 = di1c(1);
K2 = di1c(2);
A1 = [-1 4 2*K2 2*K1 K1*K2];
A2 = A1*d11c;
%A3 = [-2 2 2 K2 K1];
B1 = [2 -2 -2 A2(1) -K2 A2(2) -K1 A2(3:5)];

A1p = [-1 0 0 4 0 2*K2 2*K1 K1*K2];
A5 = A1p*d1jc(1);
A6 = A1p*d1jc(2);
A7 = -[zeros(1,4) -2*prod(St(:,1)) 2*St(2,1)^2 2*St(1,1)^2 K2*St(1,1)^2+K1*St(2,1)^2];
A8 = -[zeros(1,4) -2*prod(St(:,2)) 2*St(2,2)^2 2*St(1,2)^2 K2*St(1,2)^2+K1*St(2,2)^2];
A9 = -[0 (-2*St(2,1)) (-2*St(1,1)) (4*(St(1,1)+St(2,1))) 0 2*K2*St(1,1) 2*K1*St(2,1) 0];
A10 = -[0 (-2*St(2,2)) (-2*St(1,2)) (4*(St(1,2)+St(2,2))) 0 2*K2*St(1,2) 2*K1*St(2,2) 0];
B2 = A5+A7+A9;
B3 = A6+A8+A10;

Bs = {B1',B2',B3'}; %These are the coefficients
Bn = [6 , 12 , 12 ];

ndim = 9; % There is only 8 solutions, but when multiplying the
% equations with det(H) we get one more solution??

%%

I = [ 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 5 , 5 , 5 , 5 , 5 , 5 , 5 , 5 , 5 , 5 , 6 , 6 , 6 , 6 , 6 , 6 , 6 , 6 , 6 , 6 , 7 , 7 , 7 , 7 , 7 , 7 , 7 , 7 , 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 , 9 , 9 , 9 , 9 , 9 , 9 , 9 , 9 , 10 , 10 , 10 , 10 , 10 , 10 , 10 , 10 , 11 , 11 , 11 , 11 , 11 , 11 , 11 , 11 , 12 , 12 , 12 , 12 , 12 , 12 , 12 , 12 , 13 , 13 , 13 , 13 , 13 , 13 , 13 , 13 , 14 , 14 , 14 , 14 , 14 , 14 , 14 , 14 , 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 , 16 , 16 , 16 , 16 , 16 , 16 , 16 , 16 , 17 , 17 , 17 , 17 , 17 , 17 , 17 , 17 , 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 , 19 , 19 , 19 , 19 , 19 , 19 , 19 , 19 , 20 , 20 , 20 , 20 , 20 , 20 , 20 , 20 , 21 , 21 , 21 , 21 , 21 , 21 , 21 , 21 , 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 , 23 , 23 , 23 , 23 , 23 , 23 , 23 , 23 , 24 , 24 , 24 , 24 , 24 , 24 , 24 , 24 , 25 , 25 , 25 , 25 , 25 , 25 , 25 , 25 , 26 , 26 , 26 , 26 , 26 , 26 , 26 , 26 , 27 , 27 , 27 , 27 , 27 , 27 , 27 , 27 , 28 , 28 , 28 , 28 , 28 , 28 , 28 , 28 , 29 , 29 , 29 , 29 , 29 , 29 , 29 , 29 , 30 , 30 , 30 , 30 , 30 , 30 , 30 , 30 ]';
J = [ 3 , 5 , 6 , 10 , 15 , 16 , 17 , 25 , 26 , 33 , 12 , 15 , 16 , 21 , 24 , 25 , 26 , 32 , 33 , 37 , 5 , 7 , 8 , 12 , 18 , 19 , 20 , 28 , 29 , 35 , 15 , 18 , 19 , 22 , 27 , 28 , 29 , 34 , 35 , 38 , 16 , 19 , 20 , 23 , 28 , 29 , 30 , 35 , 36 , 39 , 25 , 28 , 29 , 31 , 34 , 35 , 36 , 38 , 39 , 40 , 1 , 2 , 3 , 5 , 12 , 15 , 16 , 25 , 9 , 11 , 12 , 15 , 22 , 24 , 25 , 32 , 10 , 12 , 13 , 16 , 23 , 25 , 26 , 33 , 21 , 22 , 23 , 25 , 31 , 32 , 33 , 37 , 2 , 4 , 5 , 7 , 15 , 18 , 19 , 28 , 11 , 14 , 15 , 18 , 24 , 27 , 28 , 34 , 3 , 5 , 6 , 8 , 16 , 19 , 20 , 29 , 12 , 15 , 16 , 19 , 25 , 28 , 29 , 35 , 22 , 24 , 25 , 28 , 32 , 34 , 35 , 38 , 13 , 16 , 17 , 20 , 26 , 29 , 30 , 36 , 23 , 25 , 26 , 29 , 33 , 35 , 36 , 39 , 31 , 32 , 33 , 35 , 37 , 38 , 39 , 40 , 1 , 2 , 3 , 5 , 12 , 15 , 16 , 25 , 9 , 11 , 12 , 15 , 22 , 24 , 25 , 32 , 10 , 12 , 13 , 16 , 23 , 25 , 26 , 33 , 21 , 22 , 23 , 25 , 31 , 32 , 33 , 37 , 2 , 4 , 5 , 7 , 15 , 18 , 19 , 28 , 11 , 14 , 15 , 18 , 24 , 27 , 28 , 34 , 3 , 5 , 6 , 8 , 16 , 19 , 20 , 29 , 12 , 15 , 16 , 19 , 25 , 28 , 29 , 35 , 22 , 24 , 25 , 28 , 32 , 34 , 35 , 38 , 13 , 16 , 17 , 20 , 26 , 29 , 30 , 36 , 23 , 25 , 26 , 29 , 33 , 35 , 36 , 39 , 31 , 32 , 33 , 35 , 37 , 38 , 39 , 40 ]';
Cnnz = 252;
iss = [ 1 , 61 , 157 ; 60 , 156 , 252 ];
Cm = 30;
Cn = 40;

Cv = zeros(Cnnz,1);
for kk = 1:length(Bs);
    Cv(iss(1,kk):iss(2,kk))=repmat(Bs{kk},Bn(kk),1);
end

C = sparse(I,J,Cv,Cm,Cn);

B = [ 28 , 33 , 34 , 35 , 36 , 37 , 38 , 39 , 40 ];
R = [ 15 , 23 , 24 , 25 , 26 , 31 , 32 ];
E = [ 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 , 16 , 17 , 18 , 19 , 20 , 21 , 22 , 27 , 29 , 30 ];
Eful = E;
Eful(8)=[];

[qq,~]=qr(C(:,Eful));
nrank = 23; % I hard coded the rank here.

CCny = qq(:,(nrank+1):end)'*C(:,[R B]);
%Rny = [ 1 , 2 , 3 , 4 , 5 , 6 , 7 ];
%Bny = [ 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 , 16 ];

brep = [-(CCny(:,1:length(R))\CCny(:,(length(R)+1):end))' eye(ndim)];
Bxny = [1:7 9 13];
actionM = brep(:,Bxny);

[vv,~]=eig(full(actionM'));
tmp = brep'*vv;
% Fix scale so mu*1 = 1
tmp = tmp./(ones(size(brep,2),1)*tmp(end,:));

xsols = tmp(13:15,:);

%% Backtrack
% Go from solutions x to H,b to R,S, for each of the solutions
Rc = cell(1,9);
Sc = cell(1,9);


ok = false(9,1);

for kk=1:size(xsols,2);
    asol = xsols(:,kk);
    if isreal(asol),
        btmp = [asol(2); asol(3)];
        Htmp = [K1+2*asol(2) asol(1);asol(1) K2+2*asol(3)];
        if all(eig(Htmp)>0),
            Ltmp = chol(inv(Htmp));
            Rtmp = (Ltmp')\Rtt;
            Stmp = Ltmp*(Stt+btmp*ones(1,3));
            Rc{kk}=Rtmp;
            Sc{kk}=Stmp;
            ok(kk)=true;
        end;
    end
end
Rc = Rc(ok);
Sc = Sc(ok);


