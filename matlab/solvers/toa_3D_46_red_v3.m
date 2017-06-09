function [sols,data0,sols_short] = toa_3D_46_red_v3(d)
%[sols] = toa_3D_64_red(d)
% Solver for minimal case TOA node calibration
%  6 receivers and 4 senders
%  measurments according to
%    d(i,j) = sqrt(sum( (x(:,i)-y(:,j)).^2 ))
%  Goal of solver is to reconstruct x and y from d.
%  There are 38 solutions counted with multiplicity and
%  complex solutions.
% In:
%    d - 6x4 matrix with distances
%    settings - settings for solver use
%               load toa_3D_46_settings
%    xt -
%    yt -
% Out:
%    sols  - cell array with all real and feasible solvers
%    stats - statistics from solver

D = 3;
m = 6;
n = 4;
d2 = d.^2;
cl = compactionmatrix(6);
cr = compactionmatrix(4);

B = cl*d2*cr'/(-2);

R = B';
S = eye(3);

xt = [zeros(D,1) R];
yt =[zeros(D,1) S];

if 0,
    linvt = r(:,2:6)/R;
    l = inv(linvt)';
    hgt = inv(l'*l);
    bgt = inv(l)*s(:,1);
    zzgt = [hgt(1,1);hgt(2,2);hgt(3,3);hgt(1,2);hgt(1,3);hgt(2,3);bgt];
end

da = d2(1,1);
db = d2(1,:)*cr';
dc = cl*d2(:,1);

H{1} = [1 0 0;0 0 0;0 0 0];
H{2} = [0 0 0;0 1 0;0 0 0];
H{3} = [0 0 0;0 0 0;0 0 1];
H{4} = [0 1 0;1 0 0;0 0 0];
H{5} = [0 0 1;0 0 0;1 0 0];
H{6} = [0 0 0;0 0 1;0 1 0];
b{1} = [1;0;0];
b{2} = [0;1;0];
b{3} = [0;0;1];

%% 5 linear equations of type C
% R(:,k)'*H*R(:,k) -2*b'*R(:,k) = dc(k);

AA = zeros(5,9);
bb = zeros(5,1);
for ii = 1:5;
    for jj = 1:6;
        AA(ii,jj) = R(:,ii)'*H{jj}*R(:,ii);
    end
    for jj = 1:3;
        AA(ii,jj+6) = -2*b{jj}'*R(:,ii);
    end
    bb(ii)=dc(ii);
end

zz0 = AA\bb;
[u,ss,v]=svd(AA);
zzb = v(:,6:9);

if 0,
    %% zz = zz0+zzb*xx;
    xxgt = zzb\(zzgt-zz0);
    xxgt(5)=det(hgt);
end

%%

data0 = [R(:); zzb(:); zz0(:); db(:); da(:)];
sols_short = solver_toa_46_v3(data0);

% find real solutions
sols_real = find_realsol(sols_short);
n_real    = size(sols_real,2);


%% recover C and b
sols_full = repmat(zz0,1,n_real)+zzb*sols_real(1:4,:);
% Each column in sols full contains 9 elements, the six elements
% of the 3x3 symmetric matrix H and the 3 elements of vector b

err =  inf;
kk = 1;
ok = 0;
x = cell(1,0);
y = cell(1,0);
err = zeros(1,0);
%keyboard;

for i = 1:size(sols_full,2);
    
    HH = zeros(3,3);
    for k = 1:6;
        HH = HH + H{k}*sols_full(k,i);
    end
    
    BB = zeros(3,1);
    for k = 1:3;
        BB = BB + b{k}*sols_full(6+k,i);
    end
    
    CC{i} = HH;
    if min(eig(HH))>0,
        l    = chol(inv(CC{i}));
        ok   = 1;
        L{kk} = l;
        BBs{kk} = BB;        
        x{kk} = inv(l')*(xt);
        y{kk} = l*(yt+repmat(BB,1,n));
        dd = sqrt( sum( (kron(ones(1,n),x{kk}) - kron(y{kk},ones(1,m)) ).^2 , 1 ) );
        dd = reshape(dd,m,n);
        % minibundle
        
        err(kk) = norm(d - dd);
        C{kk} = CC{i};
        solsvec{kk} = sols_full(:,i);
        kk = kk+1;
        
        
    else % if not positive-definite, throw away
        %         disp(' C is not positive definite')
        l    = 0;
        %continue;
    end
end

%% Store output in sols
%keyboard;
if ok
    sols.CC = CC;
    sols.L = L;
    sols.b = BBs;
    sols.x = x;
    sols.y = y;
    sols.n_real = n_real;
    sols.n_good = length(L);
    sols.err = err;
else
    sols = [];
    stats.n_good = 0;
    %     stats = [];
end

