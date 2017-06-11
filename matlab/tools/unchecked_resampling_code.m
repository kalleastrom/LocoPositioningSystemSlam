f0 = dcalc(5,:);
f1 = dcalcgt(5,:);

if 0,
    tmp = find(isnan(f1));
    f1(tmp) = 2.2*ones(1,length(tmp));
end

thresh = 0.001;
a = 2;
%tt = [-10 10];
%tt = tt(1):tt(end);

x = 1:length(f0);
xmid = x(100+1:end-100);
xmid = x(10+1:end-10);
f0t = interp1d(f0,xmid,a);
D = diag([1 1/1000 1/1000]); % obs! what is D?

%
z0 =[0;length(f1)/length(f0);1];

%
[f1t,f1td] = interp1d_with_derivative(f1,z0(2)*xmid+z0(1),a); % obs! Where is the amplitude?
f1t = z0(3)*f1t; % add the amplitude
f1td = z0(3)*f1td;
J(:,1) = f1td';
J(:,2) = (f1td.*xmid)';
J(:,3) = f1t'; % obsobs!
res = -(f0t'-f1t');
%selrows = 1:250;
selrows = 1:length(res);
selrows = 1:550;
%selrows = [1:250 630:length(res)];
dz = -D*((J(selrows,:)*D)\res(selrows));
%dz = -D*((J(selrows,:)*D)\res(selrows));
maxnorm = 20;
if abs(dz(1))>maxnorm,
    dz = dz*(maxnorm/abs(dz(1)));
end
znew = z0+dz;
[f1tnew,~] = interp1d_with_derivative(f1,znew(2)*xmid+znew(1),a); % obs! add paranthesis: znew(2)*(xmid+znew(1))?
f1tnew = znew(3)*f1tnew; % add the amplitude
[norm(res) norm(res+J*dz) norm(f0t-f1tnew)]
% if norm(resnew) > norm(res) d� minska steget.
z0 = znew;
figure(2); clf;
plot(f0t); hold on; plot(f1t);

if 0,
    dz = [0;0;1];
    litet = 0.00001;
    z1 = z0+dz*litet;
    [f1t,f1td] = interp1d_with_derivative(f1,z0(2)*xmid+z0(1),a); % obs! Where is the amplitude?
    f1t = z0(3)*f1t; % add the amplitude
    f1td = z0(3)*f1td;
    J(:,1) = f1td';
    J(:,2) = (f1td.*xmid)';
    J(:,3) = f1t'; % obsobs!
    res0 = -(f0t'-f1t');
    [f1t,f1td] = interp1d_with_derivative(f1,z1(2)*xmid+z1(1),a); % obs! Where is the amplitude?
    f1t = z1(3)*f1t; % add the amplitude
    f1td = z1(3)*f1td;
    J(:,1) = f1td';
    J(:,2) = (f1td.*xmid)';
    J(:,3) = f1t'; % obsobs!
    res1 = -(f0t'-f1t');
    [(res1-res0)/litet J*dz]
end


%%
clear dcalcgtfix
x = 1:length(f0);
xall = x(1:end);
for k = 1:6;
    f1 = dcalcgt(k,:);
    %[f1t,f1td] = interp1d_with_derivative(f1,z0(2)*xall+z0(1),a); % obs! Where is the amplitude?
    [f1t,f1td,okpos]=interp1d_with_derivative_careful(f1,z0(2)*xall+z0(1),a);
    dcalcgtfix(k,:) = f1t;
end

gts = zeros(3,size(stmp,2));
for k = 1:3;
    f1 = data.GTs(k,:);
    %[f1t,f1td] = interp1d_with_derivative(f1,z0(2)*xall+z0(1),a); % obs! Where is the amplitude?
    [f1t,f1td,okpos]=interp1d_with_derivative_careful(f1,z0(2)*xall+z0(1),a);
    gts(k,:) = f1t;
end
%
if 0,
    figure(10); clf;
    plot(dcalcgtfix');
end
gtsok = repmat(okpos',3,1);
%st = 36; sl = 555; gtsok(:,st:sl)=ones(3,(sl-st+1));
data.GTs_resamp = gts;
data.GTs_resamp_ok = gtsok;
data.sopt = stmp;
data.ropt = rtmp;
data.z0 = z0;

if 0,
    cd ../data/data
    save data5 data
    cd ../../matlab
end

figure(11);
subplot(2,1,1);
imagesc(min(dcalc,3));
subplot(2,1,2);
imagesc(min(dcalcgtfix,3));

%%
mm = mean(dcalc(:)-dcalcgtfix(:));
figure(13); clf;
plot(dcalcgtfix(:),dcalc(:)-dcalcgtfix(:)-mm,'.');

figure(14); clf;
plot(dcalc');
hold on
plot(dcalcgtfix')
plot(d')



