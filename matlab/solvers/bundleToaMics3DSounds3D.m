function [xopt,yopt]=bundleToaMics2DSounds2D(D,I,J,xt,yt,lambda, debug, weights, maxIter )
%Bundling TOA measurements with paris of microphones located on rigid constellations, with
%constnat, unknown distances between the paired microphones. ALl rigs has the same
%unknown inter distance. Futhermore, Sounds are in 3D and Microphones in
%an affine subspace of dimension 2, here assumed to be the x-y plane. Thus,
%before bundling, one has to rotate-translate the cosntellation of
%mics/sounds so that the microphones has coordinates z=0.
%Mic1, xt(:,1)is assumed to be on a rig with Mic2, xt(:,2)
%Mic3, xt(:,3)is assumed to be on a rig with Mic4, xt(:,4)
%And so on...
%What is mimized is the sum of squared residuals.
%Stop chriterion is max number of iterations, or when
%
%INPUT
%   D       #measurements x 1   Vector with the TOA measurements in
%   I       #measurements x 1   Vector with the indices for D, correspinding to
%                               which microphone each measurement is associated
%                               with
%   J       #measurements x 1   Vector with the indices for D, correspinding to
%                               which sound each measurement is associated
%                               with
%   xt      2 x #mics           The initial values for mic positions
%   yt      3 x #sounds         The initial values for sound positions
%   lambda  scalar              Parameter in the LM. If set to zero, LM is
%                               not used, but regular Gauss_Newton with step control (default=0)
%   debug   scalar              Boolean flag for debug info. (default=0)
%   weights #measurements x 1   Weights for the squared residuals. What is
%                               minimized is the sum of
%                               weight_i*(residual_i)^2. A reweighted least
%                               squares. (default = no weights, i.e. the
%                               same as ones). It is assumed all weights
%                               are positive.
%   maxIter scalar              Maximum number of iterations. This, or when 
%                               sum of squared residuals stalls, is stop criterion. (default=10) 
%OUTPUT
%   xt      2 x #mics           The optimized values for mic positions
%   yt      3 x #sounds         The optimized values for sound positions
%
%

if nargin <6
    lambda=0;
end


if nargin<7,
  debug = 0;
end;

if nargin < 9
    maxIter=10;
end

if size(xt,1) ~= 3
    error('Wrong dimesnion of receivers. Read the help for the bundleadjuster')
end

if size(yt,1) ~= 3
    error('Wrong dimesnion of senders. Read the help for the bundleadjuster')
end


for kkk = 1:maxIter;
    %kkk
    if ~exist('weights','var') %we have unweighted least sqares. Messier than setting default weights=ones, but faster
        [res,jac]=calcresandjac(D,I,J,xt,yt);
    else %we have weighted least squares
        [res,jac]=calcresandjac(D,I,J,xt,yt,weights);
    end
    %dz = -(jac\res);
    
    
    if lambda ~=0
        dz = -(jac'*jac+lambda*eye(size(jac,2)))\(jac'*res);
    else
        %         [u,s,v]=svd(full(jac),'econ');
        %         u = u(:,1:(end-2)); %taking pseudoinverse
        %         s = s(1:(end-2),1:(end-2));
        %         v = v(:,1:(end-2));
        %         dz = -v/s*u'*res;
        dz=-pinv(full(jac'*jac))*jac'*res;
    end
    
    

    [xtn,ytn]=updatexy(xt,yt,dz);
    if ~exist('weights','var') %messier than setting default weights= ones, but faster 
        [res2]=calcres(D,I,J,xtn,ytn);
    else
        [res2]=calcres(D,I,J,xtn,ytn,weights);
    end
    
    aa = [norm(res) norm(res+jac*dz) norm(res2)];
    bb = aa;
    bb=bb-bb(2);
    bb = bb/bb(1);
    cc = norm(jac*dz)/norm(res);
    % Check that the error actually gets smaller
    if norm(res)<norm(res2),
        % bad
        % check that the update is sufficiently big
        % otherwise it is maybe just numerical limitations
        if cc>1e-4,
            % OK then it is probably just the linearization that
            % is not good enough for this large step size, decrease
            kkkk = 1;
            while (kkkk<50) && (norm(res)<norm(res2)),
                dz = dz/2;
                [xtn,ytn]=updatexy(xt,yt,dz);
                if ~exist('weights','var') %messier than setting default weights= ones, but faster
                    [res2]=calcres(D,I,J,xtn,ytn);
                else
                    [res2]=calcres(D,I,J,xtn,ytn,weights);
                end
                kkkk = kkkk+1;
            end
        end
    end
    if debug,
        aa = [norm(res) norm(res+jac*dz) norm(res2)];
        bb = aa;
        bb=bb-bb(2);
        bb = bb/bb(1);
        cc = norm(jac*dz)/norm(res);
        keyboard;
        aa
        bb
        cc
    end;
    if norm(res2)<norm(res)
        xt = xtn;
        yt = ytn;
    else
        %disp([ '  stalled after ' num2str(kkk) ' iterations']);
        break
    end
end;

xopt = xt;
yopt = yt;

end

function [res,jac]=calcresandjac(D,I,J,x,y,weights)



nbrMeas = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J); %Just the differences of the positions r_i-s_j, for all combinations
Vt = V';           
dd = sqrt(sum(V.^2,1))'; %The "current" distances 
idd = 1./dd;             %inverse distances
%res = dd-D;              %residual distances, which part of what is to be minimized (sum ov squares)
                        
%% THis part is for the TOA-measurements residuals
%The unknowns are ordered as r_1x, r_1y, r_2x, r_2_y, ..., r_mx,r_my, s_1x, s1_y, s_1_z, ...., s_nx,s_ny, s_nz and
% wehre r nad s are receivers (mics) and senders(sounds)

%Each of the residuals is affected by exactly one x-coord of a receover and
%one x-coord of a sender. same for the y-coord and z-coord, ofc. Here, we
%split their partial påverkan 
II = (1:length(I))';     %Indices of mic indices.OR, rahter, indices of the equations
JJrecX = (I-1)*3+1;         %Indices of the unknowns which has somethingto do with receiver x-coords (Input: residual/equation number. Output: which of the unkonwn x coordinate or the receivers that affects that residual)
JJrecY = (I-1)*3+2;         %Indices of the unknowns which has somethingto do with receiver y-coords
JJrecZ = (I-1)*3+3;         %Indices of the unknowns which has somethingto do with receiver y-coords
JJSendX = (J-1)*3+1+3*m;      %Indices of the unknowns which has somethingto do with sender x-coords
JJSendY = (J-1)*3+2+3*m;    %Indices of the unknowns which has somethingto do with sender y-coords
JJSendZ = (J-1)*3+3+3*m;    %Indices of the unknowns which has somethingto do with sender y-coords

VVRecX = idd.*Vt(:,1); %  THe values of derivatives of  receivers x-pos.
VVRecY = idd.*Vt(:,2); %  THe values of derivatives of  receivers y-pos.
VVRecZ = idd.*Vt(:,3); %  THe values of derivatives of  receivers z-pos.
%VV3 = idd.*Vt(:,3); %  THe values of derivatives of  receivers z-pos.   WROONG! RECEOVERS DO NOT HAVE THE POSITIONS!
VVSendX = -idd.*Vt(:,1);    %The values of derivatives of senders x-pos.
VVSendY = -idd.*Vt(:,2);    %The values of derivatives of senders y-pos.
VVSendZ = -idd.*Vt(:,3);    %The values of derivatives of senders y-pos.
%VVSendZ = -idd.*Vt(:,3);    %The values of derivatives of senders z-pos.

%% THis part is for the rig distance equations residuals
%there has to be a nicer, Kalle-esque, way of doing the index IIExtra and JJExtra for the sparse
%jacobian. But, for now, this works. Feel free to change.

%each equation is of the form
%||r_i-2-r_i-1||-||r_i-r_i+1||=residual_current i=3, 5, 7...  
%so 8 unknowns in each equation
% nExtraEqs=(m/2 -1); 
% IIExtra=zeros(nExtraEqs*8,1);
% indExtra= ((length(I) + 1):  (length(I) + nExtraEqs)); %these are the indices of the new equations, uniquely
% for kk=1:nExtraEqs
%      IIExtra((1:8) +(kk-1)*8)=repmat(indExtra(kk),8,1) ;    %there is trwelve non-zero entries in each of the extra equations (rows) for the jacobian000
% end
% 
% 
% JJExtra=zeros(8*nExtraEqs,1);
% 
% for kk=1: nExtraEqs  %current extra equation
%     JJExtra((1:8) + (kk-1)*8)=((1:8) +(kk-1)*4);
% end
% 
% VVExtra=zeros(nExtraEqs*8,1);
% ind=1:2:m-1;
% VRecDiff=x(:,ind)-x(:,ind+1);
% rigDists=sqrt(sum(VRecDiff.^2,1))';
% for kk=1:nExtraEqs
%     VVExtra((1:8) + (kk-1)*8)=[ VRecDiff(:,kk)/rigDists(kk) ;  -VRecDiff(:,kk)/rigDists(kk) ; -VRecDiff(:,kk+1)/rigDists(kk+1) ; VRecDiff(:,kk+1)/rigDists(kk+1) ] ;
% end
% 
% VVExtra=VVExtra*beta; %lambda is paramter to wieight up residuals of the micrig distance equations.
% 
%             %residual distances, which part of what is to be minimized (sum ov squares)

%% And here the jacobian and residuals are made

jac = sparse([II;II;II;II;II;II],[JJrecX;JJrecY;JJrecZ;JJSendX;JJSendY;JJSendZ],[VVRecX;VVRecY;VVRecZ;VVSendX;VVSendY;VVSendZ],nbrMeas,3*m+3*n);

res = dd-D ;

if nargin >=6 %If we have weights
    R=spdiags(sqrt(weights),0,nbrMeas,nbrMeas); %so that R*R=W the weighting matrix, W which has weights on the diagonal and weight_i is multiplicated with res_i.^1
    jac=R*jac;
    res=R*res;
end
end

function res=calcres(D,I,J,x,y,weights)

V = x(:,I)-y(:,J); %Just the differences of the positions r_i-s_j, for all combinations
dd = sqrt(sum(V.^2,1))'; %The "current" distances

res = dd-D ;

if nargin >=6 %If we have weights
    nbrMeas=length(D);
    R=spdiags(sqrt(weights),0,nbrMeas,nbrMeas); %so that R*R=W the weighting matrix, W which has weights on the diagonal and weight_i is multiplicated with res_i.^1
    res=R*res;
    
end

end

function [xny,yny]=updatexy(x,y,dz)

m = size(x,2);
n = size(y,2);
dz1 = dz(1:(3*m));
dz2 = dz((3*m+1):end);
xny = x + reshape(dz1,3,m);
yny = y + reshape(dz2,3,n);
end


