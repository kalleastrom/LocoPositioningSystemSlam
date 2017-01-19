function [xopt,yopt,resVec]=bundleL1IteratedToaMics2DSounds2D(D, I,J, rec, send, lambda, maxIterOuter, maxIterInner)
%Works using bundleToaMics2DSounds3D. See help file there for details of
%vaiables
%
% INPUT
% maxIterOuter  Maximum iterations for the outer loops, where in each loop
%               you call the function  bundleToaMics2DSounds3D.
% maxIterInner  Maximum nuber of iterations inside each call of bundleToaMics2DSounds3D
% OUTPUT
%
%WARNING: Only has a chance to work with noisy measurements. Even then,
%problems with weights --> inf might occur.

resPre=calcres(D,I,J,rec,send);
resSumAbsPre=sum(abs(resPre));
xopt=rec;
yopt=send;
for ii=1:maxIterOuter
    
    
    
    weights=1./abs(resPre);
    maxLimit=min(weights)*10^4;
    weights(weights> (maxLimit))=maxLimit; %see to it that we do nota hve wieghts approaching zero
    
    
    %maxminRel=max(weights)/min(weights); %for warning debug only
    if sum(~isfinite(weights))>0  %for debug
       keyboard 
    end
    
    weights=weights/median(weights);  %for stability
    lambda=resSumAbsPre/100; %/10; %to go more towards gradient descent for big residual. Choice of increasing function of residual sum is kind of arbitrary. SHould perhaps be adaptive inside the l2-biundler instead
    [recPost,sendPost]=bundleToaMics2DSounds2D(D,I,J,xopt,yopt,lambda,0,weights,maxIterInner);
    
    
    resPost=calcres(D,I,J,recPost,sendPost);
    
    resSumAbsPost=sum(abs(resPost));
   
    
    if resSumAbsPre < resSumAbsPost
        %baaayyd!
       disp(['L1-opt stalled after' num2str(ii) ' iterations'])
       break
    else
        xopt=recPost;
        yopt=sendPost;
        resVec=recPost;
    end
    
     resPre=resPost;
     resSumAbsPre=resSumAbsPost;
end


end


function res=calcres(D,I,J,rec,send)

V = rec(:,I)-send(:,J); %Just the differences of the positions r_i-s_j, for all combinations
dd = sqrt(sum(V.^2,1))'; %The "current" distances

res = dd-D ;

if nargin >=6 %If we have weights
    nbrMeas=length(D);
    R=spdiags(sqrt(weights),0,nbrMeas,nbrMeas); %so that R*R=W the weighting matrix, W which has weights on the diagonal and weight_i is multiplicated with res_i.^1
    res=R*res;
    
end

end