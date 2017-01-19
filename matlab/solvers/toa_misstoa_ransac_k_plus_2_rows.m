function [sol,maxnrinl]=toa_misstoa_ransac_k_plus_2_rows(d,sys);
%

if nargin<2,
    sys.ransac_threshold = 0.01;
    sys.ransac_k = 5000;
    sys.rk = 3;
end;

rankplus2 = sys.rk+2;
maxnrinl = 0;
for iii = 1:sys.ransac_k;
    
    d2 = d.^2;
    inl = isfinite(d2);
    [m,n]=size(d2);
    tmprows = randperm(m);
    tmprows = tmprows(1:(sys.rk+2)); % Have to have rank + 2 rows to work with
    okcol = find(all(inl(tmprows,:)));
    
    B = d2(tmprows,okcol);
    
    ntmp = size(B,2);
    tmp2 = randperm(ntmp);
    if ntmp>(sys.rk+2),
        tmp21 = tmp2(1:(sys.rk+1));   % Use rank + 1 columns to esimate basis
        tmp22 = tmp2((sys.rk+2):end); % Use remaining columns for checking
        
        cl = compactionmatrix( sys.rk+2 );
        cr1 = compactionmatrix( sys.rk+1 );
        cr2 = compactionmatrix(size(tmp22,2));
        cr = compactionmatrix(size(tmp2,2));
        
        Btmp = cl*B(:,tmp2)*cr';
        
        B1 = Btmp(:,1:(sys.rk));     % Use rank + 1 columns to esimate basis
        B2 = Btmp(:,(sys.rk+1):end); % Use remaining columns for checking
        
        [u,s,v]=svd(B1);
        u_last = u(:,sys.rk+1);
                
        Imiss = isnan(d);
        okind = find(abs(u_last'*B2) < sys.ransac_threshold );
        inlim = zeros(size(d));
        inlim = inlim-Imiss;
        inlim(tmprows,okcol([tmp21 tmp22(okind)]))=ones(sys.rk+2,sys.rk+1+length(okind));
        nrinl = sys.rk+1+length(okind);
 
        if nrinl>maxnrinl,
            % Wow. We have found a new maximum nr of inliers
            maxnrinl = nrinl;
            % What data should be save
            % inlmatrix (only those on rows,cols are relevant??
            % Low rank approximation
            sol.rows = tmprows;
            sol.cols = okcol([tmp21 tmp22(okind)]);
            sol.row1 = sol.rows(1);
            sol.col1 = sol.cols(1);
            sol.inlmatrix = inlim;
            B = d2(sol.rows,sol.cols);
            [cl,dl]=compactionmatrix(size(B,1));
            [cr,dr]=compactionmatrix(size(B,2));
            Bhat = dl*B*dr';
            Btilde = cl*B*cr';
            [u,s,v]=svd(Btilde);
            s((sys.rk+1):end,:)=zeros(size(s,1)-sys.rk,size(s,2));
            Btilde = u*s*v';
            Bhat(2:end,2:end)=Btilde;
            sol.Bhat = Bhat;
            sol.dl = dl;
            sol.dr = dr;
        end
    end
end

% Should we bundle??
