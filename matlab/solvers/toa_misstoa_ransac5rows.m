function [sol,maxnrinl]=toa_misstoa_ransac5rows(d,sys);
%

if nargin<2,
    sys.ransac_threshold = 0.01;
    sys.ransac_k = 5000;
end;

maxnrinl = 0;
for iii = 1:sys.ransac_k;
    
    
    %%
    d2 = d.^2;
    inl = isfinite(d2);
    [m,n]=size(d2);
    tmprows = randperm(m);
    tmprows = tmprows(1:5);
    okcol = find(all(inl(tmprows,:)));
    
    %%
    B = d2(tmprows,okcol);
    %I = Ioutl(tmprows,okcol);
    
    ntmp = size(B,2);
    tmp2 = randperm(ntmp);
    if ntmp>5,
        tmp21 = tmp2(1:4);
        tmp22 = tmp2(5:end);
        
        cl = compactionmatrix(5);
        cr1 = compactionmatrix(4);
        cr2 = compactionmatrix(size(tmp22,2));
        cr = compactionmatrix(size(tmp2,2));
        
        Btmp = cl*B(:,tmp2)*cr';
        
        B1 = Btmp(:,1:3);
        B2 = Btmp(:,4:end);
        
        [u,s,v]=svd(B1);
        u4 = u(:,4);
        
        if 0,
            abs(u4'*B2)
            I(:,tmp2)
        end;
        
        Imiss = isnan(d);
        okind = find(abs(u4'*B2) < sys.ransac_threshold );
        inlim = zeros(size(d));
        inlim = inlim-Imiss;
        inlim(tmprows,okcol([tmp21 tmp22(okind)]))=ones(5,4+length(okind));
        nrinl = 4+length(okind);
        
        
        
        
        if nrinl>maxnrinl,
            
            if 0,
                okind = find(abs(u4'*B2) < sys.ransac_threshold );
                inlim = zeros(size(d));
                inlim = inlim-Imiss;
                inlim(tmprows,okcol([tmp21 tmp22(okind)]))=ones(5,4+length(okind));
                inlimgt = zeros(size(d));
                inlimgt = inlimgt-Imiss;
                inlimgt = inlimgt+Iinl;
                
                figure(1);
                subplot(2,1,1);
                colormap(gray);
                imagesc(inlim);
                subplot(2,1,2);
                colormap(gray);
                imagesc(inlimgt);
                
            end
            
            
            maxnrinl = nrinl;
            % What data should be save
            % rows
            % cols
            % row1
            % col1
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
            s(4:end,:)=zeros(size(s,1)-3,size(s,2));
            Btilde = u*s*v';
            Bhat(2:end,2:end)=Btilde;
            sol.Bhat = Bhat;
            sol.dl = dl;
            sol.dr = dr;
        end
    end
end

% Should we bundle??
