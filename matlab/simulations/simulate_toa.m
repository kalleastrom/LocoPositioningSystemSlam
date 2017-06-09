function [d,x,y] = simulate_toa_1(m,n,options)
%[d,m,n] = simulate_toa_1(m,n,option)

if nargin<3,
    options.twod=1;
    options.int = 0;
    options.normalise=1;
    options.linear=0;
    options.linearmotion=0;
end

if options.twod,
    % flat
    if options.int,
        x = round(10*randn(2,m));
        y = round(10*randn(2,n));
        if options.normalise,
            x(1,1)=0;
            x(2,1)=0;
            x(2,2)=0;
        end
        if options.linearmotion,
            x(2,1)=0;
            x(2,2)=0;
            x(2,3)=0;
        end
    else
        x = randn(2,m);
        y = randn(2,n);
        if options.normalise,
            x(1,1)=0;
            x(2,1)=0;
            x(2,2)=0;
        end
        if options.linearmotion,
            x(2,:)=zeros(1,size(x,2));
        end
    end;
    d = sqrt( sum( (kron(ones(1,n),x) - kron(y,ones(1,m)) ).^2 , 1 ) );
    d = reshape(d,m,n);
else
    % 3D
    if options.int,
        x = round(1*randn(3,m));
        y = round(1*randn(3,n));
        if options.normalise,
            x(1,1)=0;
            x(2,1)=0;
            x(3,1)=0;
            x(2,2)=0;
            x(3,2)=0;
            x(3,3)=0;
            if x(1,2)<0,
                x(1,2)=-x(1,2);
            end
            if x(2,3)<0,
                x(2,3)=-x(2,3);
            end
            if x(3,4)<0,
                x(3,4)=-x(3,4);
            end
        end
    else
        x = randn(3,m);
        y = randn(3,n);
        if options.normalise,
            x(1,1)=0;
            x(2,1)=0;
            x(3,1)=0;
            x(2,2)=0;
            x(3,2)=0;
            x(3,3)=0;
            if x(1,2)<0,
                x(1,2)=-x(1,2);
            end
            if x(2,3)<0,
                x(2,3)=-x(2,3);
            end
            if m>3,
                if x(3,4)<0,
                    x(3,4)=-x(3,4);
                end
            end
            
        end
    end;
    if options.planarmotion,
        x(3,:)=zeros(1,size(x,2));
    end
    d = sqrt( sum( (kron(ones(1,n),x) - kron(y,ones(1,m)) ).^2 , 1 ) );
    d = reshape(d,m,n);
    
end