function adfgls=adfgls(y,d,k,m,c)
%ADFGLS  Augmented Dickey-Fuller GLS Test (ERS 1996)
%
% ADF-GLS test as specified by Elliott, Rothenberg and Stock.
% The function performs the GLS detrending (demeaning) for an 
% univariate time series under local-to-unity framewok:
%
%  y"(t)-=y(t)-psi'*z(t)
%
% The transformed data is then used to run a standard ADF regression 
% to calculate the test.
%   
%       y"(t) = c + d*t + a*y(t-1)+ b_1*(y"(t-1)-y"(t-2))
%                                 + b_2*(y"(t-3)-y"(t-3))
%                                 + ... 
%                                 + b_k*(y"(t-k)-y"(t-k-1))
%                                 + e(t),
%
% The function allows for detrending via OLS as well as performing
% the standard ADF test.
%
%==========================================================================
%               adfgls(Y ,D ,K ,M , C)
%==========================================================================
%  INPUTS
%
%   Y          - Univariate time series, specified as a single vector.
%
%   D          - Deterministic component for the series:
%                   1 ----> no deterministic component
%                   2 ----> intercept
%                   3 ----> intercept and time trend
%                Note that detrending will only work under models 2 or 3.
%
%   K          - Number of lags. Minimum is 1, otherwise no regression
%                  can be performed.
%
%   M          - Method used to perform the test:
%                   1 ----> GLS detrending
%                   2 ----> OLS detrending
%                   3 ----> standard ADF test
%
%   C        -  Local-to-unity component where C is the point where 
%               detrending will be performed. The value must be negative.
%               Note that this parameter won't be necessary when M equals
%               2 or 3 and can be set to any number to allow the function
%               to run.
%  
% Note: The function is highly modifiable and can be tunned to return
%       additional outputs.
%==========================================================================
%==========================================================================
% EXAMPLE:
%   
%    % Test the CPI series from Nelson and Ploser (1982) for a unit root 
%       using and intercept and time trend  based on GLS detrending. 
%       c=-13.5 will be used to perform the detrending step:
%       
%      load Data_NelsonPlosser
%       CPI = DataTable.CPI;
%       c=-13.5;
%       test=adfgls(CPI,3,1,1,c)
%      
%
%==========================================================================
%
% References:
%   Elliott, G., Rothenberg, T. J., & James, H. (1996). Stock,
%       1996,“Efficient tests for an autoregressive unit root,”.
%       Econometrica, 64(4), 813-836.
%
%==========================================================================



y=y(:);
T=size(y,1);
x=[];

switch m
    
    case 1
        switch d
            case 2
               z=ones(T,1);
               zgls=[1; (ones(T-1,1)-1-(c/T) )];
            case 3
               z=[ones(T,1) (1:T)'];
               zgls=[ 1 1; (z(2:end,1)-(1+c/T)*z(1:end-1,1))...
                   (z(2:end,2)-(1+c/T)*z(1:end-1,2)) ];
            otherwise
                error('Unexpected option ')
        end
        ygls=[ y(1); y(2:end)-(1+c/T)*y(1:end-1)];
        psi=(zgls'*zgls)^-1*(zgls'*ygls);
        yadf=y-sum(psi'.*z,2);
        index=1;        
    case 2
        switch d
            case 2
               z=ones(T,1);
            case 3
               z=[ones(T,1) (1:T)'];
            otherwise
                error('Unexpected option ')
        end

                psi=(z'*z)^-1*(z'*y);
                yadf=y-sum(psi'.*z,2);
                index=1;
        
    case 3
        yadf=y;
        switch d
            case 1
            case 2
                x= ones(T-k,1);
            case 3
                x=[ ones(T-k,1) (1:T-k)' ];
            otherwise
                error('Unexpected option ')
        end
             index=d;
      otherwise
        error('Unexpected option ')  
       
end

  x=[ x yadf(k:end-1)];
    switch k
        case 1
        otherwise 
         for i=1:k-1
             x=[ x ( yadf(k+1-i:end-i)-yadf(k-i:end-1-i) )];
         end
    end


    phi=(x'*x)^-1*(x'*yadf(k+1:end));
    e=yadf(k+1:end)-x*phi;
    s2=(e'*e)/(size(x,1)-size(x,2));

    V=s2*(x'*x)^-1;
    adfgls=(phi(index)-1)/V(index,index)^0.5;
  
end
        
            
           
            

        
        













