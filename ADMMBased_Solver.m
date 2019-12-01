function [X,Conv,Nu,L1,Fn,out]= ADMMBased_Solver(Y, F,W_ori,lambda_lr,lambda_sp,mu)


norm_y = sqrt(mean(Y(:).^2));
Y = Y./norm_y;
F = F./norm_y;
MaxIter =1500;
Conv = inf(1,MaxIter);
tao = 1e-5;
[L, K] = size(Y);
m = size(F, 2);

Finv = (F'*F +  3*eye(m))^-1;
W = Finv*F'*Y;

%Initialization of auxiliary matrices V1,V2,V3,V4
V1 = F*W; 
V2 = W;
V3 = W;
V4 = W;

%Initialization of Lagranfe Multipliers
D1 = V1*0;
D2 = V2*0;
D3 = V3*0;
D4 = V4*0;

%current iteration number
i = 1;

%primal residual 
res_p = inf;

%dual residual
res_d = inf;

% ---------------------------------------------
 lambda1 = 0;lambda2 = 3;
%---------------------------------------------
%  ADMM iterations
%---------------------------------------------
% Configurate Waitbar
h = waitbar(0,'1','Name','LRSpILT DOSY','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0);
fprintf('Iteration\t\tPrimal Residual\t\tDual Residual\t\tConvergent error\n');

while (i <= MaxIter)
    
% Check for Cancel button press
    if getappdata(h,'canceling')
        disp('LRSpILT manually aborts in the process !');
        break
    end
%
    if mod(i, 10) == 1
        V10 = V1;
        V20 = V2;
        V30 = V3;
        V40 = V4;
    end

    if i == 800
        lambda1 =lambda_lr; lambda2 = lambda_sp;  
    end
    
    Wprev = W;
    
    W =  Finv*(F'*(V1 + D1) + (V2 + D2) + (V3 + D3) + (V4 + D4));
    V1 = 1/(1+mu)*(Y + mu*(F*W - D1));
    
    [u s v] = svd(W-D2,'econ');
    ds = diag(s);
    V2 = u*diag(max(abs(ds) - (lambda1/mu)*(1./(abs(ds)+ eps)),zeros(size(ds))))*v';%*(1./(abs(ds)+ eps))
    V3 = sign(W-D3).*max( abs(W-D3) - ...
                     (lambda2/mu),0);%*(1./abs(W-D3 + eps))
                 
    V4 = max(W - D4, 0);
    
    %update D
    D1 = D1 - F*W + V1;
    D2 = D2 - W + V2;
    D3 = D3 - W + V3;
    D4 = D4 - W + V4;
    
    if mod(i, 10) == 1
%         primal residual
        res_p = norm([V1; V2; V3; V4] - [F*W; W; W; W], 'fro');
%         dual residual
        res_d = norm([V1; V2; V3; V4] - [V10; V20; V30; V40], 'fro');
     
        if res_p > 10*res_d %%
            mu = mu*2;
            D1 = D1/2;
            D2 = D2/2;
            D3 = D3/2;
            D4 = D4/2;
        elseif res_d > 10*res_p
            mu = mu/2;
            D1 = D1*2;
            D2 = D2*2;
            D3 = D3*2;
            D4 = D4*2;
       end
    end
    

    Conv(i) = norm(W-Wprev,'fro')/norm(W,'fro');
    [Nu(i),L1(i),Fn(i),out(i)] = Obj( lambda_lr,lambda_sp,F,W,Y );
    if  (Conv(i)<tao && i>0.2*MaxIter)
        disp('-------------- LRSpILT converges !-----------------');
        break ; 
    end 
    waitbar(i/MaxIter,h,sprintf('LRSpILT processing...%0.1f%% ',i/MaxIter*100));
    fprintf('\t%d\t\t\t\t%0.5f\t\t\t\t%0.5f\t\t\t\t%0.8f\n',i,res_p,res_d,Conv(i));
    i = i + 1;
end
if i > MaxIter
     disp('       LRSpILT iteration reaches maximum and exists !');
end    
delete(h);

 
X = W.*(W>=0);


end
