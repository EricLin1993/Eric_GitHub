function [Nu,L1,Fn,out ] = Obj( lambdalr,lambdasp,K,X,Y )
%
   [~,s,~]=svd(X,'econ');
   Nu = sum(diag(s));
   L1 = norm(X(:),1);
   Fn = norm(K*X-Y,'fro');
   out = 1/2*Fn^2 + lambdalr*Nu + lambdasp*L1;

end

