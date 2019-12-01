function [X2] = rvdt(hf,X)
% RecoVe the Size of DataseT by deleting the noise region. 
%  INPUT
% X:  the solution
% hf:  the status
%  OUTPUT
% X2: the solution with orginal size
X2=zeros(size(X,1),size(hf,2));
ss=0;
for i=1:size(hf,2)
    if(hf(i)>0.5);
        ss=ss+1;
      X2(:,i)=X(:,ss);
    end
end
end

