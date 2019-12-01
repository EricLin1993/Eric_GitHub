function [rdt,hf,PL] = rsdt(cdt,thre)
% Reduce the Size of DataseT by deleting the noise region. 
%    INPUT
% cdt: the aligned dataset
% level: the level of noise in the unit of standard deviation of noise
%    OUTPUT
% rdt£º the dataset without defined noise region
sd=std(cdt(1,:));
hf=zeros(1,size(cdt,2));
sm=0;
pl=0;
np = 0;
rdt=cdt;
for i=1:size(cdt,2)
    if(cdt(1,i)>=thre)%
        sm=sm+1;
        rdt(:,sm)=cdt(:,i);
        hf(i)=1; % status
    end
end
rdt=rdt(:,1:sm);
  for i =1:length(hf)
       if hf(i)>0
           pl = pl+1;
           if i<=length(hf)-1
             if hf(i)~=hf(i+1)
                 np = np+1;
                 PL(np)= pl;
                 pl = 0;
             end
           else
               np = np+1;
               PL(np)= pl;
               pl = 0;
           end
       end   
  end
  
% plot the result to check
plot(cdt(1,:),'k');
hold on
plot(hf*max(cdt(1,:))/50,'r');
hold off
title(['Spectra above the threshold  ',num2str(thre)]);
ymax=max(max(cdt(:,:)));
xmax=size(cdt,2);
axis([1 xmax ymax*-0.05 ymax*1.05]);
end

