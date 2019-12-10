%% QGC Sample reconstruction by LRSpILT with two steps regularized parameters setting, same as shown in our paper:
clc,clear all,close all

load ('QGC.mat');% sample containning Quinine,Geranoil and Camphene and other peaks have been removed for best observation
 
 % initialization 
dt=real(NmrData.SPECTRA);% processed experimental data

 if size(dt,1)>size(dt,2)
     dt = dt.';  
 end  

g=100*NmrData.Gzlvl; % gradient values
BD=NmrData.DELTAOriginal; % diffusion time
LD=NmrData.deltaOriginal; % diffusion encoding time
cs=NmrData.Specscale;     % chemical shift
gamma = 4257.7;

 g2 = (2*pi*gamma*g*LD).^2*(BD-LD/3)*1e4;
 g2 = g2*1e-10;
 difCoef=(1:0.1:20);
 K = exp(-g2.'*difCoef); 
 cdt = dt;
 

%--------------------Thresholding for time conservation -------------------
thr = 0.0138;
cdt = cdt/max(cdt(:));
 % for view
 figure,
 plot(cs,cdt(1,:));set(gca,'Xdir','reverse');title('Spectrum');xlim([2,14]);
[rdt,hf,PL] = rsdt(cdt,thr);
rdt =rdt/max(rdt(:))*2;
 % ------------------------------  LRSpILT  ---------------------------------
   X_ori = 0;
   mu = 1;
    for k = 1:size(K,2)
        norfa(k) = norm(K(:,k),2);
        K(:,k) =  K(:,k)/ norfa(k);
    end    
       lambdalr =0.005;
       lambdasp =0.0005;
       
       [X_SPLRA,Conv,Nu,L1,Fn,out] = ADMMBased_Solver_2steps(rdt,K,X_ori,lambdalr,lambdasp,mu);
      
   for k = 1:size(K,2)
        X_SPLRA(k,:) =X_SPLRA(k,:)/ norfa(k);
    end 
     [DOSY_X] = rvdt(hf,X_SPLRA);
      DOSY_X =  DOSY_X/max( DOSY_X(:));
      Dif_Proj = sum(DOSY_X,2);
 
    
 % --------------------------------  Display -------------------------------------
%%
      figure,mesh(cs,difCoef,DOSY_X);xlabel('Chemical Shift/ppm');ylabel('Diffusion Coefficient/10^-^1^0m^2s^-^1');
      title('qgc LRSp');set(gca,'Ydir','reverse');
      figure,plot(Conv);xlabel('Iterations');title('Convergence Curve');
      figure,
      subplot(2,2,1),plot(Nu);xlabel('Iterations');title('Nuclear Norm');
      subplot(2,2,2),plot(L1);xlabel('Iterations');title('L1 Norm');
      subplot(2,2,3),plot(Fn);xlabel('Iterations');title('Frobenius Norm');
      subplot(2,2,4),plot(out);xlabel('Iterations');title('Loss Function');
%------------------ By observation selection -----------------------------
 
   Comp_pos = [4.7,7.3,10.1];
   Comp_num = length(Comp_pos);
    for i =1: Comp_num
      Comp(i,:) = DOSY_X(find(round(difCoef,3)==Comp_pos(i)),:) ; 
    end
%% ---------------------------------------------------------------------------
      figure,
      ax1 = axes('position',[0.05 0.7 0.8 0.3]);
      plot(ax1,cs,cdt(1,:));set(gca,'Xdir','reverse');axis off;
      xlim([2,14]);
      ax2 = axes('position',[0.05 0.23 0.8 0.5]);
      contour(ax2,cs,difCoef,DOSY_X,40);xlabel('Chemical Shift/ppm');ylabel('Diffusion Coefficient/10^-^1^0m^2s^-^1');
      set(ax2,'Ydir','reverse','Xdir','reverse'); 
      set(ax2,'YTick', Comp_pos);
      xlim([2,14]);ylim([1,15]);
      for i = 1:Comp_num
           line(ax2,get(ax2,'xlim'),[Comp_pos(i),Comp_pos(i)],'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
      end  
      ax3 = axes('position',[0.855 0.23 0.1 0.5]);
      plot(ax3,Dif_Proj,difCoef,'LineWidth',1);axis off;set(gca,'Ydir','reverse');
       ylim([1,15]);
      
 %% ----------------------- Extract Components -------------------------
     
figure, 
ax1 = axes('position',[0.1 0.15 0.8 0.65]);
for i = 1: Comp_num
   plot(ax1,cs,Comp(i,:)/max(Comp(i,:))+1.1*(Comp_num-i),'LineWidth',0.3,'color','blue');hold on;
end
hold off;ylim([-0.1,Comp_num*1.1]);xlim([2,14]);
set(gca,'Xdir','reverse','XMinorTick','on','YTick',[]);
xlabel('Chemical Shift/ppm');
ax2 = axes('position',[0.1 0.81 0.8 0.2]);
plot(ax2,cs,cdt(1,:),'color','blue');set(gca,'Xdir','reverse');xlim([2,14]);
axis off;


   