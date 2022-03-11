%calcultions for aquaculture


%% total demand 1960-2000
clc
clear
load WD_NONAGR_SETTORI.txt
aqua_static= WD_NONAGR_SETTORI(:,1);
%%
load hist_values.mat
area= hist_values(1,:);

%load evap_lo.txt
load p_lo.txt
load t_lo.txt
t_lo= t_lo(:,2);
dT_add = [-2,-1,1:5];
dP_add = [0.8,0.9,1.1:0.1:1.6];

 % loss evaporation
lat = 20;
dd = datenum(1961,1,1):datenum(2000,12,31);
dt = datevec(dd);
tdates = [dt(:,3) dt(:,2) dt(:,1)];
date = datetime(1961,1,1):datetime(2000,12,31);
date = yyyymmdd(date)';
date = date';
size=56;
seep_coeff= 8.5; %mm/day
demand= zeros(56,365);
for i=1:7 
    for j=1:8
        %evap loss
        T= t_lo+repmat(dT_add(i),length(t_lo),1);
      evap_lo = [date' ET_Thornthwaite(T,tdates,lat)];
        loss_evap= area.*evap_lo(:,2)*10; % multiply 10 for mm hectare/d to m3/d
        
%seepage loss
    loss_seep= area.*seep_coeff*10; % multiply 10 for mm hectare/d to m3/d
    loss_seep=repmat(loss_seep,length(t_lo),1);
    
%greenwater
gw= p_lo(:,2)*dP_add(j)*area*10;

wd_dist= loss_evap+loss_seep-gw;
wd_aqua= sum(wd_dist')'; wd_aqua(wd_aqua< 0) = 0;
 wd_aqua=reshape(wd_aqua(1:14600),365,40);
 aqua(j,:)= mean(wd_aqua');
    end
 demand(8*i-7 :8*i,:)=aqua;
end
wd_aquacult= demand; %m3/d
save('wd_aqua.mat','wd_aquacult');
%%

figure; bar(wd_aquacult); ylabel('Aquaculture water demand (m3/s)'); xlabel('preciitation scenarios');
hold on; wd_aquacult_56= reshape(wd_aquacult,size,1); ylim([0 110]);
plot(repmat(mean(wd_aquacult_56),8,1)); 
%wd_aquacult= reshape(wd_aquacult,size,1); 
figure; plot(wd_aquacult); %m3/d 
hold on;
plot(repmat(mean(wd_aquacult_56),8,1)); legend('aqua water demand','mean');
xlabel('56 sclimate scenarios'); ylabel('Aquaculture water demand (m3/d)');
wd_aquacult_56=wd_aquacult_56*(3600*24);

figure; plot(aqua_static);
ylabel('Aquaculture water demand (m3/s)'); xlabel('Scenarios by standard of water use');
ylim([0 110]);
hold on;
plot(repmat(mean(aqua_static),100,1));
legend('aqua water demand','mean');
save('wd_aquacult','wd_aquacult');
save('wd_aqua', 'wd_aquacult_56');

%% test
a=p_lo(1:14600,2);
a=reshape(a,[365,40]);
b=mean(a');