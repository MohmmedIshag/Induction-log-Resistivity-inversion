%A hybrid inversion algorithm to obtain the resistivity of the uninvaded zone based on the array induction log using the folwing formation model
    % y = a+2(b-a)/1+exp(c*x)^d
% c value is equal to 1 

clc;close all;clear all;

% Load data (load your data on text filt )
load DD.txt;
data=DD;%replace the "DD with your uploaded dat"
xdata1 =[0.254 0.508 0.762 1.524 2.286 3.046] ;
resistivity = data(:,3:8); % specify the  resistivity data column
shale= data(:,2);specify the shale volume data column it will be used for constarin
depth= data(:,1);specify the depth column
cutoff=40;% cutoff value of shale volume
rr=resistivity;

[m,n]=size(rr);
ydata=zeros(m,n);

xdata=xdata1;
% c=1;
d=1;
Rt=zeros(m,1);
Rxo=zeros(m,1);

for i=1:m
    for j=1:n
        ydata(i,j)=1./rr(i,j);%conductivity mS/m
         if shale(i)<cutoff
             lb=[ydata(i,6)/3 ydata(i,1)/3  0.1]';%set lower bound constrain 
             ub=[ydata(i,6)*3 ydata(i,1)*3  100]';%set upper bound constrain
             
             options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'MaxIter',200,'TolFun',1e-6,'InitDamping',0.02);
             A = [];
             b = [];
             Aeq = [];
             beq = [];
             
             fun = @(x,xdata) x(1)+(2*(x(2)-x(1)))./(1+exp(xdata.^x(3)));% Define the objective function
             
             x0 = [ ydata(i,6),  ydata(i,1), d];% Initial guess
               
             [xoptimal,resnorm,residual,output,flag,jacobian] = lsqcurvefit(fun,x0,xdata,ydata(i,:),lb,ub,options);
             Rt=1/xoptimal(1); %transform to the resistivity of uninvaded zone
             Rxo=1/xoptimal(2);%transform to the resistivity fluzed zone
            
             d=xoptimal(3);
             result(i,:)=[Rt Rxo  d];
         else
             Rxo=resistivity(i,1);
             Rt=resistivity(i,6);

             d=1;
             result(i,:)=[Rxo Rt  d];
         end
    end
end
   you can save your data result on your computer as excell file
