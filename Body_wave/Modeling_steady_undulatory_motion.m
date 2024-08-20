% reference：John H. Long Jr 《Flapping flexible fishPeriodic and secular body reconfigurations in swimming lamprey, Petromyzon marinus》 
% by Bingxing Chen 20181121
%   x/body_x：the axial direction  y/body_y：the lateral direction (not rotated)
% clc
% clear
% close all
disp('________________________________________________________________')
disp('Generating plots for bodywave')
load bd_data.mat

n_points=size(body_x,2);%fish body points
% DF=1.2;timestep = 1/30;num_per_cycle=round(1/DF/timestep)+1;
% NumFr=36-num_per_cycle;

for strt_frm=1:NumFr% different cycles
xData =body_x(strt_frm:strt_frm+1*num_per_cycle,:)';%% row: body points line:time series
yData = body_y(strt_frm:strt_frm+1*num_per_cycle,:)';
fishL =[];sz=size(xData);%%sz 求取两个值 sz(1)行：表示鱼体的位置 体长的数据 有多少个鱼体的点 sz(2)列：时间
for il=1:sz(2)%表示
fishL(il) =sum(sqrt(diff(xData(:,il)).^2+diff(yData(:,il)).^2));
end
fishLength=mean(fishL); %the body length  
for i_actual= 1:n_points%%body points
% U.xactual(i_actual,strt_frm)=30*(xData(i_actual,end)-xData(i_actual,1))/fishLength/ size(xData,2);%%velocity of x direction
% U.yactual(i_actual,strt_frm)=30*(yData(i_actual,end)-yData(i_actual,1))/fishLength/ size(yData,2);%%velocity of y direction
% the velocity also can use the actual velocity by the calibration for the actual fish length 
U.xactual(i_actual,strt_frm)=FPS*(xData(i_actual,end)-xData(i_actual,1))/ size(xData,2);%%velocity of x direction
U.yactual(i_actual,strt_frm)=FPS*(yData(i_actual,end)-yData(i_actual,1))/ size(yData,2);%%velocity of y direction
end
% drawnow;
% % pause
end
% 所有点的平均速度
U.Ave_x_actual=sum(U.xactual')/(size(U.xactual,2));
U.Ave_y_actual=sum(U.yactual')/(size(U.yactual,2));%the x/y components of the average velocity for individual body points
% 所有点的平均速度的平均值
U.Ave_xcom=sum(U.Ave_x_actual)/(size(U.Ave_x_actual,2));
U.Ave_ycom=sum(U.Ave_y_actual)/(size(U.Ave_y_actual,2));% the x/y components of the composite average velocity.
% 所有点的平均速度的平均值x y分量以及合量
U.Ave_com=[U.Ave_xcom  U.Ave_ycom sqrt(U.Ave_xcom^2+U.Ave_ycom^2)];% the value of the composite average velocity.
% %  calculate reconfiguration errors
U_actual=[];
for i = 1:size(U.xactual,1)
   for j = 1:size(U.xactual,2)    
    RE.x(i,j)=(U.xactual(i,j)-U.Ave_xcom)/U.Ave_com(3);%x:reconfiguration errors
    RE.y(i,j)=(U.yactual(i,j)-U.Ave_ycom)/U.Ave_com(3);%y:reconfiguration errors
    RE.Square(i,j)=sqrt(RE.x(i,j)^2+RE.y(i,j)^2);
    U_actual(i,j)=sqrt(U.xactual(i,j)^2+U.yactual(i,j)^2);
    RE.Square_UI(i,j)=RE.Square(i,j)^2*U.Ave_com(3)/U_actual(i,j);
end
end
%unsteadliness index
UI=sqrt(sum(sum(RE.Square_UI))/(size(U_actual,2)*n_points-1));%
disp( '************The calculations*************');
fprintf(1,'The composite average velocity (mm/s) is: %3.4f\n',U.Ave_com(3))
fprintf(1,'unsteadliness index is: %1.4f\n',UI)

% fprintf(1,'unsteadliness index is: %1.4f\n',U.Ave_com(3)/51/3.9)




%% the data are rotated so that the composite average velocity vectors points in the positive direction.
% theta1=atan(U.Ave_ycom/U.Ave_xcom );%%运动方向。。即运动方向与垂直成theat角度 tempx实际是图片中的y坐标。。tempy是像素中的x
theta1=atan2(U.Ave_ycom,U.Ave_xcom);
% hqx = [body_x(1,:)',body_y(end,:)'];
% hqy = [body_y(1,:)',body_x(end,:)'];
% [xD, yD] = prepareCurveData( hqx, hqy );
% % Set up fittype and options.
% ft = fittype( 'poly1' );
% % Fit model to data.
% [fitresult, gof] = fit( xD, yD, ft );
% theta1 = atan(fitresult.p1);%%运动方向。。即运动方向与垂直成theat角度 tempx实际是图片中的y坐标。。tempy是像素中的x
X = cos(theta1)*body_x'+sin(theta1)*body_y';%%the axial direction
Y = -sin(theta1)*body_x'+cos(theta1)*body_y';%the lateral direction
X=X';Y=Y';

% str_frm_Interval=strt_frm:num_per_cycle+strt_frm;
str_frm_num=1;% the number of cycle
str_frm_Interval=str_frm_num:num_per_cycle+NumFr;%%NumFr
u_actualbody=X(str_frm_Interval,:);
v_actualbody=Y(str_frm_Interval,:);
H=[];H0=[];
for i=str_frm_Interval%one cycle
% H0=[1 i*timestep (i*timestep)^2/2  sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep) cos(2*DF*2*pi*i*timestep) sin(3*DF*2*pi*i*timestep) cos(3*DF*2*pi*i*timestep) sin(4*DF*2*pi*i*timestep) cos(4*DF*2*pi*i*timestep)  ];
H0=[i*timestep (i*timestep)^2/2  sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep) cos(2*DF*2*pi*i*timestep)  ];
% H0=[1 i*timestep (i*timestep)^2/2  sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep) cos(2*DF*2*pi*i*timestep)  ];
H=[H;H0;];% 
end
% Equation:H*S=u_actualbody/v_actualbody
%% calculate the coefficient 
% u0=sum(u_actualbody)/size(u_actualbody,1);%the average value
% v0=sum(v_actualbody)/size(v_actualbody,1);
u0=u_actualbody(1,:);%the initial value
v0=v_actualbody(1,:);%初始值
%u0=0;%the initial value
%v0=0;%初始值

S_u=inv(H'*H)*H'*(u_actualbody-u0);%least square fitting method
S_v=inv(H'*H)*H'*(v_actualbody-v0);%least square fitting method
%lateral direction
Lateral.An=sqrt(S_v(3,:).^2+S_v(4,:).^2);%% the periodic amplitude term of Fourier series
% Lateral.phi_n=asin(S_v(4,:)./sqrt(S_v(3,:).^2+S_v(4,:).^2));
Lateral.phi_n=atan2(S_v(3,:),S_v(4,:));%%% the periodic phase term of Fourier series for the 1st harmonic
Lateral.An1st=sqrt(S_v(5,:).^2+S_v(6,:).^2);%% the periodic phase term of Fourier series for the 1st harmonic
Lateral.phi_n1st=atan2(S_v(5,:),S_v(6,:));%%% the periodic amplitude term of Fourier series for the 1st harmonic
% figure
% plot(1:n_points,unwrap(Lateral.phi_n),'r')%%phase changes
Lateralphi_n111=unwrap(Lateral.phi_n); 
% 展开向量 P 中的弧度相位角。每当连续相位角之间的跳跃大于或等于 π 弧度时，unwrap 
% 就会通过增加 ±2π 的整数倍来平移相位角，直到跳跃小于 π。如果 P 是矩阵，unwrap 将
% 按列运算。如果 P 是多维数组，unwrap 将对大小大于 1 的第一个维度进行运算。
% fliplr  because in the mideline analysis, the head point is 30, the tail 
% point is 1
fprintf(1,'The phase change of Lateral is: %1.4f\n',Lateralphi_n111(n_points)-Lateralphi_n111(1))
fprintf(1,'波数: %1.4f     ',abs(Lateralphi_n111(n_points)-Lateralphi_n111(1))/(2*pi))
fprintf(1,'UI: %1.4f\n',UI)

data_save.boshu = abs(Lateralphi_n111(n_points)-Lateralphi_n111(1))/(2*pi);
data_save.UI = UI;
% B = fliplr(A) 返回 A，围绕垂直轴按左右方向翻转其各列。
Lateral.har_phi_n1st=fliplr(Lateral.phi_n1st);
Lateral.har_An1st=fliplr(Lateral.An1st);%the value with respect to the value of first harmonic
Lateral.fundamental_phi_n1=fliplr(unwrap(Lateral.phi_n))+2*pi;%%unwrap

data_save.jipinxiangwei = Lateral.fundamental_phi_n1;


Lateral.fundamental_An1=fliplr(Lateral.An);
Lateral.vel=fliplr(S_v(1,:));%% the secular velocity term of Fourier series
Lateral.acc=fliplr(S_v(2,:));%% the secular acceleration of Fourier series
% anxial direction
Axial.An=sqrt(S_u(3,:).^2+S_u(4,:).^2);
Axial.phi_n=atan2(S_u(3,:),S_u(4,:));
Axial.An1st=sqrt(S_u(5,:).^2+S_u(6,:).^2);%% the periodic phase term of Fourier series for the 1st harmonic
Axial.phi_n1st=atan2(S_u(5,:),S_u(6,:));%%% the periodic amplitude term of Fourier series for the 1st harmonic
Axial.har_phi_n1st=fliplr(Axial.phi_n1st);
Axial.har_An1st=fliplr(Axial.An1st);%the value with respect to the value of first harmonic
Axial.fundamental_phi_n1=fliplr(Axial.phi_n);
Axial.fundamental_An1=fliplr(Axial.An);
Axial.vel=fliplr(S_u(1,:));%% the secular velocity term of Fourier series
Axial.acc=fliplr(S_u(2,:));%% the secular acceleration of Fourier series

%% plot 
% Reconfiguration
u_model=[];v_model=[];
u_model=H*S_u+u0;v_model=H*S_v+v0;%Reconfiguration body wave

u_average=sum(u_actualbody)/size(u_actualbody,1);
v_average=sum(v_actualbody)/size(v_actualbody,1);
u_SST=sum((u_actualbody-u_average).^2);
u_SSE=sum((u_model-u_actualbody).^2);
R_square_U_determination=1-u_SSE./u_SST;%拟合系数

v_average=sum(v_actualbody)/size(v_actualbody,1);
v_average=sum(v_actualbody)/size(v_actualbody,1);
v_SST=sum((v_actualbody-v_average).^2);
v_SSE=sum((v_model-v_actualbody).^2);
R_square_V_determination=1-v_SSE./v_SST;
% fliplr  because in the mideline analysis, the head point is 30, the tail point is 1
R_square_U_determination=fliplr(R_square_U_determination);%% 
R_square_V_determination=fliplr(R_square_V_determination);%% 

% check the fitness of all points for a single cycle period
Harmonic_structure_plot
% figure
% plot(1:n_points,unwrap(Lateral.phi_n),'r')%%鱼体波的图
% xlabel('body points)') ; ylabel('unwarp Phase changes') ;
% figure
% plot(1:n_points,Lateral.phi_n,'r')%%鱼体波的图
% xlabel('body points)') ; ylabel('Phase changes') ;
% figure
% plot(str_frm_Interval,u_model,'k',str_frm_Interval,u_actualbody,'b' )
% xlabel('time)') ;
% title('u model and u_actualbody')
% figure
% plot(str_frm_Interval,v_model,'k',str_frm_Interval,v_actualbody,'b' )
% xlabel('time)') ;
% title('v model and v_actualbody')
% figure
% plot(str_frm_Interval,v_model(:,1),'k',str_frm_Interval,v_actualbody(:,1),'b' )
% xlabel('time)') ;
% title('The v model and v_actualbody of tail point')
% body wave: check the fitness of a body wave  (raw data(rotated)   and reconfiguration data)
% figure
% plot(u_actualbody,v_actualbody,'r',u_model,v_model,'k')%%鱼体波的图 the all fish body;
% title('body wave model and actualbody')
% % body wave: check the fitness of a body wave  (raw data(rotated)   and reconfiguration data) at certain time 
% 
% num_time=1;%different time
% figure
% plot(u_actualbody(:,num_time),v_actualbody(:,num_time),'r',u_model(:,num_time),v_model(:,num_time),'k')%%鱼体波的图
% title('body wave model and actualbody of tail point')
% raw data /body wave 
body_xDatawave=fliplr(body_x);
body_yDatawave=fliplr(body_y);
for iiii=1:size(body_yDatawave,1)
% body_yDatawave(iiii,:)=body_yDatawave(iiii,:)+ (iiii-1)*5*ones(1,30);%%偏移一段距离
body_yDatawave(iiii,:)=body_yDatawave(iiii,:)+ (iiii-1)*5*ones(1,n_points);%%偏移一段距离
end
figure
plot(body_xDatawave(str_frm_Interval,:)',body_yDatawave(str_frm_Interval,:)','k-')%%鱼体波的图
title('Actual body wave')
% raw data /body wave 
u_modelDatawave=fliplr(u_model);
v_modelDatawave=fliplr(v_model);
for iiii=1:size(v_modelDatawave,1)
% v_modelDatawave(iiii,:)=v_modelDatawave(iiii,:)+ (iiii-1)*5*ones(1,30);%%偏移一段距离
v_modelDatawave(iiii,:)=v_modelDatawave(iiii,:)+ (iiii-1)*5*ones(1,n_points);%%偏移一段距离

end
figure
plot(u_modelDatawave(str_frm_Interval,:)',v_modelDatawave(str_frm_Interval,:)','k-')%%鱼体波的图
title('displaced model body wave')
% theta = atan(U.Ave_ycom/U.Ave_xcom );%%运动方向。。即运动方向与垂直成theat角度 tempx实际是图片中的y坐标。。tempy是像素中的x
u_modelDatawave=fliplr(u_model);
v_modelDatawave=fliplr(v_model);
u_model_rotated=cos(-theta1)*u_modelDatawave+sin(-theta1)*v_modelDatawave;%%the axial direction
v_model_rotated = -sin(-theta1)*u_modelDatawave+cos(-theta1)*v_modelDatawave;%the lateral direction
A_rotated = (max(v_model_rotated (:,n_points))-min(v_model_rotated (:,n_points)));%%鱼尾用Y数据  1表示鱼尾
% A_rotated = (max(v_model_rotated (:,30))-min(v_model_rotated (:,30)))%%鱼尾用Y数据  1表示鱼尾


for iiii=1:size(v_model_rotated ,1)
% v_model_rotated (iiii,:)=v_model_rotated (iiii,:)+ (iiii-1)*5*ones(1,30);%%偏移一段距离
v_model_rotated (iiii,:)=v_model_rotated (iiii,:)+ (iiii-1)*5*ones(1,n_points);%%偏移一段距离
end
figure
plot(u_model_rotated(str_frm_Interval,:)',v_model_rotated (str_frm_Interval,:)','k-')%%鱼体波的图
title('Rotated model body wave')
Velocity_vector1=[50; 200];
Velocity_vector2=[cos(-theta1) sin(-theta1); -sin(-theta1) cos(-theta1)]*[U.Ave_com(3);0]+Velocity_vector1;
line([Velocity_vector1(1) Velocity_vector2(1)],[Velocity_vector1(2) Velocity_vector2(2) ],'Color','magenta','LineStyle','-','LineWidth',2)
% annotation('arrow',Velocity_vector1/norm(Velocity_vector1),Velocity_vector2/norm(Velocity_vector2),'Color','r') 
text(Velocity_vector1(1),Velocity_vector1(2),'Initial','FontSize',15,'Color','black');
text(Velocity_vector1(1),Velocity_vector1(2),'Final','FontSize',15,'Color','black');
text(Velocity_vector1(1),Velocity_vector1(2),'Average Velocity','FontSize',15,'Color','black');
% axis off


%  one cycle frames to verify 
% u_actualbody=body_x(strt_frm:num_per_cycle+strt_frm,:);
% v_actualbody=body_y (strt_frm:num_per_cycle+strt_frm,:);
% strt_frm=1;
% H=[];H0=[];
% for i=strt_frm:num_per_cycle+strt_frm%one cycle
% H0=[1 i*timestep (i*timestep)^2/2  sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep) cos(2*DF*2*pi*i*timestep) ];
% H=[H;H0;];% 
% end
% % H*S=u_actualbody/v_actualbody
% % calculate the coefficient 
% S_u=inv(H'*H)*H'*u_actualbody;%least square fitting method
% S_v=inv(H'*H)*H'*v_actualbody;%least square fitting method
% % Reconfiguration
% u_model=[];v_model=[];
% u_model=H*S_u;v_model=H*S_v;
% 
% % body wave: check the fitness of a body wave  (raw data(rotated)   and reconfiguration data)
% figure
% plot(u_actualbody,v_actualbody,'r',u_model,v_model,'k')%%鱼体波的图 the all fish body
% % body wave: check the fitness of a body wave  (raw data(rotated)   and reconfiguration data) at certain time 
% num_time=1;%different 
% figure
% plot(u_actualbody(:,num_time),v_actualbody(:,num_time),'r',u_model(:,num_time),v_model(:,num_time),'k')%%鱼体波的图
% 
% % check the fitness of all points for a single cycle period
% time_variance=(strt_frm:num_per_cycle+strt_frm)*timestep;%the variance of time
% figure
% plot(time_variance,u_model,'k',time_variance,u_actualbody,'b' )
% figure
% plot(time_variance,v_model,'k',time_variance,v_actualbody,'b' )
% 
% % raw data /body wave 
% body_xDatawave=body_x;
% body_yDatawave=body_y;
% for iiii=1:size(body_yDatawave,1)
% body_yDatawave(iiii,:)=body_yDatawave(iiii,:)+ (iiii-1)*5*ones(1,30);%%偏移一段距离
% end
% figure
% plot(body_xDatawave(1:1:66,:)',body_yDatawave(1:1:66,:)','k-')%%鱼体波的图
% 
% % raw data /body wave 
% u_modelDatawave=u_model;
% v_modelDatawave=v_model;
% for iiii=1:size(v_modelDatawave,1)
% v_modelDatawave(iiii,:)=v_modelDatawave(iiii,:)+ (iiii-1)*5*ones(1,30);%%偏移一段距离
% end
% figure
% plot(u_modelDatawave(1:1:19,:)',v_modelDatawave(1:1:19,:)','k-')%%鱼体波的图


%% 保存图像数据

save_name = save_path_root+data_name+"_data.mat";
save(save_name,'data_save')
% disp(data_save)
