% close all
%% plot figures about harmonic structure

%% plot the the square of R in the lateral direction
fig2 = figure;

subplot(4,1,1)
% bar(R_square_U_determination(str_frm_num,:),'FaceColor',[.69 .69 .69])
x=0:n_points; y=0:0.1:0.5;
plot([0 0],[min(y) max(y)],'k');hold on;
plot([n_points n_points],[min(y) max(y)],'r');hold on;
for i=1:length(y)
    if y(i)~=0
        plot([-0.5,0],[y(i),y(i)],'k'); hold on
        plot([n_points,n_points+0.5],[y(i),y(i)],'r'); hold on
        b=text(n_points+2,y(i),num2str(abs(y(i)-0.5)),'FontSize',15,'Color','red');
        b1=text(-2,y(i),num2str(y(i)+0.5),'FontSize',15,'Color','black');
        set(b,'HorizontalAlignment','center')
        set(b1,'HorizontalAlignment','center')
    end
end
plot([-0.5,0],[0,0],'k'); hold on
b1=text(-2,0,num2str(0.5),'FontSize',15,'Color','black');
set(b1,'HorizontalAlignment','center')
plot([n_points,n_points+0.5],[0,0],'r'); hold on
b1=text(n_points+2,0,num2str(0.5),'FontSize',15,'Color','red');
set(b1,'HorizontalAlignment','center')

% plot the secular velocity term of Fourier series
for i_square=1:n_points
    R_square_square=[];
    %The coordinate of square, and the area of square is 1 its center is [0 0]
    height=R_square_U_determination(str_frm_num,i_square);
    if height>0.5|| height==0.5%if the square is sides parallel, and it indicaes the velocity is positive.. and diagonal side means negative
        R_square_square=[0   0           1        1 0 
                         0  height-0.5 height-0.5 0 0];
    end
     if height<0.5%if the square is sides parallel, and it indicaes the velocity is positive.. and diagonal side means negative
    R_square_square=[0    0            1         1 0 
                    0.5 0.5-height 0.5-height 0.5 0.5];
    end
    R_square_square=R_square_square+[(i_square-1)*ones(1,5);zeros(1,5)];
plot(R_square_square(1,:),R_square_square(2,:),'k');
if height>0.5|| height==0.5
fill(R_square_square(1,:),R_square_square(2,:),[.69 .69 .69]);
end
if height<0.5
fill(R_square_square(1,:),R_square_square(2,:),'r');
end
hold on
end
% axis([-2 32 0 0.5]) % axis([xmin xmax ymin ymax]) 
axis off
title(data_name_simple)


%% plot the the square of R in the axial direction
% figure
subplot(4,1,4)

% bar(R_square_U_determination(str_frm_num,:),'FaceColor',[.69 .69 .69])
x=0:n_points; y=0:0.1:0.5;
plot([0 0],[min(y) max(y)],'k');hold on;
plot([n_points n_points],[min(y) max(y)],'r');hold on;
for i=1:length(y)
    if y(i)~=0
        plot([-0.5,0],[y(i),y(i)],'k'); hold on
        plot([n_points,n_points+0.5],[y(i),y(i)],'r'); hold on
        b=text(n_points+2,y(i),num2str(abs(y(i))),'FontSize',15,'Color','red');
        b1=text(-2,y(i),num2str(abs(y(i)-1)),'FontSize',15,'Color','black');
        set(b,'HorizontalAlignment','center')
        set(b1,'HorizontalAlignment','center')
    end
end
plot([-0.5,0],[0,0],'k'); hold on
b1=text(-2,0,num2str(1),'FontSize',15,'Color','black');
set(b1,'HorizontalAlignment','center')
plot([n_points,n_points+0.5],[0,0],'r'); hold on
b1=text(n_points+2,0,num2str(0),'FontSize',15,'Color','red');
set(b1,'HorizontalAlignment','center')

% plot the secular velocity term of Fourier series
for i_square=1:n_points
    R_square_square=[];
    %The coordinate of square, and the area of square is 1 its center is [0 0]
    height=R_square_V_determination(str_frm_num,i_square);
    if height>0.5%if the square is sides parallel, and it indicaes the velocity is positive.. and diagonal side means negative
     R_square_square=[0    0            1         1 0 
                     0.5 height-0.5 height-0.5 0.5 0.5];
    end
     if height<0.5|| height==0.5%if the square is sides parallel, and it indicaes the velocity is positive.. and diagonal side means negative
R_square_square=[0   0           1        1 0 
                 0  height height 0 0];
    end
    R_square_square=R_square_square+[(i_square-1)*ones(1,5);zeros(1,5)];
plot(R_square_square(1,:),R_square_square(2,:),'k');
if height>0.5|| height==0.5
fill(R_square_square(1,:),R_square_square(2,:),[.69 .69 .69]);
end
if height<0.5
fill(R_square_square(1,:),R_square_square(2,:),'r');
end
hold on
end
% axis([-2 32 0 0.5]) % axis([xmin xmax ymin ymax]) 
axis off

save_name = save_path_root+data_name+"_fig2.fig";
saveas(fig2,save_name,"fig")
%% figure2
% subplot(3,1,2) %Lateral direction
fig3 = figure;
% subplot(4,1,[2 3])
% %1.2HZ
aaaa=0.01;
area_acc_triangle=0.06*aaaa;%%scaling the area
area_vel_square=0.08*aaaa;%%scaling the area
area_fundamental_Circle=0.2*0.8*0.5;%%scaling the area0.38

% area_acc_triangle=0.15;%%scaling the area
% area_vel_square=0.015;%%scaling the area
% area_fundamental_Circle=0.45;%%scaling the area0.38


% plot the secular acceleration term of Fourier series
acc_high=2;% the height ofain the graphy
Tans_matrix=[];M1=[];
for i_acc=1:n_points
%The coordinate of Triangle, and the area of square is 1, its center is[0 0]
    acc_triangle0=2/3/sqrt(3)*[-sqrt(3)/2  0  sqrt(3)/2 -sqrt(3)/2  
                                -0.5       1  -0.5     -0.5 ];
    acc_triangle=area_acc_triangle*(abs(Lateral.acc(i_acc)))^(1/1)*acc_triangle0;%different areas of vel in the Lateral direction..
    if Lateral.acc(i_acc)<0%if the point of triangle is upward, and it indicaes the accleration is positive.. and downward direction means negative
        M1=[cos(pi)  sin(pi) ;
         -sin(pi) cos(pi)   ;];%rotation Matrix   
      acc_triangle=M1*acc_triangle;% rotation, and be  downward direction
    end
    Tans_matrix=[i_acc*ones(1,4);acc_high*ones(1,4)];%Translation Matrix...the triangle would correspond to each body point
    acc_triangle=acc_triangle+Tans_matrix; 
plot(acc_triangle(1,:),acc_triangle(2,:),'k');
fill(acc_triangle(1,:),acc_triangle(2,:),'k')
hold on
end
Tans_matrix=[];M1=[];
% plot the secular velocity term of Fourier series
vel_high=4;% the height of velocity in the graphy
for i_vel=1:n_points
    %The coordinate of square, and the area of square is 1 its center is [0 0]
    vel_square0=[-0.5  0.5  0.5  -0.5 -0.5 
                -0.5 -0.5  0.5  0.5  -0.5 ];
    vel_square=area_vel_square*(abs(Lateral.vel(i_vel)))^(1/1)*vel_square0;%different areas of vel in the Lateral direction..
    if Lateral.vel(i_vel)<0%if the square is sides parallel, and it indicaes the velocity is positive.. and diagonal side means negative
        M1=[cos(pi/4)  sin(pi/4) ;
         -sin(pi/4) cos(pi/4)   ;];%rotation Matrix   
      vel_square=M1*vel_square;% rotation, and be  diagonal side square
    end
    Tans_matrix=[i_vel*ones(1,5);vel_high*ones(1,5)]; %Translation Matrix... the square would correspond to each body point
    vel_square=vel_square+Tans_matrix;
plot(vel_square(1,:),vel_square(2,:),'k');
fill(vel_square(1,:),vel_square(2,:),'k')
hold on
end


% plot the fundamental frequency term of Fourier series
fundamental_high=7;% the height of fundamental wave in the graphy
First_harmonic_high=10;% the height of fundamental wave in the graphy

% Color_fre=hsv;% colormap hsv 
% Color_num1=[(1:32)/32*pi  ];%%the corresponding radians of colors   0-pi
% Color_num2=[(1:32)/32*pi-pi];%%the corresponding radians of colors  -pi-0
Color_fre=[1 0 0;ones(44,1) (((1:44))/45)'  zeros(44,1);
           (1-((1:42))/43)' ones(42,1) zeros(42,1);
           zeros(42,1) ones(42,1) (((1:42))/43)' ;0 1 1;
           zeros(42,1) (1-((1:42))/43)'  ones(42,1)  ; 
           (((1:42))/43)' zeros(42,1)  ones(42,1);
          ones(42,1)   zeros(42,1)  (1-((1:42))/43)'];% colormap hsv 
Color_num1=[(1:128)/128*pi  ];%%the corresponding radians of colors   0-pi
Color_num2=[(1:128)/128*pi-pi];%%the corresponding radians of colors  -pi-0


for i_fundamental=1:n_points
    index=[];fundamental_Circle=[];index_First_harmonic=[];
    %The coordinate of cycle, and the area of cycleis 1 its center is [0 0]
    theta=0:pi/100:(2*pi+pi/100);  %%draw a closed circle so that the radian should be larger than 2*pi
    Circle1=sqrt(1/pi)*cos(theta);   Circle2=sqrt(1/pi)*sin(theta);  
    fundamental_Circle=area_fundamental_Circle*(abs(Lateral.fundamental_An1(i_fundamental)))^(1/1)*[Circle1;Circle2 ];%different areas of vel in the Lateral direction..  
    Tans_matrix_fundamental=[i_fundamental*ones(1,size(Circle1,2));fundamental_high*ones(1,size(Circle1,2))]; %Translation Matrix... the square would correspond to each body point
    fundamental_Circle=fundamental_Circle+Tans_matrix_fundamental;
% select the color 
if Lateral.fundamental_phi_n1(i_fundamental)>0
[Asort,index]=min(abs((Color_num1(:)-Lateral.fundamental_phi_n1(i_fundamental))));% sort the min value 
end
if Lateral.fundamental_phi_n1(i_fundamental)<0
[Asort,index]=min(abs((abs(Color_num2(:))-abs(Lateral.fundamental_phi_n1(i_fundamental)))));% sort the min value 
index=index+128;
end
selected_color =Color_fre(index,:);%
% plot(fundamental_Circle(1,:),fundamental_Circle(2,:),'k');
fill(fundamental_Circle(1,:),fundamental_Circle(2,:),selected_color,'LineStyle','none')
hold on
% First_harmonic
First_harmonic_Circle=area_fundamental_Circle*(abs(Lateral.An1st(i_fundamental)))^(1/1)*[Circle1;Circle2 ];%different areas of vel in the Lateral direction..  
Tans_matrix_First_harmonic=[i_fundamental*ones(1,size(Circle1,2));First_harmonic_high*ones(1,size(Circle1,2))]; %Translation Matrix... the square would correspond to each body point
First_harmonic_Circle= First_harmonic_Circle+Tans_matrix_First_harmonic;
if Lateral.phi_n1st(i_fundamental)>0
[Asort,index_First_harmonic]=min(abs((Color_num1(:)-Lateral.phi_n1st(i_fundamental))));% sort the min value 
end
if Lateral.phi_n1st(i_fundamental)<0
[Asort,index_First_harmonic]=min(abs((abs(Color_num2(:))-abs(Lateral.phi_n1st(i_fundamental)))));% sort the min value 
index_First_harmonic=index_First_harmonic+128;
end
selected_color_First_harmonic =Color_fre(index_First_harmonic,:);%
% plot(First_harmonic_Circle(1,:),First_harmonic_Circle(2,:),'k');
fill(First_harmonic_Circle(1,:),First_harmonic_Circle(2,:),selected_color_First_harmonic,'LineStyle','none')

hold on
end


%% in the axial direction
Lateral1=Lateral;%for convient
Lateral=Axial;%for convient

% area_acc_triangle=0.6;
% area_vel_square=0.5;
% area_fundamental_Circle=0.1;
% area_acc_triangle=0.06;%%scaling the area
% area_vel_square=0.008;%%scaling the area
% area_fundamental_Circle=0.38;%%scaling the area

% plot the secular acceleration term of Fourier series
acc_high=-2;% the height ofain the graphy
Tans_matrix=[];M1=[];
for i_acc=1:n_points
%The coordinate of Triangle, and the area of square is 1, its center is[0 0]
    acc_triangle0=2/3/sqrt(3)*[-sqrt(3)/2  0  sqrt(3)/2 -sqrt(3)/2  
                                -0.5       1  -0.5     -0.5 ];
    acc_triangle=area_acc_triangle*(abs(Lateral.acc(i_acc)))^(1/1)*acc_triangle0;%different areas of vel in the Lateral direction..
    if Lateral.acc(i_acc)<0%if the point of triangle is upward, and it indicaes the accleration is positive.. and downward direction means negative
        M1=[cos(pi)  sin(pi) ;
         -sin(pi) cos(pi)   ;];%rotation Matrix   
      acc_triangle=M1*acc_triangle;% rotation, and be  downward direction
    end
    Tans_matrix=[i_acc*ones(1,4);acc_high*ones(1,4)];%Translation Matrix...the triangle would correspond to each body point
    acc_triangle=acc_triangle+Tans_matrix; 
plot(acc_triangle(1,:),acc_triangle(2,:),'k');
fill(acc_triangle(1,:),acc_triangle(2,:),'k')
hold on
end
Tans_matrix=[];M1=[];
% plot the secular velocity term of Fourier series
vel_high=-4;% the height of velocity in the graphy
for i_vel=1:n_points
    %The coordinate of square, and the area of square is 1 its center is [0 0]
    vel_square0=[-0.5  0.5  0.5  -0.5 -0.5 
                -0.5 -0.5  0.5  0.5  -0.5 ];
    vel_square=area_vel_square*(abs(Lateral.vel(i_vel)))^(1/1)*vel_square0;%different areas of vel in the Lateral direction..
    if Lateral.vel(i_vel)<0%if the square is sides parallel, and it indicaes the velocity is positive.. and diagonal side means negative
        M1=[cos(pi/4)  sin(pi/4) ;
         -sin(pi/4) cos(pi/4)   ;];%rotation Matrix   
      vel_square=M1*vel_square;% rotation, and be  diagonal side square
    end
    Tans_matrix=[i_vel*ones(1,5);vel_high*ones(1,5)]; %Translation Matrix... the square would correspond to each body point
    vel_square=vel_square+Tans_matrix;
plot(vel_square(1,:),vel_square(2,:),'k');
fill(vel_square(1,:),vel_square(2,:),'k')
hold on
end



% plot the fundamental frequency term of Fourier series
fundamental_high=-7;% the height of fundamental wave in the graphy
First_harmonic_high=-10;% the height of fundamental wave in the graphy

% Color_fre=hsv;% colormap hsv 
% Color_num1=[(1:32)/32*pi  ];%%the corresponding radians of colors   0-pi
% Color_num2=[(1:32)/32*pi-pi];%%the corresponding radians of colors  -pi-0

for i_fundamental=1:n_points
    index=[];fundamental_Circle=[];index_First_harmonic=[];
    %The coordinate of cycle, and the area of cycleis 1 its center is [0 0]
    theta=0:pi/100:(2*pi+pi/100);  %%draw a closed circle so that the radian should be larger than 2*pi
    Circle1=sqrt(1/pi)*cos(theta);   Circle2=sqrt(1/pi)*sin(theta);  
    fundamental_Circle=area_fundamental_Circle*(abs(Lateral.fundamental_An1(i_fundamental)))^(1/1)*[Circle1;Circle2 ];%different areas of vel in the Lateral direction..  
    Tans_matrix_fundamental=[i_fundamental*ones(1,size(Circle1,2));fundamental_high*ones(1,size(Circle1,2))]; %Translation Matrix... the square would correspond to each body point
    fundamental_Circle=fundamental_Circle+Tans_matrix_fundamental;
% select the color 
if Lateral.fundamental_phi_n1(i_fundamental)>0
[Asort,index]=min(abs((Color_num1(:)-Lateral.fundamental_phi_n1(i_fundamental))));% sort the min value 
end
if Lateral.fundamental_phi_n1(i_fundamental)<0
[Asort,index]=min(abs((abs(Color_num2(:))-abs(Lateral.fundamental_phi_n1(i_fundamental)))));% sort the min value 
index=index+128;
end
selected_color =Color_fre(index,:);%
% plot(fundamental_Circle(1,:),fundamental_Circle(2,:),'k');
fill(fundamental_Circle(1,:),fundamental_Circle(2,:),selected_color,'LineStyle','none')
hold on
% First_harmonic
First_harmonic_Circle=area_fundamental_Circle*(abs(Lateral.An1st(i_fundamental)))^(1/1)*[Circle1;Circle2 ];%different areas of vel in the Lateral direction..  
Tans_matrix_First_harmonic=[i_fundamental*ones(1,size(Circle1,2));First_harmonic_high*ones(1,size(Circle1,2))]; %Translation Matrix... the square would correspond to each body point
First_harmonic_Circle= First_harmonic_Circle+Tans_matrix_First_harmonic;
if Lateral.phi_n1st(i_fundamental)>0
[Asort,index_First_harmonic]=min(abs((Color_num1(:)-Lateral.phi_n1st(i_fundamental))));% sort the min value 
end
if Lateral.phi_n1st(i_fundamental)<0
[Asort,index_First_harmonic]=min(abs((abs(Color_num2(:))-abs(Lateral.phi_n1st(i_fundamental)))));% sort the min value 
index_First_harmonic=index_First_harmonic+128;
end
selected_color_First_harmonic =Color_fre(index_First_harmonic,:);%
% plot(First_harmonic_Circle(1,:),First_harmonic_Circle(2,:),'k');
fill(First_harmonic_Circle(1,:),First_harmonic_Circle(2,:),selected_color_First_harmonic,'LineStyle','none')

hold on
end

axis([-1 n_points+5 -11 11]) % axis([xmin xmax ymin ymax]) 
axis off
x=-1:n_points+1; y=-11:11;
plot([0 0],[min(y) max(y)],'k',[min(x) max(x)],[0 0],'k');hold on;
for i=1:length(x)-1
    if x(i)~=0
        plot([x(i),x(i)],[0,0.1],'k'); hold on
      if mod(i-2,3)==0
        a=text(x(i),-0.6,num2str(x(i)));
        set(a,'HorizontalAlignment','center','FontSize',16)
      end
    end
%     if y(i)~=0
%         plot([0,0.1],[y(i),y(i)],'k'); hold on
%         b=text(-0.4,y(i),num2str(y(i)));
%         set(b,'HorizontalAlignment','center')
%     end
end
text(-2,5,'Lateral','fontsize',16,'rotation',90) 
text(-2,-5,'Axial','fontsize',16,'rotation',90) 

title(data_name_simple)

save_name = save_path_root+data_name+"_fig3.fig";
saveas(fig3,save_name,"fig")

%%
figure; %%phase radians
hold on
axis equal
% grid on
R1 = 1.5;
R2 = 3.5;
% t = linspace(0, 2*pi, 500);
% % Circles
% plot(R1*cos(t), R1*sin(t), 'k', ...
%     R2*cos(t), R2*sin(t), 'k');
% colormap(h_figure,hsv)
%% Quadrant II
% Color_fre=hsv;% colormap hsv 
for i_phase=1:256
t = linspace(i_phase*pi/128, (i_phase+1)*pi/128, 128);
rt = t(end:-1:1);
fill([R1*cos(t) R2*cos(rt)], [R1*sin(t) R2*sin(rt)], Color_fre(i_phase,:),'LineStyle','none')
end
axis off
% text(4,0,'0','fontsize',50) 
% text(0,4,'0','fontsize',50) 

Lateral=Lateral1;%for convient



