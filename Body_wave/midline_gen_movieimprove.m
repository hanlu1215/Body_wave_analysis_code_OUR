% midline generator
% generate midlines
% Modified by Jiang 2017/7/6
clear; close all;
Threshhold=0.45;%设置二值化阈值 determine the threshold value
n_points=30;%选取鱼体波中心线数据点的个数 determine the body points 
body_x = [];
body_y = [];
num_orign=73;%%要分析图片的开始值   determine the starting number of frame
num_end=135;%%要分析图片的末端值   determine the end number of frame
num_mid=104;%%图片中间跨度线  Sometimes it can be analyzed two times
n_frames=num_end-num_orign+1;% the number of frames to analyze
n_middle=num_mid-num_orign+1;
M(n_frames) = struct('cdata',[],'colormap',[]);
%% maybe first time
for i =1:n_frames
clf 
if i>num_mid %Sometimes it can be analyzed two times acoording to the conditions
%   Threshhold=input('please enter Threshhold:');%重新设置二值化阈值  0.4
  Threshhold=0.42;%重新设置二值化阈值  0.4
end
    I0=imread(strcat(num2str(i+num_orign-1),'.jpg'));%读取一帧图像
    % imshow(I)
    I = rgb2gray(I0);%转化为灰度图像
%     BW = imbinarize(I);
    BW = im2bw(I,Threshhold);%图像二值化处理 
    %显示图像
    imshow(BW)
% obtain automatically the analyzed area (pixel)  
if i<num_mid-num_orign
position =[    213  0  327  329];
p1 = position(1); p2 = position(2);p3 = position(3);  p4 = position(4);
else
  %%用鼠标选取目标区域  obtain manually the analyzed area 
    h = imrect;
    position = round(wait(h));%[X,Y,W,H]
    p1 = position(1);
    p2 = position(2);
    p3 = position(1)+position(3);
    p4 = position(2)+position(4);
end

    %清除选定区域以外的噪声 clear the noise
    BW(1:p2,:) = ones(size(BW(1:p2,:)));
    BW(p4:480,:) = ones(size(BW(p4:480,:)));
    BW(:,1:p1) = ones(size(BW(:,1:p1)));
    BW(:,p3:640) = ones(size(0*BW(:,p3:640)));
    %边缘检测 image edge detection
    BW2 = bwmorph(BW,'Majority',1);
    BWE = edge(BW2,'sobel');
%     subplot(2,1,1)
%     imshow(BW2)
%     subplot(2,1,2)
%     imshow(BWE)
    imshow(I0); 
    hold on
    
    MIDLINES = [];
    for jh = 2:479
        jv = find(BWE(jh,:)==1);
        if isempty(jv)
            continue
        end
        X = jh;
        Y = round((max(jv)+min(jv))/2);
        MIDLINES = [MIDLINES;X,Y];
        
    end
    plot(MIDLINES(:,2),MIDLINES(:,1),'ro','MarkerSize',2)
     
      temp_x = [];temp_y = [];
      %手动选择关键点
%      for num = 1:n_points
%           zoom
%          %pause
%          [tempx,tempy] = ginput(1); %鼠标获取鱼体上的10个点
%          temp_x = [temp_x;tempx]; 
%          temp_y = [temp_y;tempy];
%      end

    %自动选择关键点  obtain automatically the body points
        delta=(max(MIDLINES(:,1))-min(MIDLINES(:,1)))/(n_points-1);
        temp_x=min(MIDLINES(:,1)):delta:max(MIDLINES(:,1));
        temp_y=interp1(MIDLINES(:,1),MIDLINES(:,2),temp_x);
        plot(temp_y,temp_x,'g*','MarkerSize',5)
        
str = input('obtain the point manually: Y or N?\n','s');%%obtain automatically the body points
if str=='Y'
     temp_x = [];
      temp_y = [];
     for num = 1:n_points
          zoom
         %pause
         [tempy,tempx] = ginput(1); %鼠标获取鱼体上的30个点
         temp_x = [temp_x;tempx]; 
         temp_y = [temp_y;tempy];
     end
        pause(1)
      plot(temp_y,temp_x,'g*','MarkerSize',5)   
      temp_x=temp_x';
       temp_y= temp_y';
end
     body_x = [body_x;temp_x];%10*30的矩阵，存放鱼体上10个点共30帧的x坐标
     body_y = [body_y;temp_y];%10*30的矩阵，存放鱼体上10个点共30帧的y坐标
     M(i) = getframe(gcf);
     pause(0.1)
end
save('bd_data.mat','body_x','body_y');
movie(M) %Play the video
% v = VideoWriter('newfile.avi','Uncompressed AVI'); 
v = VideoWriter('newfile.avi','MPEG-4' );
open(v)
writeVideo(v,M)
close(v)