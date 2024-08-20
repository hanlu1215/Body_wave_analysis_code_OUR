 I=imread(strcat(num2str(119),'.jpg'));%读取一帧图像
%  Background=imread('1761.jpg');
%     I=imsubtract(I,Background);
%     I=imabsdiff(I,Background);
%     K=imadd(imabsdiff(I,Background),50);
    imshow(I)
%     return
%     I = rgb2gray(K);%转化为灰度图像
% %     BW = imbinarize(I);
%     BW = ~im2bw(I,Threshhold);%图像二值化处理 
%     %显示图像
%     imshow(BW)
    %用鼠标选取目标区域
    h = imrect;
    position = round(wait(h));%[X,Y,W,H]
   
    p1 = position(1);
    p2 = position(2);
    p3 = position(1)+position(3);
    p4 = position(2)+position(4);
    p=[ p1   p2  p3   p4 ]
    
% % p0=p;
distance=sqrt((p0(1)-p(1))^2+(p0(2)-p(2))^2)
velocity=distance/(num_end-num_orign)*30/fish_length0