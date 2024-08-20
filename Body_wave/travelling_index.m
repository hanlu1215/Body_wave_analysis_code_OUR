hll
% disp('________________________________________________________________')
% % 直接从excel读取，没有经过尺寸转换，为像素距离，经过本节转换，转为mm单位
% % data_name = '20240701（橡皮筋月牙形尾鳍）20fps_k2l4ying_1hz_data_xy';
% % data_name = '20240701（橡皮筋月牙形尾鳍）30fps_k2l4ying_2hz_data_xy';
% % data_name = '20240701（橡皮筋月牙形尾鳍）60fps_k2l4ying_4hz_data_xy';
% % data_name = '20240701（橡皮筋月牙形尾鳍）60fps_k2l4ying_6hz_data_xy';
% % data_name = '20240701（橡皮筋月牙形尾鳍）60fps_k2l4ruan_6hz_data_xy';
% % data_name = '20240718（弹簧钢不均匀刚度月牙形）_19fps_K22222_2hz_data_xy';
% data_name = '20240718（弹簧钢不均匀刚度月牙形）_19fps_K22224_1hz_data_xy';
% 
% %%
% excel_filename = [data_name,'.xlsx'];
% data_msg = readmatrix(excel_filename,'Sheet','config');
% scale = 1000/data_msg(2);
% % scale = 1000/3138.32;
% % scale = 1000/2815.000710479484;
% data_fps = data_msg(4);
% charIndex_hz = regexp(data_name, 'hz','ignorecase');
% DF = str2double(data_name(charIndex_hz-1)) ;   
% data_save.f = DF;
% charIndex_fps = regexp(data_name, 'fps','ignorecase');
% if isempty(charIndex_fps)
%     FPS = 60/3;% 数据帧率
% else
%     FPS = str2double(data_name(charIndex_fps-2:charIndex_fps-1));
% end
% 
% body_x = readmatrix(excel_filename,'Sheet','x');
% body_y = readmatrix(excel_filename,'Sheet','y');
% % body_x = body_x';body_y = body_y';
% body_x = fliplr(-body_x)*scale;body_y = fliplr(-body_y)*scale;
% timestep = 1./FPS;% 相邻曲线的时间间隔
% num_per_cycle=round(1/DF/timestep)+1;% 每个周期的曲线数
load bd_data_1.2HZ.mat
DF=1.2;timestep = 1/30;num_per_cycle=round(1/DF/timestep)+1;

Y = body_y(:,1:num_per_cycle);
Y = Y-mean(mean(Y));
%% 计算传播指数
fig = figure;
fig.Position = [2854 675 560 420];
plot(Y)
Z  = zeros(size(Y));
for mm = 1:1:size(Y,1)
    Z(mm,:) = hilbert(Y(mm,:));
end
figure
nnn = 5;
plot(Y(nnn,:),'ko');hold on
plot(real(Z(nnn,:)),'r');
plot(imag(Z(nnn,:)),'b');

R = Z*conj(Z')/size(Y,2);
[w,~] = eig(R);
W  = [real(w),imag(w)];
a = 1/cond(W);
disp("传播指数： "+num2str(a))