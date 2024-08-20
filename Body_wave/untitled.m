%% 直接从excel读取，没有经过尺寸转换，为像素距离，经过本节转换，转为mm单位
hll
% data_names = {'20240701（橡皮筋月牙形尾鳍）20fps_k2l4ying_1hz_data_xy';...
%     '20240701（橡皮筋月牙形尾鳍）30fps_k2l4ying_2hz_data_xy';...
%     '20240701（橡皮筋月牙形尾鳍）60fps_k2l4ying_4hz_data_xy';...
%     '20240701（橡皮筋月牙形尾鳍）60fps_k2l4ying_6hz_data_xy'};
% 指定文件扩展名
fileExtension = '*.xlsx'; % 这里以读取所有 .txt 文件为例
% 获取当前文件夹内所有指定扩展名的文件信息
fileInfo = dir(fileExtension);
% 提取文件名字
data_names = {fileInfo.name};
% data_names = {'20240716（弹簧钢不均匀刚度月牙形）15FPS_k62_1hz_data_xy';...
%     '20240701（橡皮筋月牙形尾鳍）15FPS_k2l4ying_1hz_data_xy'};
for ii = 10:1:size(data_names,2)
data_name = data_names{ii}(1:end-5);
Calculate_velocity
clearvars -except 'data_names';
close all;clc
end

%%
for ii = 1:1:size(data_names,1)
data_name =  data_names{ii};
load(['./data/',data_name,'_data.mat']);
disp(data_save)
data.ui(ii) = data_save.UI;
data.st(ii) = data_save.st;
data.boshu(ii) = data_save.boshu;
end


% %% 直接从excel读取，没有经过尺寸转换，为像素距离，经过本节转换，转为mm单位
% % data_name = '20240701（橡皮筋月牙形尾鳍）20fps_k2l4ying_1hz_data_xy';
% % data_name = '20240701（橡皮筋月牙形尾鳍）30fps_k2l4ying_2hz_data_xy';
% % data_name = '20240701（橡皮筋月牙形尾鳍）60fps_k2l4ying_4hz_data_xy';
% data_name = '20240701（橡皮筋月牙形尾鳍）60fps_k2l4ying_6hz_data_xy';
% 
% charIndex_hz = regexp(data_name, 'hz','ignorecase');
% DF = str2double(data_name(charIndex_hz-1)) ;   
% charIndex_fps = regexp(data_name, 'fps','ignorecase');
% if isempty(charIndex_fps)
%     FPS = 30;% 数据帧率
% else
%     FPS = str2double(data_name(charIndex_fps-2:charIndex_fps-1));
% end
% scale = 1000/3138.32;
% excel_filename = [data_name,'.xlsx'];
% body_x = readmatrix(excel_filename,'Sheet','x');
% body_y = readmatrix(excel_filename,'Sheet','y');
% body_x = body_x';body_y = body_y';
% body_x = fliplr(-body_x)*scale;body_y = fliplr(-body_y)*scale;
% 
% body_x = body_x - min(min(body_x));
% body_y = body_y - min(min(body_y));
% 
% fig = figure;
% x = body_x(2,:); y = body_y(2,:);
% plot(x,y,"+k-")
% hold on
% body_x = body_x(:,3:38);body_y = body_y(:,3:38);
% x = body_x(2,:); y = body_y(2,:);
% plot(x,y,"ro:")
% legend("full","part",'Location','northeast')
% fig.Position = [2983 770 560 420];
% axis equal
%% 
f = [1;2;4;6];
ui = [0.1924,0.1584;0.1234,0.1100;0.0673,0.0612;0.0432,0.0388];
ui(:,2) = data.ui;
boshu = [0.5931,0.8225;0.7075,0.9467;1.0665,1.0186;1.1498,1.0491];
boshu(:,2) = data.boshu;
st = [3.0748,1.9525;1.0043,0.6414;0.7752,0.4161;0.7855,0.4023;];
st(:,2) = data.st;


fig = figure;
h = plot(f,ui,'LineWidth',5);
set(h(1),'LineStyle', '-', 'Color', 'k', 'MarkerSize', 8, 'MarkerEdgeColor', 'k','LineWidth',3,Marker='square');
set(h(2),'LineStyle', ':', 'Color', 'r', 'MarkerSize', 8, 'MarkerEdgeColor', 'r','LineWidth',3,Marker='o');
xlabel('\itf');
ylabel('UI');
fig.Position = [2987 248 560 420];
legend("full","part")
set(gca,'FontSize',20);


fig = figure;

h = plot(f,1./boshu,'LineWidth',5);
set(h(1),'LineStyle', '-', 'Color', 'k', 'MarkerSize', 8, 'MarkerEdgeColor', 'k','LineWidth',3,Marker='square');
set(h(2),'LineStyle', ':', 'Color', 'r', 'MarkerSize', 8, 'MarkerEdgeColor', 'r','LineWidth',3,Marker='o');
yline([1.23,1.29])
axis([1,6,0,2])
xlabel('\itf');
ylabel('\lambda');
hold on
plot([1,1],[1.0485,0.9922],'bo','HandleVisibility','off')
fig.Position = [2988 -262 560 420];
legend("full","part",'Location','northeast')
set(gca,'FontSize',20);


fig = figure;
h = plot(f,st);
set(h(1),'LineStyle', '-', 'Color', 'k', 'MarkerSize', 8, 'MarkerEdgeColor', 'k','LineWidth',3,Marker='square');
set(h(2),'LineStyle', ':', 'Color', 'r', 'MarkerSize', 8, 'MarkerEdgeColor', 'r','LineWidth',3,Marker='o');
xlabel('\itf');
ylabel('st');
yline(0.45,'--')
yline(0.55,'--')
ylim([0,1])
fig.Position = [2989 -772 560 420];
legend("full","part",'Location','northeast')
set(gca,'FontSize',20);



%% duibi

f = 2:1:6;

st_0 = [0.8448696	0.6574154	0.5269663	0.5349771	0.6696099];
st_1 = [0.8847	0.5576	0.45	0.4363	0.4650];
a_0 = [0.1195908	0.1076185	0.0881109	0.0787526	0.0826291];
a_1 = [0.0972	0.0849	0.0699	0.0585	0.0567];
v_0 = [0.2830990	0.4910982	0.6688164	0.7360375	0.7403934];
v_1 = [0.2196905	0.4565632	0.6324066	0.6702330	0.7320484];

% d_a = (a_1-a_0)./a_0*100;
% d_v = (v_1-v_0)./st_0*100;
% d_st = (st_1-st_0)./st_0*100;
d_a = (a_1)./a_0*100;
d_v = (v_1)./v_0*100;


fig = figure;
subplot(4,1,1)
plot(f,a_0,'r')
hold on
plot(f,a_1,'k')
legend("before","now")
xlabel("f")
ylabel("A")
set(gca,'FontSize',20);

subplot(4,1,2)
plot(f,st_0,'r')
hold on
plot(f,st_1,'k')
legend("before","now")
xlabel("f")
ylabel("ST")
set(gca,'FontSize',20);

subplot(4,1,3)
plot(f,v_0,'r')
hold on
plot(f,v_1,'k')
legend("before","now")
xlabel("f")
ylabel("V")
set(gca,'FontSize',20);

subplot(4,1,4)
plot(f,d_a,'r')
hold on
plot(f,d_v,'k')
legend("d_a","d_v")
xlabel("f")
ylabel("%")
set(gca,'FontSize',20);

fig.Position = [667 217 1370 1114];