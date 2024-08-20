hll
% save_path_root = "./result_202408/";
save_path_root = "./result/";
% 获取当前文件夹内所有指定扩展名的文件信息
fileInfo = dir(save_path_root+"20240718（弹簧钢不均匀刚度月牙形）_"+"*.mat");
% 提取文件名字
data_names = {fileInfo.name};
% sort_my = 1:1:numel(data_names);
sort_my = [5 3 4 1 2 6];

fig = figure;hold on
my_line_style = get_my_style(10);
my_marker_size = 8;
my_LineWidth = 3;
for iii =1:1:numel(sort_my)
    ii = sort_my(iii);
    disp({ii,data_names{ii}})
    mat_filename = save_path_root+data_names{ii};
    load(mat_filename,"data_save")
    data1 = data_save;
%     subplot(nn,1,ii)
    PLOT_data = (data_save.jipinxiangwei-data_save.jipinxiangwei(1))*180/pi;
    plot(PLOT_data,my_line_style{ii},'Markersize',my_marker_size,'LineWidth',my_LineWidth,DisplayName=data_save.data_name_simple)
    set(gca,'FontSize',24);
end

%%
fig.Position = [251 116 1969 1197];
xlabel("bl")
ylabel("deg")
legend Location northwest

