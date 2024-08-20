hll
%%
% [file_path,data_name] = selectFolderAndGetLastNLevels(3);
% % Excel文件名
% excel_filename = [data_name,'.xlsx'];
% filename = [file_path,'/3/'];

% folderPath = uigetdir;
folderPath = 'F:\SHARE_FOR_LZQ\游动曲线\20240701（橡皮筋月牙形尾鳍）\k2l4ruan';
for exp = [4,6]
    file_path = [folderPath,'\',num2str(exp),'hz'];
    filename = [file_path,'\3\'];
    data_name = getLastNLevels(file_path, 3);
    excel_filename = [data_name,'_data_xy.xlsx'];

    disp(filename)
    disp(excel_filename)
    save_to_excel(filename,excel_filename)
end


%%自定义函数
function save_to_excel(filename,excel_filename)
% scale = 1000/3138.32;
scale = 1;
body_x = readmatrix([filename,'save_datax.csv']);
body_y = readmatrix([filename,'save_datay.csv']);

body_x = body_x*scale;
body_y = body_y*scale;
% body_x = fliplr(-body_x);body_y = fliplr(-body_y); %左右对称
% body_x = body_x(:,5:34);
% body_y = body_y(:,5:34);
body_x = body_x - min(min(body_x));
body_y = body_y - min(min(body_y));

%% excel
% 生成标题行
headers = arrayfun(@(x) sprintf('Line_%d', x), 1:size(body_y,1), 'UniformOutput', false);

% 准备数据
var1 = body_x'; %
var2 = body_y'; %

% 将矩阵转换为表格并添加标题
T1 = array2table(var1, 'VariableNames', headers);
T2 = array2table(var2, 'VariableNames', headers);

% 将表格变量写入Excel的不同工作表中
writetable(T1, excel_filename, 'Sheet', 'x','WriteMode','overwrite'); % 将T1写入Sheet1
writetable(T2, excel_filename, 'Sheet', 'y','WriteMode','overwrite'); % 将T2写入Sheet2
disp('save fine!')
end

function [folderPath,resultString] = selectFolderAndGetLastNLevels(n)
% 让用户选择一个文件夹
folderPath = uigetdir;

% 检查用户是否取消了选择
if folderPath ~= 0
    % 获取路径的最后n级组成的字符串
    resultString = getLastNLevels(folderPath, n);

    % 打印结果字符串
    disp(['最后 ', num2str(n), ' 级路径组成的字符串是: ', resultString]);
else
    disp('用户取消了选择');
    resultString = '';
end
end

function resultString = getLastNLevels(path, n)
% 将路径分割为各级
pathLevels = strsplit(path, filesep);

% 取最后n级
numLevels = length(pathLevels);
if n > numLevels
    n = numLevels;
end
lastNLevels = pathLevels(numLevels - n + 1 : numLevels);

% 用下划线连接最后n级
resultString = strjoin(lastNLevels, '_');
end

