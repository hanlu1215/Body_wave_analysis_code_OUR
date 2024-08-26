% reference��John H. Long Jr ��Flapping flexible fishPeriodic and secular body reconfigurations in swimming lamprey, Petromyzon marinus��
% by Bingxing Chen 20181121
%   x/body_x��the axial direction  y/body_y��the lateral direction (not rotated)
% clc
% clear
% close all
hll
disp('________________________________________________________________')
disp('Calculate the velocity and amplitude of fish body')
data_path_root = "./data_202408/";
save_path_root = "./result_202408/";
mkdir_han(save_path_root);
% ֱ�Ӵ�excel��ȡ��û�о����ߴ�ת����Ϊ���ؾ��룬��������ת����תΪmm��λ
data_name = '20240816(���ɸ������ζ�β)_30fps_k33334_4hz_data_xy';
%%
excel_filename = data_path_root+data_name+".xlsx";
data_msg = readmatrix(excel_filename,'Sheet','config');%fps, lin_pix, dn, data_fps
scale = 1000/data_msg(2);
% scale = 1000/3138.32;
% scale = 1000/2815.000710479484;
data_fps = data_msg(4);
charIndex_hz = regexp(data_name, 'hz','ignorecase');
DF = str2double(data_name(charIndex_hz-1)) ;   
data_save.f = DF;
FPS = data_fps;
% charIndex_fps = regexp(data_name, 'fps','ignorecase');
% if isempty(charIndex_fps)
%     FPS = 60/3;% ����֡��
% else
%     FPS = str2double(data_name(charIndex_fps-2:charIndex_fps-1));
% end
body_x = readmatrix(excel_filename,'Sheet','x');
body_y = readmatrix(excel_filename,'Sheet','y');
% ͷ���ں����ݽ�ȡ��
% body_x = body_x(:,1:end-3);body_y = body_y(:,1:end-3);
body_x = body_x';body_y = body_y';
body_x = fliplr(-body_x)*scale;body_y = fliplr(-body_y)*scale;
% body_x = fliplr(-body_x);body_y = fliplr(-body_y); %���ҶԳ�
%% ���ݳ�������
% ��ȡ�ļ��ļ�Ҫ��Ϣ����Ϊͼ��ı���
charIndex = regexpi(data_name, 'k');
data_name_simple = strrep(data_name(charIndex(1):end-8), '_', '-');
n_points = size(body_y,2);
% 51 pixel represents 100mm
% Fish_length_calibration=51.5*3.9;
Fish_length_calibration = 550;%���峤�ȣ���λ����

% ����������ֵ��
body_x = body_x - min(min(body_x));
body_y = body_y - min(min(body_y));
%% �������е�
figure(100)
for nn = 1:1:size(body_y,1)
plot(body_x(nn,:),body_y(nn,:)+nn*10);
hold on
end
hold off
%%
% DF=2;timestep = 1/30;num_per_cycle=round(1/DF/timestep)+1;
timestep = 1./FPS;% �������ߵ�ʱ����
num_per_cycle=round(1/DF/timestep)+1;% ÿ�����ڵ�������
NumFr=size(body_x,1)-num_per_cycle;% ���ٸ�����
% NumFr=1;
kkkkk = 0;%
U.x=[];
for strt_frm=1:NumFr% different cycles
    disp(['strt_frm: ',num2str(strt_frm)])
%     ȡ����n�����ڵ�����
    Data.one.x =body_x(strt_frm:strt_frm+1*num_per_cycle,:)';%% row: body points line:time series
    Data.one.y = body_y(strt_frm:strt_frm+1*num_per_cycle,:)';
    fishL =[];
    sz=size(Data.one.x);%%sz ��ȡ����ֵ sz(1)�У���ʾ�����λ�� �峤������ �ж��ٸ�����ĵ� sz(2)�У�ʱ��
%     ʹ���������ݣ��������峤��
    for il=1:sz(2)
        fishL(il) =sum(sqrt(diff(Data.one.x(:,il)).^2+diff(Data.one.y(:,il)).^2));
    end
    fishLength=mean(fishL); %the body length

    FittData.one.H=[];FittData.H0=[];
    %% obtain the rotation angle through iterations ͨ�����������ת�Ƕ�
    Data.one.X =Data.one.x;Data.one.Y =Data.one.y;%x,yΪԭʼ��
    iteration.flag=0;%
    iteration.eps=1e-2;%iteration precision ��������
    while(1)
        FittData.one.H=[];
        for i=1:num_per_cycle+1%one cycle
            % FittData.one.H0=[1 i*timestep (i*timestep)^2/2  sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep) cos(2*DF*2*pi*i*timestep) sin(3*DF*2*pi*i*timestep) cos(3*DF*2*pi*i*timestep) sin(4*DF*2*pi*i*timestep) cos(4*DF*2*pi*i*timestep)  ];
            FittData.one.H0=[1 i*timestep sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) ];% the preliminary Fourier transforms.. and it just needs the velocity term and fundamental harmonic��������Ҷ�任..��ֻ��Ҫ�ٶ���ͻ�г��
            FittData.one.H=[FittData.one.H;FittData.one.H0;];%
        end
%         ��С������Ϸ�,(9)ʽ��C����г������ϵ������
        Data.one.xInverse=Data.one.X';Data.one.yInverse=Data.one.Y';
        FittData.one.S_u=inv(FittData.one.H'*FittData.one.H)*FittData.one.H'*(Data.one.xInverse);%least square fitting method
        FittData.one.S_v=inv(FittData.one.H'*FittData.one.H)*FittData.one.H'*(Data.one.yInverse);%least square fitting method

        FittData.one.u_model=[];FittData.one.v_model=[];
        FittData.one.u_model=FittData.one.H*FittData.one.S_u;% Hc*C
        FittData.one.v_model=FittData.one.H*FittData.one.S_v;
%         �õ����ڲ���������
        FittData.one.u_model_Periodic=FittData.one.H(:,3:end)*FittData.one.S_u(3:end,:)+FittData.one.H(:,1)*FittData.one.S_u(1,:);
        FittData.one.v_model_Periodic=FittData.one.H(:,3:end)*FittData.one.S_v(3:end,:)+FittData.one.H(:,1)*FittData.one.S_v(1,:);
        FittData.P=polyfit(FittData.one.u_model_Periodic,FittData.one.v_model_Periodic,1);% ��ϵõ�����������
        FittData.f =polyval(FittData.P,FittData.one.u_model_Periodic);% �������ݵ����������ϵ���������ģ���ϵĵ�����
        % % plot to verify
        figure(101)
        plot(Data.one.xInverse,Data.one.yInverse,'r',FittData.one.u_model,FittData.one.v_model,'k')%%���岨��ͼ the all fish body;
        title('body wave model and actualbody')
        axis equal
        figure(102)%% the 'eight' of figure
        plot(FittData.one.u_model_Periodic(1:end,n_points)',FittData.one.v_model_Periodic(1:end,n_points)','r') %ת��֮ǰ�����岨����
        title('tranjectories of body points')
        if strt_frm == 1
            fig1 = figure(1);
        else
        figure(strt_frm)
        end
        for iiii = 1:size(Data.one.xInverse,1)%%˵�����ж����� Ҳ����  ��ͬʱ�������ֵ
            plot(FittData.one.u_model_Periodic(iiii,:)',FittData.one.v_model_Periodic(iiii,:)','r') %ת��֮ǰ�����岨����
            title({'periodic components of midlines....',data_name_simple})
            hold on
        end
        hold on
        plot(FittData.one.u_model_Periodic,FittData.f,'k')% the fitting curve ��ϵõ���������
        hold off
        axis equal
        xlim([min(min(FittData.one.u_model_Periodic)),max(max(FittData.one.u_model_Periodic))])
        ylim([min(min(FittData.one.v_model_Periodic)),max(max(FittData.one.v_model_Periodic))])

%         pause(0.5)
        kkkkk = kkkkk+1;
        disp(kkkkk);
%         ͨ������ϵ������ߵ�б���жϣ��Ƿ��Ѿ�ƽֱ
        if abs(atan(FittData.P(1)))<iteration.eps%%
            iteration.flag=1;
            break;%break the while
        end
%         �õ���ת�Ƕȣ�����תԭʼ����
        FittData.theta(strt_frm)=atan(FittData.P(1));% ����Ӧ��Ҫ��ת�ĽǶ�
        Data.one.X = cos(FittData.theta(strt_frm))*Data.one.x+sin(FittData.theta(strt_frm))*Data.one.y;%%%the axial direction  Rotate raw data
        Data.one.Y= -sin(FittData.theta(strt_frm))*Data.one.x+cos(FittData.theta(strt_frm))*Data.one.y;%the axial direction    Rotate raw data������ת
    end
    % after rotating the raw data
    %% Fourier fitting in order to calculate the velocity and tail amplitude����Ҷ����Լ����ٶȺ�β�����
    FittData.one.H=[];
    for i=1:num_per_cycle+1%one cycle
        % %new
        FittData.one.H0=[1 i*timestep sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep) cos(2*DF*2*pi*i*timestep) sin(3*DF*2*pi*i*timestep) cos(3*DF*2*pi*i*timestep) sin(4*DF*2*pi*i*timestep) cos(4*DF*2*pi*i*timestep)  ];
        % FittData.one.H0=[1 i*timestep sin(DF*2*pi*i*timestep) cos(DF*2*pi*i*timestep) sin(2*DF*2*pi*i*timestep)  ];
        FittData.one.H=[FittData.one.H;FittData.one.H0;];%
    end
    Data.one.xInverse=Data.one.X';Data.one.yInverse=Data.one.Y';
    FittData.one.S_u=inv(FittData.one.H'*FittData.one.H)*FittData.one.H'*(Data.one.xInverse);%least square fitting method��С������Ϸ�,(9)ʽ��C����г������ϵ������
    FittData.one.S_v=inv(FittData.one.H'*FittData.one.H)*FittData.one.H'*(Data.one.yInverse);%least square fitting method
    FittData.one.u_model=[];FittData.one.v_model=[];
    FittData.one.u_model=FittData.one.H*FittData.one.S_u;FittData.one.v_model=FittData.one.H*FittData.one.S_v;
    FittData.one.u_model=FittData.one.u_model';
    FittData.one.v_model=FittData.one.v_model';

    U.calculateFit(:,strt_frm)=FittData.one.S_u(2,:);%%the velocity term of the fitting curve in the axial direction (u). �������������ٶ���  it is the velocity of fish body
    % verify
    U.calculate(:,strt_frm)=FPS*abs(FittData.one.u_model(:,end)-FittData.one.u_model(:,1))/size(FittData.one.u_model,2);%%���е���ٶ�
    U.head(strt_frm)=FPS*abs(FittData.one.u_model(end,end)-FittData.one.u_model(end,1))/size(FittData.one.u_model,2);%%the velocity of headͷ���ٶ�
    %% but the tail amplitude should be recalculated. because the dundamental amplitude is not the entire amplitude..
    % ��Ӧ����β���������Ϊ��������������������
    A.cycle(strt_frm) =(max(FittData.one.v_model(1,:))-min(FittData.one.v_model(1,:)));%%��β��Y����  1��ʾ��β

    % for i_actual= 1:n_points%%body points
    % % U.x(i_actual)=30*(Data.one.x(i_actual,end)-Data.one.x(i_actual,1))/fishLength/ size(Data.one.x,2);%%velocity of x direction
    % % U.y(i_actual)=30*(FittData.one.v_model(i_actual,end)-FittData.one.v_model(i_actual,1))/fishLength/ size(FittData.one.v_model,2);%%velocity of y direction
    % % the velocity also can use the actual velocity by the calibration for the actual fish length
    % % % % U.one.x(i_actual)=30*(FittData.one.u_model(i_actual,end)-FittData.one.u_model(i_actual,1))/ size(FittData.one.u_model,2);%%x velocities of all body points in one cycle
    % % % % U.one.y(i_actual)=30*(FittData.one.v_model(i_actual,end)-FittData.one.v_model(i_actual,1))/ size(FittData.one.v_model,2);%%y velocity of all body points in one cycle
    % %record velocities of all body points and all cycles
    % U.all.x(i_actual,strt_frm)=U.one.x(i_actual);%
    % U.all.y(i_actual,strt_frm)=U.one.y(i_actual);%
    % end
    %the x/y components of the average velocity for individual body points
    % % U.one.mean=[mean(U.one.x(:)) mean(U.one.y(:))];% the value of the composite average velocity in one cycle.
    % % U.one.xy=sqrt(U.one.x.^2+U.one.y.^2);
    % % RE.one.x=(U.one.x-mean(U.one.x(:)))/norm(U.one.mean);RE.one.y=(U.one.y-mean(U.one.y(:)))/norm(U.one.mean);
    % % RE.one.Square=sqrt(RE.one.x.^2+RE.one.y.^2);RE.all.Square(strt_frm,:)=RE.one.Square;
    % % % The unsteadliness index in one cycle
    % % UI.one=RE.one.Square.^2*norm(U.one.mean)./U.one.xy;UI.all(strt_frm,:)=UI.one;%The unsteadliness index of all cycles
    % % %calculate theta by the composite average velocity in one cycle.
    % % theta.one.U=atan2(mean(U.one.y(:)),mean(U.one.x(:)));theta.all.U(strt_frm)=theta.one.U;
    % %
end
%%
% U�����ٶȣ�����������ٶȣ�AΪβ���ڷ�
% UI.mean=mean(UI.all(:));%The means of unsteadliness index
U.headmean=mean(U.head);
U.head_tailmean=mean(mean(U.calculate));
U.mean=mean(mean(U.calculateFit));
Stride_length=U.mean/DF;
% U.meanstandDeviation=sqrt(sum(((U.calculateFit-U.mean)/U.mean).^2)/size(U.calculateFit,2));
U.meanstandDeviation=sqrt(mean(sum(((U.calculateFit-U.mean)/U.mean).^2)/size(U.calculateFit,2)));
A.mean=mean(A.cycle);
A.standDeviation=sqrt(sum(((A.cycle-A.mean)/A.mean).^2)/size(A.cycle,2));
% The unsteadliness index in one cycle
St.mean=DF*A.mean/U.mean;
St.meanhead_tail=DF*A.mean/U.head_tailmean;

disp( '************The calculations*************');
fprintf(1,'U.mean is (the velocity term of the fitting curve in the axial direction (u)): %1.4f\n',U.mean/Fish_length_calibration)%%we use this data
fprintf(1,'The stand Deviation of U.mean of 30 body points in this trial(many cycles)is: %3.4f\n',U.meanstandDeviation)
fprintf(1,'A.mean  is: %3.4f\n',A.mean/Fish_length_calibration)
fprintf(1,'The stand Deviation of A.mean of 30 body points in this trial(many cycles) is: %3.4f\n',A.standDeviation)
fprintf(1,'St.mean is: %3.4f\n',St.mean)
fprintf(1,'Stride length is: %3.4f\n',Stride_length/Fish_length_calibration)
disp( '**********comparision to verify***************');
fprintf(1,'U.headmean (just the velocity of head point) (pixel/s) is: %3.4f\n',U.headmean/Fish_length_calibration)
fprintf(1,'St.meanhead_tail is: %3.4f\n',St.meanhead_tail)
fprintf(1,'U.head_tailmean (pixel/s) is: %3.4f\n',U.head_tailmean/Fish_length_calibration)
data_body=[DF U.mean/Fish_length_calibration A.mean/Fish_length_calibration St.mean  Stride_length/Fish_length_calibration];

data_save.v = U.mean/1000;
data_save.A = A.mean/1000;
data_save.st = St.mean;
data_save.data_name = data_name;
data_save.data_name_simple = data_name_simple;
save_name = save_path_root+data_name+"_fig1.fig";
saveas(fig1,save_name,"fig")

fprintf(1,'\nv: %1.4fm/s     ',U.mean/1000)%%we use this data
fprintf(1,'A: %3.4fm     ',A.mean/1000)
fprintf(1,'St: %3.4f\n',St.mean)
% ��ȡ���б�����
vars = who;
% Ҫ�ų��ı�������
exclude_type = 'matlab.ui.Figure';
% ����Ҫ����ı����б��ų�ָ�����͵ı���
vars_to_save = {};
for i = 1:length(vars)
    var_name = vars{i};
    if ~strcmp(class(eval(var_name)), exclude_type)
        vars_to_save{end+1} = var_name;
    end
end
% �������
save('bd_data.mat', vars_to_save{:});

%%
figure(1)
for iiii = 1:4%%˵�����ж����� Ҳ����  ��ͬʱ�������ֵ
    plot(FittData.one.u_model_Periodic(iiii,:)',FittData.one.v_model_Periodic(iiii,:)') %ת��֮ǰ�����岨����
    title({'periodic components of midlines....',data_name_simple})
    hold on
end
hold off
axis equal
xlim([min(min(FittData.one.u_model_Periodic)),max(max(FittData.one.u_model_Periodic))])
ylim([min(min(FittData.one.v_model_Periodic)),max(max(FittData.one.v_model_Periodic))])
Modeling_steady_undulatory_motion