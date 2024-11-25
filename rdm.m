%%  radarDataCu+be: 数据立方体 (numFrame x numRx*numTx x numChirp x numSampling)
%   numFrame:帧数
%   numRX*numTx：通道数，因为现在就做一发多收，所以就等于numRx--------numChannel
%   numChirp：chirp数
%   numSampling：chirp的采样点数-------numPoint

%%  准备所需的参数
c = 3e8;%光速
S = bandwidth/timeChrip;%频率斜率
lambda = c / frequency;
range_res = c / (2 * bandwidth);%距离分辨率
time_Chp = timeChrip + timeChripGap;
doppler_res = lambda / (2 *time_Chp * double(numChrip));%速度分辨率
d = lambda /2;%天线距离
angle_res = lambda / 2 / d;%角度分辨率
Q = 180;%角度FFT
fs = numChrip / timeChrip;%采样频率
%CFAR参数
train_range = 1;
train_doppler = 1;
guard_range = 1;
guard_doppler = 1;
threshold = 1.9;%门限阈值
%换个格式
numChrip = double(numChrip);
numPoint = double(numPoint);
%%  测距-- 通过计算每个采样点的相位，取平均值来得到每个chirp的相位，假设初始相位为0，以此获得相位差，然后得到时间t算出距离 2d = ct
%   上面的思路不对，应该是每个采样点对应一个频率点f，然后由θ = 2*pi*f*t得到相位，然后平均得到chirp的相位
%   仍然要改思路，现在就只有接收数据，得到中频信号之后对中频信号做2dfft，这是得到幅度谱的方法。
%   对距离维对接收信号进行fft，频谱中峰值频率就是表示🔺f，即可R =
%   C*🔺f/2/S,也*-就是差频出来的中频信号，现在不用将发送信号和接收信号做差得到，直接做fft得到即可，这又叫差频、拍频 
rdm_orign = zeros(size(radarDataCube));

target_range = zeros(numFrame);
target_spec = zeros(numFrame);
target_angle = zeros(numFrame);



%% RDM
for channel = 1:numChannel
    for frame = 1:numFrame

        range_dopplerfft = fftshift(fft2(squeeze(radarDataCube(frame,channel,:,:))),1);%2dfft后取一个通道一帧数据
        rdm_orign(frame,channel,:,:) = range_dopplerfft;
        RDM_dB = 20 * log10(abs(range_dopplerfft));
        RDM_AMP = abs(range_dopplerfft);
        %返回2dfft之后第三维和第四维的长度作为坐标，也就是chirp（速度）和point（距离）
        % RDM画图
%         
%         dopplerAxis = linspace(-numChrip/2, numChrip/2 - 1, numChrip); % Doppler 轴
%         rangeAxis = linspace(0, numPoint-1, numPoint); % Range 轴
%         [dopplerAxis,rangeAxis] = meshgrid(dopplerAxis,rangeAxis);
%         dopplerAxis = transpose(dopplerAxis);
%         rangeAxis = transpose(rangeAxis);
%         mesh(dopplerAxis,rangeAxis,RDM_AMP);
%         hold on
%         xlabel('doppler');
%         ylabel('range');
%         zlabel('Amp(db)');
%         title(['Frame ' num2str(frame) ', Channel ' num2str(channel)]);
%         colorbar;
%         hold off
%         view(45,45);
%         pause(0.000001)
    end
end
%% 测距测速测角（CFAR）rdm,train_range,train_doppler,guard_range,guard_doppler,numPoint,numChirp
for frame = 86:numFrame
    count = 1;
    axis_x_frame = zeros(numChrip*numPoint,1);
    axis_y_frame = zeros(numChrip*numPoint,1);
    spec_frame = zeros(numChrip*numPoint,1);
    RDM_Frame = squeeze(rdm_orign(frame,channel,:,:));
%    [target,noise_lever] = cfar_2d(RDM_Frame,train_range,train_doppler,guard_range,guard_doppler,numPoint,numChrip,2);
   [target,noise_lever] = cfar_2d_cross(RDM_Frame,train_range,train_doppler,guard_range,guard_doppler,numPoint,numChrip,threshold);
   for rdm_target_xlabel = 1:numChrip
       for rdm_target_ylabel = 1:numPoint
           if target(rdm_target_xlabel,rdm_target_ylabel) == 1
                phase_delta = angle(squeeze(rdm_orign(frame,1,rdm_target_xlabel,rdm_target_ylabel))) - angle(squeeze(rdm_orign(frame,2,rdm_target_xlabel,rdm_target_ylabel)));%求两通道相位差(单位;弧度)
                 % 控制相位
                [p,q] = find(phase_delta < -pi);
                phase_delta(p,q) = phase_delta(p,q) + 2 * pi;
                [p,q] = find(phase_delta > pi);
                phase_delta(p,q) = phase_delta(p,q) - 2 * pi;
                % 测距测速测角
                target_angle(frame) = asin(phase_delta * lambda  / (2 * pi * d));
                target_range(frame) = (rdm_target_ylabel - 1) * range_res;
                target_spec(frame) = (rdm_target_xlabel - numChrip / 2 - 1) * doppler_res;
                
                axis_x = target_range(frame) * cos(target_angle(frame));
                axis_y= target_range(frame) * sin(target_angle(frame));
%                 plot(axis_x,  axis_y, 'ro', 'MarkerSize', 10, 'DisplayName', '目标点'); % 绘制目标点
%                 hold on
                axis_x_frame(count) = axis_x;
                axis_y_frame(count) = axis_y;
                spec_frame(count) = (rdm_target_xlabel - numChrip / 2 - 1) * doppler_res;
                count = count + 1;
           end
       end
   end
    target_axis = [axis_x_frame, axis_y_frame];%得到CFAR出来的所有目标的直角坐标系坐标
    spec_frame= spec_frame(1:count-1);
    spec_frame = spec_frame' - spec_frame;%算出每个速度之间的差值
    target_axis(target_axis==0) = [];
    target_axis=reshape(target_axis,[],2);
    [target_clustered,isnoise] = DBSCAN(target_axis,spec_frame,0.56,1,2);%DBSCAN聚类
    
    plot(tergatTrajectory(2,frame,1), tergatTrajectory(2,frame,2), 'bs', 'MarkerSize', 10, 'DisplayName', '参考点'); % 绘制参考
    hold on
    plot(tergatTrajectory(1,frame,1), tergatTrajectory(1,frame,2), 'gs', 'MarkerSize', 10, 'DisplayName', '参考点'); % 绘制参考
    hold on
    PlotClusterinResult(target_axis,target_clustered);
    hold off
    xlabel('X 坐标 (m)');
    ylabel('Y 坐标 (m)');
    title(['Frame ' num2str(frame)]);
    xlim([0, 20]); % x 轴范围
    ylim([-10, 10]); % y 轴范围
    pause(0.01)
    

    
end


%% 目标点聚类
% 思路应该是每一帧聚类一次，然后把所有点画在一张图里做对比
% 对每一帧，做一次DBSCAN聚类，一次K-means聚类
% 不用k-means，重复做了一些工作，分出各个簇之后，算各个簇的平均值。不过需要思考一下幅值也就是能量的影响，因为大的物体返回来的能量就大嘛
% 首先是DBSCAN聚类，前期准备需要CFAR出来的目标点坐标



%% DBSCAN 聚类算法
% 输入CFAR_Data是二维数组，第一列是多普勒维，第二列是距离维
function [C_index,isnoise] = DBSCAN(CFAR_Data,Spec_Data,spec_epsilon,epsilon,MinPts)
    cfar_point_nums = size(CFAR_Data,1);%目标点数
    distance = pdist2(CFAR_Data, CFAR_Data);%所有点之间的距离
    C = 0;%簇数
    C_index = zeros(cfar_point_nums,1);%簇的索引，也就是目标点的索引

    visited = false(cfar_point_nums,1);% 布尔数组
    isnoise = false(cfar_point_nums,1);
    
    for i = 1:cfar_point_nums
        if ~visited(i)
        visited(i) = true;
         Neighbors = find(distance(i,:) <= epsilon); % 找到第 i 个点的邻居,特别注意这里Neighbors是存的邻居点的索引
         Neighbors_spec = find(abs(Spec_Data(i,:)) <= spec_epsilon);% 速度邻域半径
         Neighbors = intersect(Neighbors,Neighbors_spec); % 求符合两个条件的点
         
            if numel(Neighbors) < MinPts
                % 标记为噪声点
                isnoise(i) = true;
            else
                % 启动新簇
                C = C + 1;
                
                C_index(i) = C;%标记索引为i的点的所属簇
                k = 1;
                while true
                  j = Neighbors(k);
                    if ~visited(j)
                    visited(j) = true;
                    Neighbors2 = find(distance(j) <= epsilon);
                        if numel(Neighbors2) >= MinPts
                            Neighbors = [Neighbors Neighbors2];
                        end
                    end
                    if C_index(j) == 0
                        C_index(j) = C;%更新新邻居点的所属簇
                    end
                        k = k + 1;
                        
                    if k > numel(Neighbors)%所有邻居遍历完，且没有新的邻居点了
                    break;
                    end
                end
            end
        end
    end


end


%% 测距测速测角(最大值)
% for frame = 1:numFrame
%     if(max(abs(squeeze(rdm_orign(frame,1,:,:)))) < max(abs(squeeze(rdm_orign(frame,2,:,:)))))
%         [~,linearIndex] = max(abs(squeeze(rdm_orign(frame,2,:,:))),[], 'all','linear');%用linear模式获取线性索引，不然可能存在多个相同的最大值，返回一组向量
%         [rdm_max_xlabel, rdm_max_ylabel] = ind2sub(size(squeeze(rdm_orign(frame,2,:,:))), linearIndex);
%     else
%         [~,linearIndex] = max(abs(squeeze(rdm_orign(frame,1,:,:))),[], 'all','linear');
%         [rdm_max_xlabel, rdm_max_ylabel] = ind2sub(size(squeeze(rdm_orign(frame,1,:,:))), linearIndex);
%     end
%     
%     phase_delta = angle(squeeze(rdm_orign(frame,1,rdm_max_xlabel,rdm_max_ylabel))) - angle(squeeze(rdm_orign(frame,2,rdm_max_xlabel,rdm_max_ylabel)));%求两通道相位差(单位;弧度)
%     % 控制相位
%     [p,q] = find(phase_delta < -pi);
%     phase_delta(p,q) = phase_delta(p,q) + 2 * pi;
%     [p,q] = find(phase_delta > pi);
%     phase_delta(p,q) = phase_delta(p,q) - 2 * pi;
%     % 测距测速测角
%     target_angle(frame) = asin(phase_delta * lambda  / (2 * pi * d));
%     target_range(frame) = (rdm_max_ylabel - 1) * range_res;
%     target_spec(frame) = (rdm_max_xlabel - numChrip / 2 - 1) * doppler_res;
%     
%     % 平面坐标画图
%     axis_x(numFrame, numChannel) = target_range(frame) * cos(target_angle(frame));
%     axis_y(numFrame, numChannel) = target_range(frame) * sin(target_angle(frame));
%     plot(axis_x(numFrame, numChannel),  axis_y(numFrame, numChannel), 'ro', 'MarkerSize', 10, 'DisplayName', '目标点'); % 绘制目标点
%     hold on
%     plot(tergatTrajectory(2,frame,1), tergatTrajectory(2,frame,2), 'bs', 'MarkerSize', 10, 'DisplayName', '参考点'); % 绘制参考
%     hold on
%     plot(tergatTrajectory(1,frame,1), tergatTrajectory(1,frame,2), 'gs', 'MarkerSize', 10, 'DisplayName', '参考点'); % 绘制参考
%     hold off
%     
% 
%     xlabel('X 坐标 (m)');
%     ylabel('Y 坐标 (m)');
%     title(['Frame ' num2str(frame)]);
% 
%     xlim([0, 20]); % x 轴范围
%     ylim([-10, 10]); % y 轴范围
% 
% pause(0.000001)
% end
%% CFAR 十字形
% 在一个rdm上找出目标。
% rdm为[chirp,point]维度的二维数组，是一帧rdm的幅值数组，行代表doppler，列代表range。
function [target, noise_lever] = cfar_2d_cross(rdm,train_range,train_doppler,guard_range,guard_doppler,numPoint,numChirp,threshold)
% 构建相对于检测单元的训练单元和保护单元
    target = zeros(numChirp, numPoint);
    noise_lever = zeros(numChirp, numPoint);

for i = (guard_doppler + train_doppler + 1):(numChirp - (guard_doppler + train_doppler))
    for j = (train_range + guard_range + 1):(numPoint - (train_range + guard_range))
        
        % 训练单元内信号
        signal_train_left = sum(abs(rdm((i-(guard_doppler + train_doppler)):(i-guard_doppler-1),j)));
        signal_train_right = sum(abs(rdm((i+guard_doppler+1):(i+guard_doppler+train_doppler),j)));
        signal_train_up = sum(abs(rdm(i,(j-(guard_range + train_range)):(j - guard_range-1))));
        signal_train_dowm = sum(abs(rdm(i,(j+guard_range+1):(j+guard_range + train_range))));
        signal_train = signal_train_left + signal_train_right + signal_train_up + signal_train_dowm;

        % 估计噪声水平
        noise_lever(i,j) = signal_train/(2 * train_range + 2 * train_doppler);
        % 得到目标点
        if(rdm(i,j)/noise_lever(i,j) > threshold)
            target(i,j) = 1;
        end
    end
end
end

%% CFAR
% 在一个rdm上找出目标。
% rdm为[chirp,point]维度的二维数组，是一帧rdm的幅值数组，行代表doppler，列代表range。
function [target, noise_lever] = cfar_2d(rdm,train_range,train_doppler,guard_range,guard_doppler,numPoint,numChirp,threshold)
% 构建相对于检测单元的训练单元和保护单元
    target = zeros(numChirp, numPoint);
    noise_lever = zeros(numChirp, numPoint);

for i = (guard_doppler + train_doppler):(numChirp - (guard_doppler + train_doppler))
    for j = (train_range + guard_range):(numPoint - (train_range + guard_range))
        % 总信号
        signal_sum = rdm((i-(guard_doppler + train_doppler) + 1):(i+(guard_doppler + train_doppler)),(j-(guard_range + train_range)+1):(j+(guard_range + train_range)));
        % 保护单元内信号
        signal_guard = rdm((i-guard_doppler+1):(i+guard_doppler),(j-guard_range+1):(j+guard_range));
        % 训练单元内信号
        signal_train = sum(sum(signal_sum)) - sum(sum(signal_guard));
        % 估计噪声水平
        noise_lever(i,j) = signal_train/((guard_doppler + train_doppler)*(guard_range + train_range) - guard_doppler * guard_range);
        % 得到目标点
        if(rdm(i,j)/noise_lever(i,j) > threshold)
            target(i,j) = 1;
        end
    end
end

end
%% 测角错误尝试

% 用3dfft来测角
% radar_anglefft = fftshift(fft(radar_dopplerfft,[], 2));
% for frame = 1:numFrame
%     radar_anglefft_frame = squeeze(radarDataCube(frame,:,:,:));
%     radar_anglefft = fftn(radar_anglefft_frame);
%     [angle_max,linearIndex] = max(radar_anglefft(:));
%     [row,col,angle_index(frame)] = ind2sub(size(radar_anglefft),linearIndex);
%     
% end

% 用angle来直接算两个chirp之间的相位差
% for frame = 1:numFrame
% 
%         rx1_data = squeeze(radarDataCube(frame, 1, :, :));  % 第一个接收天线的数据
%         rx2_data = squeeze(radarDataCube(frame, 2, :, :));  % 第二个接收天线的数据
%         phase_rx1 = angle(rx1_data);  % 第一个天线的相位
%         phase_rx2 = angle(rx2_data);  % 第二个天线的相位
%         phase_diff = phase_rx2 - phase_rx1;
%         target_theta(frame) = asin(mean(phase_diff(:)) * lambda *  2 / pi / d);  % 相位差
% 
% end


% 用互相关来计算延时
% for frame = 1:numFrame
%     for chirp = 1:numChrip
%         rx1_data = squeeze(radarDataCube(frame, 1,chirp, :));  % 第一个接收天线的数据
%         rx2_data = squeeze(radarDataCube(frame, 2,chirp, :));  % 第二个接收天线的数据
%         [corr,lag] = xcorr(rx1_data,rx2_data);
%         [peaks,peak_positions] = max(abs(corr));
%         %计算相位差
%         timeDelay = lag(peak_positions);
%         delta_phi(chirp) = asin(timeDelay * c / d);
% %       delta_phi(chirp)=(peak_position / length(rx1_data)) * 2 * pi;%计算相位差
%         
% 
% 
%        
%     end
%      target_theta(frame) = asin(mean(delta_phi(:)) * lambda *  2 / pi / d);  % 相位差
% end
%% 聚类后画图
function PlotClusterinResult(X, IDX)
    k=max(IDX);
    Colors=hsv(k);
    Legends = {};
    for i=0:k
        Xi=X(IDX==i,:);
        if i~=0
            Style = 'x';
            MarkerSize = 8;
            Color = Colors(i,:);
            Legends{end+1} = ['Cluster #' num2str(i)];
        else
            Style = 'o';
            MarkerSize = 6;
            Color = [0 0 0];
            if ~isempty(Xi)
                Legends{end+1} = 'Noise';
            end
        end
        if ~isempty(Xi)
            plot(Xi(:,1),Xi(:,2),Style,'MarkerSize',MarkerSize,'Color',Color);
        end
        hold on;
    end
    hold off;
    axis equal;
    grid on;
    
    legend('Location', 'NorthEastOutside');
end
