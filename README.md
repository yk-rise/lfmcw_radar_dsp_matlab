# lfmcw_radar_dsp_matlab
# 简介
这是学习调频连续波雷达信号处理所建立的库，储存从最开始ADC采集到的数据一步一步得到点云数据、聚类算法等处理的matlab仿真程序，README里也会留存学习笔记方便以后复习以及初学者学习。
# 测距、测速和测角的公式推导

## 测量距离
总的思路就是：分别测量发出信号和反弹之后的返回信号的相位，通过 θ=2Πft 可以知道发出信号与返回信号的时间差，若再知道信号波的速度c，也就可以知道反弹物体和雷达之间的距离d了，即
$$
2d=ct---------------------------（1）
$$
此处是2d是因为这个时间t包含了过去和反弹的两端距离。
<font color=yellow>但是</font>这种方法具有很大的限制，他要求发出信号和返回信号之间的相位变化一定要在一个周期内，而对于雷达信号的频率来说，一个周期内信号传播的距离也就几厘米。并不能实际使用。
#### 实际数据处理
另外实际上，在实际信号处理中对chirp信号做1DFFT，其中频域图的峰值点所对应的频率即为中频频率。

![Pasted image 20241102151723](https://github.com/user-attachments/assets/f38751fb-f344-4d02-b6de-8fdfbc6e630b)

获得中频频率之后，就可以计算出距离R

![Pasted image 20241102151903](https://github.com/user-attachments/assets/d1e500f8-1f83-4d9c-834b-08584b2de2a9)


还可以用==距离分辨率==来算，更方便。
### 测距范围
由上图可以确定出混频得到的中频频率，
$$
 f = B/T*t=S*t=2*d*S/c  ----------------（2）
$$
T为一个chirp信号的持续时间，B为带宽，t为接收信号的延迟时间。
### 距离分辨率
由公式（2）可以推得距离分辨率与频率分辨率的关系为
$$
🔺f =2*🔺d*S/c  ----------------（2）
$$
### 速度分辨率(不考虑物体速度)
速度分辨率可以从FFT来推到，即
$$
	🔺f = fs/N1
	
$$
$$
N = fs*T
$$

其中fs为采样频率，N1为fft点数，N为采样点数，T为时宽即chirp持续时间。
则距离分辨率为 
$$
🔺d = fs*c/(N1*2*S)= fs*c*T/(N1*2*B)=N*c/(2*N1*B)
$$
当N = N1时，即采样点数等于fft点数有
$$
🔺d = c/(2*B)
$$

==PS:补充一点，快时间维度的意思就是距离维度，慢时间维度就是速度（多普勒）维度==
## 测量速度
混频器输出的信号为中频信号，它的初始相位就是发送和接收信号的相位差

![Pasted image 20241015140318](https://github.com/user-attachments/assets/bc427723-067e-46bf-87d3-8334c1520acf)

![Pasted image 20241015140329](https://github.com/user-attachments/assets/2beea0ec-784c-41c3-aa5f-0b4ec6b6b1e5)

🔺┏时往返的时延。
所以两个连续的chirp信号的相位差可以用来估计物体速度v，设两chirp信号时间间隔为TC

![Pasted image 20241015140543](https://github.com/user-attachments/assets/7e19e0dd-0b60-4368-93ec-b27d40c6f2bc)

![Pasted image 20241015140612](https://github.com/user-attachments/assets/4254c150-5e7f-4744-b9cc-5beaef931236)

#### 实际数据处理
而在实际处理中，往往是对每帧的接收数据做2DFFT,也就是对一个chirp信号的采样点维度和chirp数维度，分别能得到距离信息和速度信息。

![Pasted image 20241102152529](https://github.com/user-attachments/assets/f663d452-511f-4e65-963e-bad6071619af)


然后得到目标点里可以得到多普勒频率，之后可以算出速度。

![Pasted image 20241102152622](https://github.com/user-attachments/assets/b92924df-e3f0-4866-951a-9aba0cfe403f)

更方便的是用==速度分辨率==来计算。

### 最大测速范围
最大可以测量的相对速度是

![Pasted image 20241015140730](https://github.com/user-attachments/assets/7adea987-8daf-41aa-a1d3-177c2bc5f847)


## 调皮连续波内容
和距离测量也是同样的原理，只不过雷达在移动，或者说被测物体在移动（运动是相对的）。传输的距离就变成 2d+2vt了，而由于移动状态，返回信号会产生一个多普勒频移，入射角为Θ的信号，会产生cosΘ缩放的多普勒频移。如图

![Pasted image 20240625095348](https://github.com/user-attachments/assets/667458ee-cd72-48f6-a66a-da53b26c7e00)

根据

![Pasted image 20240625095412](https://github.com/user-attachments/assets/46e6d038-7dc9-48e0-b063-fe1f47742777)

通过算出多普勒频移可以得到物体运动的速度。




> [!NOTE] 调皮连续波的知乎文章中的问题
> ![Pasted image 20240624203029](https://github.com/user-attachments/assets/0bd0e06d-bc22-4d8c-86d5-1dab29731ae7)
> 传播的距离为什么是![Pasted image 20240624203109](https://github.com/user-attachments/assets/35c92529-bed0-4ff1-9489-cd080832cd76)
呢
> ![Pasted image 20240624203356](https://github.com/user-attachments/assets/b424c8ae-dd7c-4d7c-bc07-5c77a65e570c)
> A:此处忽略了一半的位移，应该是雷达不动，物体动，

# 测角

角度的测量需要至少2个RX天线，利用物体到每个天线的不同距离导致2D-FFT峰值的相位变化来估计到达角，

![Pasted image 20241015141015](https://github.com/user-attachments/assets/82889f2c-94e1-447a-931a-e635a3d461fd)
相位的变化在数学上可以推导出下式：

![Pasted image 20241015141323](https://github.com/user-attachments/assets/f9262aaa-0df1-4957-a1b2-a3d8bc171741)
  而根据图12中的几何关系，可以得到：
  
![Pasted image 20241015141311](https://github.com/user-attachments/assets/f9219769-bd52-4aa4-b524-ecf4167bcf9d)
  那么，
  
![Pasted image 20241015141302](https://github.com/user-attachments/assets/4d03738d-4816-410b-969a-96d824814542)
  同时，角度的准确测量也离不开∣ Δ ω ∣ < π，即
  
![Pasted image 20241015141246](https://github.com/user-attachments/assets/c63fd90f-4378-4f0d-ab64-095e76e01672)

### 最大角度测量范围
即，两个间隔为d的天线可提供的最大视角为

![Pasted image 20241015141239](https://github.com/user-attachments/assets/902a888c-a559-4902-b5ae-5e25e162d359)
当两个天线之间的间隔d=λ/2，会导致±90°的最大角视场。

==用测距的原理来测角，通过相位的改变，来测远端天线的距离，通过与天线间距离d的比值来得到角度==

#### 实际处理
实际处理中，通过取最大值来得到目标点之后，对比目标在同一帧里两个天线的接收信号的相位，计算出相位差，可以测得角度。

==注意==，要控制相位差在[-pi,pi]之间。需要加一个判断相位差的语句。

### 获取啁啾信号
通过二次瞬时相位复信号来分解为啁啾信号

![Pasted image 20240625164930](https://github.com/user-attachments/assets/fcc22f0c-2eb3-49c5-bbf1-438164ef3807)

一般表达式为

![Pasted image 20240626151921](https://github.com/user-attachments/assets/4579339b-563c-4f36-89ff-8a6ed0caf403)

# FMCW测量单目标距离
## 几何法
几何法与之前的测距原理相同，都是找发出信号和返回信号的时间差来算频率差，然后得到距离。

## 直观法

![Pasted image 20240626152358](https://github.com/user-attachments/assets/53eaaf77-8e14-4506-b0d6-010f4647460d)
 
## DSP
没看出来和之前有啥区别。



> [!NOTE] 问题3
> 测量速度的时候，由于返回信号的幅度很小，那就很容易被淹没在噪声里，如何能够确定返回信号是哪个？
> A：一个方法![Pasted image 20240626152125](https://github.com/user-attachments/assets/7b3bf72b-85ff-44c1-8dbd-55b3996ee4cc)
> 0 ┏0 fb三角形与0 B Tc三角形相似，所以![Pasted image 20240626152250](https://github.com/user-attachments/assets/6ec7f899-f1d3-424b-92f2-f6c2dd3ca478)
> Tc 、B 、 已知，只用知道┏0就可以知道fb
> 


参考文献：
[1] https://blog.csdn.net/Dandan2530/article/details/125282773?ops_request_misc=%257B%2522request%255Fid%2522%253A%252247749982-D2DD-4529-862B-45276932986F%2522%252C%2522scm%2522%253A%252220140713.130102334..%2522%257D&request_id=47749982-D2DD-4529-862B-45276932986F&biz_id=0&utm_medium=distribute.pc_search_result.none-task-blog-2~all~top_positive~default-1-125282773-null-null.142^v100^pc_search_result_base6&utm_term=fmcw%E9%9B%B7%E8%BE%BE%E6%B5%8B%E8%B7%9D%E5%8E%9F%E7%90%86&spm=1018.2226.3001.4187

[2] https://blog.csdn.net/yaozekun/article/details/132889415?spm=1001.2014.3001.5502

[3]  [FMCW雷达距离多普勒(RDM)处理方法中距离分辨率和速度分辨率的推导_dopplerfft算法-CSDN博客](https://blog.csdn.net/qq_41248471/article/details/104276739)

[4] [参考网站](https://blog.csdn.net/CUGzhazhadong/article/details/119541284)

# CFAR
## CFAR简介
CFAR（Constant False Alarm Rate）算法是一种常用的目标检测和跟踪算法。它的主要作用是在背景噪声中检测出目标信号，同时保证误检概率不变。
CFAR算法的基本思想是，对于每个雷达测量的数据点，以该点为中心，建立一个检测窗口，在该窗口内计算信号功率的平均值和方差，并将该窗口划分为若干个子窗口。然后，根据期望的误检概率和背景噪声的统计特性，计算出每个子窗口的阈值，用于判断该窗口内是否存在目标信号。

---

CFAR算法是目标检测的核心算法，由流程图可以看出，在进行2维fft之后，进行CFAR算法检测
![Pasted image 20240629145021](https://github.com/user-attachments/assets/b10b564d-be60-4829-89af-5a32d1251045)


---

CFAR算法的目标是设置阈值足够高以将虚警限制在可容忍的范围内，但又足够低以允许目标检测。（虚警是噪声被误认为是目标）也就是说CFAR是一个动态的实时调整判断阈值的算法，目的是能够尽可能的提高目标判断的准确度。

---
## CFAR窗
![Pasted image 20240629153928](https://github.com/user-attachments/assets/9150aba0-f7b8-4b7d-998e-f8df83c05291)


## 目标检测
根据CFAR窗，可以计算出门限，如果信号功率超过门限，则认为目标存在，反之则不存在，这里用<font color=yellow>xy师兄的测试图来做说明</font>。

## CA-CFAR
示意图如下
![Pasted image 20240629153712](https://github.com/user-attachments/assets/c3de088c-c450-46ff-a2ca-9012a85a01a1)


  通过示意图可以看出来，CFAR的思路就是，先根据前置窗和后置窗算出噪声平均功率，然后与门限因子相乘，得到门限值，之后和被测单元对比，被测单元的功率值大于门限则检测目标存在，反之则不存在。
$$
其中Z_{CA}是噪声功率，公式如下：
$$
$$
Z_{CA}  = 0110frac{1}{N_{train} }\sum_{i=1}^{N_{train} }X_{i}  
$$
$$
其中T_{CA}是门限因子，计算公式如下

$$
$$
T_{CA}  = {N_{train} }(P_{fa}^{-\frac{1}{{N_{train} }} } -1)
$$
## GOCA
![Pasted image 20240629160434](https://github.com/user-attachments/assets/eb248d62-4932-48a1-846b-1be7421c8e01)


示意图表示，分别将前置窗和后置窗的平均功率计算出来，然后比较取大，再与门限因子相乘得到门限。

好处：可以避免杂波边缘的虚警。杂波下 CA-CFAR 反应较晚，因为它采用整个训练单元范围的平均值，而 GOCA-CFAR 则采用最多一半的训练单元。

> [!NOTE] 问题
> 为什么能够避免杂波边缘的虚警？

## 2D CFAR
### 原理


2D CFAR的检测对象通常是距离-多普勒图像（RDM），其仿真代码链接为：https://link.zhihu.com/?target=https%3A//github.com/tooth2/2D-CFAR
示意图如下
![Pasted image 20240629161412](https://github.com/user-attachments/assets/e4030e20-fb86-4765-8729-ecb3b8d1a73d)


#### 实际仿真log
算出RDM之后，观察RDM的3D图，可以发现，在doppler维，目标大概占6格，range维度大概占20格，由此可得保护单元的大小6 * 20。由于担心把相近的目标列进训练单元。所以训练单元就往外扩一圈，也就是8 * 22。
![Pasted image 20241105202725](https://github.com/user-attachments/assets/ef99de41-f131-4bef-b228-5c9979d7760e)

![Pasted image 20241105202529](https://github.com/user-attachments/assets/8b254d02-00b2-4ba1-8c56-8628eac16a29)

但是仍然不好
train_range = 1;
train_doppler = 1;
guard_range = 10;
guard_doppler = 3;
![Pasted image 20241106104306](https://github.com/user-attachments/assets/ce8c9c59-52de-4e92-b4f2-5ef5f277f8db)


train_range = 1;
train_doppler = 2;
guard_range = 6;
guard_doppler = 12;
![Pasted image 20241106104127](https://github.com/user-attachments/assets/9975cc8d-71cc-4055-b654-b3dc3c0956df)

扩大训练单元
train_range = 4;
train_doppler = 8;
guard_range = 10;
guard_doppler = 3;
![Pasted image 20241106104443](https://github.com/user-attachments/assets/d5c3130b-039c-4063-98e9-f2bac4e2ab3a)


> [!NOTE] 问题
> 如果用正常2d_cfar，在多目标监测时，会可能把相近的目标纳入训练单元。

# 参考文献
【1】Heijne, Rainier. Comparing Detection Algorithms for Short Range Radar based on the use-case of the Cobotic-65

【2】James J. Jen. A STUDY OF CFAR IMPLEMENTATION COST AND PERFORMANCE TRADEOFFS IN HETEROGENEOUS ENVIRONMENTS

【3】MATTHIAS KRONAUGE. Fast Two-Dimensional CFAR Procedure

【4】[https://la.mathworks.com/help/phased/ug/constant-false-alarm-rate-cfar-detection.html](https://link.zhihu.com/?target=https%3A//la.mathworks.com/help/phased/ug/constant-false-alarm-rate-cfar-detection.html)

【5】Mark Richards,_Fundamentals of Radar Signal Processing_, McGraw Hill, 2005

【6】[https://zhuanlan.zhihu.com/p/652220176](https://zhuanlan.zhihu.com/p/652220176)

【7】https://zhuanlan.zhihu.com/p/508870274
