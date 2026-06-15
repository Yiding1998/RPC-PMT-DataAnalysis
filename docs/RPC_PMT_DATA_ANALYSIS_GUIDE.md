# `RPC_PMT_DataAnalysis.C` 实验数据分析说明

## 1. 文档目的

本文面向使用示波器波形研究 RPC（Resistive Plate Chamber）幅度、电荷和时间性能的实验数据分析人员，说明 [`RPC_PMT_DataAnalysis.C`](../RPC_PMT_DataAnalysis.C) 的输入约定、物理量定义、计算公式、时间游走修正、ROOT 输出和结果解释。

本文描述当前代码的实际行为。每次分析应记录数据条件、气体、气压、间隙、放大倍数、示波器设置、触发方式和 Git commit。

## 2. 程序完成的工作

程序按高压点自动执行：

1. 识别数字命名的高压子目录。
2. 读取三列 CSV：时间、RPC、PMT。
3. 完成单位换算、极性翻转和基线扣除。
4. 提取幅度、上升时间、FWHM、电荷和定时时间。
5. 计算 `PMT - RPC` 时间差。
6. 对每个高压点的分布做高斯拟合。
7. 生成参数二维相关图。
8. 建立幅度/电荷时间游走曲线并做 TA/TQ 修正。
9. 用独立验证事件比较修正前后的时间分辨率。
10. 保存事件树、拟合树、修正汇总树和 Canvas。

主调用关系：

```text
RPC_PMT_DataAnalysis()
└── analyzeDataSerial()
    ├── GetVoltages()
    ├── processVoltageSerial()
    │   ├── readCSVData()
    │   ├── baselineCorrection()
    │   ├── calculateSignalParams()
    │   └── calculateTimingDiff()
    ├── gaussianFit()
    ├── create2DPlotsOptimized()
    └── performTA_TQ_Correction()
        └── iterativeGaussianFit()
```

## 3. 运行环境与命令

依赖 CERN ROOT 和支持 C++17 的编译器：

```bash
root-config --version
root-config --cflags --libs
```

交互运行：

```bash
root -l RPC_PMT_DataAnalysis.C
```

批处理运行：

```bash
root -l -b -q 'RPC_PMT_DataAnalysis.C()'
```

输出 ROOT 文件位于运行命令所在目录。

## 4. 运行前配置

主宏顶部参数：

```cpp
const double IMPEDANCE = 50.0;
const double RPC_AMPLIFICATION = 10.0;
const double PMT_AMPLIFICATION = 1.0;
const double BASELINE_WINDOW = 10.0;
const double PMT_CFD_FRACTION = 0.5;
const double RPC_CFD_FRACTION = 0.2;
const double RPC_THRESHOLD = 10;

const bool FIT_RPC = false;
const bool FIT_PMT = true;

const string csvFolder =
    "/ustcfs/STCFUser/yzhao/LowPressure/ExpDatas/1090um/5kpa";
const char* outputFileName = "1090um_5kPa";
```

| 参数 | 当前值 | 含义 |
|---|---:|---|
| `IMPEDANCE` | 50 Ω | 示波器/读出终端阻抗 |
| `RPC_AMPLIFICATION` | 10 | RPC 电压增益，电荷换算时扣除 |
| `PMT_AMPLIFICATION` | 1 | PMT 电压增益 |
| `BASELINE_WINDOW` | 10 ns | 从波形起点开始的基线窗口 |
| `PMT_CFD_FRACTION` | 0.5 | PMT 参考时间的峰值比例 |
| `RPC_CFD_FRACTION` | 0.2 | RPC CFD 时间的峰值比例 |
| `RPC_THRESHOLD` | 10 mV | RPC 固定阈值 |
| `FIT_RPC` | `false` | RPC 是否做单脉冲高斯拟合 |
| `FIT_PMT` | `true` | PMT 是否做单脉冲高斯拟合 |
| `csvFolder` | 路径 | 数据根目录 |
| `outputFileName` | 字符串 | 输出文件名，不含 `.root` |

`RPC_AMPLIFICATION` 只影响电荷，不影响示波器端幅度和定时时间。若信号未经过放大，应设为 1。

## 5. 输入目录和 CSV

### 5.1 目录结构

```text
数据根目录/
├── 1000/
│   ├── 000001.csv
│   ├── 000002.csv
│   └── ...
├── 1050/
│   └── ...
└── ...
```

只有纯数字目录被解释为高压（V）。程序枚举纯数字文件名的 `.csv` 并按数值排序，允许中间缺号。

### 5.2 三列定义

CSV 无标题行：

```text
time_s, rpc_V, pmt_V
```

读取后：

\[
t_{\mathrm{ns}}=t_{\mathrm{s}}\times10^9,
\]

\[
V_{\mathrm{mV}}=-V_{\mathrm{V}}\times10^3.
\]

负号用于把示波器中的负脉冲翻转成正脉冲。

### 5.3 有效性要求

一个事件必须：

- 至少有 3 个有效采样点；
- 时间、RPC、PMT 数组长度一致；
- 时间有限且严格递增；
- 每行三列可以解析。

无效文件被跳过，不会终止同一高压的其他文件。

## 6. 波形预处理

### 6.1 基线

从波形第一个时间点开始取 10 ns 窗口。若窗口内有 \(N_b\) 点：

\[
B=\frac{1}{N_b}\sum_{i=0}^{N_b-1}V_i,
\qquad
V_i'=V_i-B.
\]

RPC 和 PMT 分别计算基线。必须确认前 10 ns 没有真实脉冲，否则幅度和电荷会被低估。

### 6.2 峰值搜索

只在：

\[
[\max(-14\ \mathrm{ns},t_{\min}),\min(14\ \mathrm{ns},t_{\max})]
\]

内寻找最大值和最小值，取绝对值较大者：

\[
A=\max(|V_{\max}|,|V_{\min}|).
\]

## 7. 单事件物理量

### 7.1 幅度

\[
A=|V(t_{\mathrm{peak}})|.
\]

单位 mV。事件树分支为 `amplitude_RPC`、`amplitude_PMT`。幅度没有除以前端增益。

### 7.2 阈值交点的线性插值

相邻点 \((t_1,V_1)\)、\((t_2,V_2)\) 包围阈值 \(V_{th}\) 时：

\[
t_{cross}=t_1+
\frac{V_{th}-|V_1|}{|V_2|-|V_1|}(t_2-t_1).
\]

若两点幅度近似相同，取两点时间中值。

### 7.3 10%–90% 上升时间

\[
V_{10}=0.1A,\qquad V_{90}=0.9A,
\]

\[
t_{rise}=|t_{90}-t_{10}|.
\]

单位 ns，对应 `rise_time_RPC`、`rise_time_PMT`。

### 7.4 FWHM

\[
V_{50}=0.5A,
\qquad
\mathrm{FWHM}=|t_{50,R}-t_{50,L}|.
\]

单位 ns，对应 `fwhm_RPC`、`fwhm_PMT`。

### 7.5 电荷

积分窗口由峰值 1% 的左右边界确定：

\[
V_{1\%}=0.01A.
\]

相邻采样点用梯形积分：

\[
\Delta Q_i=
\frac{V_i+V_{i+1}}{2R}(t_{i+1}-t_i),
\]

其中电压换成 V、时间换成 s、\(R=50\ \Omega\)。总电荷换算成等效电子数：

\[
N_e=\frac{|\sum_i\Delta Q_i|}{eG},
\]

\[
e=1.602176634\times10^{-19}\ \mathrm{C}.
\]

\(G\) 为对应通道放大倍数。分支 `charge_RPC`、`charge_PMT` 的单位是等效电子数，不是库仑。

### 7.6 RPC 20% CFD

\[
V_{CFD,RPC}=0.2A_{RPC}.
\]

在峰值左侧寻找上升沿交点，得到 \(t_{RPC,CFD}\)。这是离线恒比交点，不是硬件 CFD 的延迟-反相-过零算法。

### 7.7 RPC 10 mV 固定阈值

程序从 -9 ns 附近向后扫描，寻找第一个达到 10 mV 的点并插值：

\[
t_{RPC,th}=t(V=10\ \mathrm{mV}).
\]

若 RPC 幅度不高于 10 mV，`timing_threshold` 为 NaN，`threshold_valid=false`，该事件不参与阈值统计和修正。

### 7.8 PMT 参考时间

时间差使用 PMT 50% 上升沿：

\[
V_{PMT,50}=0.5A_{PMT}.
\]

当 `FIT_PMT=false` 时使用采样点插值。当 `FIT_PMT=true` 时拟合：

\[
V(t)=A\exp\left[-\frac{(t-\mu)^2}{2\sigma^2}\right].
\]

高斯上升沿达到比例 \(f\) 的时间：

\[
t_f=\mu-\sigma\sqrt{-2\ln f}.
\]

因此：

\[
t_{PMT,50}=\mu-\sigma\sqrt{-2\ln0.5}.
\]

程序进行初始拟合和两次缩小范围的重拟合。最终拟合失败、参数非有限或 \(\sigma\le0\) 时回退到采样点插值。

## 8. 时间差定义

代码不取绝对值：

\[
\Delta t_{CFD}=t_{PMT,50}-t_{RPC,CFD},
\]

\[
\Delta t_{th}=t_{PMT,50}-t_{RPC,th}.
\]

对应 `timing_CFD`、`timing_threshold`。均值包含电缆、电子学、触发和通道偏置；均值不为零不代表错误，时间性能主要看分布宽度。

## 9. 每个高压的事件树

每个高压创建 `tree_<HV>V`，例如 `tree_1000V`。

| 分支 | 单位/含义 |
|---|---|
| `amplitude_RPC` | RPC 幅度，mV |
| `amplitude_PMT` | PMT 幅度，mV |
| `rise_time_RPC` | RPC 10%–90% 上升时间，ns |
| `rise_time_PMT` | PMT 上升时间，ns |
| `fwhm_RPC` | RPC FWHM，ns |
| `fwhm_PMT` | PMT FWHM，ns |
| `charge_RPC` | RPC 等效电子数 |
| `charge_PMT` | PMT 等效电子数 |
| `timing_CFD` | PMT 50% - RPC 20%，ns |
| `timing_threshold` | PMT 50% - RPC 10 mV，ns/NaN |
| `threshold_valid` | 固定阈值时间是否有效 |

检查示例：

```cpp
TFile f("1090um_5kPa.root");
auto t = (TTree*)f.Get("tree_1000V");
t->Print();
t->Scan("amplitude_RPC:charge_RPC:timing_CFD:timing_threshold:threshold_valid",
        "", "", 10);
```

## 10. 一维分布和统计拟合

| ROOT 目录 | 对象名 | 默认范围 |
|---|---|---|
| `amplitude` | `h_amp_<HV>V` | 0–1500 mV，1.5 mV/bin |
| `rise_time` | `h_rt_<HV>V` | 0–3 ns，0.003 ns/bin |
| `fwhm` | `h_fwhm_<HV>V` | 0.5–4 ns，约 0.008 ns/bin |
| `charge` | `h_charge_<HV>V` | 0–\(3\times10^7\) e，\(2\times10^4\) e/bin |
| `timing_const` | `h_timing_CFD_<HV>V` | -5–10 ns，约 0.008 ns/bin |
| `timing_thresh` | `h_timing_thresh_<HV>V` | -5–10 ns，约 0.008 ns/bin |

阈值直方图只填有效事件。分布拟合模型：

\[
F(x)=C\exp\left[-\frac{(x-\mu)^2}{2\sigma^2}\right].
\]

初始均值取最高 bin 中心，初始宽度通常取 RMS；初次拟合后在更新的 \(\mu\pm2\sigma\) 内重复拟合两次。

### 10.1 `tree_fit_parameters`

每个高压一行：

| 分支 | 含义 |
|---|---|
| `voltage` | 高压，V |
| `mean_amp`, `sigma_amp` | 幅度均值和高斯宽度 |
| `mean_rt`, `sigma_rt` | 上升时间均值和宽度 |
| `mean_fwhm`, `sigma_fwhm` | FWHM 均值和宽度 |
| `mean_charge`, `sigma_charge` | 电荷均值和宽度 |
| `mean_timing_CFD`, `sigma_timing_CFD` | CFD 时间差均值和宽度，ns |
| `mean_timing_thresh`, `sigma_timing_thresh` | 阈值时间差均值和宽度，ns |

`sigma_*` 是分布宽度，不是均值误差。

## 11. 普通二维相关图

程序创建：

```text
Amplitude_vs_RiseTime
Amplitude_vs_FWHM
Amplitude_vs_Charge
Amplitude_vs_Timing
RiseTime_vs_FWHM
RiseTime_vs_Charge
RiseTime_vs_Timing
FWHM_vs_Charge
FWHM_vs_Timing
Charge_vs_Timing
```

其中 `Timing` 指 `timing_CFD`。这些图用于检查幅度-电荷关系、波形形状变化、多事件簇、饱和，以及 CFD 是否仍依赖幅度或电荷。

## 12. TA/TQ 时间游走修正

### 12.1 时间游走模型

定时可能满足：

\[
\Delta t=F(A)+\epsilon
\quad\text{或}\quad
\Delta t=F(Q)+\epsilon.
\]

### 12.2 事件筛选

当前实际启用的筛选只有 RPC 上升时间：

\[
\mu_{RT}-5\sigma_{RT}<RT<\mu_{RT}+5\sigma_{RT}.
\]

代码虽计算幅度和 FWHM 的 ±5σ 条件，但未启用。不要把当前结果解释为三重筛选。

### 12.3 训练集和验证集

筛选后的事件按向量索引划分：

- 偶数索引：训练时间游走曲线；
- 奇数索引：评估修正前后时间分辨率。

这避免在拟合修正函数的同一批事件上报告性能。索引来自排序并筛选后的事件序列，不等同于原始 CSV 编号。

### 12.4 Profile 和拟合

训练集建立：

1. CFD vs 幅度；
2. CFD vs 电荷；
3. 阈值时间 vs 幅度；
4. 阈值时间 vs 电荷。

二维直方图做 `ProfileX`，拟合每个 X bin 的平均时间。每个 Profile 至少需要 8 个有事件的 bin；任意一个不足时，整个高压点跳过修正。

当前函数：

\[
F(x)=p_0+\frac{p_1}{\sqrt{x}}+\frac{p_2}{x},
\]

其中 \(x=A\) 或 \(Q\)。不要外推到拟合数据范围之外。

### 12.5 修正定义

\[
\Delta t_{TA}=\Delta t-F(A),
\]

\[
\Delta t_{TQ}=\Delta t-F(Q).
\]

减去完整函数也会减去常数项，因此修正后分布通常在 0 附近，性能评价关注 \(\sigma\)。

## 13. 修正后时间分辨率

验证集生成六种分布：

- CFD 修正前、TA 后、TQ 后；
- 固定阈值修正前、TA 后、TQ 后。

每个分布用迭代高斯拟合得到 \(\sigma\) 和 ROOT 参数误差。

### 13.1 `tree_fitTAC_parameters`

| 分支 | 单位 | 含义 |
|---|---|---|
| `voltage` | V | 高压 |
| `sigma_cfd_before` | ns | CFD 修正前宽度 |
| `sigmaErr_cfd_before` | ns | 拟合误差 |
| `sigma_cfd_afterTA` | ns | CFD TA 后宽度 |
| `sigmaErr_cfd_afterTA` | ns | 拟合误差 |
| `sigma_cfd_afterTQ` | ns | CFD TQ 后宽度 |
| `sigmaErr_cfd_afterTQ` | ns | 拟合误差 |
| `sigma_thresh_before` | ns | 阈值修正前宽度 |
| `sigmaErr_thresh_before` | ns | 拟合误差 |
| `sigma_thresh_afterTA` | ns | 阈值 TA 后宽度 |
| `sigmaErr_thresh_afterTA` | ns | 拟合误差 |
| `sigma_thresh_afterTQ` | ns | 阈值 TQ 后宽度 |
| `sigmaErr_thresh_afterTQ` | ns | 拟合误差 |

终端汇总使用 ns；高压依赖 Canvas 乘以 1000，纵轴为 ps。

## 14. 时间分辨率的物理含义

程序得到 RPC-PMT 符合时间差宽度：

\[
\sigma_{meas}^2\approx
\sigma_{RPC}^2+
\sigma_{PMT}^2+
\sigma_{electronics}^2+
\sigma_{clock}^2.
\]

因此输出不能自动等同于 RPC 本征分辨率。若参考贡献独立测量且可视为独立高斯涨落：

\[
\sigma_{RPC}=
\sqrt{\sigma_{meas}^2-
\sigma_{PMT}^2-
\sigma_{electronics}^2-
\sigma_{clock}^2}.
\]

只有根号内为正、测量条件一致时，方差扣除才有意义。

## 15. 校准图和汇总图

每个成功修正的高压创建 `Cali_<HV>V/`：

- `cF1_<HV>V`：幅度、电荷、CFD 时间、阈值时间的一维分布；
- `cF2_<HV>V`：四种时间游走二维图和红色拟合曲线；
- `cF3_<HV>V`：3×2 修正前、TA 后、TQ 后对比。

根目录汇总：

- `c_res_cfd`：CFD 三条分辨率曲线；
- `c_res_th`：固定阈值三条曲线；
- `c_res_combined`：六条曲线合并。

如果某高压因 Profile bin 不足或拟合失败而跳过，汇总可能以 0 表示缺失项。必须同时检查终端日志和 `Cali_<HV>V` 内容，不能把 0 当作真实分辨率。

## 16. ROOT 输出结构

```text
1090um_5kPa.root
├── tree_<HV>V
├── tree_fit_parameters
├── tree_fitTAC_parameters
├── amplitude/
├── rise_time/
├── fwhm/
├── charge/
├── timing_const/
├── timing_thresh/
├── Amplitude_vs_RiseTime/
├── Amplitude_vs_FWHM/
├── Amplitude_vs_Charge/
├── Amplitude_vs_Timing/
├── RiseTime_vs_FWHM/
├── RiseTime_vs_Charge/
├── RiseTime_vs_Timing/
├── FWHM_vs_Charge/
├── FWHM_vs_Timing/
├── Charge_vs_Timing/
├── Cali_<HV>V/
├── c_res_cfd
├── c_res_th
└── c_res_combined
```

## 17. 推荐分析流程

1. 随机检查原始波形，确认极性、基线、触发位置、削顶、振铃和多脉冲。
2. 核对数据目录、输出名、阻抗、增益、CFD 比例和固定阈值。
3. 先用少量事件运行，确认事件数、幅度、电荷和时间数量级。
4. 检查一维直方图是否有大量 underflow/overflow。
5. 查看 `cF2_<HV>V`，确认修正曲线跟随 Profile 主趋势。
6. 查看 `cF3_<HV>V`，同时检查核心宽度、长尾、双峰和拟合范围。
7. 从树中程序化提取结果，不建议人工抄图。

```cpp
TFile f("1090um_5kPa.root");
auto fit = (TTree*)f.Get("tree_fit_parameters");
fit->Scan("voltage:mean_amp:mean_charge:sigma_timing_CFD");

auto tac = (TTree*)f.Get("tree_fitTAC_parameters");
tac->Scan("voltage:sigma_cfd_before:sigma_cfd_afterTA:sigma_cfd_afterTQ");
```

## 18. 质量检查清单

### 幅度和电荷

- 是否出现饱和平台；
- 随高压的趋势是否合理；
- 幅度-电荷图是否多簇；
- 电荷量级是否与增益一致。

### 波形形状

- 上升时间、FWHM 是否有大量 0；
- 是否超出默认直方图范围；
- 上升时间 ±5σ 筛选是否过严或过松。

### 定时

- `threshold_valid` 比例是否合理；
- PMT 参考时间是否稳定；
- CFD 是否仍依赖幅度/电荷；
- 修正线是否边界发散；
- 修正后变窄是否伴随更严重长尾。

### 统计

- 事件数是否足够；
- Profile 是否至少有 8 个有效 bin；
- 高斯拟合范围、参数误差和 χ² 是否合理；
- 不同高压的有效事件比例是否接近；
- 训练集和验证集统计量是否充足。

## 19. 常见问题

### 找不到电压目录

`1000V` 不会被识别，应命名为 `1000`。

### 处理事件数少于 CSV 数

可能是文件名非纯数字、扩展名不是小写 `.csv`、采样点少于 3、时间不递增或行解析失败。

### 阈值直方图条目少

RPC 幅度必须高于 10 mV；低幅度事件 `threshold_valid=false`。

### 电荷数量级异常

依次检查 CSV 电压单位、50 Ω 终端、增益、基线窗口和脉冲是否完整返回基线。

### 修正结果为 0 或缺少校准图

检查终端是否出现：

```text
有效Profile bin不足，跳过修正
时间游走拟合失败，跳过修正
```

### PMT 拟合很慢

`FIT_PMT=true` 每事件执行多次高斯拟合。可用小样本比较 `true/false` 后决定是否关闭；关闭时改用线性插值。

## 20. 当前限制

1. 峰值搜索固定为 -14 至 14 ns。
2. 固定阈值从 -9 ns 附近搜索，前置噪声可能提前过阈。
3. 基线只扣常数，不处理斜率和周期噪声。
4. 1% 峰值积分边界可能受噪声影响。
5. 普通直方图范围固定，范围外事件进入 underflow/overflow。
6. 实际分布可能非高斯或有长尾。
7. TA/TQ 当前只启用上升时间 ±5σ 筛选。
8. 四个 Profile 必须同时满足有效 bin 要求。
9. 奇偶划分确定但不是随机抽样；采集慢漂移可能使两组不完全等价。
10. 输出是 RPC-PMT 符合宽度，未扣除参考贡献。
11. `sigmaErr_*` 只表示拟合参数误差，不含系统误差。

## 21. 实验记录模板

```text
数据目录：
输出文件：
RPC 结构/气隙：
气体与比例：
气压、温度：
高压点：
示波器型号、采样率、时间窗口：
RPC/PMT 通道：
终端阻抗：
RPC/PMT 放大倍数：
基线窗口：
RPC CFD 比例：
RPC 固定阈值：
PMT CFD 比例：
是否拟合 PMT：
代码 Git commit：
ROOT 版本：
备注：
```

```bash
git rev-parse HEAD
root-config --version
```

记录代码版本和运行条件是保证结果可重复、可比较的必要步骤。
