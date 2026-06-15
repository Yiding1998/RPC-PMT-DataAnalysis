# RPC_PMT_DataAnalysis.C 代码修改位置与前后对照

## 1. 文档目的

本文档单独记录本轮代码检查后实际完成的修改，重点说明
**RPC_PMT_DataAnalysis.C** 中每一处修改的位置、修改前后的行为、修改原因、
影响范围和对应 Git 提交。

统计范围为提交 **4e9d497** 之后至 **8d9614f**，涉及：

- RPC_PMT_DataAnalysis.C
- RPC_PMT_DataAnalysis_Selected925.C
- RPC_WaveformSelector.C

代码改动统计为 182 行新增、66 行删除。后两个文件的修改主要用于保持筛选流程、
批处理流程与主分析程序的数据约定一致。

## 2. 提交总览

| 优先级 | 提交 | 修改主题 | 主要文件 |
|---|---|---|---|
| 配置 | d4c81c8 | 更新当前实验数据目录和输出标签 | RPC_PMT_DataAnalysis.C |
| 高 | c04c229 | 无效阈值时间改用 NaN，并显式记录有效性 | RPC_PMT_DataAnalysis.C |
| 高 | a6825d5 | 修正电荷积分越界和非等间隔积分误差 | RPC_PMT_DataAnalysis.C |
| 高 | ce43773 | 校验 CSV 列数、长度、有限值和时间单调性 | 三个分析文件 |
| 高 | 646e6a1 | 固定阈值时间独立于百分比阈值搜索 | RPC_PMT_DataAnalysis.C |
| 高 | 967c85a | 检查 PMT 拟合状态并提供插值回退 | RPC_PMT_DataAnalysis.C |
| 优化 | 61755c8 | 拟合样本与分辨率评估样本分离 | RPC_PMT_DataAnalysis.C |
| 优化 | c345598 | 简化 time-walk 模型并检查拟合质量 | RPC_PMT_DataAnalysis.C |
| 优化 | 8d11af9 | 按实际存在的 CSV 文件执行批处理 | 三个分析文件 |
| 优化 | 0f21876 | 在输出树中保存电压值 | RPC_PMT_DataAnalysis.C |
| 优化 | 9911dee | 修正循环中 ROOT 对象重名和内存泄漏 | RPC_PMT_DataAnalysis.C |

## 3. 逐项修改明细

### 3.1 更新实验目录和结果标签

**提交：** d4c81c8

**修改位置：** main() 开头的数据目录和 RPC_PMT_DataAnalysis(...) 调用参数。

**修改前：**

~~~cpp
std::string data_path =
    "/ustcfs/STCFUser/yzhao/LowPressure/ExpDatas_2606/40_50_10/3kpa/";
RPC_PMT_DataAnalysis(data_path, "590um_3kPa");
~~~

**修改后：**

~~~cpp
std::string data_path =
    "/ustcfs/STCFUser/yzhao/LowPressure/ExpDatas/1090um/5kpa/";
RPC_PMT_DataAnalysis(data_path, "1090um_5kPa");
~~~

**原因：** 将程序入口切换到当前 1090 um、5 kPa 实验数据，避免结果文件沿用旧实验标签。

**影响范围：** 只影响直接运行宏时读取的数据目录和输出文件命名，不改变分析算法。

### 3.2 正确表达未找到阈值交点的事件

**提交：** c04c229

**修改位置：**

- 头文件区新增 limits
- t_threshold、t_threshold_prev 初始化
- ROOT 输出树分支
- 时间直方图、Profile、time-walk 修正和分辨率统计的填充条件

**修改前的问题：**

未找到阈值交点时，时间变量使用 1e6 作为哨兵值。后续代码再通过
t_threshold 小于 1e2 判断结果是否有效。这会把“无测量结果”和一个数值混在一起，
并且依赖人为选择的范围上限。

**修改后：**

~~~cpp
double t_threshold = std::numeric_limits<double>::quiet_NaN();
double t_threshold_prev = std::numeric_limits<double>::quiet_NaN();
bool threshold_valid = false;

tree->Branch("threshold_valid", &threshold_valid, "threshold_valid/O");
~~~

后续统计统一使用：

~~~cpp
if (std::isfinite(t_threshold)) {
  // 填充时间统计量
}
~~~

**原因：** NaN 是浮点数据中表达“未定义/未测得”的标准方式；额外的布尔分支便于
ROOT 用户直接筛选有效事件。

**影响范围：** 无效事件不再以极大时间值污染直方图和 Profile；输出树新增
threshold_valid 分支；下游分析应使用该分支或 isfinite 选择事件。

### 3.3 修正电荷积分边界及积分公式

**提交：** a6825d5

**修改位置：** 单事件波形电荷积分循环。

**修改前：**

~~~cpp
for (int i = 0; i < sig_size; i++) {
  charge += sig[i] * (time[1] - time[0]);
}
~~~

该实现按采样点数量执行积分，但区间数量应为 sig_size - 1；同时所有区间强制使用
首个采样间隔，无法正确处理非严格等间隔时间轴。

**修改后：**

~~~cpp
for (int i = 0; i + 1 < sig_size; ++i) {
  const double dt = time[i + 1] - time[i];
  charge += 0.5 * (sig[i] + sig[i + 1]) * dt;
}
~~~

**原因：** 使用每个相邻采样点之间的实际 dt，并采用梯形积分，消除末端越界风险
以及非等间隔采样带来的系统误差。

**影响范围：** 所有基于积分电荷的直方图、二维关系、time-walk 拟合输入和输出树
电荷分支。数值可能与旧版本略有不同，这是积分定义被修正后的预期变化。

### 3.4 增加 CSV 数据完整性校验

**提交：** ce43773

**修改位置：** 三个分析文件的 CSV 读取之后。

**新增校验：**

~~~cpp
if (csv_data.size() < 3) {
  continue;
}

if (time.size() != rpc_signal.size() ||
    time.size() != pmt_signal.size()) {
  continue;
}
~~~

同时逐点检查：

- 时间和信号值必须是有限浮点数。
- 时间数组必须严格递增。
- 采样点数量必须满足后续基线、搜索和积分操作。

**原因：** 原代码默认输入文件始终完整。空列、截断列、NaN、Inf 或乱序时间均可能
造成越界、错误积分或无意义拟合。

**影响范围：** 异常 CSV 文件现在会被报告并跳过，不会终止整个批处理，也不会进入
ROOT 统计结果。

### 3.5 将固定阈值搜索与百分比阈值搜索解耦

**提交：** 646e6a1

**修改位置：** RPC 前沿时间提取逻辑。

**修改前的问题：** 固定阈值和多个百分比阈值在同一搜索过程中处理，并共享提前结束
条件。某个百分比阈值先满足后，循环可能提前退出，导致固定阈值交点没有被完整搜索。

**修改后：**

- 固定阈值从 -9 ns 起独立扫描。
- 单独保存固定阈值前一个采样点索引 idx_threshold_prev。
- 百分比阈值不再控制固定阈值搜索的生命周期。
- 固定阈值插值只使用它自身找到的相邻采样点。

**原因：** 固定阈值时间是后续 time-walk 修正的核心量，不能依赖另一组可选阈值的
搜索状态。

**影响范围：** t_threshold、固定阈值时间直方图、二维电荷-时间关系以及修正后的
时间分辨率。

### 3.6 检查 PMT 高斯拟合并增加回退算法

**提交：** 967c85a

**修改位置：** PMT 脉冲定时高斯拟合。

**修改前的问题：**

- 第二次拟合窗口使用方差值而不是标准差，缺少平方根转换。
- 未检查 TGraph::Fit 返回状态。
- 即使拟合失败或 sigma 非法，仍直接使用拟合均值作为 PMT 时间。

**修改后：**

~~~cpp
const int fit_status = g_pmt->Fit(f_pmt, "QNR");
const bool fit_valid =
    fit_status == 0 &&
    std::isfinite(f_pmt->GetParameter(1)) &&
    std::isfinite(f_pmt->GetParameter(2)) &&
    f_pmt->GetParameter(2) > 0.;
~~~

第二次拟合窗口使用标准差；当拟合无效时，程序使用峰值附近采样点插值得到 PMT
参考时间。

**原因：** PMT 时间是 RPC-PMT 时间差的参考端。一次失败拟合会直接制造异常时间差，
并进一步扭曲 time-walk 模型。

**影响范围：** PMT 到达时间、原始时间差、拟合样本和最终时间分辨率。正常拟合事件
保持原算法，只有失败或非法拟合进入回退路径。

### 3.7 分离 time-walk 拟合样本和评估样本

**提交：** 61755c8

**修改位置：** 文件遍历、二维拟合样本填充和修正后分辨率统计。

**修改后规则：**

- 偶数序号事件用于构建电荷-时间二维分布和拟合修正函数。
- 奇数序号事件只用于填充修正后的时间直方图并评估分辨率。

**原因：** 原实现用同一批事件拟合修正函数并评价修正效果，会低估泛化误差。简单的
训练/验证拆分使最终分辨率更接近独立数据上的表现。

**影响范围：** 拟合样本数量约减半；分辨率统计只包含验证样本；输出树仍可保存全部
通过质量检查的事件。

### 3.8 简化 time-walk 修正模型并检查拟合质量

**提交：** c345598

**修改位置：** 二维 Profile 的函数定义、拟合条件和结果使用。

**修改前：** 使用含高阶多项式及多个反比项的 10 参数函数。对于有限事件数，该模型
参数相关性强，容易过拟合或产生不稳定外推。

**修改后：**

~~~cpp
TF1 *fit_walk = new TF1(
    "fit_walk", "[0] + [1] / sqrt(x) + [2] / x",
    charge_min, charge_max);
~~~

拟合前要求有效 Profile bin 数不少于 8，并检查拟合状态。只有拟合成功时才应用
time-walk 修正。

**原因：** 三参数模型能够描述常见前沿触发的电荷依赖，同时显著降低小样本过拟合和
病态参数的风险。

**影响范围：** 修正曲线形状、修正后时间分布和报告的时间分辨率。低统计量数据不会
再强制套用不可靠曲线。

### 3.9 按实际存在的 CSV 文件批处理

**提交：** 8d11af9

**修改位置：** 三个分析文件的输入文件发现、编号提取和排序逻辑。

**修改前的问题：** 代码按 0.csv 到某个计数上限拼接路径。只要编号中间缺失，就可能
提前 break，导致后面实际存在的文件完全不处理。

**修改后：**

- 遍历目录中实际存在且文件名主干为数字的 .csv 文件。
- 使用 stoll 将主干转换为数字并按数值排序，避免字典序把 10.csv 排在 2.csv 前面。
- 单个文件读取失败时使用 continue，不终止后续文件。
- 波形筛选程序使用实际文件名主干生成关联结果。

**原因：** 实验数据目录可能存在删选后编号空洞、拷贝不完整或非 CSV 辅助文件。
批处理范围应由真实文件清单决定。

**影响范围：** 文件处理顺序稳定；编号不连续时后续数据不再被遗漏。

### 3.10 在输出 ROOT 树中记录实验电压

**提交：** 0f21876

**修改位置：** 输出变量定义、ROOT TTree 分支创建和每事件 Fill() 前赋值。

**修改后：**

~~~cpp
double tac_voltage = 0.;
tree->Branch("voltage", &tac_voltage, "voltage/D");
~~~

每次写入事件前，将当前目录或运行配置解析出的电压赋给 tac_voltage。

**原因：** 原输出树缺少运行电压，下游合并不同电压的数据时需要依赖文件名，容易
丢失实验条件。

**影响范围：** ROOT 树新增 voltage 分支；旧分析字段不变。

### 3.11 修正 ROOT 对象重名和循环内内存泄漏

**提交：** 9911dee

**修改位置：** 每事件创建 PMT TGraph 和 TF1 的代码。

**修改前的问题：**

- 循环内反复创建同名 TF1，ROOT 全局对象目录可能发生名称冲突或替换警告。
- 动态创建的 TGraph、TF1 未在事件结束后释放，长批处理内存持续增长。

**修改后：**

- 使用原子递增编号生成唯一拟合函数名。
- 直接以函数指针调用拟合，减少依赖 ROOT 名称查找。
- 事件处理结束后显式释放 f_pmt 和 g_pmt。

**原因：** 长时间批处理必须保证每事件临时对象及时释放，并避免 ROOT 名称注册表
产生非确定行为。

**影响范围：** 不改变物理算法；降低大数据量运行时的内存占用和 ROOT 对象冲突风险。

## 4. 按文件归纳

### RPC_PMT_DataAnalysis.C

承担全部核心修改：输入校验、文件枚举、电荷积分、阈值时间提取、PMT 拟合回退、
NaN 无效值、训练/评估拆分、time-walk 模型、电压分支和 ROOT 对象生命周期。

### RPC_PMT_DataAnalysis_Selected925.C

同步修改输入文件枚举、编号排序、CSV 数据有效性检查和单文件失败后的继续处理，
保证选定数据分析版本不会再次遗漏编号不连续的文件。

### RPC_WaveformSelector.C

同步修改输入文件枚举、CSV 校验和实际文件名使用方式，保证筛选结果与原始 CSV
编号一一对应。

## 5. 验证方式

~~~bash
git diff --check
g++ -fsyntax-only RPC_PMT_DataAnalysis.C $(root-config --cflags)
g++ -fsyntax-only RPC_PMT_DataAnalysis_Selected925.C $(root-config --cflags)
g++ -fsyntax-only RPC_WaveformSelector.C $(root-config --cflags)
~~~

逐提交检查命令：

~~~bash
git log --oneline --decorate
git show COMMIT_ID
git diff 4e9d497..8d9614f -- RPC_PMT_DataAnalysis.C
~~~

语法检查确认 C++/ROOT 接口层面的可编译性；最终物理结果仍应结合已知标定数据检查
电荷分布、原始时间分布、修正曲线和验证集分辨率是否符合实验预期。

## 6. 检查中记录但未擅自修改的问题

以下问题不属于已经提交的修改，按要求仅记录：

1. 百分比阈值交点的索引搜索方向需要结合实际脉冲极性和时间顺序再次核对。
2. 部分事件处理中存在提前 continue 路径，临时 ROOT 对象是否全部释放仍可进一步审计。
3. 编译器可能报告个别已计算但未使用的局部变量，可在确认无诊断用途后清理。
4. 固定阈值搜索从 -9 ns 开始是实验经验参数，噪声较大时可能提前命中，应以波形样本验证。
5. 偶数/奇数拆分是确定性的；若采集顺序存在时间漂移，训练集和验证集仍可能有系统差异。

这些项目尚未修改，避免在缺少实验条件确认时改变物理分析定义。

## 7. 查看历史版本

~~~bash
git show c04c229
git show a6825d5
git show ce43773
git show 646e6a1
git show 967c85a
git show 61755c8
git show c345598
git show 8d11af9
git show 0f21876
git show 9911dee

git show 4e9d497:RPC_PMT_DataAnalysis.C
git show HEAD:RPC_PMT_DataAnalysis.C
~~~
