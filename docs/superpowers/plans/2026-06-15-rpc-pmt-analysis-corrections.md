# RPC-PMT Analysis Corrections Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复十项已确认的正确性、统计方法和性能问题，并保持 ROOT 宏入口兼容。

**Architecture:** 优先在两个分析宏中同步最小修复，筛选器只修改共享的输入和文件遍历行为。每项通过针对性静态或 ROOT 检查验证后单独提交。

**Tech Stack:** C++17、CERN ROOT 6.26、Git。

---

### Task 1: 无效过阈事件
- [ ] 增加 `threshold_valid` 分支和显式有限值判定。
- [ ] 验证低幅度事件不进入阈值统计。
- [ ] 提交 `fix: 过滤无效过阈事件`。

### Task 2: 电荷积分边界
- [ ] 将右侧循环限制为 `i + 1 < size` 并使用梯形积分。
- [ ] 验证末端脉冲不越界。
- [ ] 提交 `fix: 修复电荷积分边界`。

### Task 3: CSV 输入校验
- [ ] 校验至少三个点、长度一致、时间递增。
- [ ] 验证空和短波形被拒绝。
- [ ] 提交 `fix: 校验CSV波形输入`。

### Task 4: 固定阈值搜索
- [ ] 在主宏中使用独立的 10 mV 上升沿搜索。
- [ ] 验证高幅度波形仍能得到交点。
- [ ] 提交 `fix: 独立计算固定阈值时间`。

### Task 5: PMT 高斯拟合
- [ ] 修正拟合边界公式并检查状态。
- [ ] 失败时保留插值结果。
- [ ] 提交 `fix: 修正PMT高斯定时拟合`。

### Task 6: 训练验证拆分
- [ ] 偶数事件拟合、奇数事件评估。
- [ ] 验证两集合不重叠。
- [ ] 提交 `fix: 分离时间修正训练与评估数据`。

### Task 7: 修正模型质量
- [ ] 使用三参数低阶模型。
- [ ] 检查 Profile 有效 bin 和拟合状态。
- [ ] 提交 `refactor: 简化时间游走拟合模型`。

### Task 8: 实际文件枚举
- [ ] 枚举数字命名 CSV 并排序。
- [ ] 缺号或单文件失败时继续。
- [ ] 提交 `fix: 按实际CSV文件列表处理数据`。

### Task 9: 汇总树电压分支
- [ ] 增加并填充 `voltage`。
- [ ] ROOT 检查分支存在。
- [ ] 提交 `feat: 为时间修正汇总保存电压`。

### Task 10: ROOT 对象生命周期
- [ ] 逐事件拟合对象使用唯一名称并释放。
- [ ] 运行语法检查和代表性加载检查。
- [ ] 提交 `refactor: 管理逐事件ROOT拟合对象`。
