# RPC-PMT Data Analysis Guide Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 编写面向实验数据分析人员的 `RPC_PMT_DataAnalysis.C` 详细中文手册。

**Architecture:** 以实际数据处理顺序组织文档；所有参数、公式、分支和 ROOT 对象名称从当前代码核对。README 只提供入口，不重复手册内容。

**Tech Stack:** Markdown、C++17 ROOT 宏、CERN ROOT。

---

### Task 1: 核对实现
- [ ] 读取主宏的参数、信号算法、修正流程和输出定义。
- [ ] 记录公式中的单位、符号和有效性条件。

### Task 2: 编写实验分析手册
- [ ] 新增 `docs/RPC_PMT_DATA_ANALYSIS_GUIDE.md`。
- [ ] 覆盖输入、公式、运行、输出解释和限制。

### Task 3: 增加入口并验证
- [ ] 在 `README.md` 增加链接。
- [ ] 检查对象名、参数值、Markdown 格式和 Git 差异。
- [ ] 提交并推送 GitHub。
