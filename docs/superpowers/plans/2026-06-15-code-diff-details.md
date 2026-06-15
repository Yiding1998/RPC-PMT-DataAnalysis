# RPC-PMT Code Diff Details Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 编写逐提交、逐函数的代码修改位置和行为对照文档。

**Architecture:** 以 Git 历史为事实来源，从十个独立修复提交提取差异；文档使用稳定的函数名和 commit，而不是容易漂移的绝对行号。README 提供入口。

**Tech Stack:** Git、Markdown、C++ ROOT 宏。

---

### Task 1: 提取差异事实
- [ ] 获取 `4e9d497..8d9614f` 的提交、文件统计和每个修复提交的 diff。
- [ ] 将每个 diff 映射到当前函数和代码区域。

### Task 2: 编写明细文档
- [ ] 新增 `docs/RPC_PMT_CODE_DIFF_DETAILS_2026-06-15.md`。
- [ ] 每项写明 commit、文件、函数、修改前后代码、原因、影响和验证。
- [ ] 单列未修改问题，避免与已完成修改混淆。

### Task 3: 增加入口和验证
- [ ] 在 `README.md` 增加文档链接。
- [ ] 核对十个 commit、三个文件统计和关键代码片段。
- [ ] 运行 `git diff --check`，提交并推送 GitHub。
