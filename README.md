# goldbach_dyna

筛法 ↔ 符号动力学视角下的 Goldbach 有限层模型（**mvp1**）。

本仓库整理阶段性实验、数值核对与一篇短文草稿。它**不声称**证明哥德巴赫猜想，也不声称改进 Chen / Pintz 等既有结果。

## 现状（2026-11-02）

| 层级 | 状态 |
| --- | --- |
| 奇异级数 / Haar 柱比、窗口参数 \(u\)、有限层例外集 | 已整理进短文与数值 |
| CML 同步流形（里程计层） | 数值与草稿表述 |
| 精确 \(u_*\)（周期 210，\(k\le 4\)） | 主值约 40 位；D4 大尺度相对 Haar 误差 \(<10^{-3}\) |
| \((H_4)\) | **强数值支持 / 部分**（未标已证） |
| \((H_\infty)\)、\(k\ge 5\) 地址共轭、奇偶障碍、minor arcs | **仍开放** |
| \(W_2^*W_3^*W_5^*W_7\) | **条件性命名**（论文1 身份仍待核） |

主值（十进制，与 `mvp1/u_star/REPORT.md` 一致）：

```text
u_* = 1.4855245007561368721556607962061545768951…
```

MS 词：`RLLLLLC`（候选集中最接近七位转述值者；九候选中有六个四舍五入到 `1.4855245`）。

## 目录

当前主交付在 [`mvp1/`](mvp1/)：

- `goldbach_paper.tex` — 合并短文草稿（CJK）
- `goldbach_props.tex`、`goldbach_G2.tex`–`goldbach_G4.tex` — 分阶段笔记
- `starter.py`、`validate_stage.py`、`stage2.py`、`stage3_d4.py`、`stage4_d5.py` — 实验脚本
- `results/` — V*/D* 数值与图
- `u_star/` — 高精度求根、柱字典、D4 重跑与报告
- `REPORT*.md`、`STARTER.md` — 阶段报告

更细的入口说明见 [`mvp1/README.md`](mvp1/README.md)。

## 快速复现（示意）

在装有 Python 3 与虚拟环境的本机目录中（以 `mvp1/` 为工作根）：

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt   # 若仓库提供；否则按脚本依赖安装
.venv/bin/python starter.py                 # 五步起始实验 → results/starter.json
.venv/bin/python validate_stage.py          # 阶段校验
# 精确 u_* / D4 重跑（产物只写 u_star/）：
.venv/bin/python -B u_star/roots.py
.venv/bin/python -B u_star/experiment.py
.venv/bin/python -B u_star/verify.py
```

具体成功判据与禁止事项以各 `spec*.txt` / `u_star/REPORT.md` 为准。编译短文需本地 LaTeX（CJK）。

## 明确不能越过的边界

- 零熵 clopen 可观测量仍受 Möbius 正交 / 奇偶型障碍约束；动力学改写**不能**单独绕过奇偶障碍。
- 周期 210 符号不能区分 \(n\) 与 \(n+210\) 上不同的更高层 gcd 标签，故 \(k\ge 5\) 地址字典在本设定下结构上不可建成。
- 数值稳定与柱匹配 ≠ 形式共轭证明；\(W_p\) 身份与 paper1 原文对照仍缺。

## 许可与联系

研究草稿仓库；暂无单独 LICENSE。维护者 GitHub：[@maris205](https://github.com/maris205)。
