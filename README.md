# DM-Transport

太阳暗物质输运研究项目，参考 DaMaSCUS-SUN-EVAP 的太阳背景与微观散射，研究暗物质占据数、密度和蒸发通量。

项目处于初步开发阶段，已具备有限体积网格、保守源投影、冻结 AGSS09 太阳背景，以及固定 MVP（`m_chi=0.1 GeV`、SD constant-contact、proton-only）的 63 靶直接散射率与完整三维局域单碰撞。T03 局域物理科学验证已通过；输运求解器尚未实现。

## 构建与测试

需要 CMake 3.16+、C++14 编译器和 Python 3.9+。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build -L fast --output-on-failure
```

T03 局域碰撞的科学验证需手动运行大样本；默认读取冻结的状态、种子和验收阈值，并将报告写入 `Output/Result/validation/validation_report.json`：

```sh
python3 Code/python/run_t03_validation.py
```

`--smoke` 只检查采样器与报告链路，报告会标为 `physics_validation=not_evaluated`。已按标准契约完成全量运行：`physics_validation=pass`，G0 物理门按条件关闭；本地报告位于 `Output/Result/validation/validation_report.json`。全量报告只验收选定的局域碰撞观测量；没有独立 legacy artifact，`reference_parity=not_evaluated`，不能称为与 DaMaSCUS 数值复现。自定义或过期构建只作诊断。

太阳背景测试可单独运行 `ctest --test-dir build -R solar_background --output-on-failure`；输入来自仓库内的 [AGSS09 太阳表](Code/data/solar/model_agss09.dat) 与 [核数据表](Code/data/solar/Nuclear_Data.txt)，固定来源和许可证见 [physics provenance](Code/provenance/reference_physics.json)。`sd_rate_diagnostic` 输出各靶的截面、直接散射率、占比和总率。参数见 [MVP 配置](Code/configs/benchmark/mvp.json)。[基线工具](Code/python/run_baseline.py) 用于读取本项目内保存的 DaMaSCUS-SUN-EVAP 基线产物，使用 `--help` 查看参数。参考仓库严格只读，不在其中修改、编译或运行；所需实现须移植到本项目后再构建。`fast` 测试覆盖 unit、contract 和小型物理 smoke；全量科学运行与可用时的外部 parity 由 [T03 契约](Code/configs/validation/t03_oracle_contract.json) 约束，不能由默认 CI 的通过代替。

## 文档

| 文档 | 内容 |
| --- | --- |
| [README](README.md) | 项目入口与构建说明 |
| [设计方案](Note/Proposal.md) | 科学方法、工程设计与源码依据 |
| [任务计划](Note/Task_Plan.md) | 任务依赖、进度与验收标准 |
| [变更记录](CHANGELOG.md) | 简短变更摘要 |

- 只维护上述四份文档，原位更新；跨文档使用链接，每次修改检查重复内容、过期结论和引用。
- [已有原稿](Note/archive/Proposal_2026-09-14_v1.md) 保持只读，不新增逐轮副本。代码与配置放在 `Code/`，生成产物放在 `Output/`，由 `.gitignore` 排除。
- 公开文档使用仓库相对路径或公开来源链接，不包含个人目录、私有附件地址或凭据。
