# DM-Transport

太阳暗物质输运研究项目，参考 DaMaSCUS-SUN-EVAP 的太阳背景与微观散射，研究暗物质占据数、密度和蒸发通量。

项目处于初步开发阶段，已具备有限体积网格、所属单元保守源投影、Maxwell 热浴平均相对速度，以及项目内 AGSS09 太阳背景查询（温度、半径内质量、逃逸速度和 63 个同位素靶的数密度）；完整散射链、输运求解器及物理验证尚未完成。

## 构建与测试

需要 CMake 3.16+、C++14 编译器和 Python 3.9+。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

太阳背景测试可单独运行 `ctest --test-dir build -R solar_background --output-on-failure`；输入来自仓库内的 [AGSS09 太阳表](Code/data/solar/model_agss09.dat) 与 [核数据表](Code/data/solar/Nuclear_Data.txt)，固定来源和许可证见 [physics provenance](Code/provenance/reference_physics.json)。参数见 [MVP 配置](Code/configs/benchmark/mvp.json)。[基线工具](Code/python/run_baseline.py) 用于读取本项目内保存的 DaMaSCUS-SUN-EVAP 基线产物，使用 `--help` 查看参数。参考仓库严格只读，不在其中修改或编译；所需实现须移植到本项目后再构建。默认测试覆盖上述实现与项目约定，尚不代表输运物理验收通过。

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
