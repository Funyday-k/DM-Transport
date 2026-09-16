# SolarDM-Transport

## 从 DaMaSCUS-SUN 轨迹模拟器到高光学厚度太阳暗物质相空间输运求解器

**项目工作名：** `SolarDM-Transport`
**核心目标：** 从 DaMaSCUS-SUN 已有的微观散射物理出发，建立保留轨道相位的多尺度动理学输运框架，求解暗物质的绝对相空间占据数、径向数密度、蒸发通量和太阳外弱束缚分布。直接输运求解为主线，完整轨迹和受控加速轨迹提供数值验证。

**维护说明（2026-09-16）：** 本稿保留原有 75 节科学路线，依据 DaMaSCUS-SUN-EVAP 源码、数值方法审阅与已核验文献，修正接口、捕获归一化、有限体积网格、边界和物理时间等约定。科学推导、工程设计与源码依据统一维护在本文；项目固定文档规则见 [README](../README.md). 原始文本逐字保存在 [初稿存档](archive/Proposal_2026-09-14_v1.md)。可执行科学任务、依赖与 G0–G6 验收以 [Task_Plan.md](Task_Plan.md) 为准；源码事实见 [源码依据附录](Proposal.md#code-reference)；后端、缓存、核结构及性能任务见 [工程与性能设计附录](Proposal.md#performance-design)。文中的 Milestone / PR 保留为路线说明，不另行改变 T / G 科学任务编号与顺序。

当前执行状态及证据统一见 Task_Plan。DaMaSCUS-SUN-EVAP 是严格只读素材库，全部新增代码和构建留在 DM-Transport。本稿区分已有能力与待实现能力，任何物理正确性、收敛性或性能收益都需对应验证后才能成立。V1 先闭合 CPU reference 的绝对占据与 MC 对照；optimized 后端复用同一物理和数据契约，分别验收。

设计将无碰撞假蒸发、完整单元平均和源投影列为 V1 门槛；将真实返回时间落实为 T17 的延迟边界契约；将 AP、确定性 kinetic–diffusion、DDMC 列为 T18 的三个候选。WE、伴随采样及宏态约化按证据逐步启用，不要求 MVP 同时实现所有算法。新增意见及文献结论在对应章节原位更新，不新建独立审阅或性能文档。

---

# 1. 项目为什么需要重构

当前 DaMaSCUS-SUN 的计算范式是

$$
\text{initial condition}
\rightarrow
\text{free propagation}
\rightarrow
\text{scatter}
\rightarrow
\text{free propagation}
\rightarrow
\cdots
\rightarrow
\text{escape/capture}.
$$

这对于较小截面非常自然。

但是当

$$
\sigma_{\chi p}\sim10^{-30}\ {\rm cm^2}
$$

甚至更大时，会同时发生两个问题：

1. 平均自由程快速下降；
2. 单位物理时间内散射次数迅速增加。

DaMaSCUS-SUN 当前 `Propagate_Freely()` 使用 RK45 推进轨迹，并在太阳内部根据

$$
\Gamma_{\rm scat}(r,v)
$$

限制积分步长，大致要求

$$
\Delta\tau=\int_t^{t+\Delta t}\Gamma_{\rm scat}(r(s),v(s))\,ds\lesssim\tau_{\max}.
$$

当前 normal mode 的光学深度步上限为 0.05，Capture mode 为 0.10（`Simulation_Trajectory.cpp`）；局域 rate 近似不变时对应 Delta t 约为 tau_max/Gamma。然后 `Simulate()` 每遇到散射调用 `Scatter()` 更新速度，继续传播。高截面可能同时缩小步长并增加散射事件数，成本还受轨道与自适应控制影响。

这意味着：

$$
T_{\rm CPU}
\propto
N_{\rm trajectory}
\times
N_{\rm scattering/trajectory}
$$

描述事件数主导时的粗略成本关系，不是对所有模型的精确标度；它说明逐事件轨迹在高光学厚度区域可能成为瓶颈。

我们的目标是把它改成

$$
\boxed{
\text{microscopic physics}
\rightarrow
\text{transport operator}
\rightarrow
\text{sparse linear solve}
}
$$

即：

$$
\boxed{
Q(m_\chi,\sigma)
\rightarrow
N(r,v,\mu)
}
$$

而不是：

$$
\boxed{
10^5\text{ particles}
\times
10^6\text{ scatterings}
}
$$

---

# 2. 新项目最重要的思想

我们真正需要计算的是粒子在相空间中的 **occupation measure**。

定义状态

$$
x=(r,v,\mu),
\qquad
\mu\equiv
\frac{\mathbf v\cdot\mathbf r}{vr}.
$$

其中：

* \(r\)：太阳中心距离；
* \(v\)：DM 速度；
* \(\mu\)：径向速度方向，\(\mu>0\) 向外，\(\mu<0\) 向内。

对于一个被捕获粒子，其未来运动只取决于当前

$$
(r,v,\mu)
$$

以及太阳背景。

因此 microscopic dynamics 本身是 Markovian 的。

与只保留 \((E,L)\) 的轨道平均方法相比，本项目保留轨道上的位置和径向运动方向，目标覆盖 \(t_{\rm scat}\lesssim T_{\rm orbit}\) 乃至更短的情形。该区别是待验证的物理适用范围，不表示 master equation 或非热分布本身是新方法。

<a id="related-work"></a>

## 2.1 直接相关文献与贡献边界

下表依据第一手正文及截至 2026-09-14 可核验的版本，列出与当前设计最直接的比较；不是穷尽检索或新颖性证明。

| 文献 | 已有方法与重合 | 本项目应比较的区别 |
| --- | --- | --- |
| [Liang 等，2016，JCAP 09 018，§3.1](https://arxiv.org/html/1606.02157#S3.SS1) | 用 MC 构建 E/L 散射矩阵并求非热束缚分布；明确完整状态还需轨道相位，以小截面/长自由程支持平均掉相位 | 检验相位分辨的开放源汇系统、内部解束缚后再散射与实际表面逃逸；该文 SI 模型不直接验证本项目 SD 基准 |
| [Blennow 等，EPJC 78 (2018) 386，修订 v3](https://arxiv.org/html/1802.06880v3) | 非弹性 DM 的轨道散射矩阵；式(10)已演化带捕获源的绝对粒子数，式(20)假定相互作用之间经过许多轨道 | 比较轨道平均失效的范围；绝对占据、有限年龄和平衡失效均不能单独宣称为本项目首次提出 |
| [Busoni 等，2017，JCAP 10 037，§2.3–2.4](https://arxiv.org/pdf/1703.07784) | SI/SD 与速度/动量相关相互作用、热运动及薄厚过渡；用 isothermal/LTE 插值分布和角平均消光因子，略去额外多次散射和各向异性修正 | 在匹配模型下比较近似误差；本项目拟直接计算非热分布与再散射，不是首个光学厚蒸发处理 |
| [Wen，2026，NPB 1027 117471，v3 §2](https://arxiv.org/html/2509.26192v3#S2) | 在预设的太阳介质诱导吸引势中，用 SI O1 散射和轨道平均 E/L 方程求归一化非热形状，再单独处理总粒子数 | 额外平滑势垒与普通引力下的多次再散射抑制是不同机制；不能把该文直接当作 phase-resolved 高 opacity 再散射验证 |

版本需要固定：[Blennow v3](https://arxiv.org/abs/1802.06880) 为 2019-05-08 的程序错误修订，部分重要结论改变；[Wen v3](https://arxiv.org/abs/2509.26192) 为 2026-05-06，已发表于 [Nuclear Physics B](https://doi.org/10.1016/j.nuclphysb.2026.117471)。对 Wen 的比较以 v3 明确模型为准，不能沿用早期摘要中将势垒与多次散射混写的表述。

项目贡献目标是：在统一微观物理、源、边界和真实时间约定下，连接相位分辨动力学与受控扩散极限，保留近逃逸非热尾并耦合弱束缚外域驻留。各部分及其组合需要用实现、误差预算和文献比较证明；当前不作“首次”或“据我们所知首次”的结论。跨领域算法的可迁移依据与限制集中于附录 A.10。

---

# 3. 新项目的准确数学定位

论文和代码中最好不要仅仅把它称作“Markov chain”。

更准确的名称应该是：

> **Continuous-time phase-space transport solver for captured solar dark matter**

或者：

> **A kinetic / linear-Boltzmann solver constructed from DaMaSCUS microscopic scattering physics**

Markov generator 是它数值离散以后得到的对象。

完整动力学为

$$
\boxed{
\frac{\partial f}{\partial t}
+
\mathcal L_{\rm grav} f
=
\mathcal C[f]
+
S
}
$$

其中：

* \(\mathcal L_{\rm grav}\)：引力下 collisionless streaming；
* \(\mathcal C\)：DM–solar target scattering；
* \(S\)：capture source。

如果先不考虑 DM-DM 相互作用和 annihilation，那么这是一个**线性问题**。

这非常重要。

因为这意味着最终可以写成

$$
\boxed{
\frac{d\mathbf N}{dt}
=
\mathbf C
+
Q^T\mathbf N
}
$$

其中

$$
N_\alpha
$$

就是处在离散相空间状态 \(\alpha\) 中的真实粒子数。

本文统一令 \(Q_{\alpha\beta}\) 表示 \(\alpha\to\beta\) 的速率，\(\mathbf N\) 和 \(\mathbf C\) 为列向量，前向演化始终使用 \(Q^T\)。\(N_\alpha=\int_{\mathrm{cell}\ \alpha}Jf\,dr\,dv\,d\mu\)，不是节点上的裸 \(f\)。完整物理时间演化与删去太阳外等待的内部稳态算符需分开，详见 §25、§29–33。

连续微观过程是引力确定性传播与随机碰撞组成的 Markov 过程；有限单元上的 CTMC 则是一种数值近似。把多个微观状态合成宏态后，入口位置、速度、历史和驻留年龄可能继续影响下一次出口，不能由微观 Markov 性推出宏态也严格 Markov。

T20 若用短轨迹构造宏观响应，应保存联合首次出口核 \(K_{AB}(\tau)\)、入口通量分布及必要的历史标签。只有验证驻留时间近似指数、出口概率和入口历史相容后，才能用 \(q_{AB}=p_{AB}/\langle\tau_A\rangle\) 近似 CTMC。满足半马尔可夫闭合且入口统计正确时，嵌入链与平均驻留可重建稳态占据；完整瞬态、寿命分布仍需等待核。太阳外确定性返回是这一问题的可解析测试例，不能仅靠增加样本数消除错误闭合。

---

# 4. 为什么它能够绕过几百万次散射

考虑只有 scattering 的离散状态转移矩阵 \(P\)。

一次碰撞：

$$
P.
$$

两次碰撞：

$$
P^2.
$$

三次碰撞：

$$
P^3.
$$

完整轨迹 MC 实际上是在逐个执行

$$
P\rightarrow P\rightarrow P\rightarrow\cdots
$$

直到粒子逃逸。

先用吸收事件链说明求和思想：若 P 是 transient 子矩阵且其谱半径小于 1，预期访问次数包含

$$
I+P+P^2+P^3+\cdots
$$

即

$$
\boxed{
(I-P)^{-1}.
}
$$

所以如果一个粒子平均需要

$$
10^6
$$

次散射才逃逸，trajectory MC 就必须真的做大约 \(10^6\) 次事件。

而 transport solver 解的是

$$
\boxed{
A\mathbf N=\mathbf C
}
$$

通过线性方程求解汇总多次跃迁的贡献；实际成本仍取决于网格、稀疏结构、条件数和预处理。

实际程序绝对不要计算显式矩阵逆。

应该求：

$$
\boxed{
-Q^T\mathbf N=\mathbf C.
}
$$

这里 \((I-P)^{-1}\) 只给事件链的访问次数；无吸收的归一化碰撞矩阵具有单位特征值，不能直接求该逆。物理 occupation 还需要各态驻留时间以及 streaming / boundary。实际构造连续时间生成元，再求 \(-Q^T N=C\)；这一步不会自动保证有限稳态存在或消除重尾。

---

# 5. 新项目与 DaMaSCUS-SUN 的关系

**不要把 DaMaSCUS-SUN 删除。**

它有三个角色：

### 第一，microscopic physics 素材库

目前已有：

* `Solar_Model`
* `Total_DM_Scattering_Rate`
* `DM_Scattering_Rate_Nucleus`
* `Sample_Target`
* `Sample_Target_Velocity`
* obscura 的 `Sample_Scattering_Angle_Nucleus/Electron`
* `New_DM_Velocity`
* `Scatter`
* solar potential
* solar temperature
* nuclear abundances

这些源码用于审计和最小移植，不作为本项目链接库，也不在原仓库修改或编译。参考提交、源码位置、许可证和依赖版本随移植记录；可重复性范围与验证门统一见 Task_Plan。

### 第二，ground-truth trajectory benchmark

T01 已保存 `10^{-36}` 和 `10^{-34}\,\mathrm{cm^2}` 的固定种子小样本。这些输出用于初步运行回归，不能作为独立的单碰撞 golden 或已收敛的 transport benchmark。需要新增 trajectory 统计时，应在 DM-Transport 内移植必要的受控入口，或由用户另行提供新的只读基线产物。

### 第三，capture/source 算法依据

参考 Capture mode 已在首次散射后负能量时停止。T10 在 DM-Transport 内按固定实现和基线移植首捕获状态、入射物理率、计数与误差接口，形成 freshly captured source；不回写参考仓库。

---

# 6. 不建议复制整个代码，新仓库应该这样产生

第一阶段不 fork、不修改也不在 DaMaSCUS-SUN-EVAP 中编译。T02 从固定提交审计并移植串行微观物理的最小闭包，删除 MPI、trajectory、snapshot 和参数扫描依赖；每个移植文件记录来源位置、许可证及有意差异。项目不通过外部源码路径、`add_subdirectory` 或参考仓库 build 产物建立隐式耦合。

```text
DaMaSCUS-SUN-EVAP 固定提交（只读源码事实）
              + T01 固定输出
                       │ 审计、最小移植
                       ▼
DM-Transport: transport_reference_physics
                       │
                       ▼
              transport solver
```

这意味着两边不再共享同一个编译 target；可信度由固定来源和 T03 回归维持。若参考算法以后升级，先建立新的来源版本与基线，再把差异作为独立变更移植，不能直接改动素材库或无记录地复制公式。

---

# 7. 第一阶段必须进行的项目内移植

初始基线的 `Simulation_Trajectory.cpp` 同时包含：

* propagator；
* 调用 `Solar_Model` 的 scattering rate；
* target sampling；
* thermal target velocity 与 obscura 角采样调用；
* collision；
* trajectory control。

T02 只在 DM-Transport 中拆出并移植散射链；参考文件保持原样。传播器与 trajectory control 继续作为源码依据，不进入本阶段的移植闭包。

<a id="solar-background"></a>

## 7.1 `SolarBackground`

`SolarBackground` 是冻结的 AGSS09 reference 实现，不是通用恒星背景；在出现经过验证的第二种背景模型前，不提前引入通用接口或可配置表结构。调用方显式传入项目数据路径，并以 cgs 使用以下核心接口：

```cpp
double temperature_K(double r_cm);
double mass_enclosed_g(double r_cm);
double escape_speed_cm_s(double r_cm);
double number_density_cm3(std::size_t target_index, double r_cm);
```

固定 AGSS09 表包含 1968 个数值行；项目内加载器补中心与光球合成端点，并将 29 个丰度列按核数据映射为 63 个同位素靶。各靶数密度先在每个径向节点由质量密度、丰度、同位素比例和核质量计算，再作 Steffen 单调三次插值；逃逸速度平方在节点用 30 点 Gauss–Legendre 积分构造，再作同类插值。外域延拓和非法输入的异常行为属于接口契约，调用端不依赖隐式当前目录或参考仓库路径。

这一层只描述太阳。

不得包含 trajectory logic。

---

## 7.2 `ScatteringPhysics`

参考实现中的实际链条是：

```text
Solar_Model rate
→ Sample_Target
→ Sample_Target_Velocity
→ obscura Sample_Scattering_Angle_Nucleus/Electron
→ New_DM_Velocity
→ Scatter 更新速度
```

MVP 固定为 `m_chi=0.1 GeV`、SD constant-contact 和 proton-only coupling（`a_n=0`），但仍保留全部 63 个太阳同位素靶。直接率沿用固定 obscura 归一化：对 `J_A>0`，

$$
\sigma_A^{\rm SD}
=\sigma_p\frac{4}{3}
\frac{\mu_{\chi A}^2}{\mu_{\chi p}^2}
\frac{J_A+1}{J_A}
\langle S_p\rangle_A^2,
\qquad
\Gamma_A(r,v)=n_A(r)\sigma_A^{\rm SD}\langle v_{\rm rel}\rangle_A,
$$

自旋为零的靶取 `sigma_A=0`，总率按固定靶顺序求和 `Gamma_total=sum_A Gamma_A`。直接接口返回全部 63 个分靶记录，包括零率靶，并同时返回总率；本阶段不引入 rate 插值或缓存。

T02c 先按 `Gamma_A/Gamma_total` 从上述 63 靶中选择核素。选定靶后，恒截面模型的事件条件热速度不是未加权 Maxwell 分布，而是

$$
p(\mathbf u_A\mid A,\mathbf v_\chi,\mathrm{collision})
=\frac{f_{\rm MB}(\mathbf u_A;T,m_A)
\lvert\mathbf v_\chi-\mathbf u_A\rvert}
{\langle v_{\rm rel}\rangle_A}.
$$

项目内实现保留固定参考的 Romano–Walsh 混合 proposal、拒绝步骤、方向基和 `std::mt19937` 随机调用顺序，同时只在 API 边界使用 cgs。稳定接口为：

```cpp
struct SdProtonModel {
    double dark_matter_mass_GeV;
    double proton_cross_section_cm2;
};

SdScatteringRates direct_sd_proton_scattering_rates(
    const SolarBackground&, const SdProtonModel&,
    double radius_cm, double dm_speed_cm_s);

using CartesianVelocityCmS = std::array<double, 3>;

std::size_t sample_sd_proton_target_index(
    const SdScatteringRates&, std::mt19937& rng);

CartesianVelocityCmS sample_collision_conditioned_target_velocity_cm_s(
    double temperature_K, double target_mass_GeV,
    const CartesianVelocityCmS& dm_velocity_cm_s,
    std::mt19937& rng);

struct CollisionSample {
    std::size_t target_index;
    CartesianVelocityCmS target_velocity_cm_s;
    CartesianVelocityCmS outgoing_dm_velocity_cm_s;
};

CollisionSample sample_sd_proton_collision(
    const SolarBackground&, const SdProtonModel&, double radius_cm,
    const CartesianVelocityCmS& incoming_dm_velocity_cm_s,
    std::mt19937& rng);
```

T02d 已按上述顺序组合完整局域碰撞。设入射 DM 和靶速度为 `v_chi,u_A`，质量为 `m_chi,m_A`，CM 速度为 `V_CM`，各向同性出射单位方向为 `n`，则

$$
\mathbf v_\chi'=\mathbf V_{\rm CM}
+\frac{m_A}{m_\chi+m_A}
\lvert\mathbf v_\chi-\mathbf u_A\rvert\mathbf n,
\qquad
\mathbf V_{\rm CM}=\frac{m_\chi\mathbf v_\chi+m_A\mathbf u_A}{m_\chi+m_A}.
$$

实现另提供确定性 `elastic_outgoing_dm_velocity_cm_s()`，使固定角极限、动量/能量守恒和旋转协变不依赖随机抽样即可验证。完整接口只返回靶索引、碰撞条件靶速度和出射 DM 速度，不含分箱、束缚/逃逸判定、trajectory 状态或 kernel 逻辑。

固定 legacy sampler 要求 `|v_chi|>0`；项目测试以小正速度验证解析单侧极限，不把它冒充零速 reference parity。相同 seed 的逐事件复现仅承诺相同 executable 与标准库；跨工具链只检查统计相容。

本版没有独立 `Sample_Momentum_Transfer()`。主碰撞接口以真实三维速度向量为输入/输出，transport 的 `(v,mu)` 投影由其适配层负责。靶速度、能量交换、动量转移可作为后续诊断；如输出 q，应由实际碰撞前后动量差计算。文件字段按 Task_Plan 使用 `r_cm`、`v_cm_s`、`rate_s_inv` 等单位名。

若保留 `cos_scattering_angle` 诊断，它必须标明是旧实现围绕入射实验室 DM 速度轴抽取的变量，不是一般运动靶标下的相对速度散射角。旧角采样和 rate 热平均并不自动适用于一般速度相关或各向异性相互作用。low-mass 分支是显式配置开关，不是由低质量自动启用。移植先保持旧抽样映射，任何物理修正单独建立回归基线。

---

# 8. 最重要的原则：旧 trajectory simulator 保持只读

旧 trajectory simulator 只提供固定源码事实和 T01 基线，不再承担新模块客户端。DM-Transport 在自己的命名空间中维护最小移植物理层，并以固定来源记录、固定向量和分布门检查一致性：

```text
只读 Trajectory_Simulator ──► T01 reference outputs
                                      │
                                      ▼ parity
DM-Transport ScatteringPhysics ──► TransportOperator
```

T03 只接受按 [oracle contract](../Code/configs/validation/t03_oracle_contract.json) 在 DM-Transport 外独立生成并冻结的 legacy artifact；当前物理实现不得参与期望值生成。用户规则禁止本项目构建或运行参考仓库，因此 artifact 缺失时报告仍保留 `reference_parity`，状态明确为 `not_evaluated`，不能用源码派生或解析检查冒充。`physics_validation` 是独立必过门，覆盖运动学、角分布、旋转对称性和热浴检查；具体 G0 条件由 Task_Plan 维护。

T03 在固定半径关闭 streaming，直接检验完整局域碰撞链。若入射速度从温度为当地 `T` 的 Maxwell 分布抽样，单次碰撞的弱平衡残差为

$$
R_\varphi=\mathbb E_{\mathbf v\sim f_{\rm MB}(T)}\!\left[\frac{\Gamma_{\rm tot}(r,|\mathbf v|)}{\Gamma_{\rm ref}(r)}\left(\varphi(\mathbf v')-\varphi(\mathbf v)\right)\right]\simeq0,
\quad \Gamma_{\rm ref}(r)=\Gamma_{\rm tot}(r,\sqrt{2k_BT/m_\chi}).
$$

测试函数限于预先固定的六个速度区间指标与无量纲能量 `E_chi/(k_B T)`；`Gamma` 权重不可省略，因为无条件 Maxwell 入射样本不是碰撞事件的入射分布。这些弱式检查是对选定观测量的验收，不能单独证明完整分布的平衡或详细平衡。配套检验低温 DM 加热、高温 DM 冷却、完整链在旋转后的分布相容，以及以入射实验室 DM 方向为轴的出射 CM 角余弦均匀性。另报告入射速度靠近当地逃逸速度时的 `P(|v'|>v_esc)`、分靶贡献和数值域外溢；它只是一次局域碰撞的尾部诊断，不是蒸发率。每项在运行前固定状态、种子、样本预算和验收区间，保存独立批次及失败信息；详细报告契约见 T03 配置。

冻结全量 T03 运行已通过这些选定的局域物理门，证据与 G0 条件见 [Task_Plan §8](Task_Plan.md#current-status)。与 legacy 的数值 parity 仍因缺少合规独立 artifact 而未评估；来源和差异审计不能替代数值比较，选定弱式通过也不证明完整分布或详细平衡。

---

# 9. 第一个强制回归测试

在真正开始 Markov solver 前，先完成：

```text
old Scatter()
vs
DM-Transport ScatteringPhysics::sample_collision()
```

对于固定

$$
(r,v,m_\chi,\sigma)
$$

生成例如

$$
10^6
$$

次单散射。

比较：

$$
P(v_{\rm out}),
$$

$$
P(\cos\theta),
$$

$$
\langle\Delta E\rangle,
$$

$$
P(A),
$$

以及

$$
\Gamma_{\rm scat}.
$$

要求在预先规定的统计与数值误差内一致；样本预算依尾部事件精度决定，以上数量仅为示例。直接移植尽可能保留固定种子逐事件映射；物理修正引起的随机映射变化单独记录并做分布回归。

**这一关没通过，不允许继续。**

---

# 10. 相空间变量：第一版使用 \((r,v,\mu)\)

第一版建议使用

$$
\boxed{
x=(r,v,\mu).
}
$$

原因不是它在理论上最漂亮，而是它：

* 是规则网格；
* 容易与 scattering sampler 对接；
* 容易重构；
* 容易检查；
* 不需要处理复杂的 allowed \(E-L\) domain。

定义：

$$
\mu=
\frac{v_r}{v}.
$$

于是：

$$
v_r=v\mu,
$$

$$
v_\perp=v\sqrt{1-\mu^2}.
$$

这一约化用于静态球对称背景下的径向统计。MVP 不由 `(r,v,mu)` 重建相对太阳运动方向的角各向异性信号；扩展相互作用和观测量前需重新检查对称性。

---

# 11. 引力 streaming 的精确方程

令

$$
g(r)=\frac{GM(r)}{r^2}.
$$

则 collisionless dynamics 满足

$$
\dot r=v\mu,
$$

$$
\dot v=-g(r)\mu,
$$

$$
\boxed{
\dot\mu=
(1-\mu^2)
\left[
\frac vr-\frac{g(r)}v
\right].
}
$$

因此 Boltzmann equation 可以写成

$$
\frac{\partial f}{\partial t}
+
v\mu\frac{\partial f}{\partial r}
-
g\mu\frac{\partial f}{\partial v}
+
(1-\mu^2)
\left(
\frac vr-\frac gv
\right)
\frac{\partial f}{\partial\mu}
=
\mathcal C[f]+S.
$$

上式在 \(r>0,v>0\) 的常规坐标区域使用；中心穿越、零速与 \(\mu=\pm1\) 的离散需要守恒极限或坐标对称处理，不在奇点直接代入 \(v/r\) 或 \(g/v\)。

---

# 12. 一开始就采用 conservative finite-volume formulation

不要简单在三个方向上写 finite difference。

球对称下相空间体积元是

$$
d\Gamma
=
8\pi^2r^2v^2
\,dr\,dv\,d\mu.
$$

定义

$$
J=8\pi^2r^2v^2.
$$

应该离散：

$$
\boxed{
\frac{\partial(Jf)}{\partial t}
+
\frac{\partial(J\dot r f)}{\partial r}
+
\frac{\partial(J\dot v f)}{\partial v}
+
\frac{\partial(J\dot\mu f)}{\partial\mu}
=
J\mathcal C[f]+JS.
}
$$

而不是直接离散 \(f\)。

这样能够：

* 保证 particle number conservation；
* 构造 positive transition rates；
* 让最终离散矩阵真正具有 Markov generator / M-matrix 的性质。

实现直接演化单元占据数：

$$
N_{ijk}=\int_{\mathrm{cell}_{ijk}}Jf\,dr\,dv\,d\mu,\qquad
\Delta\Gamma_{ijk}=\frac{8\pi^2}{9}
(r_{i+1/2}^3-r_{i-1/2}^3)
(v_{j+1/2}^3-v_{j-1/2}^3)
(\mu_{k+1/2}-\mu_{k-1/2}).
$$

共享面通量、正性和边界闭合要分别验证；采用有限体积写法本身不保证实现已满足这些性质。

对穿越 \(E=(v^2-v_{\rm esc}^2(r))/2=0\) 的单元，在 \(x=r^3,\ y=v^3\) 坐标中定义 \(y_{\rm esc}(x)=v_{\rm esc}^3(x^{1/3})\)。束缚测度比例是

\[
f_{\rm bound}=\frac{\int_{x_-}^{x_+}\operatorname{clip}(y_{\rm esc}(x)-y_-,0,y_+-y_-)\,dx}{(x_+-x_-)(y_+-y_-)},\qquad f_{\rm unbound}=1-f_{\rm bound}.
\]

\(\mu\) 因子约去。穿越格仍是一个全局 state；子格积分在 \(x\) 上细分，再对阈值两侧的条件 \(y\) 区间取正权节点。测度比例只用于连续单元占据的分解；一个已知 \((r,v,\mu)\) 的碰撞输出按所属 FV 单元一次沉积，不按比例拆分。最小跨语言 schema 固定 faces、`r/v/mu` 展平顺序、单位、求积阶数、阈值与点沉积约定及版本；生产数据格式留待后续任务。

<a id="streaming-gate"></a>

## 12.1 T07 必须单独排除数值假蒸发

无碰撞真实轨道保持比能量 \(E=v^2/2+\Phi(r)\) 和比角动量 \(L=rv\sqrt{1-\mu^2}\)。普通一阶 upwind 即使严格守恒粒子数，仍可能将分布数值扩散到 \(E=0\) 的另一侧。该误差可能远大于真实稀有蒸发，必须在 T07 检验，不能等到高光学厚度 T18。

采用 \(E_0=GM_\odot/R_\odot>0\)、\(L_0=R_\odot\sqrt{E_0}\)，构造 \(E/E_0=-10^{-1},-10^{-2},-10^{-3},-10^{-4}\) 的可达束缚轨道组，覆盖径向、近切向及中心穿越。关闭碰撞、使用 `kepler_return`，对声明的多个轨道周期和共同物理观测时窗比较独立精确/高精度轨道与网格传播；过长周期可用经验证的外域解析段，不能用即时回流的时间声称完成真实多周期测试。

诊断至少包含：E/L 的均值漂移与展宽、进入非束缚区域的占据以及实际外行表面逃逸的累计比例。网格量用单元内积分诊断，不能只由格心判别；初始化投影已落在 \(E\ge0\) 的质量必须单列，不归因于后续传播。近零 E、径向 L=0 时报告 \(|\Delta E|/E_0\)、\(|\Delta L|/L_0\)，仅在分母充分远离零时补充相对误差。`legacy_1au_removal` 另测物理移除，不能将其计为假蒸发。

沿 \(r,v,\mu\)、面求积及时间推进容差做收敛，预先为累计假逃逸分配绝对误差预算，并在有碰撞基准中将其影响约束到目标尾部观测量的系统误差预算内。预算未定或低于可检测精度时标记 `not_evaluated`，不能用 bulk 误差小代替。固定有限网格不承诺严格零跨越，也不能通过事后钳回 \(E<0\) 制造通过。

若细化成本不可接受或泄漏不收敛，T07 即比较沿特征线的保守重映射、阈值贴合/切割单元，或保留 r 与径向分支的 \((r,E,L,s_r)\) 网格。半拉格朗日名称本身不保证守恒、正性或不变量，候选仍通过相同门槛；详细几何约束见 §64。

---

# 13. 初始网格不要太大

首个 reference prototype 固定为 `30×40×8 = 9,600` 个状态。通过小网格验证后再考虑以下区间：

$$
N_r=30\sim40,
$$

$$
N_v=40\sim60,
$$

$$
N_\mu=8\sim12.
$$

总状态数约：

$$
2\times10^4.
$$

不要第一天就做

$$
100\times100\times50.
$$

第一版目标不是 production。

第一版目标只有：

> transport solver 与 full trajectory simulation 在可比较区域是否给出相同结果？

---

# 14. \(\mu\) 首版采用明确 cell faces 的有限体积网格

MVP 使用关于 0 对称、具有明确边界面的 \(\mu\)-bins，先从等宽划分 \([-1,1]\) 开始。单元权重来自 \(\Delta\Gamma\)，碰撞概率沉积和 angular streaming 使用相同的单元定义。

初版 \(N_\mu=8\)，随后做 \(8\to12\to16\to24\) 收敛测试。\(\mu=\pm1\) 属于边界极限，不能因解析几何退化而静默丢弃对应通道。

Gauss–Legendre 节点和权重适合未来 quadrature 或 angular-moment 路线；裸节点不是有限体积单元。若引入 GL，必须另行定义质量权重、角度面通量、碰撞沉积规则并通过相同守恒/热浴测试，不能将最近节点分箱视为等体积。

---

# 15. 速度范围必须根据质量自动调整

低质量粒子的 thermal velocity 很高。

定义核心 thermal scale：

$$
v_T=
\sqrt{\frac{2k_B T_c}{m_\chi}}.
$$

可以暂时取：

$$
v_{\max}
=
\max
\left[
1.5\,v_{\rm esc}(0),
6v_T
\right].
$$

但这一系数必须通过

$$
4v_T,\quad6v_T,\quad8v_T
$$

进行 tail convergence test。

尤其 evaporation 完全由高速尾控制，绝对不能因为 bulk distribution 看起来收敛就认为速度网格足够。

这里 `T_c` 表示热力学温度；如沿用自然单位并把 `k_B T_c` 存为能量，必须在接口中注明。`v_max` 是数值域边界：超域概率进入显式 overflow 诊断并扩域/验证闭合，不得当作物理蒸发、静默删除或钳到末格。

---

# 16. Collision operator 的构造

固定状态

$$
\alpha=(r_i,v_j,\mu_k).
$$

首先计算：

$$
\Gamma_\alpha
=
\Gamma_{\rm scat}(r_i,v_j).
$$

然后反复调用：

```cpp
sample_collision(state_alpha)
```

获得

$$
(v',\mu').
$$

统计最终进入状态 \(\beta\) 的概率：

$$
P_{\alpha\beta}.
$$

于是 collision generator 为

$$
\boxed{
Q_{\rm coll}
=
\operatorname{diag}(\Gamma)(P-I).
}
$$

即：

$$
Q^{\rm coll}_{\alpha\beta}
=
\Gamma_\alpha P_{\alpha\beta},
\qquad
\alpha\neq\beta.
$$

而 diagonal：

$$
Q^{\rm coll}_{\alpha\alpha}
=
\Gamma_\alpha
(P_{\alpha\alpha}-1).
$$

注意：

**如果一次 collision 后仍处在同一个 cell，就不能把它错误地当作离开该状态。**

归一化包括显式追踪的输出通道，零速率行直接为零。内部 \(E'>0\) 态仍在 transient 相空间内，可能再次束缚，不能在 collision builder 中直接吸收。\(v'>v_{\max}\) 是数值 overflow，必须计数并通过扩域收敛消除；碰撞代数行和应在浮点精度内闭合，MC 误差不构成漏计概率的理由。

<a id="cell-averaging"></a>

## 16.1 完整有限体积碰撞率

上述格心规则是原型近似。对单元内常数 f 的 FV 闭合，令 \(p_\beta(x)\) 为在连续入射状态 \(x=(r,v,\mu)\) 发生一次碰撞后沉积到单元 \(\beta\) 的概率，正确的单元生成元是

\[
Q^{coll}_{\alpha\beta}
=\frac{1}{\Delta\Gamma_\alpha}\int_{\Omega_\alpha}
J(x)\Gamma(x)\,[p_\beta(x)-\delta_{\alpha\beta}]\,dx.
\]

只积分 \(J\Gamma p_\beta\) 得到的是包含自事件的转移率矩阵，尚不是生成元；对角必须扣除平均事件率。定义

\[
\bar\Gamma_\alpha=\frac{\int_{\Omega_\alpha}J\Gamma\,dx}{\Delta\Gamma_\alpha},\qquad
\bar P_{\alpha\beta}=\frac{\int_{\Omega_\alpha}J\Gamma p_\beta\,dx}
{\int_{\Omega_\alpha}J\Gamma\,dx},
\quad Q^{coll}=\operatorname{diag}(\bar\Gamma)(\bar P-I).
\]

\(\bar\Gamma=0\) 的行直接置零，无需定义条件概率。\(\bar P\) 是事件率加权的单元平均，不是均匀平均 P 后再乘平均 Gamma。瞬时碰撞不改变位置，积分后仍保持径向块结构；径向传播只属于 streaming。

采用非负、已包含 J 测度并归一为总和 1 的求积权重 \(\omega_q\)，逐点累加 \(\sum_q\omega_q\Gamma(x_q)[p_\beta(x_q)-\delta_{\alpha\beta}]\)，同时保留自事件和完整去向。T05/T06 比较现有 1 点与 `r/v/μ` 的 `2×2×2` 正权求积，再按收敛需要加密；角度使用附录 A.3 的共同 FV 规则。缓存记录求积规则、输入体积/碰撞率权重及输出核压缩方式，不能混用两种离散的资产。

先在背景变化显著及 near-escape 单元验证率、联合转移、热浴稳态，再在 T16 比较绝对占据和尾部指标。MC 采样、求积、网格闭合及核压缩误差分别报告；中心近似只有在其误差低于预设预算时才可用于科学结果。

---

# 17. Collision kernel 第一版可以继续 Monte Carlo

例如每个 state：

$$
N_{\rm sample}=500
$$

开始。

先得到粗 kernel。

求解一次 distribution 后，根据 occupancy 再 adaptive refinement。

定义重要度：

$$
I_\alpha
\sim
N_\alpha
\times
\delta P_\alpha.
$$

重点增加：

* 高 occupancy states；
* near-escape states；
* 对 escape flux 敏感的 states。

而不是所有 cell 都跑 \(10^5\) 次。

零观测不等于零概率；原型必须保存各行样本数及关键稀有通道的上限，并通过独立样本重复和预算增加检查尾部。只按 occupation 加样不能可靠识别控制蒸发的稀有跃迁；加权 importance sampling 须同时验证权重、支持域和误差估计，详见 附录 A（工程与性能设计）。

---

# 18. 后续可以把 kernel MC 换成 deterministic/quasi-MC

只读参考仓库没有 `Scattering_Rates.cpp` 或独立 q 抽样 API；rate 位于 `Solar_Model.cpp`，靶选择、热靶速度和出射速度链位于 `Simulation_Trajectory.cpp`，其中角采样调用 obscura 接口。T02 在 DM-Transport 内移植这一最小调用闭包。

未来可从实际靶热速度分布、靶选择与相对运动学出发，推导匹配模型的积分变量和权重，再尝试 quadrature 或 Sobol / quasi Monte Carlo。不能假定已有独立 q 采样 API，也不能未经核验就把旧角变量当成相对速度散射角。

这是通过 reference 微观验证后的优化，不是 MVP 物理提取的前提。

---

# 19. 对 contact SD-proton，一个极大的优势是 cross section 可以分离

如果只是整体改变

$$
\sigma_p
$$

而 interaction shape 不变，那么：

$$
P_{\alpha\beta}
$$

与整体截面大小无关。

只有 rate：

$$
\Gamma_\alpha
\propto\sigma_p.
$$

因此

$$
\boxed{
Q(m_\chi,\sigma_p)
=
Q_{\rm stream}
+
\frac{\sigma_p}{\sigma_0}
Q_{\rm coll}(m_\chi,\sigma_0).
}
$$

这意味着：

**固定质量、interaction shape、太阳背景、靶组成和网格后，经验证的 normalized collision kernel 可在不同截面间复用。**

之后：

$$
10^{-38}
\rightarrow
10^{-30}\ {\rm cm^2}
$$

的扫描可以复用核，但每点的 source、求解与收敛成本仍需测量；高截面的 stiffness 不会因缓存自动消失。

当前源码未提供经模型签名保护的自动 rescale / kernel cache 接口。这是 T05 的新能力：必须验证 Gamma 随 sigma 线性、P 在误差内不变，缓存键包含全部物理与网格条件，才可启用。源 C(sigma) 单独重算。

式中可将已固定的 boundary 处理并入 streaming；显式实现则分为 \(Q=Q_{\rm stream}+Q_{\rm boundary}+(\sigma/\sigma_0)Q_{\rm coll}^{(0)}\)。Offline 构造背景、reference-rate/kernel 和 exterior 数据；online 对各截面构造 source 并求解，条件与失效规则见 [工程与性能设计附录](Proposal.md#performance-design)。

---

# 20. 但是 capture source 不能简单随 \(\sigma\) 线性缩放

高截面进入 optical-thick regime 后：

$$
C\not\propto\sigma.
$$

最终趋向 geometric saturation。

已有从 optically thin 一直到 optically thick 的太阳 capture/evaporation 工作已经专门处理过这一问题。

所以必须区分：

$$
\boxed{
Q_{\rm transport}(\sigma)
}
$$

和

$$
\boxed{
C_\alpha(\sigma).
}
$$

---

# 21. Capture source 的 MVP 实现

只读参考的 `Enable_Capture_Mode(true)` 在首次散射后 \(E<0\) 时停止；这里 E 是比能量 \(v^2/2+\Phi(r)\)。本阶段在 DM-Transport 内移植该判据并新增首捕获 `(r,v,mu)`、三维原始状态/权重、轨迹标识和计数输出。

注意旧 `E_capture_eV` 是粒子能量，不能直接代入比能量公式。Capture mode 当前只打印汇总，尚未输出完整源样本。

令首捕获条件分布为 \(\eta_\alpha\)，则：

$$
C_\alpha=R_{\rm in}p_{\rm cap}\eta_\alpha,
\qquad \sum_\alpha\eta_\alpha=1,
\qquad C_{\rm tot}=\sum_\alpha C_\alpha.
$$

\(R_{\rm in}\) 复用并独立核对 `DM_Entering_Rate` 的引力聚焦入射率，与击中 Rsun 的 halo 抽样配对；数值注入面为 1.1Rsun 不改变该抽样总体。旧 `capture_rate_raw/valid` 是概率，`Captured particle rate [1/s]` 是墙钟吞吐量，两者都不能直接作为 \(C_{\rm tot}\)。

源生成优先固定尝试数；等权且完整分类时可直接写 \(C_\alpha=R_{\rm in}n_{{\rm cap},\alpha}/n_{\rm attempts}\)。保存 attempted、captured、classified、unresolved、failure 计数及 raw/valid 口径、权重和误差，不删除未决样本后悄悄归一化。捕获概率与 \(\eta\) 来自同批样本时保留其相关性。

参考 `Data_Generation.cpp` 提供 halo/trajectory 调用依据。普通模式和 Capture mode 的光学深度推进精度不同；项目内移植需比较首次捕获分布的统计/容差收敛，不能预设同 seed 逐事件一致。本项目另实现显式已捕获初态入口，用同一 source 运行 MC benchmark；参考 `Simulate()` 默认未捕获 halo 注入，直接传入束缚态会触发保护判定。

源投影优先使用原始连续捕获样本落入所属 FV 单元的加权计数；它对经验源的单元积分已经守恒，不因“最近格”这个称呼而需要被线性平滑替换。误差来自有限样本、单元内表示和阈值交叉格的闭合。T04/T10 比较所属格计数、加密网格，以及可选的正性保守 CIC/线性沉积；后者若跨边界、\(E=0\) 或不可达域分配权重，会产生额外偏差，不能宣称天然更准确。

每次点投影须满足 \(\sum_\alpha w_\alpha=w\)、\(w_\alpha\ge0\)，不按阈值比例拆分；单元重构和积分诊断使用 §12 的阈值子格测度，并保存初始非束缚质量。原始首次捕获态的 \(E<0\) 支持不得因数值平滑变成物理非束缚源。T11 用原始连续初态与声明的单元重构初态分别运行，以隔离源表示误差；单纯保持总 C 不能证明尾部源已正确。

---

# 22. 第二阶段再消灭 capture trajectory MC

最终可以直接把 Galactic halo 当成太阳表面的 inward boundary source：

$$
S(R_\odot,v,\mu<0).
$$

然后同一个 transport solver 自动产生：

* free transit；
* reflection；
* capture；
* thermalization；
* evaporation。

即：

$$
\boxed{
\text{halo boundary flux}
\rightarrow
\text{full solar transport}
}
$$

届时甚至

$$
C
$$

都不再需要单独输入。

但是这应该属于 Version 2，不应该挡住 Version 1。

---

# 23. 最重要的 evaporation 判据

高截面区域必须放弃：

$$
E>0
\Rightarrow
\text{evaporated}
$$

这一简单判断。

因为粒子可能在：

$$
r<R_\odot
$$

被热核 kick 到

$$
E>0,
$$

但在真正离开太阳以前再次发生 scattering：

$$
E>0
\rightarrow
E'<0.
$$

因此真正的 absorbing condition 是：

$$
\boxed{
r=R_\odot,\quad
v_r>0,\quad
E\ge0.
}
$$

只有真正穿过太阳表面时才算 evaporation。

这也是为什么我们的状态必须保留位置/轨道 phase，而不能只使用 \(E,L\)。

输运物理逃逸面按 Task_Plan 取 Rsun；旧程序的终止/Kepler 匹配面为 1.1Rsun。T08/T12 用无碰撞段连接二者，并分别保存表面逃逸、旧边界到达和中间壳层占据。内部瞬间解束缚不是边界吸收。

---

# 24. 太阳表面必须存在两类完全不同的 outgoing boundary

对于

$$
r=R_\odot,\quad\mu>0
$$

分成：

## A. Unbound

$$
E\ge0.
$$

进入：

$$
\boxed{\text{ESCAPED}}
$$

absorbing state。

---

## B. Bound

$$
E<0.
$$

粒子只是暂时离开太阳。

在纯太阳势的 `kepler_return` 模型中，它沿 weakly bound orbit 运动并返回，观测半径不作为吸收面。

首个旧数据对照采用 `legacy_1au_removal`：只有远日点在 1 AU 域内的外轨道返回；超出该域者传播到首次向外穿过 1 AU 后记为独立 `OUTER_REMOVAL`。此项是模型中的外域移除，不能混同 evaporation。两种边界模型分别保存、分别验证。

旧实现以 1.1Rsun 匹配；Rsun↔1.1Rsun 的轨道映射必须保持 E/L，并守恒转移相应通量、时间及壳层占据。输出 n_inside 始终指 r<Rsun，不能把匹配面内全部区域称为太阳内部。

---

# 25. 稳态问题有一个非常有用的简化

对选定模型中确实返回的 bound outward crossing，且有限稳态存在时：

$$
E<0
$$

如果我们只关心**稳态太阳内部 distribution**，外部飞行造成的时间 delay 不改变 steady return flux。

因此可以在 interior steady-state solver 中使用：

$$
(R_\odot,v,\mu>0,E<0)
$$

直接映射为

$$
(R_\odot,v,-\mu,E<0).
$$

也就是：

$$
\boxed{
\text{outgoing bound}
\rightarrow
\text{returning bound}.
}
$$

虽然真实粒子需要在外面飞一段时间，但 steady flux 满足：

$$
F_{\rm return}=F_{\rm outward}.
$$

因此内部稳态分布不会因为忽略这段 delay 而改变。

这只构成删去外部等待的内部稳态算符 \(Q_{\rm ss}\)。其矩阵指数、特征值和 survival curve 不能解释为包含太阳外飞行的物理时间。有限年龄需要外轨道相位/飞行段或延迟返回方程；不能仅给每条轨道加一个平均等待率就认为恢复了确定性 Kepler delay。1 AU 移除通道也不能即时回流。

---

# 26. 但太阳外数密度不能忽略这个时间

这反而给了我们一个非常高效的 external-density algorithm。

一个粒子离开太阳时：

$$
E<0.
$$

定义 specific angular momentum

$$
L=R_\odot
v\sqrt{1-\mu^2}.
$$

太阳外：

$$
\Phi(r)=-\frac{GM_\odot}{r}.
$$

径向速度：

$$
\boxed{
v_r(r)=
\sqrt{
2\left(E+\frac{GM_\odot}{r}\right)
-\frac{L^2}{r^2}
}.
}
$$

粒子一次外部 excursion 在 shell \(j\) 中停留：

$$
\boxed{
\tau_j(E,L)
=
2\int_{r_j^-}^{r_j^+}
\frac{dr}{|v_r(r)|}.
}
$$

factor 2 仅适用于完整返回轨道。积分上下限须裁剪到轨道可达区间和该壳层；外域移除只取到首次外行过界的单程，因子为 1。

现有 `Compute_Bound_Kepler_Exterior_Arc` 已计算外部 dt、v²dt 与返回/移除状态；T12 复用其物理思路，参数化匹配面和任意壳边，并补纯径向、近抛物线等极限。现接口固定 1.1Rsun 和现有 histogram，不能原样当成 Rsun 上的通用 API。

---

# 27. 从 surface flux 直接得到绝对太阳外数密度

Transport solver 给出某个 bound surface channel \(a\) 的 outward flux：

$$
\Phi_a
\quad [{\rm particles/s}].
$$

那么外部 shell \(j\) 中的粒子数：

$$
\boxed{
N_j^{\rm out}
=
\sum_a
\Phi_a
\tau_{aj}.
}
$$

因此：

$$
\boxed{
n_\chi^{\rm out}(r_j)
=
\frac{1}{V_j}
\sum_a
\Phi_a\tau_{aj}.
}
$$

这是整个项目最值得实现的公式之一。

它意味着：

> weakly bound particle 不需要真的在计算机里飞到 aphelion 再飞回来。

外部数小时、数天、数月甚至更长的轨道，全部变成一次解析/一维积分。

该乘积公式适用于稳态通量及相同边界策略，所用 tau 必须区分双程返回与单程移除。有限年龄的外部占据须把历史出射通量与飞行/停留核做卷积，不能用当前通量直接乘完整 excursion 时间。

---

# 28. 这直接解决当前项目的 absolute normalization 问题

过去：

$$
N_\chi
=
C\langle\tau_{\rm evap}\rangle
$$

然后：

$$
n(r)=N_\chi g(r).
$$

问题是：

$$
\langle\tau_{\rm evap}\rangle
$$

会受到 heavy tail 强烈影响。

新方法直接求：

$$
\boxed{
-Q^T\mathbf N=\mathbf C.
}
$$

于是：

$$
N_\alpha
$$

已经是 absolute occupation。

太阳内部：

$$
\boxed{
n_\chi(r_i)
=
\frac{
\sum_{jk}N_{ijk}
}{
V_i
}.
}
$$

太阳外：

$$
\boxed{
n_\chi^{\rm out}(r_i)
=
\frac{
\sum_a\Phi_a\tau_{ai}
}{
V_i
}.
}
$$

不再需要先计算 lifetime PDF。

线性求解提供的是给定 source、物理域与离散模型的绝对幅度，不保证重尾积分有限。内部稳态加解析外部占据后才得到 N_total；若 E→0− 尾部导致外部驻留不收敛，必须报告该限制或转入有限年龄模型，不能手工 renormalize。

---

# 29. 稳态 equation

Version 1 的内部稳态采用 §25 的缩约算符：

$$
\boxed{-Q_{\rm ss}^T\mathbf N_{{\rm inside},\rm ss}=\mathbf C.}
$$

这里的相空间包括捕获后尚未最终逃逸的 \(E\ge0\) 内部状态及 recapture；\(\mathbf N_{\rm inside}\) 仍是各 `(r,v,mu)` 单元占据的列向量。外部数目由返回/移除通道的解析驻留重建。

若保留完整物理时间与外部状态，另定义 transient 生成元 \(Q_B\)：

$$
\frac{d\mathbf N}{dt}=\mathbf C+Q_B^T\mathbf N,
\qquad -Q_B^T\mathbf N_{\rm ss}=\mathbf C.
$$

两式分别属于完整模型的时间演化和稳态；不能把 \(Q_{\rm ss}\) 的时间演化解释为全域物理时间。

稳态求解前检查源可达闭合类及奇异性；有源流入无损失闭合类时不存在有限稳态。不能添加任意小漏项或裁剪负解来强行求得答案。

---

# 30. 但是绝不能默认太阳年龄足够长

这是项目科学上非常重要的一点。

如果 slowest mode 的时间尺度：

$$
\tau_{\rm slow}
\sim
\frac{1}{|\operatorname{Re}\lambda_{\rm slow}|}
$$

与

$$
t_\odot
$$

接近，那么：

$$
N_{\rm ss}
$$

不是物理答案。

有限时间解应该是：

$$
\boxed{
\mathbf N(t)
=
\int_0^t
e^{Q_B^T(t-t')}
\mathbf C(t')dt'.
}
$$

如果 capture rate 近似常数：

$$
\boxed{
\mathbf N(t)
=
\int_0^t
e^{Q_B^Ts}\mathbf C\,ds.
}
$$

所以最终程序必须同时支持：

```text
steady_state
```

和

```text
finite_age
```

两种运行模式。

以上指数作用要求 \(Q_B\) 已含外部物理时间，且示式假定 \(N(0)=0\)；非零初态另加 \(e^{Q_B^Tt}N(0)\)。\(\lambda_{\rm slow}\) 是源可激发且与观测量相关的衰减分支，非正规算符还需检查瞬态。G4/T17 完成前，不输出物理有限年龄结论；不能把 \(Q_{\rm ss}\) 代入这些公式。

<a id="return-time-kernel"></a>

## 30.1 T17 的返回时间核与外部占据

设 a 是太阳表面的束缚外行通道，b 是返回内行通道，\(F_a^{out}(t)\) 单位为粒子/秒。定义 \(K_{ba}(u)\,du\) 为离开通道 a 后，在延迟 \([u,u+du]\) 返回 b 的概率；K 的单位为秒的倒数，可含 delta。对固定背景及已验证的通道内出射分布，

\[
F_b^{in}(t)=F_b^{in,0}(t)+\sum_a\int_0^t
K_{ba}(u)F_a^{out}(t-u)\,du.
\]

对完全分辨的 \((E,L)\) Kepler 通道，\(K_{ba}(u)=R_{ba}\delta(u-T_a)\)，其中 R 是守恒返回映射，T 是从出射面到再次入射面的真实外部飞行时间。有限通道是各轨道 delay 的混合，不能仅保存平均 T；近 \(E=0^-\) 的混合还可能有很长尾。

固定卷积核另有闭合条件：单元内按**出射通量**加权的 E/L 分布必须由声明的重构规则确定，并经细分/求积验证。若这一条件分布随出射时刻 s 改变，应保留子通道或使用 \(K_{ba}(u;s)\)，不能仅因背景静态就假定聚合核时间平移不变。

此时以下占据和移除核也必须同步成为 \(H_{ja}(u;s)\)、\(M_a(u;s)\)，三者使用同一个出射时刻条件分布；积分中均取 \(s=t-u\)。归一化和外域生存恒等式对每个固定 s 成立，不能只更新 K 而保持 H/M 不变。

设 \(H_{ja}(u)\) 为出射 u 秒后仍在外域壳 j 中的概率。它是无量纲占据核，完全分辨轨道时是壳内指示函数；它不是返回概率密度 K。于是

\[
N_j^{out}(t)=N_j^{out,0}(t)+\sum_a\int_0^t
H_{ja}(u)F_a^{out}(t-u)\,du,\qquad
\tau_{aj}=\int_0^\infty H_{ja}(u)\,du.
\]

稳态且积分有限时恢复 §27 的 \(N_j^{out}=\sum_a F_a^{out}\tau_{aj}\)。初始外域粒子的相位分布必须同时产生 \(F^{in,0}\)、\(N^{out,0}\) 和初始移除流；若给定出射前史，可一致地等价处理，不能遗漏或重复加入。

`legacy_1au_removal` 另定义首次外行穿过 1 AU 的延迟密度 \(M_a(u)\)，并由它卷积得到实际移除流。对于最终返回或移除的束缚 excursion，\(\sum_b\int K_{ba}du+\int M_a du=1\)，全部外域壳的 H 之和等于尚未返回/移除的概率。若某模型允许永久外域存留，则另保留其质量。不能在太阳表面就从总粒子数中扣除仍在飞往 1 AU 的粒子。

实现的内域开放算符保留所有出射流出，返回源由上式沉积，禁止同时启用即时回流。总账目为

\[
\frac{d}{dt}(N_{inside}+N_{outside})
=C_{total}(t)-\Phi_{escape}(t)-\Phi_{outer\ removal}(t).
\]

表面非束缚逃逸按既定 Rsun 终点计数；Rsun/1.1Rsun 桥接时间与壳层只计一次。T17 对单一确定 delay、双 delay 混合、非零初始外域及延迟移除分别检验早时无提前返回、稳态极限、时间步/通道细化和上述守恒。卷积历史的截断/压缩必须报告未返回质量与误差；不能以有限内存为由删除长飞行。只有扩展成足够准确的 Markov 时间状态后才直接用矩阵指数，延迟模型本身不默认存在原维度的常系数 Q。

---

# 31. 这还能直接回答“是否达到平衡”

计算

$$
\mathbf N(t_\odot)
$$

和

$$
\mathbf N_{\rm ss}.
$$

定义：

$$
\delta_{\rm eq}
=
\frac{
\|\mathbf N(t_\odot)-\mathbf N_{\rm ss}\|
}{
\|\mathbf N_{\rm ss}\|
}.
$$

如果：

$$
\delta_{\rm eq}\ll1,
$$

才可以使用 steady state。

因此 equilibrium 不再是假设，而成为一个输出。

这也比传统单一：

$$
\tau_{\rm evap}
$$

判据更强。

此比较要求有限稳态存在，并对有限年龄解和稳态解采用相同物理域、时间定义及观测量。若稳态不存在，标记 equilibrium_error 不适用，不构造伪分母。

---

# 32. Lifetime distribution 以后也可以从矩阵反推出

如果 freshly captured distribution 为：

$$
\mathbf p_0,
$$

则 survival probability：

$$
\boxed{
S(t)
=
\mathbf 1^T
e^{Q_B^Tt}
\mathbf p_0.
}
$$

平均 residence time：

$$
\boxed{
\langle\tau\rangle
=
\mathbf 1^T
(-Q_B^T)^{-1}
\mathbf p_0.
}
$$

因此你现在 trajectory MC 中看到的：

* long tail；
* multiple timescales；
* metastable population；

以后都可以解释成 \(Q\) 的 slow eigenmodes。

这会成为项目一个很有价值的物理结果。

这些公式只对包含外部飞行时间、与目标终点一致的 \(Q_B\) 成立，且平均值公式要求有限均值。对 \(Q_{\rm ss}\) 算出的量仅是压缩外部等待后的过程。`legacy_1au_removal` 的 survival 终点包括外域移除，须与 evaporation 的分通道 first-passage 量区分。MC 的最终解束缚时刻、表面逃逸时刻及旧 1.1Rsun 终止时刻也需分别保存。

---

# 33. Heavy tail 在新框架里的物理解释

如果

$$
Q_B
$$

存在若干非常接近零且可被 source 激发的衰减模；若相关特征值为实数，可记为：

$$
0>\lambda_1>\lambda_2>\cdots,
$$

在有限维算符可对角化时，可将 survival 写为模态展开：

$$
S(t)
=
\sum_i
A_i e^{\lambda_i t}.
$$

最慢的 mode：

$$
\tau_i
=
-\lambda_i^{-1}
$$

就会产生 lifetime long tail。

如果存在多个 slow modes，就自然产生 multi-timescale evaporation。

谱分析是理解尾部的一种工具。不可对角化时可能出现带多项式因子的项；扩展网格/外域的极限可能有更复杂的尾部，必须另检谱与驻留积分收敛。

用于物理寿命解释的必须是保留外部时间的算符/延迟模型；缩约 \(Q_{\rm ss}\) 的慢模不能替代太阳外长飞行的物理尾部。有限离散系统的有限慢模也不能证明连续无限外域具有有限稳态。

一般实生成元也可能含复谱，衰减包络由 \(\operatorname{Re}\lambda_i<0\) 决定，此时不能直接用实数排序或 \(-1/\lambda_i\) 作为物理寿命。

T20 可选 QSD/Fleming–Viot 诊断考察 \(P(X_t\in dx\mid\tau_{kill}>t)\) 的晚期形状；kill 终点须按边界策略区分蒸发和外域移除，延迟模型须包含轨道年龄/相位。有限状态与大粒子数下的收敛结论，不能直接推广到无限 Kepler 外域。要检查 walker 数、谱/混合时间、域与网格及独立重复。未看到稳定条件分布可能来自混合慢或统计不足，不能据此判定连续重尾；得到 QSD 也不替代给定捕获源的绝对占据及完整寿命分布。[Asselah、Ferrari 与 Groisman](https://arxiv.org/abs/0904.3039)

---

# 34. 一个必须正视的问题：高截面不仅造成 stiffness

如果：

$$
\lambda_{\rm mfp}
=
\frac{v}{\Gamma_{\rm scat}}
$$

变得远小于 radial grid：

$$
\lambda_{\rm mfp}
\ll
\Delta r,
$$

简单的普通 upwind discretization 可能出现很强的 numerical diffusion。

即：

$$
D_{\rm numerical}
$$

可能大于真实：

$$
D_{\rm physical}.
$$

所以：

$$
\boxed{
\text{“做一个大 Markov matrix”}
}
$$

本身并不足以保证高截面结果正确。

这里 \(v/\Gamma\) 是局域 flight-length 诊断。对有靶热运动或非平凡角分布的情况，它不自动等于控制扩散的 transport mean free path；扩散系数须从已验证的 collision operator 推导/测量。

---

# 35. 因此整个项目应该分成 kinetic 和 opaque 两种 regime

定义局域 Knudsen-like parameter：

$$
\boxed{
K(r,v)
=
\frac{\lambda_{\rm mfp}(r,v)}
{H(r)}
}
$$

其中可以取：

$$
H^{-1}
=
\max\left(
\left|\frac{d\ln n}{dr}\right|,
\left|\frac{d\ln T}{dr}\right|,
\frac1r
\right).
$$

粗略分类：

$$
K\gtrsim0.1:
\quad
\text{kinetic regime},
$$

$$
10^{-2}\lesssim K\lesssim0.1:
\quad
\text{transition},
$$

$$
K\ll10^{-2}:
\quad
\text{diffusive/opaque}.
$$

这些数字只是 numerical domain decomposition 的初始值，不是固定物理界线。

\(H\) 的具体估计需记录背景量和中心处理；几何项 \(1/r\) 在中心的坐标奇点不能单独作为物理 opaque 判据。最终区域划分依据算符输运尺度、误差与收敛证据确定。

---

# 36. Version 1 不需要马上解决 asymptotic-preserving 问题

第一版应该优先在：

$$
\lambda_{\rm mfp}
\gtrsim
\Delta r
$$

或者不是极端 opaque 的 benchmark 上验证 kinetic operator。

目标：

$$
\boxed{
\text{transport result}
=
\text{trajectory result}.
}
$$

这是第一篇代码验证图最重要的内容。

---

# 37. Version 2 才进入真正的 \(10^{-30}\) production regime

T18 先在共同小模型中比较三个候选的适用范围和成本，选定生产主线及可用的交叉验证路线；不要求把三套系统都完整实现。

## 路线 A：Asymptotic-preserving kinetic scheme

保留完整

$$
(r,v,\mu)
$$

Boltzmann equation，但使用能够在

$$
\lambda_{\rm mfp}\rightarrow0
$$

时自动退化为 diffusion limit 的 numerical scheme。

优势：

* 一个统一方程；
* 光学薄厚平滑连接；
* 可在同一模型内检验扩散极限。

缺点：

* 开发复杂。

---

## 路线 B：Hybrid kinetic + diffusion

可作为 production 候选；与 AP 路线的取舍由局域输运极限、接口误差及 profile 决定。

定义：

```text
deep opaque region
        ↓
diffusion/local transport
        ↓
transition layer
        ↓
full kinetic solver
        ↓
surface evaporation layer
```

即：

$$
r<r_{\rm switch}
$$

采用 diffusion approximation。

而在 evaporation 最敏感的外层：

$$
r\gtrsim r_{\rm switch}
$$

继续完整 kinetic treatment。

因为真正决定 evaporation 的往往是 near-escape high-energy tail 和太阳外层 transport，不必在核心把每一个 microscopic mean free path 都解析出来。

Hybrid 的模型约化既可减少自由度，也可改善高 opacity 的求解成本；但 matching 需保持粒子/通量、适当能量矩和关键尾部，并用 overlap 与渐近极限验证。性能收益待实测，不以强制热化代替物理推导。

## 路线 C：轨迹 MC + diffusion random walk / DDMC

在满足受控扩散条件的厚区，用空间离散、连续时间的扩散跳跃替代大量微观散射；薄区和敏感的逃逸层继续原轨迹 MC。该思路和渐近 diffusion/transport 接口在 [Densmore 等，JCP 222 (2007) 485–503](https://doi.org/10.1016/j.jcp.2006.07.031) 的辐射问题中已有实现，但不能把光子替换成 DM 后直接沿用公式。

本项目必须从共享碰撞物理推导重力漂移、速度/能量交换、局域平衡零模和必要的慢变量。角向松弛快不保证能量已局域平衡；无法闭合到纯空间扩散时，保留能量维度或继续 kinetic。接口须给出正确入射分布、逃逸/返回概率、物理计时及能量矩，不只匹配总流量。

DDMC 的连续时钟描述扩散近似的时间过程，不自动等于原微观 first-passage law，也不替代 T17 的 Kepler 延迟。小模型通过后，将 DDMC 与普通 MC、kinetic 及已验证且可用的 WE 比较驻留、分通道首达、绝对密度和尾部，扫描切换位置/准则。其数值实现可独立于确定性 solver，但若二者共享扩散闭合或微观模块，相符不能证明这些共同假设正确。

T18 三路线共同门槛是在固定、能解析宏观背景的网格上，散射增强时趋向推导出的正确扩散解，而非要求 \(\Delta r\) 始终追随微观自由程。太阳表面和 kinetic/diffusion 接口的边界层、非负性与近逃逸尾另行收敛。若某路线在预设误差预算内没有成本收益，记录结论并停止扩展该候选，不阻塞其余路线。

---

# 38. 不要一开始直接使用传统 thermal distribution 代替 opaque region

Hybrid solver 中可以利用 local equilibrium，但必须把它作为：

$$
K\rightarrow0
$$

的 controlled limit。

不能简单预设：

$$
f_\chi
\propto
e^{-m_\chi\Phi/T_\chi}.
$$

否则我们又重新引入了项目本来就是为了检验的 thermalization assumption。

正确逻辑是：

$$
\text{kinetic operator}
\xrightarrow[K\rightarrow0]{}
\text{diffusive/local-equilibrium limit}.
$$

而不是人为强行 thermalize。

---

# 39. Solver 的第一版技术选择

Prototype 建议：

### Operator generation

C++。

原因：

* 与项目内移植的 C++ 微观闭包直接连接；
* collision sampling 高效；
* 后续可在本项目内增加受控并行。

### Sparse solve

第一版：

```text
SciPy sparse
```

足够。

输出：

```text
Matrix Market / NPZ
```

例如：

```python
from scipy.sparse.linalg import spsolve, gmres
```

这一层定义为 CPU、double precision、明确小矩阵的 reference backend，优先用可解释的简单实现核验物理与矩阵方向。Matrix Market 仅用于小测试/调试，生产缓存采用有 schema 的二进制格式；不把 GPU 或 matrix-free 作为 T02/T03 的前提。

---

# 40. Production solver

状态数达到：

$$
10^5-10^6
$$

以后建议切换：

```text
PETSc
```

使用：

* GMRES；
* BiCGSTAB；
* sparse LU 作为 benchmark；
* block preconditioner；
* later diffusion synthetic acceleration。

不要显式构造：

$$
Q^{-1}.
$$

永远只做：

$$
A x=b.
$$

optimized backend 必须与 reference 的物理、网格、source、sink 和输出契约一致。可将 streaming / boundary 稀疏项与按 r 分块的 collision 分开存储，采用 matrix-free 前向作用：

$$
y=Q^Tx=Q_{\rm stream}^Tx+
\frac{\sigma}{\sigma_0}(Q_{\rm coll}^{(0)})^Tx+
Q_{\rm boundary}^Tx.
$$

小矩阵 reference 用于检查该作用与转置，不能把附件中泛指的 \(Qx\) 混入本稿列向量前向方程。先测量装配、预处理、迭代数及内存，再选择块预处理、CPU 并行或 GPU。高 opacity 的算法与模型约化仍须通过原 G5 科学门槛。详见 [工程与性能设计附录](Proposal.md#performance-design)。

---

# 41. 一个非常重要的 matrix sanity check

完整的碰撞概率账目应满足代数守恒，非对角转移率非负：

$$
\sum_\beta Q^{\rm coll}_{\alpha\beta}=0,
\qquad Q^{\rm coll}_{\alpha\beta}\ge0\quad(\alpha\ne\beta).
$$

对 transient 子矩阵，明确的物理 sink 则满足：

$$
\sum_\beta Q_{\alpha\beta}
=-q_{\alpha\to\mathrm{escape}}
-q_{\alpha\to\mathrm{outer\ removal}}.
$$

`kepler_return` 关闭第二项。数值 velocity overflow 必须单独记录并扩域/验证闭合，不能充作上述物理 sink。严格守恒不是说输出概率已物理准确：MC 误差影响 P 的准确度，归一化及行和仍应到浮点精度。容差按行速率尺度无量纲化，并另检残差、负占据及稀有通道误差。

---

# 42. 稳态还有一个极强的 global test

在无湮灭且有限稳态存在时，所有源流入与物理 sink 流出必须闭合：

$$
C_{\rm tot}=\Phi_{\rm evap}+\Phi_{\rm outer},
\qquad
\Phi_s=\sum_\alpha N_\alpha q_{\alpha\to s}.
$$

`legacy_1au_removal` 使用两种 sink；`kepler_return` 且所有粒子最终蒸发时才化为 \(\Phi_{\rm evap}=C_{\rm tot}\)。不能把两种外域策略混在同一个基准中。

账目失配需要排查 source/sink 口径、矩阵、边界和求解误差，以及未识别的无损失闭合类。反过来，即使闭合到机器精度，也不能单独证明蒸发尾部、驻留时间或绝对密度正确；这些仍需同源 MC 与网格/核采样收敛验证。

---

# 43. 第二个极强测试：Little's law

项目内 trajectory port 或已保存基线能够给出平均 residence time benchmark 时：

$$
\boxed{
N_{\rm tot}
=
C\langle\tau\rangle.
}
$$

因此：

$$
\boxed{
\frac{N_{\rm transport}}C
=
\langle\tau\rangle_{\rm trajectory}.
}
$$

这是验证 absolute normalization 最好的测试之一。

不是只比较 normalized \(g(r)\)。

必须比较相同域及终点：\(N_{\rm inside}/C=\langle t_{\rm inside}\rangle\)，只有计入外部驻留后，\(N_{\rm total}/C\) 才对应该域内的完整寿命。旧 `bincount.txt` 是捕获后累计 dt 总和；`evaporation_times.txt` 只包含完整蒸发子样本，且主要寿命列止于最终解束缚散射，不能直接用于边界驻留的均值。

wall/step/scattering cutoff 保留的 residence prefix 不是完整寿命；数值失败的排除和完整蒸发筛选也可能带来偏差。基准优先无计算删失，否则报告验证不足或比较共同物理观测时间的有限时占据。沿用块级 jackknife 传播 capture 和 residence 的相关误差，不做事后幅度匹配。

---

# 44. 第三个测试：外部 density

选择一个已有完整固定输出或可由项目内 trajectory port 完成的 benchmark。

项目内或已保存的 Trajectory MC 统计：

$$
n_{\rm out}^{\rm MC}(r).
$$

Transport solver 使用：

$$
n_{\rm out}^{\rm op}(r)
=
\frac1{V_r}
\sum_a\Phi_a\tau_a(r).
$$

要求二者一致。

此测试要求同 source、同外域策略和公共壳边下的绝对幅度及误差一致。现有 full MC 已使用解析外轨道积分，因此还需独立 quadrature 或解析可核验轨道测试；两个客户端调用同一 helper 后相互一致，不能单独证明 helper 正确。

---

# 45. Collision-only thermalization test

固定 \(r=r_0\)，关闭 gravity/streaming，只保留温度 \(T(r_0)\) 的碰撞热浴。首版应比较积分后的 Maxwell–Boltzmann 单元占据：

$$
N_{{\rm eq},\alpha}
=\int_{\mathrm{cell}\ \alpha}Jf_{\rm MB}\,dr\,dv\,d\mu,
\qquad Q_{\rm coll}^T\mathbf N_{\rm eq}\simeq0.
$$

这里的 `Q_coll` 是 T05 构造的有限体积离散生成元，属于 T06 验收；T03 先用 §8 的弱式检验尚未分箱的完整单碰撞过程。两项不能互相代替。

不能把节点上的裸 \(f_{\rm MB}\) 代入本稿占据数方程。随速度/角度网格和核样本数提高，检查速度分布、角各向同性、能量矩及稳态残差收敛；速度域尾部不能静默截断。

详细平衡 \(N_{{\rm eq},\alpha}Q_{\alpha\beta}=N_{{\rm eq},\beta}Q_{\beta\alpha}\) 是额外测试，不能由 \(Q^TN_{\rm eq}\simeq0\) 推出。随机采样、有限体积投影和物理适用范围导致的偏差分别记录。

---

# 46. 必须建立如下 test hierarchy

## Unit tests

```text
solar_model_test
scattering_rate_test
single_collision_test
kernel_normalization_test
collision_conservation_test
streaming_conservation_test
surface_boundary_test
cross_section_scaling_test
```

## Physics integration tests

```text
thermal_bath_test
ballistic_limit_test
weak_scattering_test
trajectory_transport_test
external_orbit_test
steady_flux_test
finite_age_test
```

## Numerical convergence tests

```text
Nr convergence
Nv convergence
Nmu convergence
kernel sample convergence
vmax convergence
solver tolerance convergence
```

---

# 47. 推荐的新仓库结构

文档固定为四份；其职责与维护规则见 README。Code/include、src、tests、configs 及 Python 基线工具已开始落地；下图仍包含尚待实现的求解/应用模块，不能由目录示意推定能力已完成。

```text
DM-Transport/
├── README.md                 # 项目入口、运行说明与文档维护规则
├── CHANGELOG.md              # 简短修订记录
├── CMakeLists.txt
├── Note/
│   ├── Proposal.md           # 科学路线、工程设计与源码依据
│   ├── Task_Plan.md          # 任务、依赖、状态和验收
│   └── archive/              # 已有 v1 原稿冻结，不新增逐轮备份
├── Code/
│   ├── include/transport/    # 网格、碰撞核、传播/边界、源、外域与观测接口
│   ├── src/                  # C++ 算符及源构建
│   ├── apps/                 # build_kernel、build_source、benchmark_mc
│   ├── python/               # 稳态/有限年龄求解、收敛、图表
│   ├── tests/
│   └── configs/benchmark/
└── Output/
    ├── Result/               # 运行数据、JSON 验证/性能报告与 manifest
    └── Figure/
```

DaMaSCUS-SUN-EVAP 的固定提交只提供源码事实和已保存基线；项目内 `SolarBackground`、`ScatteringPhysics` 记录逐文件来源及有意差异。reference 与 optimized 在 DM-Transport 内共享这一移植物理层、schema 和验证输入，不为加速再复制一套散射公式。后端细节见附录 A；具体执行任务以 Task_Plan 为准。

---

# 48. 推荐的数据结构

不要让 transport code 到处使用三个 index。

定义：

```cpp
struct PhaseSpaceIndex {
    int ir;
    int iv;
    int imu;
};
```

提供：

```cpp
size_t flatten(ir, iv, imu);
PhaseSpaceIndex unflatten(size_t id);
```

然后所有 sparse matrix 统一使用：

```text
global_state_id
```

这会显著降低 bug 数量。

---

# 49. Collision kernel 文件应该独立缓存

例如：

```text
kernels/
  mchi_0p10GeV_SD/
      grid.h5
      rate_unit.h5
      transition_kernel.h5
      metadata.json
```

metadata 必须记录：

```text
DM mass
interaction
reference sigma
solar model
Nr Nv Nmu
vmax
kernel samples
random seed
git commit
date
```

缓存键还必须含耦合、靶列表、太阳表 SHA256、网格面/速度域、物理实现版本、抽样/沉积规则及单位。加载时核验全部签名并拒绝错配。仅修改不改变这些条件的 solver 才可复用核；截面复用须先满足 §19 的验证条件。

---

# 50. RNG 必须从一开始就保证 reproducibility

不要让结果依赖 MPI worker 数。

建议每个 state 使用：

$$
\text{seed}
=
H(
\mathrm{master\ seed},
\mathrm{model/grid\ fingerprint},
m_\chi,
state\_id,
sample\_id
).
$$

即使：

```text
MPI -n 20
```

换成：

```text
MPI -n 200
```

也应该得到同一 kernel。

这对于以后论文复现非常重要。

哈希或 counter-based 映射应规范、跨进程稳定，不能使用实现相关的 `std::hash` 或 MPI rank 决定逻辑样本。并行浮点归并允许声明的舍入差异；要区分同一随机样本集合、整数计数相同和最终浮点逐位相同。T02 项目内移植先保留与旧 trajectory 相容的随机调用顺序，避免把 RNG 更换混入物理回归。

---

# 51. 推荐的开发阶段

## Milestone 0 — Freeze the old simulator reference

目标：

> 固定只读 DaMaSCUS 来源和可复现的旧结果，不再改动参考仓库。

完成：

* 建立 baseline benchmark；
* 保存固定 seed 输出；
* 保存 scattering rate；
* 保存 reflection/capture ratio；
* 保存若干 trajectory。

Gate：

$$
\text{saved legacy reference}
=
\text{T01 manifest and outputs}.
$$

比较对象必须来自 T01 已保存的固定输出；不为生成新一轮“before”结果而修改或编译参考仓库。当前通过状态与尚缺的验证证据只在 [Task_Plan §8](Task_Plan.md#current-status) 维护。

以下 Milestone 保留为主题路线，不是实际顺序；实施依赖、近期 T00–T03 范围以及 G0–G6 科学验收以 Task_Plan 为准。尤其参数化 exterior 是边界与人工源闭环的前置工作；有限年龄 G4 和高 opacity G5 不因性能任务而交换含义。

---

# 52. Milestone 1 — Port microscopic physics

T02 在 DM-Transport 内实现：

```text
SolarBackground（最小太阳表接口）
ScatteringPhysics（串行 rate 与单碰撞）
Provenance manifest（来源位置、许可证、有意差异）
```

旧 trajectory simulator 保持只读；项目内 physics target 不带入 MPI/legacy 调度，并提供 T03 系统回归所需的固定输入接口。

这一阶段**不写 Markov solver**。

Gate：

* 项目内 consumer 只链接 physics target 即可查询背景、直接 rate 和单碰撞；
* 来源位置、许可证及有意差异记录完整；
* 零速极限、稳定分支和固定数值检查通过。rate、单碰撞分布与固定种子的系统 parity 由 T03 验收。

---

# 53. Milestone 2 — Collision kernel prototype

实现：

```text
PhaseSpaceGrid
CollisionKernelBuilder
```

固定：

$$
m_\chi=0.1\ {\rm GeV}
$$

先只考虑：

$$
\text{SD proton}.
$$

使用：

$$
30\times40\times8
$$

小网格。

输出：

$$
\Gamma_\alpha
$$

和：

$$
P_{\alpha\beta}.
$$

Gate：

$$
\sum_\beta P_{\alpha\beta}=1
$$

且：

$$
Q_{\rm coll}^T \mathbf N_{\rm thermal}
\simeq0
$$

在固定 thermal bath test 中成立；N_thermal 是 §45 的单元积分占据。overflow 单独追踪、扩域收敛，P 行和不能靠静默丢弃后归一化实现。

reference 先直接按 \((r,v,\mu)\) 验证。经 T06 确认局域旋转对称后，optimized builder 可尝试复用 \(K_{\rm loc}(r,v\to v',\cos\psi)\)，用解析方位角沉积映射到 \(\mu'\)；必须核验与同一有限体积质量测度的 reference 核一致。该结构优化的适用条件、误差和退化极限见 附录 A（工程与性能设计），不预设加速倍率。

---

# 54. Milestone 3 — Add gravitational streaming

实现：

```text
StreamingOperator
BoundaryOperator
```

构建：

$$
Q=
Q_{\rm stream}
+
Q_{\rm coll}.
$$

此阶段 source 可以人工设置。无碰撞不变量和假蒸发必须先通过 §12.1 的 T07 门槛。

例如在某个 state 注入：

$$
C_\alpha=1.
$$

然后观察：

$$
N(r,v,\mu).
$$

Gate：

* number conservation；
* ballistic test；
* reflective boundary test；
* escape test。

边界验收包含 bound-return、unbound escape 和配置选择的 1 AU removal；Rsun/1.1Rsun 桥接须先验证。以上 Q 可把固定 boundary 并入 streaming，实现仍分别记录各分量及 sink。

---

# 55. Milestone 4 — Real capture source

复用已存在的 Capture mode 停止链：

```text
halo entry
→ scattering
→ first E<0 state
→ stop
```

新增首捕获状态、入射物理率和完整计数输出，经 §21 归一化产生：

$$
C_\alpha.
$$

然后第一次真正求：

$$
-Q^TN=C.
$$

输出：

$$
n_\chi(r).
$$

这时项目第一次真正实现：

$$
\boxed{
\text{absolute density without complete evaporation trajectories}.
}
$$

验收 sum(C_alpha)=C_total，保存源投影误差、raw/valid 概率及 unresolved/failed 比例；配套已捕获初态 MC 入口。已有 stop 逻辑不是待从零实现的交付物。

---

# 56. Milestone 5 — Benchmark against frozen and project-local trajectories

至少选择三档：

### A. Optically thin

例如：

$$
\sigma\sim10^{-38}
$$

### B. Intermediate

例如：

$$
\sigma\sim10^{-35}
$$

### C. Highest cross section where project-local full MC remains feasible

例如：

$$
10^{-33}\sim10^{-32}
$$

具体根据 DM-Transport 内实现的运行速度决定。

比较：

$$
n(r),
$$

$$
f(v),
$$

$$
N_{\rm tot},
$$

$$
\Phi_{\rm evap},
$$

$$
n_{\rm out}(r).
$$

候选截面仅作覆盖示例，实际首点由 T01 计时与完整分类/轨迹比例选择。对照用同一 source 和边界，新增时间加权相空间/speed 统计；当前径向 dt 与 v²dt 不能重建完整速度分布。报告数值失败、计算删失及外域移除，不能只拿完整蒸发子样本代表所有捕获粒子。

<a id="rare-event-validation"></a>

## 56.1 T11/T15 的可选稀有事件轨迹分支

当普通 MC 的稀有首达/尾部采样不足时，先尝试 Weighted Ensemble（WE）小原型。按 \(E/E_0\) 或 \((r,E)\) 等进度坐标分组，在固定物理时间间隔做带权 splitting/merging，底层调用 DM-Transport 内已验证的轨迹实现。分箱只分配计算资源，不代替完整动力学状态。正确重采样保持条件期望中的加权路径分布；原理见 [Zhang、Zuckerman 与 Jasnow](https://arxiv.org/abs/0810.1963)。

T11 先准备可完整续跑的状态：坐标/速度、真实时间、捕获和终点标记、必要的散射/外域相位状态、统计锚点与随机流。分裂后子路径权重和等于父权重，未来 RNG 独立；合并按权重随机选保留路径并承接总权重，不能平均成一条不存在的轨迹。

T15 以同一 source、终点和物理时窗比较普通 MC 与 WE 的加权驻留、有限时间首达概率、分通道通量和源汇；与 transport 先比较通过 V1 契约的稳态占据。保存权重、谱系、重采样策略及删失，按独立 ensemble 重复估计区间；克隆粒子不是独立样本，普通权重 ESS 也不足以刻画谱系相关。对 walker 数、分箱、重采样间隔和种子做收敛。

实施时分开时间能力：T15 可先对照普通 MC 与 WE 的有限时间首达，并与 transport 比较已通过 V1 契约的稳态占据；任何涉及真实外域时间的 transport 瞬态/首达分布对照都须等待 T17，不能用即时回流 Q_ss 提前完成。

若使用吸收后重新注入的稳态 flux/MFPT 估计，重新注入必须匹配指定捕获源，并验证稳态及率比估计的有限样本误差；不能把 QSD 初态或一个 flux 倒数普遍当成原捕获寿命。WE 降低稀有事件采样困难，不自动减少单条路径的散射成本。先做有停止条件的 pilot；未获得已验证收益时保留普通 MC 或相同有限时窗对照，不把 WE 设为 V1 必选项。宏态 Milestoning 留 T20，遵守 §3 的入口与记忆约束。

---

# 57. Milestone 6 — External orbit module

基于已有解析外轨道能力，参数化实现：

```text
ExteriorOrbit
```

输入：

$$
E,L
$$

输出：

```text
r_peri
r_apo
outside_time
shell_residence_time[]
```

然后将 bound outgoing surface flux 转成：

$$
n_{\rm out}(r).
$$

输入还包括匹配半径、任意壳边及 outer policy，输出分离完整 period、实际 excursion time、返回/移除终态与单/双程驻留。先完成独立轨道和极限验证，再将结果用于外部密度；湮灭后处理属于扩展研究，不属于本模块通过就自动验证的结论。

---

# 58. Milestone 7 — High-opacity solver

这才真正进入：

$$
\sigma\sim10^{-31}-10^{-30}\ {\rm cm^2}.
$$

首先画：

$$
K(r)
=
\frac{\lambda_{\rm mfp}}H.
$$

明确哪些位置：

```text
kinetic
transition
diffusive
```

然后决定：

```text
AP kinetic solver
```

还是：

```text
kinetic + diffusion hybrid
```

并按 §37 评估 `trajectory + DDMC` 的有条件小原型。选择依据包括物理极限、接口/尾部误差和端到端成本，不要求三路线全部生产化。

没有这一步之前，不应该把简单 coarse-grid Markov 结果当作最终高截面物理结果。

---

# 59. Milestone 8 — Finite solar age

加入：

$$
N(t_\odot).
$$

比较：

$$
N(t_\odot)
$$

与：

$$
N_{\rm ss}.
$$

计算 slow eigenvalues。

这一步会直接回答：

> 低质量、高截面 DM 是否真的达到 evaporation/capture steady state？

G4/T17 必须先恢复外部飞行时间。可用外轨道相位/飞行段或延迟返回模型；只对瞬时回流 Q_ss 求矩阵指数或慢特征值不足以完成本里程碑。有限年龄外部数目由历史出射与停留核重建；非稳态总账目需包含域内粒子增长率。

---

# 60. 最后才加入 annihilation

不要在第一版 transport solver 中加入 annihilation。

第一阶段保持：

$$
\text{linear problem}.
$$

完成并验证：

$$
n_\chi(r)
$$

以后，再计算：

$$
\boxed{
\Gamma_A
=
\frac12
\langle\sigma_Av\rangle
\int
n_\chi^2(r)dV
}
$$

对于 s-wave。

如果以后需要让 annihilation 反过来改变 distribution，再加入 nonlinear sink：

$$
\frac{dN_\alpha}{dt}
=
C_\alpha
+
(Q^TN)_\alpha
-
A_\alpha[N].
$$

这是独立的 Version 3。

---

# 61. Version 1 最终应该输出什么

每一个：

$$
(m_\chi,\sigma_p)
$$

必须产生：

```text
phase_space_distribution.h5
radial_density.txt
velocity_distribution.h5
surface_bound_flux.h5
surface_escape_flux.h5
external_density.txt
summary.json
```

`summary.json` 至少包括：

```text
mchi
sigma
capture_rate_particles_s
capture_probability_raw_valid
source_attempted_classified_unresolved
outer_policy
physical_escape_radius
legacy_matching_radius
N_inside
N_outside
N_total
evaporation_flux
outer_removal_flux
mean_inside_residence_time
mean_total_residence_time
physical_time_model
finite_age_validation_status
solver_residual
grid_info
kernel_info
```

同时保存单位、各项误差、数值 overflow/失败/删失计数及归一化定义。V1/G3 不强制提供含外部时间的 slowest_timescale 或 equilibrium_error；G4 验证后再输出这些物理时间量，缺失时标记 not_evaluated，不能以 Q_ss 结果填充。

---

# 62. 最重要的结果图

项目完成以后，优先形成以下图。

### Figure 1

$$
n_\chi(r)
$$

比较：

```text
full trajectory
transport operator
traditional thermal theory
```

---

### Figure 2

计算时间：

$$
T_{\rm CPU}
$$

vs

$$
\sigma.
$$

Full MC 成本受散射数、轨道长度和截断策略影响；transport 受核/source、条件数与迭代数影响，都应实测。分别展示 full MC、首次 transport solve、缓存后追加截面点，记录核、source、装配、预处理、求解与外域重建成本，不预设随 sigma 的增长率或加速倍率。

---

### Figure 3

$$
n_\chi^{\rm out}(r)
$$

比较：

```text
trajectory residence histogram
surface-flux occupation reconstruction
```

---

### Figure 4

$$
S(t)
$$

以及含外部物理时间的算符/延迟模型的 slow modes；此图在 G4 验证后制作。

解释 heavy tail。

---

### Figure 5

$$
K(r)
$$

展示：

```text
optically thin
transition
optically thick
```

区域。

---

### Figure 6

二维：

$$
(m_\chi,\sigma)
$$

map：

```text
thermalized
nonthermal
evaporation dominated
capture–evaporation equilibrium
finite-age nonequilibrium
optically thick
```

这最终会非常接近一个真正的 **solar-DM transport regime map**。

---

# 63. 后续候选：Legendre 角矩表示

T20 可研究

\[
f(r,v,\mu)=\sum_{\ell=0}^{\ell_{max}}f_\ell(r,v)P_\ell(\mu).
\]

在微观核已经验证旋转协变的前提下，方位角平均满足

\[
\langle P_\ell(\mu')\rangle_\phi=P_\ell(\mu)P_\ell(c),\qquad
K_\ell(r;v'|v)=\int_{-1}^{1}K_{loc}(r,v;v',c)P_\ell(c)\,dc.
\]

这里 K_loc 沿用附录 A.3 的 `dv' dc` 密度约定；上式定义条件转移的角矩，装配 f 的前向算符仍须包含输入/输出的速度体积权重和事件率。碰撞在不同 l 之间分块对角，但每个 l 仍耦合 v 与 v'，不能在有能量交换时简化为单个衰减率。角模分解不要求解已各向同性；只有低阶截断的有效性依赖实际角结构。

streaming、引力及半空间出射/返回边界仍耦合角模。厚区 bulk 接近各向同性不保证高速逃逸群体或边界层也如此，必须对 l_max、半空间通量和尾部单独收敛。

角矩系数不是非负单元占据，截断的重构 f 可能为负，不能继续对整个 moment 矩阵强加 FV 生成元的逐项非负判据。FV/moment 转换保留质量权重；正性 filter/limiter 是需要验证的新离散，不允许事后截零掩盖误差。候选与 FV reference 比较守恒、热浴、尾部、扩散极限及相同误差下的存储和求解成本；不进入 MVP。

---

# 64. 不变量坐标作为 T07 的条件回退

首版从 \((r,v,\mu)\) 开始；若 §12.1 的假蒸发误差无法以可接受成本达标，应在 T07 提前评估 \((r,E,L,s_r)\)，而不是等验证成功后才允许换坐标。其中 \(s_r=\operatorname{sign}(v_r)\)，无碰撞时 \(\dot E=\dot L=0\)，保留 r 和径向分支即保留轨道相位，不能退化成只存 E/L 的轨道平均。

该坐标的可达域受 \(2[E-\Phi(r)]-L^2/r^2\ge0\) 约束。需要推导正确 Jacobian、体积和面通量，处理径向转向点的分支交换、中心/径向退化、阈值切割与碰撞后的守恒重映射；不能仅因不变量形式简单就宣称 streaming 已精确。

先做与普通 FV 相同的 collisionless 小基准，比较不变量漂移、假逃逸、正性及成本。若需要更换主网格，同步调整 T04/T05/T08/T10 的几何、核、边界和源投影并重跑验收；T20 只保留进一步推广，不承担补救 V1 已失败门槛的责任。

---

# 65. 不要做的事情

以下路径应明确禁止。

### 不要 1

继续把：

```text
maximum_scatterings
```

从

$$
10^6
$$

调到：

$$
10^8
$$

然后靠 HPC 硬跑。

---

### 不要 2

用 median lifetime 代替：

$$
\langle\tau\rangle.
$$

median 可以描述典型粒子，但不能给 absolute occupation。

---

### 不要 3

简单截断 lifetime tail，然后手工 renormalize。

---

### 不要 4

把：

$$
E>0
$$

直接定义成 evaporation。

高截面区域这是错误的。

---

### 不要 5

直接使用旧的：

$$
(E,L)
$$

orbit-averaged Markov model。

它依赖 long-mean-free-path assumption。

---

### 不要 6

做一个 coarse：

$$
(r,v,\mu)
$$

网格以后，不做 mean-free-path / grid-size check 就把 \(10^{-30}\) 结果当作可靠结果。

---

### 不要 7

一开始同时加入：

```text
capture
thermalization
evaporation
annihilation
gamma signal
parameter scan
```

这样几乎一定无法定位错误。

---

# 66. 第一版的最小科学目标

MVP 不需要解决整个太阳 DM 问题。

只回答一个问题：

> **给定一组 freshly captured DM source，在忽略 annihilation 时，transport operator 是否能够在不完整追踪 evaporation trajectories 的情况下，重现 full trajectory Monte Carlo 得到的绝对 phase-space occupation？**

如果答案是 yes，项目就已经成立。

这里“相同”指同一捕获源、物理域、吸收终点及驻留定义下的绝对幅度与分布在预定误差预算内一致。该目标不包含尚未验证的有限年龄、无限外域重尾或最高截面生产结论。

---

# 67. MVP 推荐 benchmark

优先：

$$
m_\chi=0.1\ {\rm GeV}
$$

SD proton。

选择一个完整 trajectory 尚可承受的截面，例如：

$$
10^{-36}
$$

或：

$$
10^{-34}\ {\rm cm^2}.
$$

具体由 T01 核验已有数据、计时、完整事件比例和失败/删失后选择；不能仅因目录里存在数据就认定该点可作 clean benchmark。

先不要碰：

$$
10^{-30}.
$$

先证明：

$$
\boxed{
\text{method is correct}.
}
$$

之后再证明：

$$
\boxed{
\text{method extends to the inaccessible regime}.
}
$$

---

# 68. 第一张 validation 图应该极其简单

横轴：

$$
r/R_\odot.
$$

纵轴：

$$
n_\chi(r)
$$

absolute units。

画两条：

```text
Frozen/project-local trajectory MC
Transport operator
```

如果二者在统计误差内一致，这比做几十张复杂 phase-space 图更重要。

---

# 69. 第二张 validation 图

画：

$$
N_\chi
$$

vs

$$
\sigma.
$$

比较：

```text
trajectory
transport
```

在 overlap 区域一致。

然后 trajectory curve 因计算成本停止。

只有通过原 G5 的 AP/kinetic–diffusion、高 opacity 收敛及所需有限年龄验证后，transport curve 才可延伸到：

$$
10^{-30}\ {\rm cm^2}.
$$

这张图会非常直观地说明新方法存在的必要性。

---

# 70. 项目的最终计算链

完整 production pipeline 应该变成：

$$
\boxed{
m_\chi,\text{ interaction model}
}
$$

↓

$$
\boxed{
\text{build normalized collision kernel}
}
$$

↓

$$
\boxed{
Q_{\rm coll}^{(0)}
}
$$

↓

对通过 §19 缩放验证、保持同一物理/网格条件的 \(\sigma\)：

$$
\boxed{
Q(\sigma)
=
Q_{\rm stream}
+
Q_{\rm boundary}
+
\frac{\sigma}{\sigma_0}Q_{\rm coll}^{(0)}
}
$$

↓

计算：

$$
\boxed{
C_\alpha(\sigma)
}
$$

↓

求解：

$$
\boxed{
-Q^T N=C
}
$$

或者：

$$
\boxed{
N(t_\odot)
}
$$

↓

得到：

$$
\boxed{
n_\chi^{\rm in}(r)
}
$$

↓

surface bound flux：

$$
\boxed{
\Phi(E,L)
}
$$

↓

analytic exterior occupation：

$$
\boxed{
n_\chi^{\rm out}(r)
}
$$

↓

最终：

$$
\boxed{
n_\chi(r),
N_\chi,
E_{\rm evap},
\Gamma_A,
\text{observable signals}.
}
$$

这里 \(Q_{\rm coll}^{(0)}\) 明确表示 reference \(\sigma_0\) 下的生成元。Offline 缓存 kernel/rate/exterior，online 独立构造 \(C(\sigma)\) 并求解。稳态与 finite-age 分别选 \(Q_{\rm ss}\) 和包含物理时间的模型；后者外部占据需历史卷积。湮灭及观测信号保留为后续扩展，性能后端不改变这些科学门槛。

---

# 71. 最终科学上能够得到什么

这个项目最终并不只是一个“DaMaSCUS acceleration”。

它可以回答：

### 1. 什么时候传统 orbit-averaged treatment 失效？

通过：

$$
t_{\rm scat}/T_{\rm orbit}.
$$

### 2. 什么时候 global thermal distribution 失效？

直接比较 kinetic solution。

### 3. 高截面下 evaporation 为什么被 rescattering 抑制？

直接看到：

$$
E>0
\rightarrow
E<0
$$

population flux。

### 4. weakly bound exterior population 到底有多少？

直接由 surface flux × orbital residence time 得到。

### 5. heavy evaporation tail 从哪里来？

通过包含外部物理时间的算符/延迟模型及其慢模与尾部收敛，不能由 Q_ss 的谱直接代替。

### 6. capture–evaporation steady state 是否在太阳年龄内建立？

通过：

$$
N(t_\odot)
\quad{\rm vs}\quad
N_{\rm ss}.
$$

### 7. 从 optically thin 到 optically thick 的 distribution 如何连续变化？

最终形成完整 regime map。

这比单纯再多计算几个 DaMaSCUS benchmark 的科学空间大很多。

---

# 72. 推荐的第一批 Git commits / PR

以下 PR 标签是原路线的主题拆分；实际工作包、前置依赖和通过标准统一以 Task_Plan 的 T00–T20 / G0–G6 为准，不据本列表改变执行顺序。

## PR-0 — Baseline

```text
Freeze trajectory benchmarks and fixed seeds.
No physics changes.
```

## PR-1 — Physics port

```text
Port the minimal SolarBackground and ScatteringPhysics closure.
Keep provenance and match frozen legacy outputs.
```

## PR-2 — Single-collision tests

```text
Add statistical regression tests for rate and post-collision distributions.
```

## PR-3 — Phase-space grid

```text
Implement PhaseSpaceGrid(r,v,mu).
No transport yet.
```

## PR-4 — Collision kernel

```text
Build Gamma_alpha and P_alpha_beta.
Write sparse operator to disk.
```

## PR-5 — Collision-only solver

```text
Verify conservation and local thermal equilibrium.
```

## PR-6 — Streaming operator

```text
Implement conservative gravitational streaming.
```

## PR-7 — Bound/unbound surface boundary

```text
Implement escape, bound-return, and configured outer-removal channels; bridge Rsun/1.1Rsun.
```

## PR-8 — Capture source

```text
Reuse existing first-post-scatter-bound stop; export source states and absolute C_alpha.
```

## PR-9 — Steady solver

```text
Solve -Q^T N = C.
Output absolute n(r).
```

## PR-10 — Trajectory validation

```text
Compare same-source absolute densities and residence times with matched domains and censoring policy.
```

## PR-11 — Exterior orbit occupation

```text
Reconstruct n_out(r) analytically from surface flux.
```

## PR-12 — High-opacity diagnostics

```text
Compute mean free path and Knudsen map.
Identify need for AP/diffusion treatment.
```

之后才进入 production high-\(\sigma\) solver。

---

# 73. 当前关键实施任务

T00–T04 已完成，G0 局域物理门按约定条件关闭；当前主线是完成 T05 reference collision kernel，再进行 T06 离散热浴验收。具体状态、证据和剩余门槛只在 [Task_Plan §8](Task_Plan.md#current-status) 维护。

---

# 74. 项目 Definition of Done

当以下条件全部成立时，第一代项目才可以认为完成：

$$
\boxed{
\text{microscopic scattering parity}
}
$$

通过；

$$
\boxed{
\text{particle conservation}
}
$$

通过；

$$
\boxed{
\text{full-trajectory overlap benchmark}
}
$$

通过；

$$
\boxed{
N/C=\langle\tau\rangle
}
$$

通过；

$$
\boxed{
\Phi_{\rm evap}+\Phi_{\rm outer}=C
}
$$

按所选边界模型稳态测试通过；kepler_return 时 Phi_outer=0。守恒闭合不替代物理准确度验证；

$$
\boxed{
n_{\rm out}^{\rm operator}
=
n_{\rm out}^{\rm trajectory}
}
$$

通过；

完成：

$$
N_r,N_v,N_\mu
$$

收敛；

完成：

$$
v_{\max}
$$

收敛；

完成：

$$
\lambda_{\rm mfp}/\Delta r
$$

diagnostic；

明确标记哪些结果：

```text
kinetic solver reliable
```

以及哪些结果：

```text
require opaque-regime treatment
```

之后才能正式把结果推进到：

$$
\sigma\sim10^{-30}\ {\rm cm^2}.
$$

上述 V1 科学闭环对应 G3，所有平均值需同 source/域/终点且无未解释的计算删失。G4 物理时间与 G5 高 opacity 验证独立通过后才能扩展相关结论。optimized 后端还须通过 reference parity 和性能报告；实现加速不替代任何原科学 gate。

---

# 75. 项目的核心原则

整个新项目始终坚持一句话：

$$
\boxed{
\textbf{Do not simulate the history if the occupation measure can be solved directly.}
}
$$

DaMaSCUS-SUN-EVAP 在本项目中的角色是只读素材：

> **提供固定的 microscopic scattering 源码事实与已保存的 trajectory-level ground truth。**

而新的 SolarDM-Transport 负责：

> **把这些局域微观过程组合成长期、绝对归一化的太阳暗物质相空间分布。**

这才是从当前轨迹模拟器向高截面、低质量、强蒸发区域真正可扩展的方法。

---

<a id="performance-design"></a>

## 附录 A：工程与性能设计

设计更新：2026-09-14。本附录结合用户提供的性能建议、数值方法与文献意见和第一手来源，定义参考实现、优化实现及 P00–P09 的技术契约。科学任务和阶段编号仍以 [Task_Plan.md](Task_Plan.md) 为准；以下优化尚未实现或测得效果。

结论：采纳分层后端、局域核复用、解析方位角沉积、离线/在线划分和独立性能任务；修正矩阵方向、预条件零模、随机数一致性和重要性采样的表述。性能优化从设计期预留接口，启用时必须通过参考实现对照。

### A.1 当前设计决策与启用条件

以下维护最终决策，相关推导集中在对应章节；不逐轮追加审阅副本。

| 主题 | 当前决策 | 验收或适用范围 |
| --- | --- | --- |
| 定位与文献 | 相位分辨的多尺度动理学框架 | §2.1 的贡献目标待验证，不作首次声明 |
| 科学阶段与后端 | 保留 G0–G6，reference/optimized 共用契约 | 硬件 P 任务不替换有限年龄与高 opacity 科学门槛 |
| streaming | T07 提前检验不变量与假蒸发 | §12.1；细化不达标时提前评估不变量坐标/保守特征线 |
| 碰撞与源投影 | 完整平均 J Gamma(P-I)，原始源所属格计数为参考 | §16.1、§21；CIC 为有正性/支持约束的可选投影 |
| 局域核与角积分 | 联合 (v',c_lab) 加解析 phi 和 FV 输入平均 | 降低角分箱噪声，保留速度/角度相关性及稀有尾验证 |
| 离线/在线与表格 | 有指纹、失效规则和独立误差预算 | 微观生成元可条件缩放，宏观多散射响应不能照搬 |
| CPU 并行/RNG | 私有可写状态、逻辑随机地址、分层复现 | T02 保留 legacy RNG；优化后端另验线程安全与统计 |
| matrix-free/预条件 | ApplyQT 为前向；块求解加宏观粗校正候选 | 处理零模，不能以加假 sink 或仅加速 matvec 代替 |
| 伴随与重要性 | 目标导向 pilot、支持完整、冻结算符、独立验证 | 无偏 Q 不保证求逆后的 N 无偏；固定源总稳态逃逸通常无灵敏度 |
| 物理时间 | T17 采用相位状态或经验证的返回/占据核 | §30.1；保留长延迟、初始外域及实际移除时间 |
| 稀有轨迹 | T11 预留续跑/权重，T15 可选 WE pilot | 保留底层动力学；检查谱系误差，不承诺减少单路径散射 |
| 高 opacity | T18 比较 AP、确定性 hybrid、DDMC 三候选 | 从 DM 碰撞/引力推导极限与接口；选择通过误差和成本门槛的路线 |
| 宏态/Milestoning | T20 保留入口通量、等待核及必要历史 | 均值或稳态可简化的条件不等于全瞬态精确 |
| QSD/角矩 | T20 的专题诊断和表示候选 | QSD 不替代绝对占据；角矩截断不自动保正 |
| GPU/I/O | 先 profile optimized CPU；选择二进制和设备布局 | 不预先排名 GPU 调度；首次/缓存与最终误差共同计费 |

---

### A.2 两套后端，共用物理契约

`reference_cpu`：CPU、float64、小网格、明确矩阵、直接调用 DM-Transport 内已验证的移植物理层。保留 `(r,v,mu)` 直接采样路线作为独立对照，不因优化后端出现而删除。

`optimized_cpu`：经过验证的局域核复用、只读表/私有缓存、MPI/OpenMP、分块或 matrix-free 算符、预条件及 continuation。GPU kernel builder 与 GPU/distributed solver 是可选的两个独立能力，不互为必需前提。

主碰撞接口面向 `sample_collision(r, velocity_in, model, rng) -> velocity_out`，总率接口面向 `total_rate(r, speed, model)`。靶编号、靶速度、动量转移、能量交换属于诊断。生产表格由同一已验证物理生成；任何抽样重写/插值都视为新数值近似，独立比较联合分布和最终观测量。

所有后端共用状态索引、FV 几何、单位、边界、source、输出 schema 和误差定义。参考后端使用小矩阵并不意味着任意高 opacity 下都正确；两个后端相符也不能替代 T18 的扩散极限检验。

### A.3 局域旋转对称核与解析角度沉积（P02）

#### A.3.1 核的对象与适用条件

对各向同性、无偏振热靶且无外部优选方向的局域环境，旋转协变要求：同时旋转输入和输出速度，散射概率不变。因此可缓存

\[
K_{loc}(r,v;v',c),\qquad
c=\frac{\boldsymbol v\cdot\boldsymbol v'}{|\boldsymbol v||\boldsymbol v'|}.
\]

K_loc 表示一次碰撞条件下的联合概率；若使用连续密度，约定相对于 `dv' dc`，全域积分为 1。从 `d³v'` 密度转换时必须包含对应 Jacobian。这里 c 是 **实验室系 DM 碰撞前后速度夹角**，不是旧 `cos_alpha`，也不是质心系入射相对速度的散射角。必须保留 v' 与 c 的联合分布，不可分别采样其边缘分布。微观样本来源可以是固定轴上的向量碰撞，再通过速度点积提取 c。

旋转对称不要求微分截面在散射角上为常数；只依赖相对散射角的前向/后向分布也可能满足。当前 MVP 限制为已验证 SD 支路，是源码物理正确性的限制；通过旋转测试不能证明旧一般散射支路的参考轴正确。

相同 `(r_i,v_j)` 的样本可以用于全部入射 mu 单元。约 `N_mu` 倍只指原来重复的微观样本调用数上限；角沉积、存储、FV 积分和求解仍有成本，实际加速需测量。相空间状态数和速度/角度耦合仍保持三维。

#### A.3.2 方位角积分

给定非零输入/输出速度及 c，在输入速度轴周围方位角均匀。对给定入射 mu：

\[
\mu'=A+B\cos\phi,\quad A=\mu c,\quad
B=\sqrt{(1-\mu^2)(1-c^2)},\quad \phi\sim U[0,2\pi).
\]

输出单元 l 的边界为 `[mu_l^-,mu_l^+)`，最后一个单元包含 +1。`B>0` 时令

\[
z_-={\rm clip}((\mu_l^--A)/B,-1,1),\quad
z_+={\rm clip}((\mu_l^+-A)/B,-1,1),
\]

则其概率为

\[
W_l(\mu,c)=\frac{\arccos z_--\arccos z_+}{\pi}.
\]

`z_+<=z_-` 直接取零。接近端点时可用下式避免相近反余弦相减，并验证数值等价：

\[
W_l=\frac{2}{\pi}\operatorname{atan2}
\left(z_+-z_-,\sqrt{(1-z_-)(1+z_-)}+\sqrt{(1-z_+)(1+z_+)}\right).
\]

`B=0` 是 `mu'=A` 的确定性分箱；不得除零。对于很小但非零 B，要检查其支持是否跨越单元面，不能一律当成 delta。浮点越界仅修正已确认的舍入误差并记录容差；v 或 v' 为零时 c 未定义，走 T03 的静止/退化速度契约，不能直接套此式。

#### A.3.3 有限体积平均与验收

对于入射 mu 单元 k 内常数 f，角度 measure 为 `dmu`，应使用

\[
\overline W_{kl}(c)=\frac1{\Delta\mu_k}
\int_{\mu_k^-}^{\mu_k^+} W_l(\mu,c)\,d\mu.
\]

入射单元积分用正权求积并做求积收敛；若 MVP reference 采用中心代表值，两后端应先用相同规则对照，再把 cell-average 作为共同的离散改进。r/v 的代表点或单元求积也必须一致；完整事件率加权公式及求积门槛见 §16.1，不能分别平均 Gamma 与 P 后相乘。对 `(v',c)` 先粗分箱再取中心会新增核压缩误差，不能把它当作纯复用。

P02 的验收包括：`sum_l W=1`、W 非负、c=±1/入射 mu=±1 的端点处理、对称反射、显式 phi 抽样与积分一致；对于精确 FV 几何平均还应满足 `Delta_mu_k Wbar_kl = Delta_mu_l Wbar_lk`，并保持各向同性占据 `Delta_mu/2`。这些几何测试之外，还要重跑热浴、near-escape 联合分布及最终绝对占据比较。

独立条件角矩检查为 `E_phi[mu']=mu*c`、`E_phi[mu'^2]=mu²c²+(1-mu²)(1-c²)/2`；从分箱中心重建的角矩另有分箱误差，不与连续条件矩混为一谈。

解析 phi 沉积是条件平均，只减少角度分箱噪声。给定 r，能量阈值是否被跨越由 v' 决定；未采到的高速 v' 事件不会凭角积分出现。共享样本使不同 mu 行的误差相关，误差重复应按独立 local-kernel 样本批次或 seed 进行，不能把展开后的各行当独立 MC。

### A.4 离线资产、在线参数点与表格化（P03）

| 阶段 | 内容 | 复用条件 |
| --- | --- | --- |
| 离线 | 太阳背景、分靶率/Gamma0、局域联合核、FV 沉积几何、streaming 模板、外域停留核 | 模型、质量、耦合/靶集合、背景、网格/速度域、单位、边界和 schema 指纹一致 |
| 在线 | C(sigma)、缩放碰撞算符、选择/更新预条件器、求解、重建观测量 | 每点重新验证来源和残差；C 不沿用 sigma 线性缩放 |

固定条件下 `Q_coll(sigma)=s Q_coll,0`，`s=sigma/sigma0`；源、外域策略或 adaptive 网格变化会使相关缓存失效。hybrid 接口随 sigma 改变时，也不能宣称整个混合算符只需乘一个系数。

表格化先从背景、分靶 rate 和目标选择概率开始。更深的表只能按经过推导的条件链保存，或直接保存联合 `(v',c)`；附件提出的 `F(q,cos(theta)|...)` 没有定义完整的热靶条件关系，不作为主 API。对 SD，简单解析采样可能比高维 CDF 更省成本，P00 决定是否值得建表。

新增表格误差单列预算：CDF 单调性/端点、速率精度、联合分布、平均能量交换、速度尾、热浴平衡与最终 N。table interpolation、角度几何积分、统计误差分开评估。表越界显式失败或回到 CPU reference；不使用未经验证外推。

### A.5 CPU 并行与随机数（P01、P07）

#### A.5.1 先处理共享可写状态

参考源码的插值查询会更新内部搜索缓存，核素表还存在惰性初始化；详见 [源码依据附录](Proposal.md#code-reference)。项目内移植不能把含可写缓存的背景对象直接放进 `omp parallel for`。

先串行完成依赖初始化；worker 使用私有物理/插值缓存或经过验证的只读表，以及私有沉积缓冲。MPI 分配质量点或稳定编号的 state groups，OpenMP 处理节点内样本块；static/dynamic 调度按负载实测选择。合并顺序固定，控制 MPI×OpenMP×BLAS 线程数，避免超额占核。

参考代码请求 `MPI_THREAD_FUNNELED`，但 T02 不移植其 MPI 调度。P01 若在本项目引入 MPI，须重新声明并检查实际线程等级；内核构建和 source 生成分别测强/弱缩放，不保证同一种分块策略适用两者。[Open MPI 文档](https://docs.open-mpi.org/en/main/man-openmpi/man3/MPI_Init_thread.3.html)

#### A.5.2 随机地址，而非线程编号

T02/T03 的项目内移植保留 legacy `std::mt19937` 及调用顺序完成回归。固定 obscura API 明确要求 `std::mt19937&`；counter engine 无法不改接口就直接替换。优化后端再引入 Philox 等候选，并记录版本和测试向量。counter 方法的原始参考见 [Salmon 等，2011](https://www.thesalmons.org/john/random123/papers/random123sc11.pdf)。

建议逻辑地址：

```text
key     = (master_seed, normalized_model_grid_fingerprint, rng_version)
counter = (local_state_id, sample_id, sampling_stage,
           rejection_attempt, draw_index, output_lane)
```

字段编码、范围、溢出策略和不同抽样步骤的 domain separation 都要固定；拒绝采样不能只分配有限几个共享随机数。normalized-kernel key 不因单独缩放 sigma 改变；source 使用对应物理参数和独立随机域。并行调度不改变 sample ID，adaptive 追加也不能复用旧 ID。

验收分三级：相同 key/counter 的原始随机整数逐位一致；同平台/算法及确定性归约的结果可复现；跨 CPU/GPU 的变换后样本和统计量满足数值/统计容差。`log/sqrt/trig`、融合运算、拒绝判定和归约次序可能不同，不承诺最终逐位一致。GPU 初期维持 double 精度基准，降精度另立误差验证。

### A.6 分量算符与预条件（P05、P06、P08）

#### A.6.1 明确主算符是 Q 的转置

继续使用行速率 Q、列占据 N。固定 s 与已冻结的核样本：

\[
\operatorname{ApplyQT}(x)=Q_{stream}^Tx+sQ_{coll,0}^Tx+Q_{boundary}^Tx,
\qquad \operatorname{ApplyA}(x)=-\operatorname{ApplyQT}(x).
\]

接口同时提供 `ApplyQ`，用于伴随问题和内积检验。reference 中的显式矩阵与三个分量在随机向量上逐项比较，并测试 `<u,Qv>=<Q^Tu,v>`、非负转移和源汇账目。一轮普通 Krylov 求解内必须冻结核；每次 matvec 重新随机抽核会改变线性算符。

streaming 用邻接 stencil，boundary 用通道映射，collision 可用径向块、局域核展开或经过验证的压缩表示。matrix-free 避免全局装配/重复拷贝，但如果仍保存所有稠密块，其碰撞存储阶数仍是 `Nr(Nv Nmu)^2`。核展开的每次应用成本也需实测，必要时缓存常用块。

PETSc shell matrix 可提供自定义算符应用；自定义数据或 sigma 变化时须通知 operator state，并处理预条件器失效。[PETSc MatCreateShell](https://petsc.org/release/manualpages/Mat/MatCreateShell/)

#### A.6.2 碰撞零模与块预条件

P05 分为基础算符/伴随接口、块预条件和后续粗校正三个可独立验收部分。P04 只依赖前两者中足以正确求解前向/伴随系统的能力，不必等待 DSA 完成；高 opacity 粗校正与 T18 的扩散推导衔接。

设 `B_i=Q_coll,0` 的第 i 个径向块。初始 block-Jacobi 候选可取实际 `A=-Q^T` 的径向对角块，或较便宜的近似

\[
M_i=-sB_i^T+D_i,\qquad
D_i=-\operatorname{diag}(Q_{stream}^T+Q_{boundary}^T)_i.
\]

这是预条件矩阵的候选，不是普适可逆公式。`B_i 1=0`，故纯碰撞块奇异；右平衡模属于 `B_i^T`。D 能否解除零模取决于实际状态连通性和流出，必须逐块检查。不可直接写 `B_i^{-1}`、默认 SPD 或使用未经证明的 CG/Cholesky。

失败时考虑受约束零模求解、micro–macro 分解、宏观输运粗空间或扩散校正。允许在预条件器内使用有说明的数值移位，但不能给真实物理 A 增加虚假漏失。时间隐式步应针对 `I-dt Q^T` 构造预条件器，与稳态问题区分。[PETSc 零空间接口](https://petsc.org/release/manualpages/Mat/MatSetNullSpace/)

独立径向块构成 block-Jacobi；块内可用 ILU/局部求解，跨块 block-ILU 则不再完全独立。预条件器的内存/构建时间、全局迭代数和跨 sigma 退化都纳入评估。[PETSc PCBJACOBI](https://petsc.org/release/manualpages/PC/PCBJACOBI/)

P05 进一步比较 DSA 类宏观扩散粗校正。若局域碰撞块具有单一归一平衡模 \(e_i\)，令 \(B_i^Te_i=0\)、\(\mathbf1^Te_i=1\)，可将向量分为 \(e_i\rho_i+g_i\)，其中 \(\mathbf1^Tg_i=0\)：micro 部分在受约束子空间求解，macro 部分由与原离散一致的跨半径粗算符修正。若存在多个守恒/慢模，粗空间必须相应扩展，不能强行只留密度。

低阶校正系数和边界来自本项目碰撞/引力的输运极限，避免套用单速中子的 \(D=v\lambda/3\)；在每次 sigma 扫描比较块预条件与粗校正的构建成本、迭代增长和实际残差。DSA 的任务是加速同一个 \(A=-Q^T\) 的求解，AP 的任务是保证离散本身具有正确极限，两者不能相互替代。薄厚交界/异质背景下的粗校正可能退化，须做局部与整体稳定性测试。[异质介质 DSA 原始研究](https://arxiv.org/abs/2001.09196)

#### A.6.3 continuation 与失效规则

同一状态空间扫描相邻 sigma 时可用上一个 N 作初猜；每点使用自己的 C 和 A，并独立检验真实残差及源汇。预条件器按迭代数、停滞或参数/网格/边界变化触发重建。若网格变化，先做保守投影且重新验证；不得直接复用向量索引。

块预条件器和 GPU matvec 处理代数成本，不能修复粗网格 kinetic 离散在强散射区的错误扩散极限。T18/P09 的 AP/hybrid/DDMC 方法同时承担物理正确性与成本削减目标，后续才讨论加速到何种硬件。

### A.7 稀有通道、自适应与伴随（P04）

#### A.7.1 选择能反映核误差的观测量

从粗核得到 N 后，联合考虑高占据区域、near-escape、surface outgoing 及 escape/rebinding 通道。每个 local state 仍保留最低覆盖；pilot 只用于分配下一批预算，独立验证批次用于评估结果，避免数据依赖的停止规则被误当固定样本无偏估计。

对于 `AN=C`、`A=-Q^T`、`J=g^TN`，定义伴随 `A^T z=g`。固定离散空间下的一阶变化为

\[
\delta J=z^T\delta C+z^T\delta Q^TN+(\delta g)^TN.
\]

若 `q_alpha,beta += epsilon`、`q_alpha,alpha -= epsilon` 且 g/C 不变，则 `delta J = epsilon N_alpha (z_beta-z_alpha)`。共享局域核同时改变多条 mu 行，需要把这些相关变化一起计算，而非逐行当独立误差。

在有限稳态存在、固定 C、唯一损失为逃逸且扰动保持粒子守恒时，`Phi_escape=C_total`；因此总逃逸流对这类核扰动的灵敏度恒为零，不适合作为唯一采样目标。优先选择 N、平均驻留、外部密度、逃逸谱/角通道，或有限年龄指标。若 g 自身含逃逸率，必须包含上式 `delta g` 项，它参与上述抵消。伴随线性化只辅助分配预算，最终结论仍需独立重采样和误差检验。

#### A.7.2 重要性采样不等于自动无偏归一化

设目标局域碰撞分布为已规范化 p(z)，proposal 为 h(z)，在目标支持上 h>0，权重 `w=p/h` 可计算。对固定样本数，`mean(w a_beta)` 是单元概率的无偏估计；它的有限样本行和未必为 1。除以 `sum w` 的自归一估计有有限样本偏差，不能一面强制归一化、一面仍声称严格无偏。[Owen，Importance sampling，第 9 章](https://artowen.su.domains/mc/Ch-var-is.pdf)

本项目可优先研究对**生成元**的逐样本保守估计。设 `a_beta(z)` 为含角积分的非负沉积权重，完整去向上 `sum_beta a_beta=1`：

\[
\widehat Q_{\alpha\beta}=\frac{\Gamma_\alpha}{n}
\sum_{k=1}^{n} w_k\big[a_\beta(z_k)-\delta_{\alpha\beta}\big].
\]

于是每个样本贡献行和为零，非对角元非负，在已知正确权重和固定 n 条件下期望给出目标生成元。对角离开率使用同批加权流出之和；不要把中间加权 histogram 当成固定 Gamma 下已精确归一的 P。数值 overflow 仍须按 T05 处理，不能变成物理 escape。

实施前推导完整条件链的 proposal 密度与 Jacobian；只给出最终态 `p/h` 符号而无法计算归一密度不够。记录权重均值、ESS、最大权重、批次方差和尾通道区间，并对照普通 MC 可验证区域。proposal 要覆盖 bulk 与尾部；自适应归一/截权若使用，必须另报偏差。无偏 Q 也不保证有限样本求解后的 N 无偏，需传播核误差。

#### A.7.3 目标导向采样的实施闭环

P04 采用 `粗参考核 → 前向 N / 伴随 z → 样本预算或可计算的 proposal → 冻结新核重求解 → 独立批次验证`。用局域样本对目标 J 的贡献及成本估计方差，而不只按 N 排序；相关的 mu 行作为一个采样组，必要时为多个预先指定观测量分配预算。

这借鉴 [CADIS/FW-CADIS](https://www.ornl.gov/file/forward-weighted-cadis-method-global-variance-reduction/display) 的前向/伴随重要度思路。原方法的源偏置与权窗并不自动提供本项目随机核的正确 likelihood；若只做样本数分配，不额外修改动力学。若改变碰撞或源 proposal，先构造带基线覆盖、归一密度可计算的条件抽样，再分别更新核权重或源权重。检查样本支持、权重极值、批次方差及最终 N 的偏差，不能仅凭伴随图形合理就宣称无偏加速。

### A.8 GPU、I/O 与性能验收（P00、P07–P09）

先测 optimized CPU。GPU 优先尝试局域独立碰撞样本，使用已验证的只读物理资产和联合抽样；不把整个 obscura 对象图迁入设备。GPU 的随机变换、拒绝采样、沉积/归约、传输和内存都计入端到端成本。目标设备、double 性能和软件支持在 P07 前调查，未指定 GPU 时不预设 CUDA 可用。

matrix-free/GPU/distributed 求解属于 P08，可独立于 P07 采用 CPU 生成的核。高 opacity 的离散方法与预条件器先通过 T18/T19；GPU 不是达到 `1e-30 cm²` 物理可靠性的门槛。

reference/debug 可保存 Matrix Market。production 按分块访问、索引宽度和元数据需要选择 HDF5、原生 CSR 二进制或 PETSc binary；`nnz≈1e7` 触发序列化预算检查，不是普适理论阈值。保存 schema、校验和、块布局和压缩方法。无需显式全局 Q 的后端，不为满足旧文件名额外生成巨大 Q。

正式 `performance_report.json` 至少包含：rate/table、kernel sampling、角沉积、source、assembly/operator setup、precondition、solve、exterior、I/O/传输的时间；samples/s、拒绝次数、GMRES iterations、matvec time、peak RAM/device memory；硬件、线程/MPI 布局、精度、RNG、cache hit、误差预算及关联 validation ID。

比较四项：项目内 full trajectory MC（首个点可含 T01 固定基线）、首次 reference transport、首次 optimized transport、缓存后新增 sigma。首次运行包含离线和预条件成本；缓存扫描仍包含 source 重算、失效重建、I/O 和求解。各方案必须达到同一统计和离散误差目标；只报告实测加速，不预先排名各优化的收益。

### A.9 实施顺序与失败回退

近期完成 T00–T06，P00 从 T01 开始记录成本；通过物理、截面缩放和旋转对称验证后开展 P02/P01，T07–T16 的参考科学路线继续推进。表格化、matrix-free、自适应、continuation 各按自己的依赖接入，GPU 等待 optimized CPU 的证据与 profile。

任何优化若不能在相同离散/边界/源下通过局部及最终观测量对照，保留 reference 路径并将该能力标记为未验证。性能任务的准确依赖、交付物和验收条件列在 [Task_Plan.md](Task_Plan.md) 的 P00–P09 表。

### A.10 跨领域方法依据与迁移限制

此处保留影响设计的第一手来源；附件中的检索站首页、聚合页或未核验的条目不作为实施依据。下列迁移判断是设计推论，尚无本项目数值验证。

| 方法与来源 | 本项目采用的思路 | 不能直接外推的结论 |
| --- | --- | --- |
| [Lemou–Mieussens，micro–macro AP](https://doi.org/10.1137/07069479X)；[Lemou–Méhats，边界层](https://arxiv.org/abs/1202.1994) | T18 推导宏观/微观分解，并同时验证体内与边界的扩散极限 | 原文模型不自动给出含引力、能量交换 DM 的闭合；AP 不自动保证尾部和正性 |
| [Densmore 等，2007，DDMC](https://doi.org/10.1016/j.jcp.2006.07.031) | T18 的连续时间扩散跳跃与 kinetic 接口候选 | 扩散计时不等于微观首达定律；依赖正确厚区闭合 |
| [Zhang–Zuckerman–Jasnow，WE](https://arxiv.org/abs/0810.1963) | T11/T15 带权路径重采样，不更换底层动力学 | 克隆不增加等量独立样本，不能预设任何 progress coordinate 都有效 |
| [Bello-Rivas–Elber，Exact Milestoning](https://pmc.ncbi.nlm.nih.gov/articles/PMC4352169/)；[Vanden-Eijnden 等，原理条件](https://pubmed.ncbi.nlm.nih.gov/19045328/) | T20 的入口 crossing 分布、短片段及驻留时间统计 | 任意宏态索引或起点平衡分布不保证精确；正确均值公式不保证完整瞬态 |
| [ORNL，FW-CADIS](https://www.ornl.gov/file/forward-weighted-cadis-method-global-variance-reduction/display) | P04 用目标/伴随信息配置预算与重要度 | 源偏置权窗的理论不是自适应随机 Q 求逆的无偏证明 |
| [Southworth 等，异质介质 DSA](https://arxiv.org/abs/2001.09196) | P05 用低阶输运/扩散校正跨块慢模 | 加速原离散不修复错误扩散极限，薄厚交界仍需验证 |
| [Hoagland 等，2021，ITMM](https://doi.org/10.1080/00295639.2021.1898309) | 已有微观/输运程序构造局部响应，再全局迭代的架构参照 | 含多次散射的 response 不同于单碰撞核，通常不能沿用线性 sigma 缩放；瞬态响应还需时间核 |
| [Asselah–Ferrari–Groisman，QSD](https://arxiv.org/abs/0904.3039) | T20 检查晚期幸存者分布 | 有限状态、大粒子数结论不直接覆盖无限外域；QSD 不替代绝对占据 |

ITMM 作为现有离线/在线与 P05 的参照，不再新增独立实施任务。角矩按 §63 作为表示候选；GPU 的 history/event 分组按 A.8 实测，不从其他粒子输运代码的排名推定本项目速度。

---

<a id="code-reference"></a>

## 附录 B：源码事实与复用接口

源码检查日期为 2026-09-14，参考快照固定为 `b5678f5b193aa567ca10715c2a6c764c9e72eec7`；下表链接均指向该只读提交。DaMaSCUS-SUN-EVAP 不接受本项目修改，也不参与后续构建；“规划处理”只描述在 DM-Transport 内的移植或独立实现。本附录记录源码事实；运行证据和任务状态见 [Task_Plan.md](Task_Plan.md)，优化契约见 [工程与性能设计附录](Proposal.md#performance-design)。

### B.1 文档和源码的优先关系

- [本项目 Proposal](Proposal.md) 是科学目标依据。v2 已同步修正本文指出的文件名和能力描述；[v1 存档](archive/Proposal_2026-09-14_v1.md) 保留原始表述，以下差异指原稿与本次检查源码的差异。
- [DaMaSCUS README](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/README.md) 说明当前 capture、evaporation、边界与输出工作流。
- [PROJECT_DOCUMENTATION](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/PROJECT_DOCUMENTATION.md) 有完整物理说明，但仍含历史逐轨迹文件、MPI 批量同步等描述；当前行为以源码和 README 为准。
- 2026-09-10 代码审查（未随参考仓库发布）针对旧提交 `543660f…`，不能把整份问题列表直接当作当前未修复缺陷。本轮只将再次定位到当前代码的问题用于规划。

### B.2 可以复用的真实接口

| 模块 | 源码证据 | 已有能力 | 规划处理 |
| --- | --- | --- | --- |
| 太阳背景 | [Solar_Model.hpp:62](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/include/Solar_Model.hpp#L62) | Mass、Temperature、Local_Escape_Speed、核/电子数密度 | T02 在项目内移植最小只读表接口；不带入轨迹状态 |
| 散射率 | [Solar_Model.cpp:386](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Solar_Model.cpp#L386) | 逐靶 rate、总 rate、直接计算/插值路径 | T02 移植 SD 恒截面直接 rate 闭包；T03 验证数值与适用范围 |
| rate 插值 | [Solar_Model.cpp:518](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Solar_Model.cpp#L518) | 在 MPI_COMM_WORLD 上建立规则 r/v 表，速度上限固定 0.75 自然单位 | T02 不移植 MPI 建表；T05 在项目内按实际速度域另建缓存 |
| 靶选择 | [Simulation_Trajectory.cpp:2460](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L2460) | 根据各靶 rate 抽样、预分配核 rate 缓存 | T02c 已移植源顺序核靶 CDF 与显式 RNG；固定 MVP 无电子通道 |
| 热靶速度 | [Simulation_Trajectory.cpp:2498](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L2498) | 碰撞条件下的热靶速度采样 | T02c 已移植恒截面条件 sampler；零速保留为独立极限验证 |
| 单碰撞 | [Simulation_Trajectory.cpp:2587](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L2587) | 选择靶、抽靶速度、调用 obscura 角采样、修改速度 | T02d 已在项目内移植完整局域流程并显式接收 RNG；不改 legacy `Scatter` |
| 引力传播 | [Simulation_Trajectory.hpp:392](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/include/Simulation_Trajectory.hpp#L392) | Free_Particle_Propagator 已有独立类 | 作为 T07 只读算法依据；在本项目独立实现和验证 |
| 首次捕获即停 | [Simulation_Trajectory.cpp:2698](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L2698) | 散射后更新 capture state，capture mode 立即终止 | T10 在本项目移植判据并新增源输出 |
| 最终事件 | [Simulation_Trajectory.hpp:286](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/include/Simulation_Trajectory.hpp#L286) | Trajectory_Result.final_event 有位置和速度 | T10 在本项目转换成 r/v/mu，并保留连续原始状态 |
| halo 初态 | [Simulation_Utilities.cpp:292](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Utilities.cpp#L292) | 引力聚焦 impact parameter，按面积抽样入射盘 | T10 在本项目移植并与入射率统一归一化及采样权重 |
| 总入射率 | [Reflection_Spectrum.cpp:66](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Reflection_Spectrum.cpp#L66) | pi R² (rho/m) [mean(u)+vesc² mean(1/u)] | T10 在本项目移植成 source 公共函数，不带入反射谱/KDE 依赖 |
| 外部轨道 | [Simulation_Trajectory.hpp:84](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/include/Simulation_Trajectory.hpp#L84) | BoundKeplerExteriorArc、返回状态及壳层贡献相关接口 | T12 在本项目移植思路并参数化匹配面、壳边和退化极限 |
| 径向占据 | [Simulation_Trajectory.hpp:205](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/include/Simulation_Trajectory.hpp#L205) | TrajectoryBincount 保存时间加权统计与生存标记 | T11 在本项目移植径向方法，新增速度/角度统计与同源启动 |
| 参考构建 | [src/CMakeLists.txt:35](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/CMakeLists.txt#L35) | lib_damascus_sun 含物理、轨迹、MPI、快照、数据与参数扫描 | 不接入本项目构建；T02 在 DM-Transport 新建最小串行 target |

### B.3 与提案不同、需要变更任务描述的地方

#### B.3.1 不存在的文件与能力

参考源码没有 `Scattering_Rates.cpp`，也没有提案式独立 `Sample_Momentum_Transfer()`。实际角采样由 obscura 的 DM 类完成；靶选择、靶速度和出射速度构造都在 `Simulation_Trajectory.cpp` 内部，T02 必须在项目内划出最小边界。

当前 rate 插值表是重建并缓存数组，未发现提案所称“同质量更改截面时自动 rescale”的接口或模型签名失效保护。截面分离仍是首版 SD 恒截面模型可验证的物理性质，但 T05 需要新建、测试并记录缩放缓存规则。

#### B.3.2 首版 SD 模型必须显式限定

当前 [New_DM_Velocity](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L2577) 围绕 DM 实验室速度取散射方向；[角采样调用](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L2604) 也传入实验室速度。对于一般各向异性模型，这不同于围绕入射相对速度轴处理碰撞。当前 [rate 实现](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Solar_Model.cpp#L400) 则把 sigma(v_DM) 提到热平均相对速度之外，不能直接推广至速度相关截面。

旧审查中的这两类问题仍可定位到当前源码。MVP 因而选 SD 恒截面且各向同性支路；它仍需独立热浴与统计验证，不能把“与旧代码一致”当作物理正确性的全部证据。一般相互作用修正单独排期。

在当前固定版本 obscura 的 SD 实现中，总截面无速度依赖，low-mass 分支使用均匀 cos 角采样；low-mass 是 [显式配置开关](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Parameter_Scan.cpp#L558)，不是因输入质量小而自动启用。T01 保存被实际链接的依赖源码版本，T03 重新核对这些性质。

配置解析另有强制覆盖规则：`DM_isospin_conserved=true` 时 obscura 会设置质子/中子相同耦合并忽略 `(1,0)`；MVP 必须显式 `false` 加 `DM_relative_couplings=(1,0)`，并核对实际模型。旧示例的测试结果不能因此重新命名为质子单独耦合基准。Configuration.cpp（参考构建目录内 `_deps/obscura-src/src/Configuration.cpp:273`）

新接口若需要 q，应从实际动量差获得；旧 cos_alpha 的参考轴并非入射相对速度轴，不能不加检查就按相对速度散射角公式反推 q。修正角坐标可能改变固定种子轨迹，需要区分逐事件回归与分布回归。

#### B.3.3 rate 支持零速，不代表 sampler 支持

热平均相对速度有零 DM 速度极限，但 [Sample_Target_Velocity](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L2506) 拒绝 `vDM<=0`。MVP 有限体积速度中心可以避开零点，接口仍应说明此限制，并设计零速极限测试；不得用随意小速度替代而不检验误差。

#### B.3.4 source-conditioned MC 需要新入口

[Simulate 初始化](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L2626) 将初态当作未捕获 halo 粒子；传播又有 uncaptured-bound 数值异常防护。直接把 T10 的负能量状态传给现入口，不等于正确“从首次捕获继续模拟”。T11 需要显式设置捕获状态、计时起点、统计锚点和物理参考能量。

#### B.3.5 三种 rate 名称必须分开

[Capture summary](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Data_Generation.cpp#L2789) 同时包含捕获概率和计算吞吐量。新输出明确命名 `capture_probability`、`incident_rate_particles_s`、`capture_rate_particles_s` 与 `simulation_throughput_s_inv`。只有物理入射率乘捕获概率才能给出本计划的绝对 C。

数值注入面为 1.1 R_sun，但当前入射样本仍按“能命中 R_sun”抽取，不能把归一化公式的面积直接改成 `(1.1 R_sun)²`。另外，[capture mode 的光学深度推进路径](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Simulation_Trajectory.cpp#L1911) 与普通模式精度控制不同；T10 验证两种模式的首捕获分布统计相容，不把固定 seed 的逐事件相同当作未经核验的前提。

#### B.3.6 当前外域是具体的物理模型

[CMake 默认域](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/CMakeLists.txt#L47) 是 1 AU，外域壳数 423；直方图布局随编译配置变化。[README Outputs](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/README.md) 明确束缚轨道超域后在外行过界时移除，贡献单程驻留而不计作蒸发事件。

因此 Proposal 的“所有束缚外行都返回”和当前默认 MC 不是同一个边界问题；`Phi_escape=C` 不能直接用于含 outer removal 的旧数据。外部占据重建也要区分单程移除和双程返回。

当前 [Kepler 测试](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/tests/test_Kepler_Return.cpp#L89) 明确拒绝退化径向轨道；输运面通道涉及角度端点/极限时，T12 必须补支持或给出受控的端点处理。现有三维返回实现已有专门测试，不把旧报告“未更新三维位置”的结论直接沿用到当前版本。

#### B.3.7 轨迹统计不是未经筛选的完整寿命分布

当前 README 约定：数值失败不进入驻留统计；部分计算截断保留已接受驻留前缀；`evaporation_times.txt` 只含完整有效解束缚事件；1 AU 移除不是蒸发。该表的 `lifetime_unbinding_sec` 与实际逃逸面交叉时间也不是同一量。

T11/T15 必须按首次捕获后的全部有效样本与相同终点构造 occupation；不能仅用完整蒸发子样本均值来归一化所有捕获粒子。可以复用 `residence_jackknife_blocks.tsv` 的块级误差思想，但不能假定旧每列字段就等价于新统计量。

### B.4 可移植的验证逻辑

| 只读测试源码 | 可移植的验证逻辑 | 新项目仍缺的部分 |
| --- | --- | --- |
| [test_Solar_Model.cpp](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/tests/test_Solar_Model.cpp) | 背景、各靶 rate、插值节点/SD 路径 | 项目内移植 parity、模型签名与截面缩放 |
| [test_Physics_Validation.cpp:251](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/tests/test_Physics_Validation.cpp#L251) | Kepler 收敛、壳层驻留、太阳背景、入射面积律、角 CDF、零速平均速率 | 完整单碰撞联合分布、transport 数值误差及稀有尾 |
| [test_Simulation_Trajectory.cpp](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/tests/test_Simulation_Trajectory.cpp) | 轨迹分类、守恒径向统计与多种数值边界 | 已捕获初态入口、时间加权相空间占据 |
| [test_Kepler_Return.cpp](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/tests/test_Kepler_Return.cpp) | 返回几何与特定退化输入拒绝 | 任意壳边/匹配面、径向与近抛物极限的输运接口 |
| [test_Data_Generation.cpp](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/tests/test_Data_Generation.cpp) | capture 分母、截断、输出契约/错误路径 | 绝对 source 单位、原始首捕获样本与 source-conditioned 误差 |
| [validation/physics_validation.py](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/validation/physics_validation.py) | 速率插值/多种子运行矩阵与报告框架 | 新 schema、绝对归一化、统一边界、完整速度分布和核 MC 误差 |

旧代码审查后的固定参考已包含 quickstart 配置、输出失败测试和 Kepler 返回测试。T01 已记录当时的实际执行结果；后续开发只消费已保存基线。若未来需要更新参考版本，先由用户单独处理素材库并提供新的固定提交与产物，本项目不在其中重跑构建。

### B.5 对开发顺序的直接影响

1. 把提案“重写 capture 停止逻辑”改为“输出首次捕获样本并补绝对率”。
2. 在真正的 transport–MC 对比前新增“从已捕获状态启动的受控 MC 入口”。
3. 外域策略、统计终点和源归一化在 T00 固定；ExteriorOrbit 前移为边界构建依赖。
4. 项目内物理移植保持最小闭包；不搬动整个旧代码树，也不复刻 MPI/snapshot 调度。
5. 单碰撞 parity 与独立物理验证并列；一般相互作用问题不能被 parity 掩盖。
6. 有限年龄必须使用包含外部飞行时间的模型；高 opacity 必须在已验证 kinetic 基础上另做扩散极限验证。

### B.6 并行与随机数的源码约束

| 源码证据 | 本轮确认的事实 | 对性能任务的影响 |
| --- | --- | --- |
| libphysica 插值查询:168（参考构建目录内 `patched_dependencies/libphysica/Numerics.cpp:168`） | 查询更新 `correlated_calls`、`jLast` 搜索缓存 | 多线程共享同一插值对象存在可写状态竞争风险；P01 使用线程私有对象/缓存或只读表 |
| obscura 核素表:176（参考构建目录内 `_deps/obscura-src/src/Target_Nucleus.cpp:176`） | 全局 `all_nuclei` 在空表时惰性赋值 | P01 在并行区外串行初始化；不能只给 RNG 加锁就认为线程安全 |
| obscura 角采样 API:87（参考构建目录内 `_deps/obscura-src/include/obscura/DM_Particle.hpp:87`） | 形参明确为 `std::mt19937&` | Philox 不能直接作为替代参数；T02 项目内移植保留旧接口，优化后端另做抽样/接口验证 |
| [MPI 初始化:23](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/main.cpp#L23) | 请求 `MPI_THREAD_FUNNELED` | P01 若在 DM-Transport 新增 MPI，调用限初始化线程并检查实际提供等级 |
| [rate MPI 建表:532](https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP/blob/b5678f5b193aa567ca10715c2a6c764c9e72eec7/src/Solar_Model.cpp#L532) | 参考实现直接操作 `MPI_COMM_WORLD` | T02 不移植；P01 的离线构建与只读使用在本项目独立实现 |

依赖文件位置相对于参考构建目录，不作为公开下载链接；T01 已保存依赖 commit、补丁及代码指纹，依赖升级时必须重新定位证据。以上为源码中可写状态和接口的检查，未执行线程竞争测试，也不代表已存在 OpenMP 后端。

### B.7 对补充建议的技术澄清

- 局域旋转协变允许将微观样本保存为 `(r,v)->(v',c_lab)`；这是一项待实现优化，当前源码没有此缓存。其适用性不等同于“CM 微分角分布必须均匀”。
- 不新增人为必需的 `Sample_Momentum_Transfer` 或未定义联合条件的 `F(q,cos(theta))` API。保留向量碰撞链，表格化后仍需验证热靶与最终态的联合相关性。
- 本项目主动力学是 `Q^T N`，matrix-free 实现需提供 `ApplyQT`，并给伴随提供 `ApplyQ`；不能照抄附件的 `y=Qx` 用于占据演化。
- 保守碰撞块有零模，径向块预条件要针对 `A=-Q^T` 检查可逆性。矩阵存储优化和线性求解收敛不能代替高 opacity 的扩散极限验证。
- CPU/GPU 复现区分原始随机位、固定平台数值结果和跨设备统计相容。旧 RNG 回归与新随机后端验收分别记录，避免同时改物理提取和随机映射。

这些结论已落实为 P00–P09 及相应验收条款；逐项建议采纳表和推导见 [工程与性能设计附录](Proposal.md#performance-design)。
