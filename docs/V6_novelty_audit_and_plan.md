# V6 准备：问题审计、文献定位、定位陈述、逐节修改计划与实验 TODO

**作者决定已并入（2026-10-07）**：标题定为 "Sharp Robust QoS
Certification from Six-Scalar Spectral Summaries"（RIS 后缀版为备选）；
§2.3 重写为 "Six-Scalar Information Model and Certification Hierarchy"
（信息映射 ℐ、类 𝒞(s)、最强一致证书 B\*(s)=inf q\* 作为问题定义；
Level 0–4 层级含"信息增⇒类缩⇒证书可紧"单调原理，并以此解释 L_quad
有时更紧=信息更多；工程场景压缩为一句+免责句）；§2.4 改为"数学→无线
含义→旧界为何松"三段（I₀/a 改称幅度；删除 "physically impossible" 与
"precisely"；结论句="same information, used more fully"）；术语统一为
six-scalar information class / budget / partial-information certificate。
最终英文见 V6_replacement_text.md（已更新）。

依据：V5 全文（paper_v5.md / V5 DOCX）、五轮实验数据、2026-10-07 系统
文献检索（检索域：Buzano 不等式族、PSD 补全、算子插值/收缩延拓、最优
恢复、鲁棒 SINR 界、RIS 鲁棒波束、CSI 误差方向分解、压缩统计最坏二次
优化）。检索来源链接见文末。

---

## A. V5 问题审计（Critical / Major / Moderate / Minor）

### Critical（不改则投稿即处劣势）

**A1. 文献定位缺位且引用倒退。** V3 曾引用 Bottazzi–Conde
《Generalized Buzano Inequality》（Filomat 2023），V4.1/V5 将其删除，
参考文献只剩 4 条。而式 (7) 的交叉项界 |dᴴℓ| ≤ min{αs, δ√(β₀²−s²)}
的第二支正是 Buzano 型修正 Cauchy–Schwarz；定理 2 的补全论证属于
Parrott/算子球收缩补全与 Grone–Johnson–Sá–Wolkowicz PSD 补全（2×2
弦图特例）的经典框架；定理 2 的"信息类下确界"在概念上是最优恢复中
的 radius of information（Micchelli–Rivlin 1977；Melkman–Micchelli
1979）。当前引用状态两头吃亏：既没给审稿人看到作者知道经典谱系
（易被"这不就是 Buzano/Parrott"一击），又因不引而显得查新未做。
**修正**：恢复并扩充参考文献（Buzano 1974/Bottazzi–Conde 2023、
Grone et al. 1984、Parrott 1978 及其现代强化版、Micchelli–Rivlin
1977），novelty 表述统一采用限定式："To the best of our knowledge,
we have not found an existing result that gives this closed form under
exactly these retained summaries."

**A2. "previous envelope"（旧共享预算界 Q_old）无出处。** 全文把
min{β₀r, αx+βy} 型旧界当作比较基线但从未归属。若不处理，审稿人会
认定 strawman。检索表明：CSI 误差按"期望方向/正交补"分解并在两支上
分别取最坏的做法在 bounded-error 鲁棒波束文献中是常见手法（如
Sun et al., NOMA 认知无线电鲁棒设计 2018 明确做平行/正交分解；
Zhang–Ng–Schober 一线的 norm-bounded worst-case SINR 体系为背景）。
**修正**：在 Related Work 明确旧界的身份——"the envelope implicitly
obtained when the three interference gains are used separately, as is
natural in the bounded-error robust beamforming line"，并区分它与本文
的联合谱兼容使用；若旧界出自作者先前备忘录，也应如实注明。

**A3. 摘要与卖点排序违反理论论文定位。** V5 摘要 ~350 词、含十余个
数字（6960、424/472、+2.12、+0.63、调用数、runtime），像 rebuttal
数据表；结论同病。按 V6 要求压缩至 180–230 词，卖点排序改为：
信息类定理 > 包络 > 边界 > 互补认证 > 经验增益 > runtime。

### Major

**A4. 六标量信息模型的动机缺失。** 无"Why a fixed-size information
interface"小节；读者会认为信息类是为定理人造的。需新增小节（V6 §2.3）
+ 信息层级框架（Level 0–4）。

**A5. 六标量无物理释义。** 需新增"Physical interpretation"小节（V6
§2.4），并说明谱兼容性=三量同源、不可同时各自取最坏——这正是旧界
保守的来源。

**A6. 边界小节缺直觉铺垫。** β₀=α 的特殊性（期望方向达到算子范数
⇒ 奇异向量等式条件强制 d ⊥ Cᴴ(v⊥) ⇒ 附加相位耦合）应先给直觉再进
定理；0.59 vs −0.21 需配"信息类切换"示意图（非退化类极限 −0.21 ↔
边界类 0.59）。

**A7. 规模实验混杂因子。** L 配置（3 用户/小区）每流功率 1/3 与
S/M（1/2）不同，接受率跨配置不可比（V5 已声明但仍并列在一张表）。
需按 F-B 重跑受控版本。

**A8. ρ 与工程 CSI 精度无对应。** ρ=0.06/0.15/0.30 的物理意义未建立
（相对半径 vs dB NMSE）。需按 F-A 补解释性实验。

**A9. Related Work 缺整节。** 仅引言带过 4 篇。V6 需独立小节 + 对照
表（见 B 部分）。

**A10. 证书术语需一句精确区分。** 增加："The analytical inequality
is exact; the reported implementation uses floating-point arithmetic with
conservative numerical safeguards rather than formally verified directed
rounding." 关键案例（翻转、边界、近阈值）可选做区间算术/高精度外向
舍入 audit（F-C），否则保持限制措辞。

**A11. runtime 主导痕迹残留。** 摘要仍有 "twice as fast"；2× 须全部
限定为 "against the original fixed-iteration reference implementation"，
主表用容差停止参考，200 次版本降为补充。

### Moderate

**A12. 防御性重复。** "not claimed as / does not establish" 在 V5 出现
约 15 处；按"每个限制只在其最合适位置讲一次"合并（S-procedure 非创
新→Related Work 一次；runtime 不普适→计时节+局限两次；有限实验不证
安全→数值范围节一次；查新待完成→Related Work/Discussion 一次）。

**A13. 缺架构图与信息层级图。** 新增 Fig. 1（certification layer 流程
图：候选→全几何/六标量两路→max→正即认证→否则 TRS 回退）与信息层级
示意（文字框图即可）。

**A14. L_quad 引入位置过晚**（§7.4 才正式出现）：放入信息层级 Level 4
元素，在 §2.3 一并引出。

### Minor

A15 标题候选与重评（见 E）；A16 Table 5 星标脚注排版；A17 图注风格
统一（"diagnostics vs certified"已统一，检查全文）；A18 参考文献条目
核对（[1][2][3][4] 已实，新增条目需 DOI 核对）；A19 关键词更新；
A20 §7 各小节编号在重组后需连续。

### 定理审稿式检查（要点版）

| 项 | Thm 1 | Prop 1 | Thm 2 | Thm 3 |
|---|---|---|---|---|
| 假设 | (3) 兼容 | 任意兼容实例 | n≥5,m≥3, β>0, β₀>α, p,γ>0, r≥0, a,I₀≥0 | β₀=α>0, β>0, I₀≥αa, n,m≥2 |
| 隐藏假设 | 无（(3) 已含 δ≤β, ν≤α） | 无 | ε∈(0,1)；S(0)<β₀² 严格性已证；满行秩由 Mε≻0 | 边界类非空条件已列 |
| 退化情形 | α=0/β=0/端点/零算子均已归位 | x>a 时 √p\|a−x\| 已处理 | α=0、r=0 单独构造；零填充 | β=0、r=0 单独构造 |
| 复/实一致 | 复 Hermitian 语境，实为特例 | 同 | 构造全实坐标（相位显式） | 相位耦合显式 |
| inf vs min | max 在紧区间达到 | — | **inf，不必达到**（‖ĥ‖→∞ 族）——文中已声明，V6 需在定理后加一句强调 | **达到**（显式构造） |
| 紧性/连续 | 连续函数紧区间 | 逐点不等式 | 1-D min 紧性；ε→0 取极限 | 1-D min 紧性 |
| 数值含义 | 闭式直接实现 | 区间下界安全 | 构造达界（机器精度核查） | 构造达界（精确） |
| 待独立复核 | 无 | 无 | 一般 m≥3 秩/象空间核对建议由独立数学复核确认（流程项，非已知缺口） | 无 |

**Proof issue / required verification（不擅自修补，如实列出）**：
(i) Thm 2 中 t_ε>0 的严格性论证依赖"等号迫使最大特征向量 ⟂ ℓ₂"——
逻辑闭合，但建议投稿前请独立第三方复核该步；(ii) 一般 m≥3 的零填充
保持六摘要与达界误差范数的逐条引理化（当前为正文叙述）；(iii) Thm 3
的 F_c 分支切换与 min 的联合取法可在行文中再明确一句。均非已知错误。

---

## B. Novelty landscape（文献定位表）

| Reference | Problem | 信息量 | 保留信息 | 界类型 | 共享球预算？ | 谱兼容？ | 信息类 sharp？ | 边界？ | 领域 | 与本文关系 |
|---|---|---|---|---|---|---|---|---|---|---|
| Zhou et al., TSP 2020 [1] | IRS 鲁棒传输设计 | 完整 | 完整矩阵 | 精确(S-proc/SDP) | 是（单球） | 否（不需要） | 否 | 否 | 无线 | 全信息设计代表；本文认证层用其同型精确求解作参考 |
| Zheng et al., TCOMM 2025 [4] | 多功能 RIS 能效鲁棒设计 | 完整 | 完整矩阵 | 精确(S-proc LMI) | 是 | 否 | 否 | 否 | 无线 | 同上；本文 SDP 基线的经典技术来源 |
| Zhang–Ng–Schober 线（如 Sun et al. 2018） | norm-bounded worst-case SINR | 完整 | 完整 | 精确/保守界 | 是 | 否 | 否 | 否 | 无线 | 平行/正交误差分解的常见出处；旧界 Q_old 的文献近亲 |
| Wang et al.（概率 SINR 约束系列） | 概率鲁棒波束 | 完整 | 分布/矩 | 保守近似 | 否 | 否 | 否 | 否 | 无线 | 对照：本文确定性最坏情形、压缩信息 |
| Rendl–Wolkowicz 1997; Fortin–Wolkowicz 2004 [2,3] | TRS 精确求解 | 完整 | 完整 | 精确 | — | — | — | hard case | 优化 | 本文参考求解器与 KKT 框架 |
| Grone–Johnson–Sá–Wolkowicz 1984 | 部分矩阵 PSD 补全 | 压缩 | 指定元素 | 可行性刻画 | — | — | — | — | 数学 | 本文 2×2 Gram 补全是其弦图特例 |
| Parrott 1978；"A strong Parrott theorem"；arXiv 2506.10485 | 给定列/子空间作用的收缩补全 | 压缩 | 指定算子值 | 集合刻画（算子球） | — | 部分 | — | Parrott 反例(三维以上) | 数学 | Thm 1 极值的经典框架；本文是其在 (α,β,β₀;x,y) 组合下的显式求值 |
| Buzano 1974；Bottazzi–Conde 2023 | 修正 Cauchy–Schwarz \|⟨Tx,y⟩\| | 压缩 | 范数标量 | 双边不等式 | — | 部分 | 否 | 否 | 数学 | 式(7)交叉项界的直接亲缘；必须引用 |
| Micchelli–Rivlin 1977；Melkman–Micchelli 1979；Foucart 2024 | 最优恢复 / radius of information | 压缩 | 有限数据泛函 | 最坏误差半径 | — | — | 是（一般框架） | — | 数学/计算 | Thm 2 概念等价框架；本文给出特定六标量 QoS 问题的闭式半径 |
| 本文 | 固定候选认证 | **压缩** | **六标量** | **安全下界+类最优** | 是 | **是** | **是（非退化类）** | **是** | 无线+数学 | — |

**分层判定**：
- Thm 1 属类别 **C**（已有数学框架的新 specialization）：收缩补全 +
  Buzano 型交叉项 + 初等凹最大化，全部经典；在"恰好保留 (α,β,β₀) 且
  误差预算 x²+y²≤r²"的组合下的三段闭式，**未检索到同型结果**。
- Thm 2：radius-of-information 框架的既有概念；本文贡献=该特定问题
  的**闭式半径求值 + 显式极值族（含 ‖ĥ‖→∞ 与边界跳变现象）**。
- Thm 3：奇异向量等式条件标准；相位耦合包络 + 0.59/−0.21 类切换
  **未检索到先例**。
- 禁用 first/novel 措辞；统一限定式表述。

**本文比已有结果多解决的一步**：在"压缩为六个标量谱摘要"这一固定
信息预算下，给出了从方向极值到 QoS 证书、再到信息类最优性与退化
边界分类的**完整闭环**——文献中各层均有先例，闭环没有。

---

## C. Revised positioning（四句话）

- **Problem**：对已固定的 beamforming+RIS 候选配置，当认证接口只保留
  六个标量谱摘要时，能对球约束最坏 QoS 给出多强的一致证书？
- **Main theoretical contribution**：显式三段谱兼容方向极值 → 安全且
  支配旧共享预算界的 QoS 包络，并在非退化信息类上 sharp；退化谱边界
  需单独的相位耦合包络（0.59 vs −0.21 的类切换）。
- **Engineering relevance**：该证书与全几何二次界无支配关系、可安全
  取 max 组合，适合作为 cheap bound → six-scalar bound → exact
  fallback 分级认证接口的中间层。
- **Limitation**：sharpness 只对所述无约束信息类成立（极值族需大
  ‖ĥ‖）；经验增益小且批次相关；不替代全信息精确求解。

---

## D. Section-by-section modification plan（V5 → V6）

| 节 | 保留 | 删除/降级 | 移动 | 新增 | 理由 |
|---|---|---|---|---|---|
| Title | — | — | — | 6 候选+分析（E 部分） | 现题过抽象 |
| Abstract | 骨架 | 全部数字细节、SDP 实现细节、runtime 细节 | 数字入正文 | 180–230 词新版（E） | A3 |
| §1 Intro | P1/P2 骨架 | 防御句 | — | 5 段式重写（E）；信息层级一句预告 | A3/A12 |
| 新 §1.x Contributions | — | 5 点列表 | — | 压缩为 3 点（E） | 用户规范 |
| §2 模型 | 2.1/2.2 全部 | — | — | **§2.3 Why a fixed-size information interface**（动机 7 场景+fixed budget 定位+层级 Level 0–4+L_quad 预引出）；**§2.4 Physical interpretation**（六标量无线释义+谱兼容=保守性来源） | A4/A5/A14 |
| §3 包络 | 全部（含证明） | — | — | — | 理论主线不动 |
| §4.1/4.2 | 全部 | — | — | Thm 2 后加"inf 不必达到"强调句；证明状态脚注（独立复核中） | 审稿检查项 |
| §4.3/§5 边界 | 全部数学 | — | — | 直觉铺垫段（等式条件⇒正交⇒相位耦合）+ 类切换示意小图；0.59/−0.21 叙述改为"不同信息类的下确界变化" | A6 |
| §5 七标量 | 全部 | — | — | 定位为 additional-information example；七标量包络=Future Work | 用户规范 |
| 新 Related Work（§2 前或后） | — | — | — | 独立小节+对照表（B 部分）+旧界归属声明 | A1/A2/A9 |
| §6 认证接口 | 6.1/6.2 主体 | 乘子迭代细节降级 | 6.3 SDP 判定规则细节→附录 | 一句 analytical vs floating-point 区分（A10） | A10/A11 |
| §7 实验 | 主体结果 | — | 电池/哈希/历史说明→附录 | 重组为 7.1–7.7（下表）；规模表换 F-B 数据；ρ-dB 对应段（F-A） | A7/A8 |
| §8 lifecycle | 全部（压半） | — | — | — | 次要 |
| §9 Discussion | 要点 | 重复防御句 | — | 三层贡献分离重写（E） | A12 |
| §10 Conclusion | — | 全部数字复述 | — | 三层式重写（E） | 用户规范 |
| References | [1–4] | — | — | +Buzano/Bottazzi–Conde、Grone et al.、Parrott(+强化版)、Micchelli–Rivlin(+Melkman)、Sun et al. 代表作 | A1 |
| 附录/补充 | — | — | 数值电池、配置哈希、v1→v2 更正史、200 次计时表 | — | A12 |

**§7 重组映射**（问题导向替代代码史导向）：
7.1 Analytical implementation checks（定理 1/构造/边界核验）→
7.2 Controlled threshold-crossing examples（谱扫描=确定性见证）→
7.3 Fixed-candidate certification study（易池+敏感性合并）→
7.4 Independent validation（预注册批，最重要）→
7.5 Complementary certificates and staged screening（联合+分级+规模）→
7.6 Independent exact cross-validation（SDP）→
7.7 Timing and implementation sensitivity（容差参考为主）。

---

## F. Experiment TODO（仅确实需重跑者）

**F-A. Realistic CSI-error interpretation**
- Purpose：给 ρ 建立工程含义；验证结论在较真实链路缩放下仍出现。
- Independent variable：CSI 误差水平（相对半径 ρ∈{0.06,0.15,0.30} 与
  NMSE −20/−15/−10 dB 的对应说明；二者只做口径映射，不硬等价）。
- Controlled variables：候选生成规则、γ、接受阈值、参考实现。
- Setup：(i) 口径段——相对球 r=ρ‖ĥ‖ 的 dB 解释；(ii) 小型现实链路池：
  2 BS、RIS 16、距离相关路损、BS–RIS Rician、RIS–UE Rician/Rayleigh、
  有界归一化 CSI 误差；每档 ρ 重算六标量与证书。
- Metrics：各档直接认证数、旧界/新界/联合覆盖、翻转数；与主池结论
  的方向一致性。
- Expected outputs：一张 ρ–dB 对照表 + 一张现实链路小表（或图）。
- Supports：结论对现实缩放稳健；ρ 的工程可读性。
- Cannot support：任何统计覆盖概率声明；增益幅度外推。

**F-B. Controlled scaling experiment（替换现 Table 6）**
- Purpose：维度增大时，标量界评估 vs 精确参考的计算伸缩。
- Independent variable：堆叠维度 dim ∈ {4, 8, 16, 32}。
- Controlled variables（全部固定）：用户/天线比、每基站总功率=1、
  每流功率=1/用户数（各档一致的用户数）、γ、噪声、信道归一化、
  候选生成规则、ρ 归一化、τ 规则。
- Pipelines：容差停止 TR 纯参考、Q_G+回退、L_quad+回退、选定级联。
- Metrics：单次界评估成本、单次参考成本、参考调用减少、联合覆盖
  （同配置内比较）；逐种子中位数与离散度。
- Expected outputs：替换 Table 6（cost-vs-dim 表或图）。
- Supports：筛选相对成本的维度趋势（限本实现）。
- Cannot support：跨配置接受率比较；复杂度阶主张。

**F-C.（可选）关键案例的高精度/区间算术 audit**
- Purpose：把"analytical vs floating-point"区分落到实处。
- Setup：仅对确定性翻转例、边界四例、验证批近阈值单元（|q\*|≤10τ）
  用 mpmath 高精度+外向舍入重算区间下界与参考值。
- Metrics：双精度结果与高精度下界的偏差；是否存在符号翻转。
- Supports：工程结论在关键案例的数值稳健性。
- Cannot support：全体验证的形式化证书（除非全面区间化，明确不做）。

**不新增**：更多随机种子批次；七标量包络推导；GLTR。

---

## 检索来源（2026-10-07）

- [Buzano 不等式与数值半径应用（Survey of trace inequalities; q-numerical radius 等）](https://ajmaa.org)、[Buzano 原始文献线索（Rend. Sem. Mat. Torino 31）](https://www.ele-math.com)
- [Grone–Johnson–Sá–Wolkowicz 1984 PSD 补全（弦图定理）](https://scholar.google.com.tr/citations?user=WmTBubUAAAAJ&hl=de)、[Wolkowicz《Matrix completion problems》手册章](https://backend.orbit.dtu.dk)
- [Sarason 算子延拓视角综述](https://library.slmath.org/books/Book33/files/sarason.pdf)、[Dritschel–Rovnyak Kreĭn 空间收缩延拓](https://uva.theopenscholar.com/files/james-rovnyak/files/p38.pdf)、[A strong Parrott theorem](https://www.researchgate.net/publication/243064171_A_strong_Parrott_theorem)、[arXiv:2506.10485（列收缩补全的精确 Parrott 型结果）](https://arxiv.org/pdf/2506.10485)
- [Micchelli–Rivlin《A Survey of Optimal Recovery》1977](https://www.semanticscholar.org/paper/A-Survey-of-Optimal-Recovery-Micchelli-Rivlin/17cc6ad5cfdbf5a0f5fd9a5c0c892a0550fa265c)、[Melkman–Micchelli 型线性算子最优估计线索](https://link.springer.com/article/10.1007/s10444-004-7616-1)、[Foucart 2024 radius of information](https://foucart.github.io/publi/OR_L1Noise_v3.pdf)、[SIAM 多保真 worst-case 学习（Chebyshev 半径=QP 平方根）](https://epubs.siam.org/doi/10.1137/24M1671025)
- [Zheng et al. 2025 TCOMM（IEEE Xplore）](https://ieeexplore.ieee.org)、[arXiv:2412.08032](https://arxiv.org/abs/2412.08032)
- [Zhou et al. 2020 TSP（arXiv:2001.07054）](https://www.arxiv.org/abs/2001.07054)
- 鲁棒 SINR/方向分解线索：[Sun et al. 2018 NOMA 认知无线电鲁棒设计（平行/正交误差分解）](https://eprints.soton.ac.uk)、[Yao et al. 2023 RIS cell-free 鲁棒波束](https://arxiv.org)、[概率 SINR 约束系列（Wang et al., NTHU）](https://www.ee.nthu.edu.tw)
