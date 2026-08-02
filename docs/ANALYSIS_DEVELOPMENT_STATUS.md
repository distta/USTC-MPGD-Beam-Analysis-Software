# BeamAnalysis 分析逻辑与后续开发基线

本文记录当前代码已经实现的分析逻辑、run346 的可复现基线，以及下一阶段仍需解决的问题。修改重建或标定逻辑前，应先按文末的回归流程保存对照结果，避免把几何变化、匹配变化和算法变化混在一次测试中。

## 1. 数据流和事件标识

分析主链为：原始数据转换 → `HitProcessor` → `ClusterBuilder` → `ClusterReconstructor` → `TrackAnalysis` → `DUTAnalysis` / `PadDUTAnalysis` / `TimeResolution`。

- `eventID` 是处理后 ROOT 文件中的事件序号，用来回读同一重建事件。
- `rawEventID` 是采集系统事件号，只用于与示波器外部 T0 对应，可能回绕或重复。
- `trackIndex` 区分同一重建事件内的多条径迹。
- 时间分析内部使用 `(eventID, trackIndex)` 作为径迹唯一键，绝不以 `rawEventID` 合并径迹。
- 同一个 `rawEventID` 对应多个重建事件或多个示波器条目时，该 ID 被标记为歧义并从外部 T0 匹配中排除。

这个区分必须保留。后续若改变 ROOT tree schema，需要同时检查 `TrackInfo.root/TrackerValidation`、`DUTInfo.root/DUTTree` 和 `TimeResolution.root/TimingTree` 的关联方式。

## 2. Hit、cluster 和电荷定义

### HitProcessor

`HitProcessor` 同时支持波形输入和 VMM 已给出的 `adc/hit_time_ns`。波形时间由探测器或脚本配置的拟合模式产生，并保存 `time`、`timeError`、上升时间和过阈宽度等信息。VMM 的 `hit_time_ns` 是相对 trigger 的时间。

### ClusterBuilder

- 相邻条带按 readout plane 分组；VMM 可用 `timeWindowNs` 限制相邻条带时间差。
- `cluster.size` 是实际参与的唯一条带数，`cluster.range` 是最大、最小 strip ID 之间的跨度。
- `cluster.charge` 的正式定义是 cluster 内所有 hit 的 ADC 幅度 `strip.amp` 之和。
- `cluster.maxAmp` 是最大单条 ADC 幅度。
- 电荷重心仍使用 `strip.charge` 加权 strip ID；若总电荷无效则退回首尾 strip 中点。
- `cluster.time` 在 builder 阶段只是最早条带时间；时间分辨率脚本会根据配置重新从选中 cluster 的条带计算时间。

## 3. TrackAnalysis

- 正式模式要求三个或更多 Tracker；此前试验性的 two-tracker 概算模式已经删除。
- 初始对齐样本要求每个 Tracker 恰好一个 `LocalHit`。
- 最终寻迹允许多击中，通过空间 gate、`chi2/ndf` 和 hit 冲突消解选择径迹。
- 性能统计使用最终选中的全部径迹，以 leave-one-out 残差和联合直线拟合协方差计算单层分辨率、pointing resolution 和角分辨率。
- `useEstimatedResolution=false` 时保留配置分辨率，不使用运行中估计值。
- 径迹 T0 写入 `Tracks.t0/hasT0`。每个命中探测器先把所选 X/Y clusters 的有效条带按 `amplitude²` 加权，然后对各探测器时间等权平均。这个 T0 用作 DUT uTPC 输入；它与 `TimeResolution.timeMethod` 的可配置时间是两条独立路径，后续应消除这处重复实现。

## 4. DUTAnalysis

### 匹配与效率

- 空间匹配采用 cluster envelope：track 预测交点落在实际最小/最大 strip 中心向外扩展 `margin` 的范围内即匹配。
- X、Y 和 2D 效率均按总事件数加权；低统计 bin 和显式排除 bin 不参与汇总。
- fake efficiency 使用不同事件的径迹和未被本事件径迹使用的 clusters 配对，作为独立假匹配概率，不从效率中扣除。
- DUT tree 保留 hit、ADC、cluster charge、cluster size、位置、残差和匹配标志等逐事件信息。
- `HitMaps/` 输出全部重建 hit 与 track 匹配 hit 的细粒度二维图；X/Y profile 合并在最后一张 profile canvas 中。

### DUT 几何对齐

几何对齐必须使用电荷重心而不是 uTPC 位置，否则时间标定产生的位置斜率会被旋转参数吸收，导致“uTPC 对齐后 CC 歪掉”。当前实现：

1. 每次试探几何都从所有 clusters 重新匹配，避免一开始选错 cluster 后永久锁死。
2. 以软匹配从宽到严逐级优化，并要求最佳候选相对次佳候选有足够 dominance。
3. 在正式最小化前扫描 `rotX` 初值。
4. `dz` 固定，因为当前束流角度跨度不足以与平移独立约束。
5. 最终性能仍使用正式 cluster-envelope 选择，不使用对齐阶段的软权重。

仍需补充一个独立的 alignment validation tree，保存每阶段候选数量、dominance 和被拒原因。目前只有最终几何及汇总，难以审计局部极小值。

## 5. TimeResolution

### 时间计算方法

`TimeResolution.config.timeMethod` 是时间算法的唯一用户入口，支持：

- `XMean`、`YMean`：对应面的条带时间按 `1/timeError²` 加权；无有效误差时普通平均。
- `XMin`、`YMin`：对应面的最早时间。
- `XYMean`：X、Y 两面全部有效条带按 `1/timeError²` 联合加权。
- `XYMin`：两面全部有效条带的最早时间。
- `XYWeighted`：两面全部有效条带按 `amplitude^p` 加权，`p` 由 `amplitudeWeightPower` 配置。

时间方法不再写进 `ClusterBuilder` 或每个 detector 的 reconstruction 配置，避免同一个 cluster 定义被空间重建和时间分析重复解释。

注意：这里的“唯一入口”只指 `TimeResolution` 脚本自身。`TrackAnalysis` 写入的 track T0 当前仍固定采用 A² 条带权重，不能由该字段配置。

### 无外部 T0

- 三个 Tracker 必须在同一个 `(eventID, trackIndex)` 上都有有效时间。
- 对三组两两时间差做单高斯拟合，由三个方差方程解出各 Tracker 分辨率。
- track time 按各 Tracker 的逆方差组合。
- DUT 分辨率由 `sigma(DUT-track)^2 - sigma(track)^2` 开方得到。

### 有外部 T0

- 外部 T0 来自 `OscilloscopeAnalysis.root/Events.referenceTime`，按无歧义 `rawEventID` 匹配。
- 每根条保存 `t_strip-T0`、幅度、上升时间、过阈宽度等到 `StripTimingTree`，可用于独立诊断。
- 对每个探测器拟合单条 `t_strip-T0` 与幅度的三次关系，逐条修正后再按所选 `timeMethod` 重算 cluster time。
- 默认结果直方图使用修正后的时间；修正参数写入 `TimeAmplitudeCorrections`，不单独输出拟合曲线 canvas。
- `TimeAmplitude/` 保留每根条的时间—幅度二维关系及 profile，用来验证 time-walk 是否真实存在。
- 外部 T0 的已知误差按逐事件 `timeError` 进入结果，并从测得宽度中扣除；不可再把它当成零分辨理想参考。
- DUT timing efficiency 的分母为“有 T0 且 track 进入有效区”的径迹，空间未匹配、无有效 DUT 时间和时间窗外事件仍保留在分母。残差拟合图上画出最终选择的时间窗。
- `maximumStripAmplitude > 0` 时，任一 cluster 条带超过阈值的事件不进入外部 T0 分辨率与 timing-efficiency 计算；负值表示禁用。

## 6. rawUTPC 重建

当前 `rawUTPC` 只在存在 track T0 且 cluster 至少三条时启用，否则退回电荷重心。对第 `i` 根条：

```text
deltaT_i = strip.time_i - trackT0 + timeOffset
           + stripTimeSlope * (stripID_i - centerStrip)
drift_i  = deltaT_i * driftVelocity
stripX_i = stripID_i * pitch
```

用 `(stripX_i, drift_i)` 做直线拟合。条带纵向误差为 `timeError_i * driftVelocity`；同时加入位于半气隙深度的电荷重心约束点，其横向误差由 `chargeCentroidXError` 配置。拟合迭代把 X 误差投影到纵向，最终取直线与读出参考面的交点作为位置。

`timeOffset`、`driftVelocity`、`stripTimeSlope` 和几何 `rotX` 存在明显相关性，不能一次全部自由扫描后只报告最小残差。正确标定顺序应是：

1. 用电荷重心锁定 DUT 几何和 cluster 匹配集合。
2. 固定几何，以 track time 标定 DUT–Tracker 时间 offset。
3. 用漂移场/气隙的物理范围约束漂移速度。
4. 扫描 `stripTimeSlope`，同时检查残差均值、`resY vs predY` 斜率、cluster-size 分段和效率。
5. 只用位置平移消除剩余常数偏移，不再改变旋转。

## 7. run346 可复现基线

2026-08-02 从 uTPC 斜率扫描阶段的 VS Code 历史配置恢复并在临时目录验证，未覆盖正式 `result/346/DUTInfo.root`。当前 `data/Micro-TPC_40deg/config.json` 已恢复为：

```json
"position": [0.24178, 8.01503, 391.63505],
"rotation": [-0.63297, 0.01572, -0.00875],
"Y pitch": 0.399,
"driftVelocity": 0.02,
"timeOffset": 8.9,
"stripTimeSlope": 1.0,
"chargeCentroidXError": 0.3
```

当前代码、现有 `TrackInfo.root` 和上述配置的临时验证结果：

| 指标 | run346 |
|---|---:|
| Y residual mean | 1.67 µm |
| Y residual sigma68 | 125.33 µm |
| X residual sigma68 | 68.36 µm |
| 2D efficiency | 96.77 ± 0.18% |
| fake-match probability | 0.047% |

原始历史扫描点在未做最后 Y 平移时为 `sigma68 = 127.71 µm`、单高斯核心 `112.24 µm`，但均值约 `-0.776 mm`。因此“约 120 µm”确实来自当时的 `stripTimeSlope=1.0` 与几何联合标定，而不是 HitProcessor 的 Fit 模式。最后只调整 detector Y 平移，把均值移到零附近，并保留宽度。

这组数值是当前软件回归基线，不等价于已经完成物理唯一标定。特别是它和后来独立 CC 对齐几何并不相同，说明 uTPC 标定与几何仍有退化；发表或跨 run 比较前必须按上一节顺序重新拆分。

## 8. 后续开发优先级

### P0：保证结果可复现

- 为每次参数扫描写出独立目录和完整配置快照，不覆盖正式 ROOT 文件。
- 增加机器可读的 run summary：输入文件指纹、git revision、配置、匹配数、均值、sigma68、Gaussian sigma、效率。
- 建立 run304 时间分辨率和 run346 空间分辨率两个自动回归样本。
- 统一文档与实际输出：ROOT 目录名、T0 误差扣除方式和幅度修正是否生效必须由测试验证。

### P1：拆开几何与 uTPC 标定

- 先冻结 CC 几何，再生成固定的高纯匹配样本；扫描过程只更新时间参数。
- 对不同 `predY` 区域、cluster size、最大条幅度和条带号分别画残差，确认多峰是否来自区域、坏道或错误匹配。
- 对 20°、40° 和不同漂移场 run 比较 `stripTimeSlope` 与漂移速度；电子学项应跨 run 接近，随漂移场变化的项应符合漂移物理。
- 加入训练/验证样本划分，禁止在同一批事件上既拟合 time-walk 又报告最终分辨率。

### P2：大 cluster 与系统误差

- 比较大 cluster 中去掉边缘低电荷条、鲁棒直线拟合和时间误差重标定，但每种算法必须同时报告效率损失。
- 使用 `StripTimingTree` 检验上升时间和过阈宽度；除非在独立验证样本优于纯幅度修正，否则不进入默认重建。
- 建立 channel-level 坏道/时间 offset 图，区分局部坏道、多区域增益和全局 time-walk。
- 将外部 T0 分辨率、tracker pointing error、cluster matching 假匹配分别作为系统误差报告。

## 9. 安全回归流程

1. 复制数据集配置到 `/tmp/<scan-name>/config.json`。
2. 复用只读的 processed data 和已有 TrackInfo；在临时目录运行 DUT/TimeResolution。
3. 保存每个扫描点的配置、终端摘要和 ROOT 指标，不覆盖 `<base>/result/<run>/`。
4. 一次只改变一类参数：几何、匹配、时间标定或重建算法。
5. 至少比较均值、sigma68、Gaussian 核心、残差—位置斜率、效率、假匹配率和各 cluster-size 分段。
6. 只有验证样本也改善后，才把参数或代码写回正式配置。
