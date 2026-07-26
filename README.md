# BeamAnalysis

BeamAnalysis 将数据转换、探测器初始化、ROOT 数据读取和分析脚本统一到一个非交互式运行流程中。正式调用格式只有一种：

```bash
./BeamAnalysis <base_dir> <run_id>
```

例如：

```bash
./BeamAnalysis /data/beam 1591
```

不提供 `run`、`convert`、`validate` 或 `display` 子命令。事件显示由配置中的 `mode` 选择。
`run_id` 按现有 run 命名约定执行严格十进制数字校验。

## 构建

依赖：

- 支持 C++17 的编译器
- CMake 3.14 或更高版本
- ROOT（包含 Core、RIO、Tree、Hist、Graf、Gpad、MathCore）
- Eigen3
- nlohmann/json 头文件

构建命令：

```bash
cmake -S . -B build
cmake --build build -j
```

生成的程序为 `build/BeamAnalysis`。也可以将其复制或链接到项目根目录后按上面的正式格式运行。

## 配置

完整分析示例见 `config/example_analysis.json`，事件显示示例见 `config/example_display.json`。使用时将所选示例保存为 `<baseDir>/config.json`。核心结构如下：

```json
{
  "mode": "analysis",
  "conversion": {
    "enabled": true,
    "overwrite": false,
    "type": "APV25SRS"
  },
  "detectors": [],
  "scripts": []
}
```

第一个命令行参数是数据集根目录，第二个参数是 run ID。程序固定使用以下结构：

```text
<baseDir>/
├── config.json
├── geometry.json
├── raw/
│   └── run<run_id>.root
├── processed/
│   └── run<run_id>.root
└── result/
    └── <run_id>/
```

`BTAPVDat` 转换器的原始输入扩展名为 `.dat`，pedestal 默认读取 `<baseDir>/raw/run<run_id>_pedestal.dat`。转换器的 channel map 等专用参数仍写在 `conversion` 中。

程序会自动创建缺失的 `raw/`、`processed/` 和 `result/<run_id>/` 目录。配置默认读取 `<baseDir>/config.json`，并兼容旧位置 `<baseDir>/raw/config.json`；VMM3a 的 geometry 默认读取 `<baseDir>/geometry.json`。需要读取的输入文件不存在时直接退出。

`mode` 支持：

- `analysis`：按数组顺序运行所有 `enabled: true` 的脚本，任一脚本失败即停止。
- `eventDisplay`：不运行脚本，读取输出目录中的 `TrackInfo.root` 并启动交互式事件显示。

`TrackAnalysis.config.runAlignment` 控制 tracker alignment；分析模式不会询问 y/n。旧字段 `performAlignment` 仍可读取，但新配置应使用 `runAlignment`。

TrackAnalysis 先使用每个 Tracker 都恰好有一个 `LocalHit` 的事件做 alignment 和可选的初始 hit 分辨率估计，再用得到的 X/Y hit 分辨率扫描全部事件，包括多击中事件，通过 gate、`chi2/ndf` 和 hit 冲突消解选出径迹。设置 `useEstimatedResolution=false` 可以保留配置中的 `resolutionX`/`resolutionY`。最终 `Performance` 会清除初始 single-hit 统计，并使用所有最终选中的径迹（包括多击中事件中的径迹）计算每层 leave-one-out 残差及合并后的等效单层残差；这些分布使用单高斯拟合的 `sigma`，再由联合直线拟合协方差计算 pointing resolution。`TrackInfo.root/Performance/Resolution` 中包含单高斯拟合图、随 z 变化的 X/Y tracking 分辨率曲线、各 Tracker/DUT 位置的 pointing 分辨率、共同单层分辨率和角度分辨率，长度单位为微米、角度单位为微弧度。

Tracker alignment 的 `alignmentIterations` 只是最大迭代数。参数更新连续稳定，或相对 loss 改善连续 `alignmentConvergencePatience` 轮低于 `alignmentRelativeLossTolerance` 时会提前结束；达到最大迭代数仍不满足条件会报告未收敛。

`TimeResolution` 默认使用 detector 配置中恰好三个 `role: "Tracker"` 的探测器，也可以在脚本配置中用 `"trackerIDs": [id1, id2, id3]` 明确指定。每个探测器使用 `XYMean` 时间，即与该 track 匹配的 X/Y clusters 所包含有效读出通道拟合时间的逆方差加权平均值；Pad DUT 同样使用匹配 cluster 内有效通道拟合时间按 `1/timeError²` 加权。没有任何有效正误差时才退回普通平均。`timingWaveform` 可显式配置时间拟合的 `mode`、`cfdFraction` 等 `HitProcessor` 参数。三组时间差只使用三个 tracker 时间都有效的共同事件样本，以单高斯拟合宽度解出各层分辨率，再按 `1/sigma_i^2` 构造加权 `trackTime`。`TimeResolution.root/TrackTimes` 保存三层时间、权重和加权 track 时间；启用 `analyzeDUTTiming` 后，`DUTTimes` 和 `DUTTimeResolution/` 保存 DUT 时间差及扣除 track reference 后的 DUT 分辨率。设置 `applyDUTTimeWalkCorrection: true` 后，程序用最大通道幅度的二次模型做两折交叉验证；只有总样本至少 0.5% 改善且两个独立验证折都变窄时才采用校正，否则自动保留原始结果。原始/校正 residual、验收结果、两折宽度和模型参数均写入 ROOT。`DUTTimeResolution/DUT_<id>/Correlations/` 保存 DUT-track residual 相对于最大/平均幅度、cluster charge、size、centroid、local X/Y 及 track predicted X/Y 的二维分布、profile 和 Pearson 系数。

`OscilloscopeAnalysis` 独立分析 LeCroy 分段 CSV。未配置目录时自动读取 `<base>/scope/<runID>/` 下的 C1–C4 默认文件名：C4 按负逻辑串行帧 `start(1) + 16-bit eventID + stop(0)` 解码，位周期为 25 ns；C1–C3 使用负脉冲 CFD 定时，并从 C4 帧起点向前选择最近的信号，忽略更早的堆积和 eventID 之后的脉冲。默认输出 `<base>/result/<runID>/OscilloscopeAnalysis.root`，其中 `OscilloscopeEvents` 保存逐 segment 的 eventID、三路时间、幅度、baseline、CFD 阈值和候选脉冲数；`Timing/` 保存三组时间差，`TimingResolution` 保存高斯宽度及由三探测器方程解出的单路时间分辨率。必要时可覆盖 `csvDirectory`、`outputFile`、`cfdFraction`、`minPulseAmplitude`、`histogramBins` 和 `maxWaveformFiles`。

`PadDUTAnalysis` 的图形输出位于 `PadDUTInfo.root/DUT_<id>/`。`Residuals/` 保存 X/Y 残差分布、二维残差，以及 `resX/resY` 分别相对于 track 预测交点 `HitX/HitY` 的全部四种相关分布；X/Y 分布均使用 100 bins，显示范围为样本 `mean ± 5 RMS`，在确定合适模型前不做高斯拟合。`HitMaps/` 分别保存全部重建 hit 和 track 匹配 hit 的二维分布。`Efficiency/` 保存 X/Y/2D 三张二维效率 map、2D 效率沿 X/Y 的两个投影和一张 margin 扫描图；`Fake/` 只保存三张 map 和两个投影，不做 margin 扫描。输出不再同时写入内容重复的 histogram、graph 和 canvas。所有二维图默认使用 `COLZ TEXT`，Fake map 和投影的显示范围固定为 0–5%。

`DUTAnalysis` 只使用 cluster-envelope 定义 X/Y/2D 效率。平均效率统一按总事件数加权，即所有有效 bin 的总命中数除以总事件数，不再计算各 bin 等权平均。`margin` 是 cluster 最小/最大边缘条带中心向外扩展的距离，不叠加 tracking sigma 或其他参数。`minEntriesPerBin` 设置一个 bin 参与汇总所需的最小事件数；低于阈值、被 `excludeXBins`/`excludeYBins` 排除或为空的 bin 不参与汇总。

对于同时具有有效 X、Y selected clusters 的平面 DUT，`DUTAnalysis` 使用两面
全部有效 ChannelHit 时间的联合均值作为 `XYMean`。VMM 的 ChannelHit 时间
已经是相对 trigger 的 `hit_time_ns`，因此
`DUTInfo.root/DUT_<id>/TimeResolution/` 直接保存 `hXYMeanTime`、单高斯
拟合画布 `cXYMeanTimeResolution` 和 `TimingSummary`。

效率不均匀度定义为有效空间 bin 效率的变异系数 `sigma(efficiency_i) / mean(efficiency_i)`，各有效 bin 等权，分别计算 X、Y 和 2D。结果写入终端汇总、效率 map 和 `EfficiencySummary` 的 `nonuniformityX/Y/2D` 分支。

核心配置示例：

```json
"efficiencyMap": {
  "minEntriesPerBin": 30,
  "margin": 0.6,
  "marginScan": { "min": 0.0, "max": 2.0, "step": 0.05 },
  "fake": { "enabled": true, "seed": 12345, "partnersPerEvent": 20 }
}
```

Cluster-envelope 使用 cluster 内实际最小/最大 strip 中心，并在两端直接增加 `margin`。Fake 分析只使用“单径迹且每个 Tracker 恰好有一个 LocalHit”的事件，结果仅作为独立假效率报告，不修正效率。

结果按 `DUTInfo.root/DUT_<id>/` 排列。`Efficiency/` 中严格只有七个 canvas：`cEventCount`、`cEffMapX`、`cEffMapY`、`cEffMap2D`、`cEff2DProjectionX`、`cEff2DProjectionY` 和 `cEffVsMargin`。两张投影图分别展示 2D 效率随 local X/Y 的变化；扫描图仅显示按总事件数计算的 X/Y/2D 曲线。`EfficiencySummary` 保存配置 margin 下的数值和误差。假效率分析对每个源事件 A 移除所有被 A 自身径迹匹配的 DUT clusters，保留其余 clusters，再与不同事件 B 的径迹匹配。Fake 分析直接复用效率分析的区域、binning、排除 bin、`margin` 和 `minEntriesPerBin`；`fake` 配置块只保留开关、随机种子和 partner 数。`Fake/` 保存 `hFakeEfficiencyMapX/Y/2D` 和 `FakeSummary`；`AnalysisConfig` 保存复现参数。

原有 `detectors`、`scripts` 和 `conversion` 内的转换器专用字段保持可用。旧的菜单入口、命令行目录覆盖选项、`input`/`output` 路径块及 `config/defaults.json` 配置档案已移除。

## 转换和输入

- `conversion` 不存在或 `enabled: false`：直接读取派生出的标准输入。
- `conversion.enabled: true`：从派生出的 raw 路径转换到 processed 路径。
- 输出已存在且 `overwrite: false`：复用已有输出。
- 输出不存在或 `overwrite: true`：调用 `conversion.type` 对应的转换器。

当前转换器类型为 `APV25SRS`、`BTAPVDat` 和 `VMM3a`。

`VMM3a` 首先从 `hits` 树构造事件。触发通道、ADC 阈值、通道时间
offset、符合窗以及触发前后事件时间窗均配置在
`conversion.eventBuilder` 中，完整示例见 `VMM/config.json`。VMM3a 输出
schema v3：

```text
Events:   evt, time, det, plane, id0, id1, chip, channel,
          waveform, adc, rawhitid, hit_time_ns
Channels: chip, channel, det, plane, id0, id1,
          pedestal, noise_sigma, gain, polarity
Metadata: schema_version
```

其中 `Events.time` 是触发的全局时间，`hit_time_ns` 是每个 hit 相对该
触发的时间。`HitProcessor` 自动识别 waveform 输入或 VMM 的
`adc/hit_time_ns` 直接输入，两种模式统一使用一个 `threshold`。VMM
聚类可用 `ClusterBuilder.timeWindowNs` 限制相邻 strip 的时间差。

VMM 原始 `time` 满足 `time = readout_time + chip_time`。事件构建必须使用
不会在芯片周期处回绕的全局 `time`；`readout_time` 缺少细时间，
`chip_time`、`bcid` 和 `tdc` 只适合电子学诊断。标准事件因此只保存触发
`time` 和相对 `hit_time_ns`，不重复保存这些组成量。

VMM geometry 使用 `(detector,fec,vmm)` 定位硬件映射，`plane` 是写入
schema v3 的输出标签，因此修改 plane 不会再破坏映射查找。plane 编号不
等同于物理 X/Y：实际测量方向由 detector 配置中的
`readoutPlanes.angle` 定义，0° 测量 local X，90° 测量 local Y。

## 输出

分析脚本的 ROOT 文件及其他结果写入 `<baseDir>/result/<run_id>/`。运行状态和最终统计只输出到终端，不生成额外的 JSON summary。

## 常见错误

- 参数数量错误：打印 usage 并返回退出码 2。
- JSON 无法读取或语法错误：在 ROOT 数据初始化之前返回非零。
- `scripts[i].type` 或 `conversion.type` 未注册：配置验证阶段报告完整字段位置。
- `DUTAnalysis`/`TimeResolution` 前没有启用 `TrackAnalysis`，且输出目录中没有 `TrackInfo.root`：配置验证失败。
- Event Display 找不到 `TrackInfo.root`：先用 analysis 配置生成径迹结果。
- 配置不存在：检查 `<baseDir>/config.json` 或兼容位置 `<baseDir>/raw/config.json`。
- 转换输入不存在：检查 `<baseDir>/raw/run<run_id>.*`。
