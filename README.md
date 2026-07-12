# BeamAnalysis

BeamAnalysis 将数据转换、探测器初始化、ROOT 数据读取和分析脚本统一到一个非交互式运行流程中。正式调用格式只有一种：

```bash
./BeamAnalysis <run_id> <base_dir>
```

例如：

```bash
./BeamAnalysis 1591 /data/beam
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
    "type": "SRS"
  },
  "detectors": [],
  "scripts": []
}
```

第二个命令行参数是数据集根目录。程序固定使用以下结构：

```text
<baseDir>/
├── config.json
├── raw/
│   └── run<run_id>.root
├── processed/
│   └── run<run_id>.root
└── result/
    └── <run_id>/
```

`BTAPVDat` 转换器的原始输入扩展名为 `.dat`，pedestal 默认读取 `<baseDir>/raw/run<run_id>_pedestal.dat`。转换器的 channel map 等专用参数仍写在 `conversion` 中。

程序会自动创建缺失的 `raw/`、`processed/` 和 `result/<run_id>/` 目录。`baseDir` 和 `config.json` 必须已经存在；需要读取的输入文件不存在时直接退出。

`mode` 支持：

- `analysis`：按数组顺序运行所有 `enabled: true` 的脚本，任一脚本失败即停止。
- `eventDisplay`：不运行脚本，读取输出目录中的 `TrackInfo.root` 并启动交互式事件显示。

`TrackAnalysis.config.runAlignment` 控制 tracker alignment；分析模式不会询问 y/n。旧字段 `performAlignment` 仍可读取，但新配置应使用 `runAlignment`。

原有 `detectors`、`scripts` 和 `conversion` 内的转换器专用字段保持可用。旧的菜单入口、命令行目录覆盖选项、`input`/`output` 路径块及 `config/defaults.json` 配置档案已移除。

## 转换和输入

- `conversion` 不存在或 `enabled: false`：直接读取派生出的标准输入。
- `conversion.enabled: true`：从派生出的 raw 路径转换到 processed 路径。
- 输出已存在且 `overwrite: false`：复用已有输出。
- 输出不存在或 `overwrite: true`：调用 `conversion.type` 对应的转换器。

当前转换器类型为 `SRS` 和 `BTAPVDat`。`config/srs.json` 与 `config/bt_apve.json` 给出了各自示例。

## 输出

分析脚本的 ROOT 文件及其他结果写入 `<baseDir>/result/<run_id>/`。运行状态和最终统计只输出到终端，不生成额外的 JSON summary。

## 常见错误

- 参数数量错误：打印 usage 并返回退出码 2。
- JSON 无法读取或语法错误：在 ROOT 数据初始化之前返回非零。
- `scripts[i].type` 或 `conversion.type` 未注册：配置验证阶段报告完整字段位置。
- `DUTAnalysis`/`TimeResolution` 前没有启用 `TrackAnalysis`，且输出目录中没有 `TrackInfo.root`：配置验证失败。
- Event Display 找不到 `TrackInfo.root`：先用 analysis 配置生成径迹结果。
- 配置不存在：检查 `<baseDir>/config.json`。
- 转换输入不存在：检查 `<baseDir>/raw/run<run_id>.*`。
