# BeamAnalysis 编译说明

## 环境要求

由于当前环境中 cmake 依赖 libssl.so.10 库，需要在 SLC7 环境中进行编译。

## 编译步骤

### 方法1：使用 SLC7 环境（推荐）

```bash
# 1. 进入 SLC7 环境
slc7

# 2. 启动 bash
bash

# 3. 进入项目目录
cd /ustcfs/STCFUser/qxhuang/workarea/BeamAnalysis

# 4. 使用构建脚本
bash BeamAnalysisBuild.sh
```

### 方法2：手动编译

```bash
# 在 SLC7 环境的 bash 中执行：
cd /ustcfs/STCFUser/qxhuang/workarea/BeamAnalysis
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

## 验证编译结果

编译成功后，可执行文件位于：
```
/ustcfs/STCFUser/qxhuang/workarea/BeamAnalysis/build/BeamAnalysis
```

## 运行程序

**注意**：必须从项目根目录运行程序，否则会出现配置文件路径错误。

```bash
cd /ustcfs/STCFUser/qxhuang/workarea/BeamAnalysis
./build/BeamAnalysis <run_number> [config_file]
```

## 新增功能

### Track 分析模式
- 仅重建 tracker cluster
- 支持粗对齐（含残差图输出）
- 支持精对齐
- 输出 Tracks TTree（供 DUT 使用）
- 输出 TrackerValidation TTree（用于 tracker 性能验证）

### DUT 分析模式
- 加载 track 信息
- 仅重建 DUT cluster
- 计算 DUT 残差
- 输出 DUTTree

## 代码修改概要

已完成以下修改：
1. ✓ 添加 ProcessMode 枚举类型
2. ✓ 修改 ProcessEvent 支持模式选择
3. ✓ 实现 RunTrackClustering（仅重建 tracker cluster）
4. ✓ 实现 RunCoarseAlignment（粗对齐 + 残差图输出）
5. ✓ 实现 RunFineAlignment（调用 AlignmentOptimizer）
6. ✓ 实现 RunTrackSelection（筛选并输出 Tracks TTree）
7. ✓ 实现 OutputTrackerValidation（输出 TrackerValidation TTree）
8. ✓ 修改 Run 方法实现 Track/DUT 分析模式选择
9. ✓ 实现 RunDUTClustering 和修改 RunDUTAnalysis 支持加载 track 信息

所有代码修改已通过静态检查，无编译错误。
