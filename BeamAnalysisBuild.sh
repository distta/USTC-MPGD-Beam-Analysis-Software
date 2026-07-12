#!/bin/bash

# BeamAnalysis 快速构建脚本
# Author: Huang Qixuan

set -e

echo "=== BeamAnalysis 项目构建 ==="

# 检查ROOT环境
if ! command -v root-config &> /dev/null; then
    echo "错误: 未找到ROOT框架，请检查ROOT环境配置"
    exit 1
fi

echo "ROOT版本: $(root-config --version)"

# 创建构建目录
if [ ! -d "build" ]; then
    mkdir build
    echo "创建构建目录: build/"
fi

cd build

# CMake配置
echo "配置项目..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# 编译
echo "编译项目..."
make -j$(nproc 2>/dev/null || echo 4)

# 检查结果
if [ -f "BeamAnalysis" ]; then
    echo "✓ 构建成功！"
else
    echo "✗ 构建失败：可执行文件未生成"
    exit 1
fi