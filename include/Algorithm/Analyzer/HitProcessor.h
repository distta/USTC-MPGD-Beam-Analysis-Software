#ifndef HIT_PROCESSOR_H
#define HIT_PROCESSOR_H

#include "AlgorithmFactory.h"
#include "IAlgorithm.h"

#include "Config.h"
#include "DataModel.h"

// 前向声明，避免循环依赖
class DetectorFrame;

/**
 * @brief 波形处理算法 - 负责RawData到ChannelHit的转换
 *
 * 功能：
 * - 处理原始波形数据（ADC采样点）
 * - 提取关键物理量（峰值、电荷、时间等）
 * - 判断信号有效性
 */
class HitProcessor : public IAlgorithm {
   public:
    HitProcessor() = default;
    virtual ~HitProcessor() = default;

    // 实现IAlgorithm接口
    std::string GetName() const override { return "HitProcessor"; }
    std::string GetVersion() const override { return "1.0.0"; }

    void LoadConfig(const json& config) override {
        m_config.loadFrom(config);
    }

    void Print() const override {
        std::cout << "[" << GetName() << " v" << GetVersion() << "]" << std::endl;
        m_config.print();
    }

    // 统一接口: 处理DetectorFrame中的所有RawData
    bool Process(DetectorFrame& frame) override;

    ChannelHit ProcessHit(const RawData& rawData);

   private:
    HitProcessorConfig m_config;

    // 私有处理方法
    ChannelHit processDirectHit(const RawData& rawData);
    ChannelHit processWaveformDefault(const RawData& rawData);
    ChannelHit processWaveformLeadingEdgeFit(const RawData& rawData);
    ChannelHit processWaveformMode1(const RawData& rawData);
};

#endif  // HIT_PROCESSOR_H
