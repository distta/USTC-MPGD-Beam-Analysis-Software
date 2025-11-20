#ifndef WAVEFORM_PROCESSOR_H
#define WAVEFORM_PROCESSOR_H

#include "AlgorithmFactory.h"
#include "IAlgorithm.h"

#include "Config.h"
#include "DataModel.h"

/**
 * @brief 波形处理算法 - 负责RawData到StripHit的转换
 *
 * 功能：
 * - 处理原始波形数据（ADC采样点）
 * - 提取关键物理量（峰值、电荷、时间等）
 * - 判断信号有效性
 */
class WaveformProcessor : public IAlgorithm {
   public:
    WaveformProcessor() = default;
    virtual ~WaveformProcessor() = default;

    // 实现IAlgorithm接口
    std::string GetName() const override { return "WaveformProcessor"; }
    std::string GetVersion() const override { return "1.0.0"; }

    void LoadConfig(const json& config) override {
        m_config.loadFrom(config);
    }

    void Print() const override {
        std::cout << "[" << GetName() << " v" << GetVersion() << "]" << std::endl;
        m_config.print();
    }

    StripHit ProcessWaveform(const RawData& rawData);

   private:
    WaveformConfig m_config;

    // 私有处理方法
    StripHit processWaveformDefault(const RawData& rawData);
    StripHit processWaveformLeadingEdgeFit(const RawData& rawData);
};

#endif  // WAVEFORM_PROCESSOR_H
