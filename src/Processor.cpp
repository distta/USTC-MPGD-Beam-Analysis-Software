#include "Processor.h"
#include "DataModel.h"
#include "TGraph.h"
#include <cstddef>
#include <vector>

bool Processor::initialize(const nlohmann::json& config) {

   if (config.contains("waveform")) {
      auto wf = config["waveform"];
      if (wf.contains("cfdFraction")) m_config.waveform.cfdFraction = wf["cfdFraction"].get<double>();
      if (wf.contains("noiseThreshold")) m_config.waveform.noiseThreshold = wf["noiseThreshold"].get<double>();
      if (wf.contains("saturationLevel")) m_config.waveform.saturationLevel = wf["saturationLevel"].get<double>();
      if (wf.contains("riseTimeStart")) m_config.waveform.riseTimeStart = wf["riseTimeStart"].get<double>();
   }

   if (config.contains("cluster")) {
      auto cl = config["cluster"];
      if (cl.contains("maxGap")) m_config.cluster.maxGap = cl["maxGap"].get<int>();
      if (cl.contains("minClusterSize")) m_config.cluster.minClusterSize = cl["minClusterSize"].get<int>();
      if (cl.contains("maxClusterSize")) m_config.cluster.maxClusterSize = cl["maxClusterSize"].get<int>();
   }

   return true;
}

std::vector<Cluster> Processor::buildCluster(std::map<StripInfo, RawData>& rawDataBuffer) {

   std::vector<Cluster> Clusters;
   Cluster aCluster;
   int stripTypeFlag = -1;
   int stripIDFlag = -1;

   for (const auto& [stripInfo, rawData] : rawDataBuffer) {

      if (stripInfo.stripType != stripTypeFlag || stripInfo.stripID > stripIDFlag + m_config.cluster.maxGap + 1) {
         aCluster.type = stripTypeFlag;
         aCluster.calculateProperties();
         if (aCluster.clusterSize > m_config.cluster.minClusterSize) {
            Clusters.push_back(aCluster);
         }
         aCluster.clear();
      }

      StripData stripData;
      stripData.stripID = stripInfo.stripID;
      stripData.type = stripInfo.stripType;

      processWaveform(rawData, stripData);

      aCluster.strips.push_back(stripData);
   }

   return Clusters;
}

void processWaveform(const RawData& rawData, StripData& stripData) {
   const std::vector<double>& waveform = rawData.adc;

   int wavePeakAmp = 0;
   int wavePeakTime = 0;
   int waveOverThStart = -1;
   int waveOverThEnd = -1;
   double waveFitTime = 0;
   TGraph SignalWaveGraph;

   for (size_t index = 0; index < waveform.size(); index++) {
      if (wavePeakAmp < waveform.at(index)) {
         wavePeakAmp = waveform.at(index);
         wavePeakTime = index;
      }

      if(waveform.at(index))

      // double corWave = 0;
      // for (int m = -2; m <= 2; m++) {
      //    if (k + m >= 0 && k + m < (*apv_q).at(j).size())
      //       corWave += (*apv_q).at(j).at(k + m) / pow(2, fabs(m));
      // }
      // corWave = corWave / 2.5;
      SignalWaveGraph.AddPoint(index, waveform.at(index));
   }
}