#pragma once

#include "DataModel.h"
#include "Detector.h"
#include <TChain.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

class Pipeline {
  public:
   Pipeline() = default;
   ~Pipeline();

   void Initialize(const std::string& configFile);

   void AddDetector(std::shared_ptr<Detector> det) {
       m_dets[det->GetID()] = det;
   }
   
   void Run(const std::string& inputFile);

   bool EventFilter(Event& event);

  private:
   std::tuple<int, int, int> ElectronicMap(int boardID, int channelID);

  private:
   std::map<int, std::shared_ptr<Detector>> m_dets;

   // ROOT 变量
   TChain* rawChain_ = nullptr;
   unsigned int apv_evt_ = 0;
   std::vector<unsigned int>* apv_id_ = nullptr;
   std::vector<unsigned int>* apv_ch_ = nullptr;
   std::vector<unsigned int>* mm_strip_ = nullptr;
   std::vector<std::vector<double>>* apv_q_ = nullptr;
};