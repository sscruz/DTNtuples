#include "DTNtupleTPGSimAnalyzer.h"
#include "TVector2.h"
#include "TF1.h"

// AM
// 1,2 -> 3 hits
// 3,4 -> 4 hits
// 5 -> 3+2
// 6 -> 3+3
// 7 -> 4+2
// 8 -> 4+3
// 9 -> 4+4


DTNtupleTPGSimAnalyzer::DTNtupleTPGSimAnalyzer(const TString & inFileName,
                                               const TString & outFileName,
                                               const TString & quality = "",
                                               Int_t index = -99,
                                               Int_t maxindex = +99
                                              ):
  m_outFile(outFileName,"RECREATE"), DTNtupleBaseAnalyzer(inFileName), quality_(quality), index_(index), maxindex_(maxindex)
{

  m_minMuPt = 20;

//  m_maxMuSegDPhi = 0.2; // Old
//  m_maxMuSegDEta = 0.3; // Old
   m_maxMuSegDPhi = 0.1;  // New
   m_maxMuSegDEta = 0.15; // New

  m_minSegHits = 4;

  m_maxSegTrigDPhi = 0.1;
  m_maxMuTrigDPhi  = 0.2;

}



DTNtupleTPGSimAnalyzer::~DTNtupleTPGSimAnalyzer()
{

}



void DTNtupleTPGSimAnalyzer::Loop()
{
  book();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEvent(jentry);   nbytes += nb;

      if (jentry % 100 == 0)
  std::cout << "[DTNtupleTPGSimAnalyzer::Loop] processed : "
      << jentry << " entries\r" << std::flush;

      fill();
    }

  std::cout << std::endl;

  endJob();
}



void DTNtupleTPGSimAnalyzer::book()
{
  m_outFile.cd();

//   std::vector<std::string> algoTag  = {"HB",      "AM", "AM+RPC"};
  std::vector<std::string> algoTag  = {"HB",      "AM"};
  std::vector<std::string> totalTag = {"matched", "total"};
  std::vector<std::string> chambTag = {"MB1",     "MB2", "MB3", "MB4"};

  for (const auto & algo : algoTag)
  {
    for (const auto & chamb : chambTag)
    {
      for (const auto & total : totalTag)
      {
      m_plots["Eff_" + chamb + "_" + algo + "_" + total] = new TH1D(("hEff_" + chamb + "_" + algo + "_" + total).c_str(),
                                                ("Efficiency for " + chamb + " " + algo + "; Sector; efficiency").c_str(),
                                                5, -2.5, +2.5);
      }
    }
  }
}



void DTNtupleTPGSimAnalyzer::fill()
{
  std::vector<std::string> chambTags = { "MB1", "MB2", "MB3", "MB4"};
  std::vector<std::string> whTags    = { "Wh.-2", "Wh.-1", "Wh.0", "Wh.+1", "Wh.+2"};
  std::vector<std::string> secTags   = { "Sec.1", "Sec.2", "Sec.3", "Sec.4", "Sec.5", "Sec.6", "Sec.7", "Sec.8","Sec.9","Sec.10","Sec.11","Sec.12","Sec.13","Sec.14"};

  for (std::size_t iGenPart = 0; iGenPart < gen_nGenParts; ++iGenPart)
  {
    if (std::abs(gen_pdgId->at(iGenPart)) != 13 || gen_pt->at(iGenPart) < m_minMuPt) continue;

    // CB this should not be a vector ...
    std::vector<std::size_t> bestSegIndex = { 999, 999, 999, 999 };
    std::vector<Int_t> bestSegNHits       = { 0, 0, 0, 0 };


    for (std::size_t iSeg = 0; iSeg < seg_nSegments; ++iSeg)
    {
      Int_t segSt    = seg_station->at(iSeg);
      Int_t segNHits = seg_phi_nHits->at(iSeg);
      
      Double_t muSegDPhi = std::abs(acos(cos(gen_phi->at(iGenPart) - seg_posGlb_phi->at(iSeg))));
      Double_t muSegDEta = std::abs(gen_eta->at(iGenPart) - seg_posGlb_eta->at(iSeg));
      

      if (muSegDPhi <  m_maxMuSegDPhi &&
          muSegDEta <  m_maxMuSegDEta &&
          segNHits  >= m_minSegHits   &&
          segNHits  >= bestSegNHits.at(segSt - 1)
          && (TMath::Abs(seg_phi_t0->at(iSeg)) < 15) // CUT ON TIME
          && ((seg_z_nHits->at(iSeg) >= 4 && (segSt < 4)) || (segSt == 4)) // SL2 requirement
         )
      {
        bestSegNHits[segSt - 1] = segNHits;
        bestSegIndex[segSt - 1] = iSeg;
      }
    }

    Int_t minQualityAM = -99;
    Int_t minQualityHB = -99;
    std::vector<Int_t> vetoedqualitiesAM; vetoedqualitiesAM.clear();
    Bool_t doQualityOR = false;
    if (quality_ == "nothreehits")
    {
      minQualityAM = 3;
      minQualityHB = 4;
    }
    else if (quality_ == "higherthanfour")
    {
      minQualityAM = 5;
      minQualityHB = 5;
    }
    else if (quality_ == "higherthanfourvetoing")
    {
      minQualityAM = 5;
      minQualityHB = 5;
      vetoedqualitiesAM.push_back(5);
      vetoedqualitiesAM.push_back(7);
    }
    else if (quality_ == "qualityOR")
    {
      doQualityOR = true;
      minQualityAM = 3;
      minQualityHB = 4;
    }
    

    // ==================== VARIABLES FOR THE HOUGH TRANSFORM BASED ALGORITHM
    for (const auto & iSeg : bestSegIndex)
    {
      if (iSeg == 999) continue;


      Int_t segWh  = seg_wheel->at(iSeg);
      Int_t segSec = seg_sector->at(iSeg);
      if (segSec == 13) segSec = 4;
      if (segSec == 14) segSec = 10;
      Int_t segSt  = seg_station->at(iSeg);
      
      std::string chambTag = chambTags.at(segSt - 1);
      std::string whTag    = whTags.at(segWh + 2);
      std::string secTag   = secTags.at(segSec - 1);

      Int_t    bestTPHB = -1;
      Double_t bestSegTrigHBDPhi = 1000;
      Double_t bestHBDPhi = 0;
      Int_t    besttrigHBBX = 0;
      for (std::size_t iTrigHB = 0; iTrigHB < ph2TpgPhiEmuHb_nTrigs; ++iTrigHB)
      {
        Int_t trigHBWh  = ph2TpgPhiEmuHb_wheel->at(iTrigHB);
        Int_t trigHBSec = ph2TpgPhiEmuHb_sector->at(iTrigHB);
        Int_t trigHBSt  = ph2TpgPhiEmuHb_station->at(iTrigHB);
        Int_t trigHBBX  = ph2TpgPhiEmuHb_BX->at(iTrigHB);

        if (segWh  == trigHBWh && segSec == trigHBSec &&  segSt  == trigHBSt)
        {
          Double_t trigGlbPhi    = trigPhiInRad(ph2TpgPhiEmuHb_phi->at(iTrigHB),trigHBSec);
          Double_t finalHBDPhi   = seg_posGlb_phi->at(iSeg) - trigGlbPhi;
          Double_t segTrigHBDPhi = abs(acos(cos(finalHBDPhi)));

          if (doQualityOR)
          {
            if ((segTrigHBDPhi < m_maxSegTrigDPhi) && (trigHBBX == 20) && (bestSegTrigHBDPhi > segTrigHBDPhi) && (ph2TpgPhiEmuHb_quality->at(iTrigHB) >= minQualityHB || ph2TpgPhiEmuHb_quality->at(iTrigHB) == -1) && (((ph2TpgPhiEmuHb_index->at(iTrigHB) < maxindex_) && (index_ < 0)) || ((ph2TpgPhiEmuHb_index->at(iTrigHB) == index_) && (index_ >= 0))) )
            {
              bestTPHB          = iTrigHB;
              besttrigHBBX      = trigHBBX;
              bestSegTrigHBDPhi = segTrigHBDPhi;
              bestHBDPhi        = TVector2::Phi_mpi_pi(finalHBDPhi);
            }
          }
          else
          {
            if ((segTrigHBDPhi < m_maxSegTrigDPhi) && (trigHBBX == 20) && (bestSegTrigHBDPhi > segTrigHBDPhi) && (ph2TpgPhiEmuHb_quality->at(iTrigHB) >= minQualityHB) && (((ph2TpgPhiEmuHb_index->at(iTrigHB) < maxindex_) && (index_ < 0)) || ((ph2TpgPhiEmuHb_index->at(iTrigHB) == index_) && (index_ >= 0))) )
            {
              bestTPHB          = iTrigHB;
              besttrigHBBX      = trigHBBX;
              bestSegTrigHBDPhi = segTrigHBDPhi;
              bestHBDPhi        = TVector2::Phi_mpi_pi(finalHBDPhi);
            }
          }
        }
      }

      if (bestTPHB > -1 && seg_phi_t0->at(iSeg) > -500)
      {
        m_plots["Eff_" + chambTag + "_HB_matched"]->Fill(segWh);
      }

      if (seg_phi_t0->at(iSeg) > -500)
      {
        m_plots["Eff_" + chambTag + "_HB_total"]->Fill(segWh);
      }


      // ==================== VARIABLES FOR THE ANALYTICAL METHOD ALGORITHM
      Int_t    bestTPAM = -1;
      Int_t    AMRPCflag= -1;
      Double_t bestSegTrigAMDPhi = 1000;
      Double_t bestAMDPhi = 0;
      Int_t    besttrigAMBX = 0;
      for (std::size_t iTrigAM = 0; iTrigAM < ph2TpgPhiEmuAm_nTrigs; ++iTrigAM)
      {
        Int_t trigAMWh  = ph2TpgPhiEmuAm_wheel->at(iTrigAM);
        Int_t trigAMSec = ph2TpgPhiEmuAm_sector->at(iTrigAM);
        Int_t trigAMSt  = ph2TpgPhiEmuAm_station->at(iTrigAM);
        Int_t trigAMBX  = ph2TpgPhiEmuAm_BX->at(iTrigAM);

        if (segWh == trigAMWh && segSec == trigAMSec && segSt  == trigAMSt)
        {
          Double_t trigGlbPhi    = trigPhiInRad(ph2TpgPhiEmuAm_phi->at(iTrigAM),trigAMSec);
          Double_t finalAMDPhi   = seg_posGlb_phi->at(iSeg) - trigGlbPhi;
          Double_t segTrigAMDPhi = abs(acos(cos(finalAMDPhi)));
          if (doQualityOR)
          {
            if ((segTrigAMDPhi < m_maxSegTrigDPhi) && (trigAMBX == 20) && (bestSegTrigAMDPhi > segTrigAMDPhi) && (ph2TpgPhiEmuAm_quality->at(iTrigAM) >= minQualityAM || ph2TpgPhiEmuAm_quality->at(iTrigAM) == -1) && (((ph2TpgPhiEmuAm_index->at(iTrigAM) < maxindex_) && (index_ < 0)) || ((ph2TpgPhiEmuAm_index->at(iTrigAM) == index_) && (index_ >= 0))))
            {
              if (vetoedqualitiesAM.size() == 0)
              {
                bestTPAM          = iTrigAM;
                besttrigAMBX      = trigAMBX;
                bestSegTrigAMDPhi = segTrigAMDPhi;
                bestAMDPhi        = TVector2::Phi_mpi_pi(finalAMDPhi);
                AMRPCflag         = ph2TpgPhiEmuAm_rpcFlag->at(iTrigAM);
              }
              else
              {
                Bool_t vetoed = false;
                for (UInt_t i = 0; i < vetoedqualitiesAM.size(); i++)
                {
                  if (ph2TpgPhiEmuAm_quality->at(iTrigAM) == vetoedqualitiesAM.at(i)) vetoed = true;
                }
                if (! vetoed) {
                  bestTPAM          = iTrigAM;
                  besttrigAMBX      = trigAMBX;
                  bestSegTrigAMDPhi = segTrigAMDPhi;
                  bestAMDPhi        = TVector2::Phi_mpi_pi(finalAMDPhi);
                  AMRPCflag         = ph2TpgPhiEmuAm_rpcFlag->at(iTrigAM);
                }
              }
            }
          }
          else
          {
            if ((segTrigAMDPhi < m_maxSegTrigDPhi) && (trigAMBX == 20) && (bestSegTrigAMDPhi > segTrigAMDPhi) && (ph2TpgPhiEmuAm_quality->at(iTrigAM) >= minQualityAM) && (((ph2TpgPhiEmuAm_index->at(iTrigAM) < maxindex_) && (index_ < 0)) || ((ph2TpgPhiEmuAm_index->at(iTrigAM) == index_) && (index_ >= 0))))
            {
              if (vetoedqualitiesAM.size() == 0)
              {
                bestTPAM          = iTrigAM;
                besttrigAMBX      = trigAMBX;
                bestSegTrigAMDPhi = segTrigAMDPhi;
                bestAMDPhi        = TVector2::Phi_mpi_pi(finalAMDPhi);
                AMRPCflag         = ph2TpgPhiEmuAm_rpcFlag->at(iTrigAM);
              }
              else
              {
                Bool_t vetoed = false;
                for (UInt_t i = 0; i < vetoedqualitiesAM.size(); i++)
                {
                  if (ph2TpgPhiEmuAm_quality->at(iTrigAM) == vetoedqualitiesAM.at(i)) vetoed = true;
                }
                if (! vetoed) {
                  bestTPAM          = iTrigAM;
                  besttrigAMBX      = trigAMBX;
                  bestSegTrigAMDPhi = segTrigAMDPhi;
                  bestAMDPhi        = TVector2::Phi_mpi_pi(finalAMDPhi);
                  AMRPCflag         = ph2TpgPhiEmuAm_rpcFlag->at(iTrigAM);
                }
              }
            }
          }
        }
      }

      if (bestTPAM > -1 && seg_phi_t0->at(iSeg) > -500)
      {
        m_plots["Eff_" + chambTag + "_AM_matched"]->Fill(segWh);
//         if (AMRPCflag > 0) m_plots["Eff_" + chambTag + "_AM+RPC_matched"]->Fill(segWh);
      }

      if (seg_phi_t0->at(iSeg) > -500)
      {
        m_plots["Eff_" + chambTag + "_AM_total"]->Fill(segWh);
//         m_plots["Eff_" + chambTag + "_AM+RPC_total"]->Fill(segWh);
      }
//       if (iSeg == 0)
//       {
//         for (std::size_t iTrigAM = 0; iTrigAM < ph2TpgPhiEmuAm_nTrigs; ++iTrigAM)
//         {
//           Int_t trigAMWh  = ph2TpgPhiEmuAm_wheel->at(iTrigAM);
//           Int_t trigAMSt  = ph2TpgPhiEmuAm_station->at(iTrigAM);
//           Int_t trigAMBX  = ph2TpgPhiEmuAm_BX->at(iTrigAM);
//           chambTag = chambTags.at(trigAMSt - 1);
//           if (trigAMBX > 5)
//           {
//             m_plots["Eff_" + chambTag + "_AM_total"]->Fill(trigAMWh);
//             m_plots["Eff_" + chambTag + "_AM+RPC_total"]->Fill(trigAMWh);
//           }
//         }
//       }

    }
  }
}



void DTNtupleTPGSimAnalyzer::endJob()
{
  m_outFile.cd();

  m_outFile.Write();
  m_outFile.Close();
}



Double_t DTNtupleTPGSimAnalyzer::trigPhiInRad(Double_t trigPhi, Int_t sector)
{
  return trigPhi / 65536. * 0.8 + TMath::Pi() / 6 * (sector - 1);
}
