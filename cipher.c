#include "PMTDisplay.h"
int main(int argc, char** argv)
{

  const bool simple = false;
  const bool debug = false;
  const string imgFormat = "gif";  // "jpg"

  // Set the RAT DB configuration
  RAT::DB* db = RAT::DB::Get();
  //db->SetServer("https://deapdb.physics.carleton.ca/deapdb");  // no longer necessary
  db = db;  // silence warning
  RAT::RunConfig* rConfig = RAT::RunConfig::GetRunConfig();
  RAT::DetectorConfig *dc = new RAT::DetectorConfig();

  /////////////////////////////////////////////
  // Check the number of command line arguments and parse them.
  if (argc < 2 || argc > 4) {
    cout << "Usage: pmtdisplay txtFileName [outFileName inFilePath]" << endl;
    return -1;
  }
  vector< string > argsVec = ParseCommandLineArguments(argc, argv);
  string txtFileName = argsVec[1];
  string outFileName = txtFileName;
  if (argc > 2)
    outFileName = argsVec[2];
  else
    cout << "outFileName not specified, setting as '" << txtFileName << "'" << endl;
  string inFilePath = "";
  if (argc > 3) {
    inFilePath = argsVec[3];
    cout << "using specified path for CAL file: " << inFilePath << endl;
  }
  else
    cout << "inFilePath not specified, using CAC default path for CAL file" << endl;

  // Open up the text file and see how many lines there are
  // i.e. how many entries
  ifstream fileD;
  string txtFile = txtFileName + ".txt";
  fileD.open(txtFile.c_str());
  Int_t nLines = 0;
  ParseFile(fileD, nLines);
  cout << "Text file has: " << nLines << " listed events" << endl;
  Int_t nEffLines = nLines;

  // Initialise new pointers to the runIDs, subRunIDs and eventIDs
  // all to zero.
  Int_t* runIDs = new Int_t[ nLines ];
  Int_t* subRunIDs = new Int_t[ nLines ];
  Int_t* evIDs = new Int_t[ nLines ];
  for (Int_t iE = 0; iE < nLines; iE++) {
    runIDs[iE] = 0;
    subRunIDs[iE] = 0;
    evIDs[iE] = 0;
  }

  // Now loop over all the entries in the text file and fill
  // in their values
  Int_t curL = 0;
  while (fileD) {
    fileD >> subRunIDs[ curL ] >> runIDs[ curL ] >> evIDs[ curL ];
    curL++;
  }
  fileD.close();

  // Temporary variables used in the loop.
  Int_t curRunID = -1;
  Int_t curSubRunID = -1;
  Int_t curEvID = -1;
  Int_t prevSubRunID = -1;
  Int_t prevRunID = -1;
  Int_t firstEv = -1;
  Int_t nEv = -1;

  // Pointers to various data structures.
  TTree* curT = NULL;
  TFile* curF = NULL;
  //TTree* curRawT = NULL;
  //TFile* curRawF = NULL;

  // Create vectors for the PMT charge
  const double DEAP_Radius = 851.0;
  const int nPMT = 255;
  const int nTimeBins = 10;    // TimeFit2 uses 400 bins of width 40 ns
  Int_t timeBin[nTimeBins+1] = { 0, 10, 20, 30, 40, 80, 200, 400, 800, 1600, 16000 };
  vector< Double_t > pmtQ(nPMT);
  vector< Double_t > pmtQ_evt(nPMT);
  vector< vector< Double_t > > pmtQ_evt_time(nTimeBins+1, vector< Double_t >(nPMT) );
  vector< vector< Double_t > > pmtQ_evt_sum(nTimeBins+1, vector< Double_t >(nPMT) );
  TVector3* PMTpos = new TVector3[nPMT];
  for (int iPMT = 0; iPMT < nPMT; iPMT++) {
    dc->GetPMTPosition(iPMT, PMTpos[iPMT]);
    PMTpos[iPMT] = DEAP_Radius * PMTpos[iPMT].Unit();
  }
  Double_t subpeak_Q;
  Double_t subpeak_time;
  Int_t subpeak_time_bin;

  // Event loop
  for (Int_t iK=0; iK<nLines; iK++) {

    curRunID = runIDs[iK];
    curSubRunID = subRunIDs[iK];
    curEvID = evIDs[iK];
    if (prevRunID != curRunID) {
      rConfig->SetCurrentRun(curRunID);
      prevRunID = curRunID;
    }
    PMTInfoUtil pmtInfo;

    if (prevSubRunID != curSubRunID || prevRunID != curRunID) {
      if (inFilePath == "") {
        // with default file paths
        curT    = GetCALTreePointer(curRunID,curSubRunID,curF);
        //curRawT = GetRAWTreePointer(curRunID,curSubRunID,curRawF);
      } else {
        // with explicit CAL file path
        curT = GetTreePointer(inFilePath,curF);
        //curRawT = GetRAWTreePointer(curRunID,curSubRunID,curRawF);
      }
      prevSubRunID = curSubRunID;

      firstEv = GetFirstEntryEventID(curT);
      nEv = curT->GetEntries();
      if (debug) {
        cout << "There are " << nEv << " entries in the current CAL file" << endl;
        cout << "  GetFirstEntryEventID returned event number " << firstEv << endl;
      }
    }

    RAT::DS::Root* ds = new RAT::DS::Root();
    curT->SetBranchAddress("ds", &ds);
    //RAT::DS::Root* dsRaw = new RAT::DS::Root();
    //curRawT->SetBranchAddress("ds", &dsRaw);

    Int_t curIndex = curEvID - firstEv;
    if (curIndex < 0) {
      cout << "Requested event " << curEvID << " < " << firstEv << " first event number in file. Skip this event" << endl;
      nEffLines--;
      continue;
    }
    if (curIndex > nEv) {
      cout << "Requested event " << curEvID << " has calculated index " << curIndex << " > " << nEv << " events in file. Skip this event" << endl;
      nEffLines--;
      continue;
    }

    curT->GetEntry(curIndex);
    //curRawT->GetEntry(curIndex);

    if (ds->GetTS(0)->GetEventID() == curEvID) {
      cout << "----------------" << endl;
      cout << "On event: " << curEvID << endl;
      Int_t nPMT = ds->GetCAL(0)->GetPMTCount();

      cout << "PMT count: " << nPMT << endl;
      cout << "QPE count: " << ds->GetCAL(0)->GetQPE() << endl;
      cout << "Fprompt: " << ds->GetCAL(0)->GetFprompt() << endl;

      Double_t timefitT0            = ds->GetCAL(0)->GetTimeFit()->GetT0();
      Double_t timefitR             = ds->GetCAL(0)->GetTimeFit()->GetPosition().Mag();
      Double_t timefitCosTheta      = ds->GetCAL(0)->GetTimeFit()->GetPosition().CosTheta();
      Double_t timefitPhi           = ds->GetCAL(0)->GetTimeFit()->GetPosition().Phi();
      Double_t mblikelihoodR        = ds->GetEV(0)->GetMBLikelihood()->GetPosition().Mag();
      Double_t mblikelihoodCosTheta = ds->GetEV(0)->GetMBLikelihood()->GetPosition().CosTheta();
      Double_t mblikelihoodPhi      = ds->GetEV(0)->GetMBLikelihood()->GetPosition().Phi();
      //Double_t shellfitR            = ds->GetEV(0)->GetShellFit()->GetPosition().Mag();
      //Double_t shellfitCosTheta     = ds->GetEV(0)->GetShellFit()->GetPosition().CosTheta();
      //Double_t shellfitPhi          = ds->GetEV(0)->GetShellFit()->GetPosition().Phi();
      //Double_t centroidR            = ds->GetEV(0)->GetCentroid()->GetPosition().Mag();
      //Double_t centroidCosTheta     = ds->GetEV(0)->GetCentroid()->GetPosition().CosTheta();
      //Double_t centroidPhi          = ds->GetEV(0)->GetCentroid()->GetPosition().Phi();

      cout << "TimeFit T0: " << timefitT0 << endl;
      cout << "TimeFit R: " << timefitR << endl;
      cout << "TimeFit CosTheta: " << timefitCosTheta << endl;
      cout << "TimeFit Phi: " << timefitPhi << endl;
      cout << "MBLikelihood R: " << mblikelihoodR << endl;
      cout << "MBLikelihood CosTheta: " << mblikelihoodCosTheta << endl;
      cout << "MBLikelihood Phi: " << mblikelihoodPhi << endl;
      //cout << "ShellFit R: " << shellfitR << endl;
      //cout << "ShellFit CosTheta: " << shellfitCosTheta << endl;
      //cout << "ShellFit Phi: " << shellfitPhi << endl;
      //cout << "Centroid R: " << centroidR << endl;
      //cout << "Centroid CosTheta: " << centroidCosTheta << endl;
      //cout << "Centroid Phi: " << centroidPhi << endl;

      TGraph* TFgraph = new TGraph();
      TFgraph->SetPoint(0, timefitPhi, timefitCosTheta);
      TFgraph->SetMarkerSize(1);
      TFgraph->SetMarkerColor(kBlue);
      TFgraph->SetMarkerStyle(20);

      TGraph* MBgraph = new TGraph();
      MBgraph->SetPoint(0, mblikelihoodPhi, mblikelihoodCosTheta);
      MBgraph->SetMarkerSize(1);
      MBgraph->SetMarkerColor(kGreen+3);
      MBgraph->SetMarkerStyle(21);

      TLegend* legend_graph = new TLegend(0.62, 0.9, 0.75, 0.99);
      legend_graph->AddEntry(TFgraph, Form("TF R = %.0f mm", timefitR), "P");
      legend_graph->AddEntry(MBgraph, Form("MB R = %.0f mm", mblikelihoodR), "P");

      for (Int_t iPMT = 0; iPMT < nPMT; iPMT++) {
        RAT::DS::PMT* pmt = ds->GetCAL(0)->GetPMT(iPMT);
        Int_t pmtID = pmt->GetID();
        if (!pmtInfo.IsPMTGood(pmt)) continue;
        pmtQ[pmtID] += pmt->GetQPE();
        pmtQ_evt[pmtID] = pmt->GetQPE();

        if (simple) continue;

        // Get pulse sub-peak times and charges for animated display
        if (debug) cout << "PMT " << pmtID << " with QPE " << pmt->GetQPE() << ", pulse count " << pmt->GetPulseCount() << endl;
        for (int ipulse = 0; ipulse < pmt->GetPulseCount(); ++ipulse) {
          RAT::DS::Pulse *pulse = pmt->GetPulse(ipulse);
          if (debug) cout << "    pulse " << ipulse << " sub-peak count " << pulse->GetSubpeakCount() << endl;
          for (int ipeak = 0; ipeak < pulse->GetSubpeakCount(); ++ipeak) {
            subpeak_Q = pulse->GetSubpeakCharge(ipeak)/pmtInfo.GetSPECharge(pmt);
            subpeak_time = pulse->GetSubpeakTime(ipeak) - timefitT0;
            //subpeak_time = pulse->GetSubpeakTimeWithTOFCorrection(ipeak) - timefitT0;   // suggested by Thomas McElroy, to test
            if (debug) cout << "        subpeak_Q " << subpeak_Q << ", subpeak_time " << subpeak_time << endl;

            // Place in correct time bin
            subpeak_time_bin = nTimeBins;
            for (int ibin = 0; ibin < nTimeBins+1; ibin++) {
              if (subpeak_time < timeBin[ibin]) {
                subpeak_time_bin = ibin;
                break;
              }
            }

            if (debug) {
              cout << "        subpeak_time_bin " << subpeak_time_bin;
              if (subpeak_time_bin == 0)
                cout << ", underflow bin for times < " << timeBin[0] << " ns" << endl;
              else if (subpeak_time_bin == nTimeBins)
                cout << ", overflow bin for times > " << timeBin[nTimeBins-1] << " ns" << endl;
              else
                cout << ", for times between " << timeBin[subpeak_time_bin-1] << " and " << timeBin[subpeak_time_bin] << " ns" << endl;
            }

            // fill arrays
            pmtQ_evt_time[subpeak_time_bin][pmtID] += subpeak_Q;
            for (int ibin = subpeak_time_bin; ibin < nTimeBins+1; ibin++) {
              pmtQ_evt_sum[ibin][pmtID] += subpeak_Q;
            }

          } // end loop over sub-peaks
        } // end loop over pulses
      } // end loop over PMTs

      // Draw this event's overall PMT map
      TString h_evt_title = Form("Run %d, Event %d", curRunID, curEvID);
      TH2D* h_evt = PMTMap(PMTpos, pmtQ_evt, h_evt_title);
      TCanvas* c_evt = DrawPMTMap(h_evt);
      TFgraph->Draw("Psame");
      MBgraph->Draw("Psame");
      legend_graph->Draw("same");
      c_evt->Update();
      c_evt->SaveAs(Form("%s_event%d.%s", outFileName.c_str(), curEvID, imgFormat.c_str()));
      delete h_evt;
      delete c_evt;

      if (simple) continue;

      // Draw this event's PMT maps vs time
      for (int ibin = 0; ibin < nTimeBins+1; ibin++) {
        TString h_sum_title = Form("Run %d, Event %d, Time < %d ns", curRunID, curEvID, timeBin[ibin]);
        TString h_time_title = h_sum_title;  // correct for underflow ibin == 0
        if (ibin > 0) {
          h_time_title = Form("Run %d, Event %d, Time %d to %d ns", curRunID, curEvID, timeBin[ibin-1], timeBin[ibin]);
        }
        if (ibin == nTimeBins) {
          h_sum_title = h_evt_title;
          h_time_title = Form("Run %d, Event %d, Time > %d ns", curRunID, curEvID, timeBin[nTimeBins-1]);  // overflow bin
        }

        TH2D* h_time = PMTMap(PMTpos, pmtQ_evt_time[ibin], h_time_title);
        TCanvas* c_time = DrawPMTMap(h_time);
        TFgraph->Draw("Psame");
        MBgraph->Draw("Psame");
        legend_graph->Draw("same");
        c_time->Update();
        c_time->SaveAs(Form("%s_event%d_time%02d.%s", outFileName.c_str(), curEvID, ibin, imgFormat.c_str()));
        delete h_time;
        delete c_time;

        TH2D* h_sum = PMTMap(PMTpos, pmtQ_evt_sum[ibin], h_sum_title);
        TCanvas* c_sum = DrawPMTMap(h_sum);
        TFgraph->Draw("Psame");
        MBgraph->Draw("Psame");
        legend_graph->Draw("same");
        c_sum->Update();
        c_sum->SaveAs(Form("%s_event%d_sum%02d.%s", outFileName.c_str(), curEvID, ibin, imgFormat.c_str()));
        delete h_sum;
        delete c_sum;
      }

      delete TFgraph;
      delete MBgraph;

    } // end if GetTS(0)->GetEventID()...
  } // end loop over events


  // Draw the averaged PMT map over all events considered
  cout << "Normalizing averaged PMT map over " << nEffLines << " processed events" << endl;
  for (int kpmt = 0; kpmt < nPMT; kpmt++) {
    pmtQ[kpmt] = pmtQ[kpmt]/nEffLines;   // effective lines: not counting events not found in root files
  }
  TH2D* h_avg = PMTMap(PMTpos, pmtQ, "Average PMT charges");
  TCanvas* c_avg = DrawPMTMap(h_avg);
  c_avg->SaveAs(Form("%s_SummedMap.%s", outFileName.c_str(), imgFormat.c_str()));
  delete h_avg;
  delete c_avg;

}

////////////////////////////////////////////////
////////////////////////////////////////////////

void ParseFile(ifstream& fFile, Int_t &nLines)
{
  string tmpStr = "";
  nLines = -1;
  while (fFile) {
    getline(fFile, tmpStr);
    nLines++;
  }

  fFile.clear();
  fFile.seekg(0, ios::beg);
}

////////////////////////////////////////////////
////////////////////////////////////////////////

vector< string > ParseCommandLineArguments(int argc, char** argv)
{
  vector< string > vecStr;
  stringstream tmpStream;
  string tmpStr = "";
  for (Int_t iArg = 0; iArg < argc; iArg++) {
    tmpStr = ""; tmpStream.clear(); tmpStream.str("");
    tmpStream << argv[iArg];
    tmpStream >> tmpStr;
    vecStr.push_back(tmpStr);
  }
  return vecStr;
}

////////////////////////////////////////////////
////////////////////////////////////////////////

TTree* GetTreePointer(string filePath, TFile* curF)
{
  if (!CheckFileExists(filePath)) {
    cout << "File does not exist: " << filePath << endl;
    return NULL;
  }
  curF = new TFile((filePath).c_str(), "READ");
  TTree* myTree = (TTree*)curF->Get("T");
  cout << "Opened file: " << (filePath).c_str() << endl;
  return myTree;
}

////////////////////////////////////////////////
////////////////////////////////////////////////

TTree* GetCALTreePointer(Int_t runID, Int_t subRunID, TFile* curF)
{
  stringstream tmpStream;
  string runIDStr = "";
  string subRunIDStr = "";

  tmpStream << setw(6) << setfill('0') << runID;
  tmpStream >> runIDStr; tmpStream.clear();

  tmpStream << setw(4) << setfill('0') << subRunID;
  tmpStream >> subRunIDStr; tmpStream.clear();

  std::string dirPrefix = "/home/hpc3832/data/cal/run" + runIDStr;
  std::string subRunPrefix = "/deap_cal_"
                             + runIDStr
                             + "_" + subRunIDStr
                             + ".root";
  return GetTreePointer(dirPrefix+subRunPrefix, curF);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
/*
TTree* GetRAWTreePointer(Int_t runID, Int_t subRunID, TFile* curF)
{
  stringstream tmpStream;
  string runIDStr = "";
  string subRunIDStr = "";

  tmpStream << setw(6) << setfill('0') << runID;
  tmpStream >> runIDStr; tmpStream.clear();

  tmpStream << setw(4) << setfill('0') << subRunID;
  tmpStream >> subRunIDStr; tmpStream.clear();

  std::string dirPrefix = "/home/hpc3832/data/raw/run" + runIDStr;
  std::string subRunPrefix = "/deap_raw_"
                             + runIDStr
                             + "_" + subRunIDStr
                             + ".root";
  return GetTreePointer(dirPrefix+subRunPrefix, curF);
}
*/
////////////////////////////////////////////////
////////////////////////////////////////////////

string IntToString(Int_t intVal)
{
  stringstream tmpStream;
  string tmpString = "";
  tmpStream << intVal;
  tmpStream >> tmpString;
  return tmpString;
}


////////////////////////////////////////////////
////////////////////////////////////////////////

Bool_t CheckFileExists(std::string fileName)
{
  ifstream f(fileName.c_str());
  if (f.good()) { return true; }
  else { return false; }
}

////////////////////////////////////////////////
////////////////////////////////////////////////

Double_t GetMaximumValue(TH1D* hist)
{
  Int_t binV = hist->GetMaximumBin();
  return hist->GetBinContent(binV);
}

////////////////////////////////////////////////
////////////////////////////////////////////////

Double_t GetMaximumValue(TH2D* hist)
{
  Int_t binV = hist->GetMaximumBin();
  return hist->GetBinContent(binV);
}

////////////////////////////////////////////////
////////////////////////////////////////////////

// Int_t GetFirstEntryEventID(TTree* myT)
// {
//   RAT::DS::Root* ds = new RAT::DS::Root();
//   myT->SetBranchAddress("ds", &ds);
//   myT->GetEntry(0);
//   Int_t evID = ds->GetTS(0)->GetEventID();
//   return evID;
}

////////////////////////////////////////////////
////////////////////////////////////////////////

TH2D* PMTMap(TVector3* PMTpos, vector< Double_t > pmtQ, TString title)
{
  TString titlesXYZ = title + ";PMT #phi; PMT Cos#theta; QPE";
  TH2D* pmtMapPtr = new TH2D(title, titlesXYZ, 600, -TMath::Pi(), TMath::Pi()*1.04, 600, -1.0, 1.);

  const double LG_radius = 95.0;
  TVector3 binvec = TVector3(1.0, 0.0, 0.0);
  TVector3 diff   = TVector3(1.0, 0.0, 0.0);

  unsigned int curpmt = 0;
  TVector3 curdiff = TVector3(1.0, 0.0, 0.0);

  for (int i=0; i<pmtMapPtr->GetNbinsX(); i++) {
    for (int j=0; j<pmtMapPtr->GetNbinsY(); j++) {

      double theta = acos(pmtMapPtr->GetYaxis()->GetBinLowEdge(j + 1));
      double phi   = pmtMapPtr->GetXaxis()->GetBinLowEdge(i + 1);

      binvec = TVector3(1.0, 0.0, 0.0);
      binvec.SetMagThetaPhi(855., theta, phi);

      // first check whether we are still in the same PMT as the previous bin
      curdiff = PMTpos[curpmt] - binvec;
      if (curdiff.Mag() < LG_radius) {
        if (curdiff.Mag() > LG_radius - 5.) {
          pmtMapPtr->SetBinContent(i + 1, j + 1, 1e6);   // PMT boundaries drawn black ie very high "QPE value"
        } else {
          pmtMapPtr->SetBinContent(i + 1, j + 1, pmtQ[curpmt]);
        }
        continue;   // we are in the same PMT as before and filled this histogram bin: now go to next bin
      }

      // if we get here, we are no longer in the same PMT, so search
      for (unsigned int kpmt = 0; kpmt < pmtQ.size(); kpmt++) {
        //if (kpmt == curpmt) continue;
        diff = PMTpos[kpmt] - binvec;
        if (diff.Mag() < LG_radius) {
          if (diff.Mag() > LG_radius - 5.) {
            pmtMapPtr->SetBinContent(i + 1, j + 1, 1e6);   // PMT boundaries drawn black ie very high "QPE value"
          } else {
            pmtMapPtr->SetBinContent(i + 1, j + 1, pmtQ[kpmt]);
          }
          curpmt = kpmt;
          break;   // only one PMT matches each histogram bin: stop looking when found
        }
      } // pmts
    } // y bins
  } // x bins

  return pmtMapPtr;
}

////////////////////////////////////////////////
////////////////////////////////////////////////

TCanvas* DrawPMTMap(TH2D* pmtMapPtr)
{
  RAT::DEAPStyle *fStyle = RAT::DEAPStyle::GetStyle();

  TCanvas* c = new TCanvas();
  c->cd();
  //c->SetLogz();
  c->SetGridx(0);
  c->SetGridy(0);
  gPad->SetRightMargin(0.15);
  pmtMapPtr->SetMinimum(0.1);
  pmtMapPtr->SetMaximum(30);
  pmtMapPtr->GetZaxis()->SetTitleOffset(0.9);
  //fStyle->chooseDSPalette(1,0);  // DEAP purple
  fStyle->chooseDSPalette(0,1);    // inverted blackbody
  pmtMapPtr->Draw("COLZ");

  fStyle->drawDEAP(0);
  c->Update();
  return c;
}

