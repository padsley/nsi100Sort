#define ProcessRuns_cxx
#include "ProcessRuns.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double sigma = 0.03; //Energy matching sigma

double m_FrequenceClock = 2.052;

void ProcessRuns::Loop()
{
    //   In a ROOT session, you can do:
    //      Root > .L ProcessRuns.C
    //      Root > ProcessRuns t
    //      Root > t.GetEntry(12); // Fill t data members with entry number 12
    //      Root > t.Show();       // Show values of entry 12
    //      Root > t.Show(16);     // Read and show values of entry 16
    //      Root > t.Loop();       // Loop on all entries
    //
    
    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0) return;
    
    ReadTimeTable();

    Long64_t nentries = fChain->GetEntriesFast();
    
    TFile *fout = new TFile("../processed/output.root","RECREATE");
    TTree *trout = new TTree("DATA","DATA");
    
    //Make branches here
    int position = -1;
    
    trout->Branch("position",&position,"position/I");
    
    int run = -1;
    trout->Branch("run",&run,"run/I");
    
    int tickN = -1;
    trout->Branch("tickN",&tickN,"tickN/I");
    
    double deltaE = -1;
    trout->Branch("deltaE",&deltaE,"deltaE/D");
    
    double plasG = -1;
    trout->Branch("plasG",&plasG,"plasG/D");
    
    double wire = -1;
    trout->Branch("wire",&wire,"wire/D");
    
    double rho = -1;
    trout->Branch("rho",&rho,"rho/D");
    
    double Brho = -1;
    trout->Branch("Brho",&Brho,"Brho/D");
    
    double Bfield = -1;
    trout->Branch("Bfield",&Bfield,"Bfield/D");
    
    double Ex = -1;
    trout->Branch("Ex",&Ex,"Ex/D");
    
    int SiliconHits = -1;
    
    trout->Branch("SiliconHits",&SiliconHits,"SiliconHits/I");
    
    int SiliconHitsPerDetector[6] = {-1,-1,-1,-1,-1,-1};
    trout->Branch("SiliconHitsPerDetector",&SiliconHitsPerDetector,"SiliconHitsPerDetector[6]/I");
    
    vector<double> SiliconEnergy;
    trout->Branch("SiliconEnergy",&SiliconEnergy);
    
    vector<int> SiliconTime;
    trout->Branch("SiliconTime",&SiliconTime);
    
    vector<int> SiliconDetector;
    trout->Branch("SiliconDetector",&SiliconDetector);
    
    vector<int> SiliconADCPChannel;
    trout->Branch("SiliconADCPChannel",&SiliconADCPChannel);
    
    vector<int> SiliconADCNChannel;
    trout->Branch("SiliconADCNChannel",&SiliconADCNChannel);
    
    vector<double> SiliconTheta;
    trout->Branch("SiliconTheta",&SiliconTheta);
    
    vector<double> SiliconPhi;
    trout->Branch("SiliconPhi",&SiliconPhi);
    
    vector<double> SiliconPath;
    trout->Branch("SiliconPath",&SiliconPath);
    
    vector<int> SiliconTDCChannel;
    trout->Branch("SiliconTDCChannel",&SiliconTDCChannel);
    
    vector<int> SiliconPStrip;
    trout->Branch("SiliconPStrip",&SiliconPStrip);
    
    vector<int> SiliconNStrip;
    trout->Branch("SiliconNStrip",&SiliconNStrip);
    
    vector<int> TimeDifference;
    trout->Branch("TimeDifference",&TimeDifference);
    
    vector<double> TimeDifference2;
    trout->Branch("TimeDifference2",&TimeDifference2);
    
    vector<double> CoincidenceQValue;
    trout->Branch("CoincidenceQValue",&CoincidenceQValue);
    
    vector<double> SiliconCMPhi;
    trout->Branch("SiliconCMPhi",&SiliconCMPhi);
    
    vector<double> SiliconCMTheta;
    trout->Branch("SiliconCMTheta",&SiliconCMTheta);
    
    int FPTime = -1;
    trout->Branch("FPTime",&FPTime,"FPTime/I");
    
    int FPOtherTime = -1;
    trout->Branch("FPOtherTime",&FPOtherTime,"FPOtherTime/I");
    
   int CurrentRun = -1;
    
    TGraph *NMR = new TGraph();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        
        if(CurrentRun != RunNumber) // Check what the run number is - if the run number is different, change the current run number and read in the NMR file
        {
            CurrentRun = RunNumber;
            ReadNMRForRun(RunNumber, NMR); //Get the B-field file for the run
        }
        
        //Set the branches to some initial value
        
        position = -1;
        deltaE = -1;
        plasG = -1;
        wire = -1;
        rho = -1;
        Brho = -1;
        Ex = -1;
        Bfield = -1;
        SiliconHits = -1;
        FPTime = -1;
        FPOtherTime = -1;
	tickN = -1;
        run = -1;
        
        //CLEAR THE VECTORS - OTHERWISE THE DISK FILLS UP WITH CRUD!!!!!!
        
        SiliconEnergy.clear();
        SiliconTime.clear();
        SiliconDetector.clear();
        
        SiliconADCPChannel.clear();
        SiliconADCNChannel.clear();
        SiliconPStrip.clear();
        SiliconNStrip.clear();
        SiliconTDCChannel.clear();
        SiliconTheta.clear();
        SiliconPhi.clear();
        SiliconPath.clear();
        SiliconCMPhi.clear();
        SiliconCMTheta.clear();
        CoincidenceQValue.clear();
        
        SiliconRecoilPhi = 0;
        SiliconRecoilTheta = 0;
        //Qvalue = 0;
        
        TimeDifference.clear();
        TimeDifference2.clear();
        
        run = RunNumber;
        
        for(int i=0;i<6;i++)SiliconHitsPerDetector[i]=-1;
        
        for(int i=0;i<tdcN;i++)//Loop over TDCs to look for good TDC hits
        {
            if(tdcList[i] == 126) //One focal-plane timing value
            {
                if(FPTime == -1)FPTime = tdcData[i];
                //else printf("Error1\n");
            }
            if(tdcList[i] == 117) //The **other** FP timing value
            {
                if(FPOtherTime == -1)FPOtherTime = tdcData[i];
                //else printf("Error2\n");
            }
        }
        
        for(int i=0;i<adcN;i++)
        {
            if(adcList[i]==193) //deltaE
            {
                deltaE = adcData[i];
            }
            if(adcList[i]==194) //Wire
            {
                wire = adcData[i];
            }
            if(adcList[i]==196) //Plastic at large Brho
            {
                plasG = adcData[i];
            }
	    if(adcList[i]==195)
	      {
		//	plasP = adcData[i];
	      }
            if(adcList[i]==192) //Position
            {
                position = adcData[i];
                rho = 0.62126 + 4.06966e-5 * position; // position->rho
                
                
                //Work out the NMR with some copy of what's in NPTool
                m_RunStart   = m_TimeTable[RunNumber].first;
                m_RunStop    = m_TimeTable[RunNumber].second;
                m_RunLength  = m_RunStop.AsDouble() - m_RunStart.AsDouble();
                //                 m_CurrentNMR = m_NMRTable[m_CurrentRunNumber]; 
                
		tickN = tick;
                
                double fAbsoluteTick = m_RunStart.AsDouble() + tickN/m_FrequenceClock;
		//cout << tickN << endl;
		//cout << fAbsoluteTick << endl;
                
                Bfield = GetBField(fAbsoluteTick, NMR); // Get the B-field
		//    cout << "Bfield: " << Bfield << endl;
                Brho = rho * GetBField(fAbsoluteTick, NMR); //Calculate Brho
                //                 cout << "Brho: " << Brho << endl;
                Ex = BrhoToEx(Brho); //Calculate Ex from Brho
            }
        }
        
        if(Brho!=-1)
        {
            for(int i=0;i<adcN;i++)
            {
                //                 cout << "L260" << endl;
                //cout << "adcList[i]: " << adcList[i] << endl;
                if(((adcList[i]-adcList[i]%16)/16)%2==0 && adcList[i]<192)//look for p-side channel in the ADCs
                {//cout << "L262" << endl;
                    int DetectorNumber = (adcList[i] - adcList[i]%32)/32; // work out the detector number
                    for(int j=0;j<adcN;j++)
                    {//cout << "L265" << endl;
                        if(((adcList[j]-adcList[j]%16)/16)%2==1)//look for n-side ADC channel
                        {//cout << "L267" << endl;
                            //                             cout << DetectorNumber << endl;
                            //                             cout << (adcList[j] - adcList[j]%32)/32 << endl << endl;
                            if(DetectorNumber==(adcList[j] - adcList[j]%32)/32)//match p and n side on same DetectorNumber
                            {//cout << "L269" << endl;
                                double EnergyPside = Calibration(adcData[i],adcList[i]);//calculate energies
                                double EnergyNside = Calibration(adcData[j],adcList[j]);
                                if(abs(EnergyPside-EnergyNside)<5*sigma && EnergyPside>0.3 && EnergyNside>0.3) //if both hits are >300 keV and the energies agree within 5*resolution parameter defined at the top
                                {
                                    //cout << "L274" << endl;
                                    for(int k=0;k<tdcN;k++)//Check that TDC is from the same detector as the ADC - then check that it's the same channel
                                    {//cout << "L276" << endl;
                                        if(tdcList[k]<96 && (DetectorNumber == (tdcList[k]-tdcList[k]%16)/16) && (tdcList[k]%16 == adcList[i]%16))//Only accept the event if it's in the same channel on the same detector and if there is a good timing value for that detector
                                        {
                                            //Store the information
                                            //                                             cout << "Silicon Fill" << endl;
                                            SiliconTime.push_back(tdcData[k]);
                                            
                                            SiliconEnergy.push_back(EnergyPside); 
                                            
                                            SiliconDetector.push_back(DetectorNumber+1);
                                            
                                            SiliconTDCChannel.push_back(tdcList[k]);
                                            
//                                             TimeDifference.push_back(tdcData[k] - FPTime - TDCOffsets[tdcList[k]]);
                                            TimeDifference.push_back(tdcData[k] - FPTime);
                                            
//                                             TimeDifference2.push_back(tdcData[k] - FPOtherTime - TDCOffsets[tdcList[k]]);
                                            TimeDifference2.push_back(tdcData[k] - FPOtherTime);

                                            SiliconADCPChannel.push_back(adcList[i]);
                                            SiliconADCNChannel.push_back(adcList[j]);
                                            
                                            SiliconPStrip.push_back(adcList[i]%16);
                                            SiliconNStrip.push_back(adcList[j]%16);
                                            
                                            TVector3 VectorToSiliconHit = SiliconFlightPath(DetectorNumber,adcList[i]%16,adcList[j]%16);
                                            
                                            SiliconTheta.push_back(VectorToSiliconHit.Theta());
                                            //                                         cout << "SiliconTheta: " << VectorToSiliconHit.Theta() << endl;
                                            SiliconPhi.push_back(VectorToSiliconHit.Phi());
                                            SiliconPath.push_back(VectorToSiliconHit.Mag());
                                            
                                            CoincidenceQValue.push_back(KinematicsToQValue(Brho,EnergyPside,938.782980,VectorToSiliconHit));
                                            
                                            SiliconCMTheta.push_back(SiliconRecoilTheta);
                                            SiliconCMPhi.push_back(SiliconRecoilPhi);
                                            
                                            SiliconHits++;
                                            SiliconHitsPerDetector[DetectorNumber]++;
                                        }
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if(SiliconTime.size()!=SiliconEnergy.size())printf("BOLLOCKS!\n SiliconEnergy: %lu \t SiliconTime: %lu\n",SiliconEnergy.size(),SiliconTime.size()); //Print error if mismatched vector sizes
        if(position>0)trout->Fill(); //If good position, fill the tree
    }
//     hPos->Write();
    trout->Write();
    fout->Close();     
}

void ReadTimeTable()
{
    cout << "\tReading time table file..." << endl;
    
    ifstream file("TimeTable.txt");
    
    // define variables
    Int_t run_midas, run_narval;
    Int_t date1, date2, time1, time2;
    TTimeStamp start, stop;
    pair<TTimeStamp, TTimeStamp> ptime;
    
    // read file
    while (file >> run_narval >> run_midas >> date1 >> time1 >> date2 >> time2) {
        start.Set(date1, time1, 0, 1, 0);
        stop.Set(date2, time2, 0, 1, 0);
        ptime.first  = start;
        ptime.second = stop;
        // fill maps
        m_TimeTable[run_midas] = ptime;
        m_NarvalMidasTable[run_narval] = run_midas;
              cout << "Narval: " << run_narval << " Midas: " << run_midas << endl;
    }
    
    // close file
    file.close();
}


void ReadNMRForRun(int RunNumber, TGraph *graphy)
{
    Double_t old_x=-10.;
    if(graphy->GetN()>0)
    {
        for(int i=0;i<graphy->GetN();i++)
        {
            graphy->RemovePoint(i);
        }
    }
    
    
    double fDelay = 3657.;
    //     fDelay = 3500.;
    
    ifstream input;
    char buffer[256];
    
    cout << "RunNumber: " << RunNumber << endl;
    
    int narvalRun = m_NarvalMidasTable.find(RunNumber)->second;
    map<int,int>::iterator it;
    int key = -1;
    for(it = m_NarvalMidasTable.begin(); it != m_NarvalMidasTable.end();it++)
    {
        if(it->second == RunNumber)
        {
            key = it->first;
        }
    }
    narvalRun = key;
    
    cout << "Narval Run Number: " << narvalRun << endl;
    
    sprintf(buffer,"../rmn/Data_%d.dat",narvalRun);
    
    input.open(buffer);
    
    if(!input.is_open())
    {
        std::cout<<"\n Cant open Rmn data file "<<buffer<<" ...\n"<<std::flush;
        std::cout<<"\n Nor fRmnFile "<<buffer<<" ...\n"<<std::flush;
        
        return ;
    }
    ////////////////Copied in
    char tmp[1000];
    input.getline(tmp,sizeof(tmp),'\n');
    //  TTimeStamp a;
    input.getline(tmp,sizeof(tmp),':');
    Int_t year,month,day,hour,minutes,seconds;
    
    char ch_tmp;
    input>>year;input>>ch_tmp;input>>month;input>>ch_tmp;input>>day;
    input.getline(tmp,sizeof(tmp),'-');
    input>>hour;input>>ch_tmp;input>>minutes;input>>ch_tmp;input>>seconds;
    
    TTimeStamp fOpenFileTime;
    fOpenFileTime.Set(year,month,day,hour,minutes,seconds, 0, 1, 0);
    //  fOpenFileTime.Print();
    Double_t time=fOpenFileTime.AsDouble(),x=1.,y;
    Int_t i=0;
    
    Double_t ymin = 10;
    Double_t ymax = -1;
    Double_t ymean = 0;
    
    input.getline(tmp,sizeof(tmp),'\n');
    input.getline(tmp,sizeof(tmp),'\n');
    for(;;){
        input>>x;input>>ch_tmp;input>>y;
        y*=1.e-7;
        // apply large field
        //     if (fIsLargeField) y += 1;
        
        if(x > old_x) {
            // TGraph
            graphy             -> SetPoint(i++, x+time+fDelay, y);
            
            // min, max
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
            // mean
            ymean += y;
        }
        else break;
        old_x=x;
        x=0.;
        
    }
    input.close();
    
    
    char saveBuffer[256];
    graphy->SetName("NMR");
    sprintf(saveBuffer,"figures/NMRGraphRun%d.root",RunNumber);
    graphy->SaveAs(saveBuffer);
    
    // set mean, min and max field values
    //   fMin = ymin;
    //   fMax = ymax;
    //   fMean = ymean / fRmnRelativeTime->GetN();
    // set min and max absolute time
    //   fTMin = fRmn->GetX()[0];
    //   fTMax = fRmn->GetX()[fRmn->GetN()-1];
    
    //   return 0;
}
