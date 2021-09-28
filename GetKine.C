#include "CREXdata.h"
#include <fstream>


//For a given detector 
void GetKine(int run,double ADCcut,double peak,string det ="up",string Det ="Upstream",double E=2.18132){

 gStyle->SetOptStat(1002211);


 //Load CREX asymmetry table
 LoadTable("ca48_fsu.dat",0);

 char arm[5]; 
 double th0; //removing th0 as an argument. I am setting th0 to be positive
 int sign; // need for lab angle calculation

 if(run < 10000) { sprintf(arm,"L"); th0 = 4.789; sign = 1; } else{ sprintf(arm,"R"); th0 = 4.771; sign = -1;} 

 double cth0 = TMath::Cos(th0*d2r);
 double sth0 = TMath::Sin(th0*d2r);

 //Histograms
 TH1F *Lab = new TH1F("Lab",Form(" %sHRS #theta_{lab} for %s Detector, Run %d",arm,Det.c_str(),run),150,2,10);
 TH1F *Qsq = new TH1F("Qsq",Form("%sHRS Q^{2} for %s Detector, Run %d",arm,Det.c_str(),run),150,0,0.07);
 TH1F *Asym = new TH1F("Asym",Form("%sHRS Asymmetry for %s Detector, Run %d",arm,Det.c_str(),run),150,1.8,3.0); 
 TH1F *thtg = new TH1F("thtg",Form("%sHRS #theta_{tg} for %s Detector, Run %d",arm,Det.c_str(),run),150,-0.06,0.06);
 TH1F *phtg = new TH1F("phtg",Form("%sHRS #phi_{tg} for %s Detector, Run %d",arm,Det.c_str(),run),150,-0.03,0.03);

 Lab->GetXaxis()->SetTitle("#theta_{lab} (deg)");
 Qsq->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
 Asym->GetXaxis()->SetTitle("Asymmetry (ppm)");
 thtg->GetXaxis()->SetTitle("#theta_{tg} (rad)");
 phtg->GetXaxis()->SetTitle("#phi_{tg} (rad)");

 TChain *T = new TChain("T");
 T->Add(Form("../Rootfiles/prex%sHRS_%d_-1*.root",arm,run));

 //Relevant data from root tree
 double thisu1, thisv1, thisu2, thisv2;
 double thisP, thisTh, thisPh, thisEvt;
 double thisADC;
 //Kinematics
 double thisCos, thisQsq, thisAsym, thisAng;
 
 //I am going keep this just for consistency with what we have been doing 
 double efact = E/(peak+0.001);

 vector<double> qsq, Angle, asym, Thtg, Phtg;


 T->SetBranchAddress(Form("%s.vdc.u1.nclust",arm),&thisu1); T->SetBranchAddress(Form("%s.vdc.v1.nclust",arm),&thisv1);
 T->SetBranchAddress(Form("%s.vdc.u2.nclust",arm),&thisu2); T->SetBranchAddress(Form("%s.vdc.v2.nclust",arm),&thisv2); 
 T->SetBranchAddress(Form("%s.gold.th",arm),&thisTh); T->SetBranchAddress(Form("%s.gold.ph",arm),&thisPh);
 T->SetBranchAddress(Form("%s.gold.p",arm),&thisP); T->SetBranchAddress("P.evtypebits",&thisEvt);
 T->SetBranchAddress(Form("P.%sQadc%s",det.c_str(),arm),&thisADC); 

 long n = T->GetEntries();
 
 //Looping over tree
 for(long i = 0; i < n; i++){
   T->GetEntry(i);

    int thisevent = (int) thisEvt;
    if(thisu1 == 1 && thisv1 == 1 && thisu2 == 1 && thisv2 == 1 && abs(thisTh)<0.08 && abs(thisPh)<0.05 && thisP > E*0.96 && thisP < E*1.002 && ((thisevent&2)==2) ) { 

      thisCos = (cth0 - sign*thisPh*sth0)/(TMath::Sqrt(1+thisTh*thisTh+thisPh*thisPh));  
      thisP *= efact;
      thisQsq = 2*E*thisP*(1-thisCos); 
      thisAng = r2d*TMath::ACos(thisCos);
      thisAsym = Interpolate(thisP*1000,thisAng,0,1);
      thisAsym /= ppm; 

        if(thisADC > ADCcut){
          qsq.push_back(thisQsq);
          Angle.push_back(thisAng);
          Thtg.push_back(thisTh);
          Phtg.push_back(thisPh);
          asym.push_back(thisAsym);
        }    
 
    }

 } 

 for(int j = 0; j < qsq.size(); j++){

    Qsq->Fill(qsq[j]);
    Lab->Fill(Angle[j]);
    thtg->Fill(Thtg[j]);
    phtg->Fill(Phtg[j]);
    Asym->Fill(asym[j]);

 }

 TCanvas *c_l = new TCanvas();
 Lab->Draw();
 c_l->SaveAs(Form("Lab%s_Run%d.pdf",Det.c_str(),run));


 TCanvas *c_q = new TCanvas();
 Qsq->Draw();
 c_q->SaveAs(Form("Qsq%s_Run%d.pdf",Det.c_str(),run));
 
 TCanvas *c_a = new TCanvas();
 Asym->Draw();
 c_a->SaveAs(Form("Asym%s_Run%d.pdf",Det.c_str(),run));

 TCanvas *c_th = new TCanvas();
 thtg->Draw();
 c_th->SaveAs(Form("thtg%s_Run%d.pdf",Det.c_str(),run));
 
 TCanvas *c_ph = new TCanvas();
 phtg->Draw();
 c_ph->SaveAs(Form("phtg%s_Run%d.pdf",Det.c_str(),run));

 gSystem->Exec(Form("pdfunite Lab%s_Run%d.pdf Qsq%s_Run%d.pdf Asym%s_Run%d.pdf thtg%s_Run%d.pdf phtg%s_Run%d.pdf Kine%sHRS_%s_Run%d.pdf",Det.c_str(),run,Det.c_str(),run,Det.c_str(),run,Det.c_str(),run,Det.c_str(),run,arm,Det.c_str(),run));



 long nentries = Qsq->GetEntries(); double val = TMath::Sqrt(nentries);
 
 std::cout << run << "   " << Lab->GetMean() << "   " << Lab->GetRMS()/val << "   " << Qsq->GetMean() << "   " << Qsq->GetRMS()/val << "   " << Asym->GetMean() << "  " << Asym->GetRMS()/val << "   " << thtg->GetMean() << "   " << thtg->GetRMS()/val << "   " << phtg->GetMean() << "   " << phtg->GetRMS()/val << std::endl;












}
