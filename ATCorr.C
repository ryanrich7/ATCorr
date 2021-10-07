#include "TLine.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "ATCorr.h"
#include <iostream>
#include <fstream>

using namespace std;

ATCorr::ATCorr(string file1, string file2){

  ifstream File1, File2;
  File1.open(file1.c_str());
  File2.open(file2.c_str());

  double rL, qL, qeL, tL, teL, pL, peL, lL, leL, aL, aeL, cL, ceL, sL, seL,vL,veL, hL, heL;
  double rR, qR, qeR, tR, teR, pR, peR, lR, leR, aR, aeR, cR, ceR, sR, seR, vR, veR, hR, heR;


  while(File1 >> rL >> lL >> leL >> qL >> qeL >> aL >> aeL >> tL >> teL >> pL >> peL >> cL >> ceL >> sL >> seL >> vL >> veL >> hL >> heL  ){
  RunL.push_back(rL); AngleL.push_back(lL); AngleerL.push_back(leL); QsqL.push_back(qL); QsqerL.push_back(qeL); asymL.push_back(aL);
  asymerL.push_back(aeL); ThtgL.push_back(tL); ThtgerL.push_back(teL); PhtgL.push_back(pL); PhtgerL.push_back(peL); 
  cosphL.push_back(cL); cospherL.push_back(ceL); sinphL.push_back(sL); sinpherL.push_back(seL);
  VermomL.push_back(vL); VermomerL.push_back(veL); HormomL.push_back(hL); HormomerL.push_back(heL); 

  }

  while(File2 >> rR >> lR >> leR >> qR >> qeR >> aR >> aeR >> tR >> teR >> pR >> peR >> cR >> ceR >> sR >> seR >> vR >> veR >> hR >> heR ){
  RunR.push_back(rR); AngleR.push_back(lR); AngleerR.push_back(leR); QsqR.push_back(qR); QsqerR.push_back(qeR); asymR.push_back(aR);
  asymerR.push_back(aeR); ThtgR.push_back(tR); ThtgerR.push_back(teR); PhtgR.push_back(pR); PhtgerR.push_back(peR);
  cosphR.push_back(cR); cospherR.push_back(ceR); sinphR.push_back(sR); sinpherR.push_back(seL);
  VermomR.push_back(vR); VermomerR.push_back(veR); HormomR.push_back(hR); HormomerR.push_back(heR);
 
  }



}

ATCorr::ATCorr(string file1, string file2, string file3, string file4){
   
  ifstream File1, File2, File3, File4;
  File1.open(file1.c_str());
  File2.open(file2.c_str());
  File3.open(file3.c_str());
  File4.open(file4.c_str());

 double rL, qL, qeL, tL, teL, pL, peL, lL, leL, aL, aeL, cL, ceL, sL, seL, vL, veL, hL, heL;
 double rL1, qL1, qeL1, tL1, teL1, pL1, peL1, lL1, leL1, aL1, aeL1, cL1, ceL1, sL1, seL1,vL1,veL1, hL1, heL1;
 double rR, qR, qeR, tR, teR, pR, peR, lR, leR, aR, aeR, cR, ceR, sR, seR,vR,veR, hR, heR;
 double rR1, qR1, qeR1, tR1, teR1, pR1, peR1, lR1, leR1, aR1, aeR1,cR1, ceR1, sR1, seR1,vR1,veR1, hR1, heR1;



  while(File1 >> rL >> lL >> leL >> qL >> qeL >> aL >> aeL >> tL >> teL >> pL >> peL  >> cL >> ceL >> sL >> seL >> vL >> veL >> hL >> heL){
  RunL.push_back(rL); AngleL.push_back(lL); AngleerL.push_back(leL); QsqL.push_back(qL); QsqerL.push_back(qeL); asymL.push_back(aL);
  asymerL.push_back(aeL); ThtgL.push_back(tL); ThtgerL.push_back(teL); PhtgL.push_back(pL); PhtgerL.push_back(peL);
  cosphL.push_back(cL); cospherL.push_back(ceL); sinphL.push_back(sL); sinpherL.push_back(seL);
  VermomL.push_back(vL); VermomerL.push_back(veL); HormomL.push_back(hL); HormomerL.push_back(heL);
  }

  while(File2 >> rR >> lR >> leR >> qR >> qeR >> aR >> aeR >> tR >> teR >> pR >> peR  >> cR >> ceR >> sR >> seR >> vR >> veR >> hR >> heR){
  RunR.push_back(rR); AngleR.push_back(lR); AngleerR.push_back(leR); QsqR.push_back(qR); QsqerR.push_back(qeR); asymR.push_back(aR);
  asymerR.push_back(aeR); ThtgR.push_back(tR); ThtgerR.push_back(teR); PhtgR.push_back(pR); PhtgerR.push_back(peR);
  cosphR.push_back(cR); cospherR.push_back(ceR); sinphR.push_back(sR); sinpherR.push_back(seR);
  VermomR.push_back(vR); VermomerR.push_back(veR); HormomR.push_back(hR); HormomerR.push_back(heR);
  }


  while(File3 >> rL1 >> lL1 >> leL1 >> qL1 >> qeL1 >> aL1 >> aeL1 >> tL1 >> teL1 >> pL1 >> peL1  >> cL1 >> ceL1 >> sL1 >> seL1 >> vL1 >> veL1 >> hL1 >> heL1){
  RunL1.push_back(rL1); AngleL1.push_back(lL1); AngleerL1.push_back(leL1); QsqL1.push_back(qL1); QsqerL1.push_back(qeL1); asymL1.push_back(aL1);
  asymerL1.push_back(aeL1); ThtgL1.push_back(tL1); ThtgerL1.push_back(teL1); PhtgL1.push_back(pL1); PhtgerL1.push_back(peL1);
  cosphL1.push_back(cL1); cospherL1.push_back(ceL1); sinphL1.push_back(sL1); sinpherL1.push_back(seL1);
  VermomL1.push_back(vL1); VermomerL1.push_back(veL1); HormomL1.push_back(hL1); HormomerL1.push_back(heL1);
  }

  while(File4 >> rR1 >> lR1 >> leR1 >> qR1 >> qeR1 >> aR1 >> aeR1 >> tR1 >> teR1 >> pR1 >> peR1  >> cR1 >> ceR1 >> sR1 >> seR1 >> vR1 >> veR1 >> hR1 >> heR1){
  RunR1.push_back(rR1); AngleR1.push_back(lR1); AngleerR1.push_back(leR1); QsqR1.push_back(qR1); QsqerR1.push_back(qeR1); asymR1.push_back(aR1);
  asymerR1.push_back(aeR1); ThtgR1.push_back(tR1); ThtgerR1.push_back(teR1); PhtgR1.push_back(pR1); PhtgerR1.push_back(peR1);
  cosphR1.push_back(cR1); cospherR1.push_back(ceR1); sinphR1.push_back(sR1); sinpherR1.push_back(seR1);
  VermomR1.push_back(vR1); VermomerR1.push_back(veR1); HormomR1.push_back(hR1); HormomerR1.push_back(heR1);
  }

}

ATCorr::~ATCorr(){

  cout << "Signing off" << endl;
}

int ATCorr::Init(){

//Doing basic computations here
  ComputeQsqAvg();ComputeQsqDiff();
  ComputePhtgAvg();ComputePhtgDiff();
  ComputeThtgAvg();ComputeThtgDiff();
  ComputeAsymAvg(); ComputeAsymDiff();
  ComputeVerAvg(); ComputeVerDiff();
  ComputeHorAvg(); ComputeHorDiff();
  ComputeCosAvg(); ComputeCosDiff();
  ComputeSinAvg(); ComputeSinDiff();

  return 1;

}

int ATCorr::Init(int i){

//Do 4 detector combo here
  ComputeAsym();
  ComputeTh();

return i;

}


void ATCorr::PlotMainDetPlots(){

  ComputeXi(1);
  //ComputeXi(1,2);

  int color[3] = {4,2,4};
  
   //This function is not the most efficient but it works 

  //Separate into different run periods
  double wien1[3] = {1,2,3}; double wien1err[3] = {0};
  double wien2[1] = {4}; double wien2err[1] = {0};
  double wien3[5] = {5,6,7,8,9}; double wien3err[5] = {0};

  //Detector combos
  double qsqavg1[3], qsqavg1err[3], qsqavg2[1], qsqavg2err[1], qsqavg3[5], qsqavg3err[5];
  double qsqdd1[3], qsqdd1err[3], qsqdd2[1], qsqdd2err[1], qsqdd3[5], qsqdd3err[5];
  double thtgdd1[3], thtgdd1err[3], thtgdd2[1], thtgdd2err[1], thtgdd3[5], thtgdd3err[5];
  double thtgavg1[3], thtgavg1err[3], thtgavg2[1], thtgavg2err[1], thtgavg3[5], thtgavg3err[5]; 
  double phtgdd1[3], phtgdd1err[3], phtgdd2[1], phtgdd2err[1], phtgdd3[5], phtgdd3err[5];
  double phtgavg1[3], phtgavg1err[3], phtgavg2[1], phtgavg2err[1], phtgavg3[5], phtgavg3err[5];
  double asymavg1[3], asymavg1err[3], asymavg2[1], asymavg2err[1], asymavg3[5], asymavg3err[5];  
  double asymdd1[3],asymdd1err[3], asymdd2[1], asymdd2err[1], asymdd3[5], asymdd3err[5];
  double veravg1[3], veravg1err[3], veravg2[1], veravg2err[1], veravg3[5], veravg3err[5];
  double verdd1[3],verdd1err[3], verdd2[1], verdd2err[1], verdd3[5], verdd3err[5];
  double horavg1[3], horavg1err[3], horavg2[1], horavg2err[1], horavg3[5], horavg3err[5];
  double hordd1[3],hordd1err[3], hordd2[1], hordd2err[1], hordd3[5], hordd3err[5];
  double cosavg1[3], cosavg1err[3], cosavg2[3], cosavg2err[3], cosavg3[5], cosavg3err[5];
  double sinavg1[3], sinavg1err[3], sinavg2[3], sinavg2err[3], sinavg3[5], sinavg3err[5];
  double cosdd1[3], cosdd1err[3], cosdd2[3], cosdd2err[3], cosdd3[5], cosdd3err[5];
  double sindd1[3], sindd1err[3], sindd2[3], sindd2err[3], sindd3[5], sindd3err[5];
  double xi1[3], xi1err[3], xi2[1], xi2err[1], xi3[5], xi3err[5];


  //Individual Detectors themselves
  double qsL1[3], qsL1er[3], qsL2[1], qsL2er[1], qsL3[5], qsL3er[5];
  double qsR1[3], qsR1er[3], qsR2[1], qsR2er[1], qsR3[5], qsR3er[5];
  double thsL1[3], thsL1er[3], thsL2[1], thsL2er[1], thsL3[5], thsL3er[5];
  double thsR1[3], thsR1er[3], thsR2[1], thsR2er[1], thsR3[5], thsR3er[5];
  double phsL1[3], phsL1er[3], phsL2[1], phsL2er[1], phsL3[5], phsL3er[5];
  double phsR1[3], phsR1er[3], phsR2[1], phsR2er[1], phsR3[5], phsR3er[5];
  double asyL1[3], asyL1er[3], asyL2[1], asyL2er[1], asyL3[5], asyL3er[5];
  double asyR1[3], asyR1er[3], asyR2[1], asyR2er[1], asyR3[5], asyR3er[5];
  double labL1[3], labL1er[3], labL2[1], labL2er[1], labL3[5], labL3er[5];
  double labR1[3], labR1er[3], labR2[1], labR2er[1], labR3[5], labR3er[5];
  double verL1[3], verL2[1], verL3[5], horL1[3], horL2[1], horL3[5];  
  double verR1[3], verR2[1], verR3[5], horR1[3], horR2[1], horR3[5];
  double verL1er[3], verL2er[1], verL3er[5], horL1er[3], horL2er[1], horL3er[5];
  double verR1er[3], verR2er[1], verR3er[5], horR1er[3], horR2er[1], horR3er[5];
  double cosL1[3], cosL2[1], cosL3[5], sinL1[3], sinL2[1], sinL3[5];  
  double cosR1[3], cosR2[1], cosR3[5], sinR1[3], sinR2[1], sinR3[5];
  double cosL1er[3], cosL2er[1], cosL3er[5], sinL1er[3], sinL2er[1], sinL3er[5];
  double cosR1er[3], cosR2er[1], cosR3er[5], sinR1er[3], sinR2er[1], sinR3er[5];


  int ii = 0; int iii = 0;
   //Separate data into Wien states
  for(int i = 0; i < 9; i++){ 
  
    if(i < 3) { qsqavg1[i] = QsqAvg[i]; qsqavg1err[i] = QsqAvgerr[i]; qsqdd1[i] = QsqDiff[i]; qsqdd1err[i] = QsqDifferr[i]; 
                thtgavg1[i] = ThtgAvg[i]; thtgavg1err[i] = ThtgAvgerr[i]; thtgdd1[i] = ThtgDiff[i]; thtgdd1err[i] = ThtgDifferr[i];
                phtgavg1[i] = PhtgAvg[i]; phtgavg1err[i] = PhtgAvgerr[i]; phtgdd1[i] = PhtgDiff[i]; phtgdd1err[i] = PhtgDifferr[i];
                asymavg1[i] = AsymAvg[i]; asymavg1err[i] = AsymAvgerr[i]; asymdd1[i] = AsymDiff[i]; asymdd1err[i] = AsymDifferr[i];  
                veravg1[i] = VerMomAvg[i]; veravg1err[i] = VerMomAvgerr[i]; verdd1[i] = VerMomDiff[i]; verdd1err[i] = VerMomDifferr[i];
                horavg1[i] = HorMomAvg[i]; horavg1err[i] = HorMomAvgerr[i]; hordd1[i] = HorMomDiff[i]; hordd1err[i] = HorMomDifferr[i];
                cosavg1[i] = CosAvg[i]; cosavg1err[i] = CosAvgerr[i]; cosdd1[i] = CosDiff[i]; cosdd1err[i] = CosDifferr[i];
                sinavg1[i] = SinAvg[i]; sinavg1err[i] = SinAvgerr[i]; sindd1[i] = SinDiff[i]; sindd1err[i] = SinDifferr[i];
                qsL1[i] = QsqL[i]; qsL1er[i] = QsqerL[i]; qsR1[i] = QsqR[i]; qsR1er[i] = QsqerR[i]; 
                thsL1[i] = ThtgL[i]; thsL1er[i] = ThtgerL[i]; thsR1[i] = ThtgR[i]; thsR1er[i] = ThtgerR[i];
                phsL1[i] = PhtgL[i]; phsL1er[i] = PhtgerL[i]; phsR1[i] = PhtgR[i]; phsR1er[i] = PhtgerR[i];
                asyL1[i] = asymL[i]; asyL1er[i] = asymerL[i]; asyR1[i] = asymR[i]; asyR1er[i] = asymerR[i];
                labL1[i] = AngleL[i]; labL1er[i] = AngleerL[i]; labR1[i] = AngleR[i]; labR1er[i] = AngleerR[i];
                verL1[i] = VermomL[i]; verL1er[i] = VermomerL[i]; verR1[i] = VermomR[i]; verR1er[i] = VermomerR[i];
                horL1[i] = HormomL[i]; horL1er[i] = HormomerL[i]; horR1[i] = HormomR[i]; horR1er[i] = HormomerR[i];
                cosL1[i] = cosphL[i]; cosL1er[i] = cospherL[i]; cosR1[i] = cosphR[i]; cosR1er[i] = cospherR[i];
                sinL1[i] = sinphL[i]; sinL1er[i] = sinpherL[i]; sinR1[i] = sinphR[i]; sinR1er[i] = sinpherR[i];
                xi1[i] = Xi[i]; xi1err[i] = Xierr[i];

              }
     if(i == 3){              
                qsqavg2[ii] = QsqAvg[i]; qsqavg2err[ii] = QsqAvgerr[i]; qsqdd2[ii] = QsqDiff[i]; qsqdd2err[ii] = QsqDifferr[i];
                thtgavg2[ii] = ThtgAvg[i]; thtgavg2err[ii] = ThtgAvgerr[i]; thtgdd2[ii] = ThtgDiff[i]; thtgdd2err[ii] = ThtgDifferr[i];
                phtgavg2[ii] = PhtgAvg[i]; phtgavg2err[ii] = PhtgAvgerr[i]; phtgdd2[ii] = PhtgDiff[i]; phtgdd2err[ii] = PhtgDifferr[i];
                asymavg2[ii] = AsymAvg[i]; asymavg2err[ii] = AsymAvgerr[i]; asymdd2[ii] = AsymDiff[i]; asymdd2err[ii] = AsymDifferr[i];
                veravg2[ii] = VerMomAvg[i]; veravg2err[ii] = VerMomAvgerr[i]; verdd2[ii] = VerMomDiff[i]; verdd2err[ii] = VerMomDifferr[i];
                horavg2[ii] = HorMomAvg[i]; horavg2err[ii] = HorMomAvgerr[i]; hordd2[ii] = HorMomDiff[i]; hordd2err[ii] = HorMomDifferr[i];
                qsL2[ii] = QsqL[i]; qsL2er[ii] = QsqerL[i]; qsR2[ii] = QsqR[i]; qsR2er[ii] = QsqerR[i];
                thsL2[ii] = ThtgL[i]; thsL2er[ii] = ThtgerL[i]; thsR2[ii] = ThtgR[i]; thsR2er[ii] = ThtgerR[i];
                phsL2[ii] = PhtgL[i]; phsL2er[ii] = PhtgerL[i]; phsR2[ii] = PhtgR[i]; phsR2er[ii] = PhtgerR[i];
                asyL2[ii] = asymL[i]; asyL2er[ii] = asymerL[i]; asyR2[ii] = asymR[i]; asyR2er[ii] = asymerR[i];
                labL2[ii] = AngleL[i]; labL2er[ii] = AngleerL[i]; labR2[ii] = AngleR[i]; labR2er[ii] = AngleerR[i]; 
                verL2[ii] = VermomL[i]; verL2er[ii] = VermomerL[i]; verR2[ii] = VermomR[i]; verR2er[ii] = VermomerR[i];
                horL2[ii] = HormomL[i]; horL2er[ii] = HormomerL[i]; horR2[ii] = HormomR[i]; horR2er[ii] = HormomerR[i];
                cosL2[ii] = cosphL[i]; cosL2er[ii] = cospherL[i]; cosR2[ii] = cosphR[i]; cosR2er[ii] = cospherR[i];
                sinL2[ii] = sinphL[i]; sinL2er[ii] = sinpherL[i]; sinR2[ii] = sinphR[i]; sinR2er[ii] = sinpherR[i];
                cosavg2[ii] = CosAvg[i]; cosavg2err[ii] = CosAvgerr[i]; cosdd2[ii] = CosDiff[i]; cosdd2err[ii] = CosDifferr[i];
                sinavg2[ii] = SinAvg[i]; sinavg2err[ii] = SinAvgerr[i]; sindd2[ii] = SinDiff[i]; sindd2err[ii] = SinDifferr[i];
                xi2[ii] = Xi[i]; xi2err[ii] = Xierr[i];                
                ii++;
               }
     if( i >= 4) {  qsqavg3[iii] = QsqAvg[i]; qsqavg3err[iii] = QsqAvgerr[i]; qsqdd3[iii] = QsqDiff[i]; qsqdd3err[iii] = QsqDifferr[i];
            thtgavg3[iii] = ThtgAvg[i]; thtgavg3err[iii] = ThtgAvgerr[i]; thtgdd3[iii] = ThtgDiff[i]; thtgdd3err[iii] = ThtgDifferr[i];
            phtgavg3[iii] = PhtgAvg[i]; phtgavg3err[iii] = PhtgAvgerr[i]; phtgdd3[iii] = PhtgDiff[i]; phtgdd3err[iii] = PhtgDifferr[i];
            asymavg3[iii] = AsymAvg[i]; asymavg3err[iii] = AsymAvgerr[i]; asymdd3[iii] = AsymDiff[i]; asymdd3err[iii] = AsymDifferr[i];
            veravg3[iii] = VerMomAvg[i]; veravg3err[iii] = VerMomAvgerr[i]; verdd3[iii] = VerMomDiff[i]; verdd3err[iii] = VerMomDifferr[i];
            horavg3[iii] = HorMomAvg[i]; horavg3err[iii] = HorMomAvgerr[i]; hordd3[iii] = HorMomDiff[i]; hordd3err[iii] = HorMomDifferr[i];
            cosavg3[iii] = CosAvg[i]; cosavg3err[iii] = CosAvgerr[i]; cosdd3[iii] = CosDiff[i]; cosdd3err[iii] = CosDifferr[i];
            sinavg3[iii] = SinAvg[i]; sinavg3err[iii] = SinAvgerr[i]; sindd3[iii] = SinDiff[i]; sindd3err[iii] = SinDifferr[i];
            qsL3[iii] = QsqL[i]; qsL3er[iii] = QsqerL[i]; qsR3[iii] = QsqR[i]; qsR3er[iii] = QsqerR[i];
            thsL3[iii] = ThtgL[i]; thsL3er[iii] = ThtgerL[i]; thsR3[iii] = ThtgR[i]; thsR3er[iii] = ThtgerR[i];
            phsL3[iii] = PhtgL[i]; phsL3er[iii] = PhtgerL[i]; phsR3[iii] = PhtgR[i]; phsR3er[iii] = PhtgerR[i];
            asyL3[iii] = asymL[i]; asyL3er[iii] = asymerL[i]; asyR3[iii] = asymR[i]; asyR3er[iii] = asymerR[i];
            labL3[iii] = AngleL[i]; labL3er[iii] = AngleerL[i]; labR3[iii] = AngleR[i]; labR3er[iii] = AngleerR[i]; 
            verL3[iii] = VermomL[i]; verL3er[iii] = VermomerL[i]; verR3[iii] = VermomR[i]; verR3er[iii] = VermomerR[i];
            horL3[iii] = HormomL[i]; horL3er[iii] = HormomerL[i]; horR3[iii] = HormomR[i]; horR3er[iii] = HormomerR[i];
            cosL3[iii] = cosphL[i]; cosL3er[ii] = cospherL[i]; cosR3[iii] = cosphR[i]; cosR3er[ii] = cospherR[i];
            sinL3[iii] = sinphL[i]; sinL3er[iii] = sinpherL[i]; sinR3[iii] = sinphR[i]; sinR3er[iii] = sinpherR[i];
            xi3[iii] = Xi[i]; xi3err[iii] = Xierr[i];
            iii++;
               }

 
  }

   Xifac = new TMultiGraph();  
   Xifac->SetTitle("#xi(Left/Right Apparatus Asymmetry); Run Number; #xi");
   Xifac->SetMinimum(-0.5); Xifac->SetMaximum(0.5);   

   Q2_L = new TMultiGraph();
   Q2_L->SetTitle("LHRS Q^{2} (GeV/c)^{2}; Run Number; Q^{2} (GeV/c)^{2}");
   Q2_L->SetMinimum(0.029); Q2_L->SetMaximum(0.032);

   Q2_R = new TMultiGraph();
   Q2_R->SetTitle("RHRS Q^{2} (GeV/c)^{2}; Run Number; Q^{2} (GeV/c)^{2}");
   Q2_R->SetMinimum(0.029); Q2_R->SetMaximum(0.032);

   Q2Av = new TMultiGraph();
   Q2Av->SetTitle("Average Q^{2} (GeV/c)^{2}; Run Number; Q^{2} (GeV/c)^{2}");
   Q2Av->SetMinimum(0.029); Q2Av->SetMaximum(0.032);

   Q2DD = new TMultiGraph();
   Q2DD->SetTitle("Q^{2} Double Difference (GeV/c)^{2}; Run Number; Q^{2} (GeV/c)^{2}");
   Q2DD->SetMinimum(-0.002); Q2DD->SetMaximum(0.002);

   Ph_L = new TMultiGraph();
   Ph_L->SetTitle("LHRS <#phi_{tg}>; Run Number; <#phi_{tg}>");
   Ph_L->SetMinimum(-0.01); Ph_L->SetMaximum(0);

   Ph_R = new TMultiGraph();
   Ph_R->SetTitle("RHRS <#phi_{tg}>; Run Number; <#phi_{tg}>");
   Ph_R->SetMinimum(0); Ph_R->SetMaximum(0.01);

   PhAv = new TMultiGraph();
   PhAv->SetTitle("Average <#phi_{tg}>; Run Number; <#phi_{tg}>");
   PhAv->SetMinimum(-5e-4); PhAv->SetMaximum(7e-4);

   PhDD = new TMultiGraph();
   PhDD->SetTitle("<#phi_{tg}> Double Difference; Run Number; <#phi_{tg}>");
   PhDD->SetMinimum(-0.01); PhDD->SetMaximum(0);

   Th_L = new TMultiGraph();
   Th_L->SetTitle("LHRS <#theta_{tg}>; Run Number; <#theta_{tg}>");
   Th_L->SetMinimum(-0.01); Th_L->SetMaximum(0.01);

   Th_R = new TMultiGraph();
   Th_R->SetTitle("RHRS <#theta_{tg}>; Run Number; <#theta_{tg}>");
   Th_R->SetMinimum(-0.01); Th_R->SetMaximum(0.01);

   ThAv = new TMultiGraph();
   ThAv->SetTitle("Average <#theta_{tg}>; Run Number; <#theta_{tg}>");
   //ThAv->SetMinimum(-0.01); ThAv->SetMaximum(0.01);

   ThDD = new TMultiGraph();
   ThDD->SetTitle("<#theta_{tg}> Double Difference; Run Number; <#theta_{tg}>");
   ThDD->SetMinimum(-0.005); ThDD->SetMaximum(0.005);

   APV_L = new TMultiGraph();
   APV_L->SetTitle("LHRS <A_{PV}> (ppm); Run Number; <A_{PV}> (ppm)");
   APV_L->SetMinimum(2.4); APV_L->SetMaximum(2.44);

   APV_R = new TMultiGraph();
   APV_R->SetTitle("RHRS <A_{PV}> (ppm); Run Number; <A_{PV}> (ppm)");
   APV_R->SetMinimum(2.4); APV_R->SetMaximum(2.44);

   APVAvg = new TMultiGraph();
   APVAvg->SetTitle("Average <A_{PV}> (ppm); Run Number; <A_{PV}> (ppm)");
   APVAvg->SetMinimum(2.4); APVAvg->SetMaximum(2.44);

   APVDD = new TMultiGraph();
   APVDD->SetTitle("<A_{PV}> Double Difference (ppm); Run Number; <A_{PV}> (ppm)");
   APVDD->SetMinimum(-0.01); APVDD->SetMaximum(0.01);

   Lab_L = new TMultiGraph();
   Lab_L->SetTitle("LHRS #theta_{lab} (deg); Run Number; #theta_{lab} (deg)");
 //  Lab_L->SetMinimum(4.5); Lab_L->SetMaximum(4.56);

   Lab_R = new TMultiGraph();
   Lab_R->SetTitle("RHRS #theta_{lab} (deg); Run Number; #theta_{lab} (deg)");
//   Lab_R->SetMinimum(4.54); Lab_R->SetMaximum(4.62);

   Ver_L = new TMultiGraph();
   Ver_L->SetTitle("LHRS sin(#theta)sin(#phi); Run Number; sin(#theta)sin(#phi)");
//   Ver_L->SetMinimum(0.05); Ver_L->SetMaximum(0.1);

   Ver_R = new TMultiGraph();
   Ver_R->SetTitle("RHRS sin(#theta)sin(#phi); Run Number; sin(#theta)sin(#phi)");
//   Ver_R->SetMinimum(0.05); Ver_R->SetMaximum(0.1);

   VerAv = new TMultiGraph();
   VerAv->SetTitle("LHRS and RHRS Average sin(#theta)sin(#phi); Run Number; sin(#theta)sin(#phi)");
//   VerAv->SetMinimum(0.05); VerAv->SetMaximum(0.1);

   VerDD = new TMultiGraph();
   VerDD->SetTitle("LHRS and RHRS sin(#theta)sin(#phi) Double Difference; Run Number; sin(#theta)sin(#phi)");
//   VerDD->SetMinimum(-7e-4); VerDD->SetMaximum(0);

   Hor_L = new TMultiGraph();
   Hor_L->SetTitle("LHRS sin(#theta)cos(#phi); Run Number; sin(#theta)cos(#phi)");
   Hor_L->SetMinimum(-0.005); Hor_L->SetMaximum(0);  
 
   Hor_R = new TMultiGraph();
   Hor_R->SetTitle("RHRS sin(#theta)cos(#phi); Run Number; sin(#theta)cos(#phi)");
   Hor_R->SetMinimum(0); Hor_R->SetMinimum(0.005);

   HorAv = new TMultiGraph();
   HorAv->SetTitle("LHRS and RHRS Average sin(#theta)cos(#phi); Run Number; sin(#theta)cos(#phi)");
   HorAv->SetMinimum(-5e-4); HorAv->SetMaximum(5e-4);

   HorDD = new TMultiGraph();
   HorDD->SetTitle("LHRS and RHRS sin(#theta)cos(#phi) Double Difference; Run Number; sin(#theta)cos(#phi)");
   HorDD->SetMinimum(-0.005); HorDD->SetMaximum(0);

   Cos_L = new TMultiGraph();
   Cos_L->SetTitle("LHRS cos(#phi); Run Number; cos(#phi)");
   Cos_L->SetMinimum(-0.06); Cos_L->SetMaximum(0);

   Cos_R = new TMultiGraph();
   Cos_R->SetTitle("RHRS cos(#phi); Run Number; cos(#phi)");
   Cos_R->SetMinimum(0); Cos_R->SetMaximum(0.06); 

   CosAv = new TMultiGraph();
   CosAv->SetTitle("LHRS and RHRS Average cos(#phi); Run Number; cos(#phi)");
   CosAv->SetMinimum(-0.006); CosAv->SetMaximum(0.006);

   CosDD = new TMultiGraph();
   CosDD->SetTitle("LHRS and RHRS cos(#phi) Double Difference; Run Number; cos(#phi)");
   CosDD->SetMinimum(-0.08); CosDD->SetMaximum(-0.02);
  

   Sin_L = new TMultiGraph();
   Sin_L->SetTitle("LHRS sin(#phi); Run Number; sin(#phi)");
//   Sin_L->SetMinimum(0.9);Sin_L->SetMaximum(1.05);

   Sin_R = new TMultiGraph();
   Sin_R->SetTitle("RHRS sin(#phi); Run Number; sin(#phi)");
//   Sin_R->SetMinimum(0.9); Sin_R->SetMaximum(1.05);

   SinAv = new TMultiGraph();
   SinAv->SetTitle("LHRS and RHRS Average sin(#phi); Run Number; sin(#phi)");
//   SinAv->SetMinimum(0.9); SinAv->SetMaximum(1.05);

   SinDD = new TMultiGraph();
   SinDD->SetTitle("LHRS and RHRS sin(#phi) Double Difference; Run Number; sin(#phi)");
//   SinDD->SetMinimum(-0.006); SinDD->SetMaximum(0.006);



     //Detector Combos
     TGraphErrors *xiwien[3], *qavwien[3], *qdwien[3], *thavwien[3], *thdwien[3], *phavwien[3], *phdwien[3], *asavwien[3], *asdwien[3];
     TGraphErrors *veravwien[3], *verdwien[3], *horavwien[3], *hordwien[3];
     TGraphErrors *sinavwien[3], *sindwien[3], *cosavwien[3], *cosdwien[3];
     //Individual detectors 
     TGraphErrors *qLwien[3], *qRwien[3], *thLwien[3], *thRwien[3], *phLwien[3], *phRwien[3], *asLwien[3], *asRwien[3], *lbLwien[3], *lbRwien[3]; 
     TGraphErrors *verLwien[3], *verRwien[4], *horLwien[3], *horRwien[3];
     TGraphErrors *sinLwien[3], *sinRwien[3], *cosLwien[3], *cosRwien[3];


     xiwien[0] = new TGraphErrors(3,wien1,xi1,wien1err,xi1err); xiwien[1] = new TGraphErrors(1,wien2,xi2,wien2err,xi2err); 
     xiwien[2] = new TGraphErrors(5,wien3,xi3,wien3err,xi3err);
   
     qavwien[0] = new TGraphErrors(3,wien1,qsqavg1,wien1err,qsqavg1err); qdwien[0] = new TGraphErrors(3,wien1,qsqdd1,wien1err,qsqdd1err);
     qavwien[1] = new TGraphErrors(1,wien2,qsqavg2,wien2err,qsqavg2err); qdwien[1] = new TGraphErrors(1,wien2,qsqdd2,wien2err,qsqdd2err); 
     qavwien[2] = new TGraphErrors(5,wien3,qsqavg3,wien3err,qsqavg3err); qdwien[2] = new TGraphErrors(5,wien3,qsqdd3,wien3err,qsqdd3err);

     thavwien[0] = new TGraphErrors(3,wien1,thtgavg1,wien1err,thtgavg1err); thdwien[0] = new TGraphErrors(3,wien1,thtgdd1,wien1err,thtgdd1err);
     thavwien[1] = new TGraphErrors(1,wien2,thtgavg2,wien2err,thtgavg2err); thdwien[1] = new TGraphErrors(1,wien2,thtgdd2,wien2err,thtgdd2err);
     thavwien[2] = new TGraphErrors(5,wien3,thtgavg3,wien3err,thtgavg3err); thdwien[2] = new TGraphErrors(5,wien3,thtgdd3,wien3err,thtgdd3err);

     phavwien[0] = new TGraphErrors(3,wien1,phtgavg1,wien1err,phtgavg1err); phdwien[0] = new TGraphErrors(3,wien1,phtgdd1,wien1err,phtgdd1err);
     phavwien[1] = new TGraphErrors(1,wien2,phtgavg2,wien2err,phtgavg2err); phdwien[1] = new TGraphErrors(1,wien2,phtgdd2,wien2err,phtgdd2err);
     phavwien[2] = new TGraphErrors(5,wien3,phtgavg3,wien3err,phtgavg3err); phdwien[2] = new TGraphErrors(5,wien3,phtgdd3,wien3err,phtgdd3err);

     asavwien[0] = new TGraphErrors(3,wien1,asymavg1,wien1err,asymavg1err); asdwien[0] = new TGraphErrors(3,wien1,asymdd1,wien1err,asymdd1err);
     asavwien[1] = new TGraphErrors(1,wien2,asymavg2,wien2err,asymavg2err); asdwien[1] = new TGraphErrors(1,wien2,asymdd2,wien2err,asymdd2err);
     asavwien[2] = new TGraphErrors(5,wien3,asymavg3,wien3err,asymavg3err); asdwien[2] = new TGraphErrors(5,wien3,asymdd3,wien3err,asymdd3err);

    veravwien[0] = new TGraphErrors(3,wien1,veravg1,wien1err,veravg1err); verdwien[0] = new TGraphErrors(3,wien1,verdd1,wien1err,verdd1err);
    veravwien[1] = new TGraphErrors(1,wien2,veravg2,wien2err,veravg2err); verdwien[1] = new TGraphErrors(1,wien2,verdd2,wien2err,verdd2err);
    veravwien[2] = new TGraphErrors(5,wien3,veravg3,wien3err,veravg3err); verdwien[2] = new TGraphErrors(5,wien3,verdd3,wien3err,verdd3err);

    horavwien[0] = new TGraphErrors(3,wien1,horavg1,wien1err,horavg1err); hordwien[0] = new TGraphErrors(3,wien1,hordd1,wien1err,hordd1err);
    horavwien[1] = new TGraphErrors(1,wien2,horavg2,wien2err,horavg2err); hordwien[1] = new TGraphErrors(1,wien2,hordd2,wien2err,hordd2err);
    horavwien[2] = new TGraphErrors(5,wien3,horavg3,wien3err,horavg3err); hordwien[2] = new TGraphErrors(5,wien3,hordd3,wien3err,hordd3err);

    cosavwien[0] = new TGraphErrors(3,wien1,cosavg1,wien1err,cosavg1err); cosdwien[0] = new TGraphErrors(3,wien1,cosdd1,wien1err,cosdd1err);
    cosavwien[1] = new TGraphErrors(1,wien2,cosavg2,wien2err,cosavg2err); cosdwien[1] = new TGraphErrors(1,wien2,cosdd2,wien2err,cosdd2err);
    cosavwien[2] = new TGraphErrors(5,wien3,cosavg3,wien3err,cosavg3err); cosdwien[2] = new TGraphErrors(5,wien3,cosdd3,wien3err,cosdd3err);

    sinavwien[0] = new TGraphErrors(3,wien1,sinavg1,wien1err,sinavg1err); sindwien[0] = new TGraphErrors(3,wien1,sindd1,wien1err,sindd1err);
    sinavwien[1] = new TGraphErrors(1,wien2,sinavg2,wien2err,sinavg2err); sindwien[1] = new TGraphErrors(1,wien2,sindd2,wien2err,sindd2err);
    sinavwien[2] = new TGraphErrors(5,wien3,sinavg3,wien3err,sinavg3err); sindwien[2] = new TGraphErrors(5,wien3,sindd3,wien3err,sindd3err);


     qLwien[0] = new TGraphErrors(3,wien1,qsL1,wien1err,qsL1er); qRwien[0] = new TGraphErrors(3,wien1,qsR1,wien1err,qsR1er);
     qLwien[1] = new TGraphErrors(1,wien2,qsL2,wien2err,qsL2er); qRwien[1] = new TGraphErrors(1,wien2,qsR2,wien2err,qsR2er);  
     qLwien[2] = new TGraphErrors(5,wien3,qsL3,wien3err,qsL3er); qRwien[2] = new TGraphErrors(5,wien3,qsR3,wien3err,qsR3er);

     thLwien[0] = new TGraphErrors(3,wien1,thsL1,wien1err,thsL1er); thRwien[0] = new TGraphErrors(3,wien1,thsR1,wien1err,thsR1er);
     thLwien[1] = new TGraphErrors(1,wien2,thsL2,wien2err,thsL2er); thRwien[1] = new TGraphErrors(1,wien2,thsR2,wien2err,thsR2er);
     thLwien[2] = new TGraphErrors(5,wien3,thsL3,wien3err,thsL3er); thRwien[2] = new TGraphErrors(5,wien3,thsR3,wien3err,thsR3er);

     phLwien[0] = new TGraphErrors(3,wien1,phsL1,wien1err,phsL1er); phRwien[0] = new TGraphErrors(3,wien1,phsR1,wien1err,phsR1er);
     phLwien[1] = new TGraphErrors(1,wien2,phsL2,wien2err,phsL2er); phRwien[1] = new TGraphErrors(1,wien2,phsR2,wien2err,phsR2er);
     phLwien[2] = new TGraphErrors(5,wien3,phsL3,wien3err,phsL3er); phRwien[2] = new TGraphErrors(5,wien3,phsR3,wien3err,phsR3er);

     asLwien[0] = new TGraphErrors(3,wien1,asyL1,wien1err,asyL1er); asRwien[0] = new TGraphErrors(3,wien1,asyR1,wien1err,asyR1er);
     asLwien[1] = new TGraphErrors(1,wien2,asyL2,wien2err,asyL2er); asRwien[1] = new TGraphErrors(1,wien2,asyR2,wien2err,asyR2er);
     asLwien[2] = new TGraphErrors(5,wien3,asyL3,wien3err,asyL3er); asRwien[2] = new TGraphErrors(5,wien3,asyR3,wien3err,asyR3er);

     lbLwien[0] = new TGraphErrors(3,wien1,labL1,wien1err,labL1er); lbRwien[0] = new TGraphErrors(3,wien1,labR1,wien1err,labR1er);
     lbLwien[1] = new TGraphErrors(1,wien2,labL2,wien2err,labL2er); lbRwien[1] = new TGraphErrors(1,wien2,labR2,wien2err,labR2er);
     lbLwien[2] = new TGraphErrors(5,wien3,labL3,wien3err,labL3er); lbRwien[2] = new TGraphErrors(5,wien3,labR3,wien3err,labR3er);

     verLwien[0] = new TGraphErrors(3,wien1,verL1,wien1err,verL1er); verRwien[0] = new TGraphErrors(3,wien1,verR1,wien1err,verR1er);
     verLwien[1] = new TGraphErrors(1,wien2,verL2,wien2err,verL2er); verRwien[1] = new TGraphErrors(1,wien2,verR2,wien2err,verR2er);
     verLwien[2] = new TGraphErrors(5,wien3,verL3,wien3err,verL3er); verRwien[2] = new TGraphErrors(5,wien3,verR3,wien3err,verR3er);

     horLwien[0] = new TGraphErrors(3,wien1,horL1,wien1err,horL1er); horRwien[0] = new TGraphErrors(3,wien1,horR1,wien1err,horR1er);
     horLwien[1] = new TGraphErrors(1,wien2,horL2,wien2err,horL2er); horRwien[1] = new TGraphErrors(1,wien2,horR2,wien2err,horR2er);
     horLwien[2] = new TGraphErrors(5,wien3,horL3,wien3err,horL3er); horRwien[2] = new TGraphErrors(5,wien3,horR3,wien3err,horR3er);

     cosLwien[0] = new TGraphErrors(3,wien1,cosL1,wien1err,cosL1er); cosRwien[0] = new TGraphErrors(3,wien1,cosR1,wien1err,cosR1er);
     cosLwien[1] = new TGraphErrors(1,wien2,cosL2,wien2err,cosL2er); cosRwien[1] = new TGraphErrors(1,wien2,cosR2,wien2err,cosR2er);
     cosLwien[2] = new TGraphErrors(5,wien3,cosL3,wien3err,cosL3er); cosRwien[2] = new TGraphErrors(5,wien3,cosR3,wien3err,cosR3er);

     sinLwien[0] = new TGraphErrors(3,wien1,sinL1,wien1err,sinL1er); sinRwien[0] = new TGraphErrors(3,wien1,sinR1,wien1err,sinR1er);
     sinLwien[1] = new TGraphErrors(1,wien2,sinL2,wien2err,sinL2er); sinRwien[1] = new TGraphErrors(1,wien2,sinR2,wien2err,sinR2er);
     sinLwien[2] = new TGraphErrors(5,wien3,sinL3,wien3err,sinL3er); sinRwien[2] = new TGraphErrors(5,wien3,sinR3,wien3err,sinR3er);


     for(int i = 0; i < 3; i++){ 
      xiwien[i]->SetMarkerStyle(21); xiwien[i]->SetMarkerColor(color[i]);
      qavwien[i]->SetMarkerStyle(21); qavwien[i]->SetMarkerColor(color[i]);
      qdwien[i]->SetMarkerStyle(21); qdwien[i]->SetMarkerColor(color[i]);
      thavwien[i]->SetMarkerStyle(21); thavwien[i]->SetMarkerColor(color[i]);
      thdwien[i]->SetMarkerStyle(21); thdwien[i]->SetMarkerColor(color[i]);
      phavwien[i]->SetMarkerStyle(21); phavwien[i]->SetMarkerColor(color[i]);
      phdwien[i]->SetMarkerStyle(21); phdwien[i]->SetMarkerColor(color[i]);
      asavwien[i]->SetMarkerStyle(21); asavwien[i]->SetMarkerColor(color[i]);
      asdwien[i]->SetMarkerStyle(21); asdwien[i]->SetMarkerColor(color[i]);
      veravwien[i]->SetMarkerStyle(21); veravwien[i]->SetMarkerColor(color[i]);
      verdwien[i]->SetMarkerStyle(21); verdwien[i]->SetMarkerColor(color[i]);
      horavwien[i]->SetMarkerStyle(21); horavwien[i]->SetMarkerColor(color[i]);
      hordwien[i]->SetMarkerStyle(21); hordwien[i]->SetMarkerColor(color[i]);
      cosavwien[i]->SetMarkerStyle(21); cosavwien[i]->SetMarkerColor(color[i]);
      cosdwien[i]->SetMarkerStyle(21); cosdwien[i]->SetMarkerColor(color[i]);
      sinavwien[i]->SetMarkerStyle(21); sinavwien[i]->SetMarkerColor(color[i]);
      sindwien[i]->SetMarkerStyle(21); sindwien[i]->SetMarkerColor(color[i]);
      qLwien[i]->SetMarkerStyle(21); qLwien[i]->SetMarkerColor(color[i]);
      qRwien[i]->SetMarkerStyle(21); qRwien[i]->SetMarkerColor(color[i]);
      thLwien[i]->SetMarkerStyle(21); thLwien[i]->SetMarkerColor(color[i]);
      thRwien[i]->SetMarkerStyle(21); thRwien[i]->SetMarkerColor(color[i]);
      phLwien[i]->SetMarkerStyle(21); phLwien[i]->SetMarkerColor(color[i]);
      phRwien[i]->SetMarkerStyle(21); phRwien[i]->SetMarkerColor(color[i]);
      asLwien[i]->SetMarkerStyle(21); asLwien[i]->SetMarkerColor(color[i]);
      asRwien[i]->SetMarkerStyle(21); asRwien[i]->SetMarkerColor(color[i]);
      lbLwien[i]->SetMarkerStyle(21); lbLwien[i]->SetMarkerColor(color[i]);
      lbRwien[i]->SetMarkerStyle(21); lbRwien[i]->SetMarkerColor(color[i]);
      verLwien[i]->SetMarkerStyle(21); verLwien[i]->SetMarkerColor(color[i]);
      verRwien[i]->SetMarkerStyle(21); verRwien[i]->SetMarkerColor(color[i]);
      horLwien[i]->SetMarkerStyle(21); horLwien[i]->SetMarkerColor(color[i]);
      horRwien[i]->SetMarkerStyle(21); horRwien[i]->SetMarkerColor(color[i]);
      cosLwien[i]->SetMarkerStyle(21); cosLwien[i]->SetMarkerColor(color[i]);
      cosRwien[i]->SetMarkerStyle(21); cosRwien[i]->SetMarkerColor(color[i]);
      sinLwien[i]->SetMarkerStyle(21); sinLwien[i]->SetMarkerColor(color[i]);
      sinRwien[i]->SetMarkerStyle(21); sinRwien[i]->SetMarkerColor(color[i]);



      Xifac->Add(xiwien[i]);
      Q2_L->Add(qLwien[i]); Q2_R->Add(qRwien[i]); 
      Lab_L->Add(lbLwien[i]); Lab_R->Add(lbRwien[i]); 
      Th_L->Add(thLwien[i]); Th_R->Add(thRwien[i]);
      Ph_L->Add(phLwien[i]); Ph_R->Add(phRwien[i]);
      APV_L->Add(asLwien[i]); APV_R->Add(asRwien[i]);
      Q2Av->Add(qavwien[i]); Q2DD->Add(qdwien[i]);
      ThAv->Add(thavwien[i]); ThDD->Add(thdwien[i]);
      PhAv->Add(phavwien[i]); PhDD->Add(phdwien[i]);
      APVAvg->Add(asavwien[i]); APVDD->Add(asdwien[i]);
      Ver_L->Add(verLwien[i]); Ver_R->Add(verRwien[i]);
      Hor_L->Add(horLwien[i]); Hor_R->Add(horRwien[i]); 
      VerAv->Add(veravwien[i]); VerDD->Add(verdwien[i]);
      HorAv->Add(horavwien[i]); HorDD->Add(hordwien[i]);
      Cos_L->Add(cosLwien[i]); Cos_R->Add(cosRwien[i]);
      Sin_L->Add(sinLwien[i]); Sin_R->Add(sinRwien[i]);
      CosAv->Add(cosavwien[i]); CosDD->Add(cosdwien[i]);
      SinAv->Add(sinavwien[i]); SinDD->Add(sindwien[i]);


     }
 
   //Covid lines -- hopefully this can removed once DB is fixed
   TLine *tq = new TLine(4.5,0.029,4.5,0.032); tq->SetLineStyle(2);
   TLine *tx = new TLine(4.5,-0.5,4.5,0.5); tx->SetLineStyle(2);
   TLine *tqd = new TLine(4.5,-0.002,4.5,0.002);tqd->SetLineStyle(2);
   TLine *tth = new TLine(4.5,-0.01,4.5,0.01); tth->SetLineStyle(2);
   TLine *tthd = new TLine(4.5,-0.005,4.5,0.005); tthd->SetLineStyle(2);
   TLine *tlab = new TLine(4.5,4.54,4.5,4.62); tlab->SetLineStyle(2);
   TLine *tlab1 = new TLine(4.5,4.5,4.5,4.56); tlab1->SetLineStyle(2);
   TLine *taL = new TLine(4.5,2.4,4.5,2.45); taL->SetLineStyle(2);
   TLine *ta = new TLine(4.5,2.40,4.5,2.44); ta->SetLineStyle(2);
   TLine *tad = new TLine(4.5,0,4.5,0.01); tad->SetLineStyle(2);
   TLine *tphL = new TLine(4.5,-0.01,4.5,0); tphL->SetLineStyle(2);
   TLine *tphR = new TLine(4.5,-0.01,4.5,0.01); tphR->SetLineStyle(2);
   TLine *tphA = new TLine(4.5,-5e-4,4.5,7e-4); tphA->SetLineStyle(2);
   TLine *tv = new TLine(4.5,0.05,4.5,0.1); tv->SetLineStyle(2);
   TLine *tvv = new TLine(4.5,-7e-4,4.5,0); tvv->SetLineStyle(2);
   TLine *thL = new TLine(4.5,-0.005,4.5,0); thL->SetLineStyle(2);
   TLine *thh = new TLine(4.5,-5e-4,4.5,5e-4); thh->SetLineStyle(2);
   TLine *thR = new TLine(4.5,0,4.5,0.005); thR->SetLineStyle(2);

   //Dithering line 
   TLine *dit = new TLine(1,-0.037,9,-0.037);


   //Legends
   TLegend *lav = new TLegend(0.55,0.7,0.9,0.9);
   TLegend *qdit = new TLegend(0.55,0.7,0.9,0.9);
 
  for(int i = 0; i < 3; i++) { if( i != 1 ) { lav->AddEntry(xiwien[i],"Wien Right","p"); qdit->AddEntry(xiwien[i],"Wien Right","p");} 
                        else { lav->AddEntry(xiwien[i], "Wien Left","p");  qdit->AddEntry(xiwien[i], "Wien Left","p"); } }

   lav->AddEntry(tq,"COVID","l");
   qdit->AddEntry(dit,"#xi = -0.037","l");
  

  //Now to draw and fit all of this
/*  TCanvas *cqL = new TCanvas();
  Q2_L->Draw("AP");
  Q2_L->Fit("pol0","W","",1,4);
  qLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fl1 = Q2_L->GetFunction("pol0");
  fl1->SetLineColor(1); fl1->SetLineWidth(1);
  TPaveText *ql1 = new TPaveText(1,0.029,4,0.0295);
  ql1->SetTextColor(kRed);
  ql1->AddText(Form("p0 = %0.5e#pm%0.1e",fl1->GetParameter(0),fl1->GetParError(0)));
  TF1 *fl2 = qLwien[2]->GetFunction("pol0");
  fl2->SetLineColor(1); fl2->SetLineWidth(1);
  TPaveText *ql2 = new TPaveText(5,0.029,9,0.0295);
  ql2->SetTextColor(kBlue);
  ql2->AddText(Form("p0 = %0.5e#pm%0.1e",fl2->GetParameter(0),fl2->GetParError(0)));
  ql2->SetTextAlign(12);
  ql1->Draw("sames");
  ql2->Draw("sames");
  tq->Draw();
  lav->Draw();
  cqL->SaveAs("Q2usLfit.png");

  TCanvas *cqR = new TCanvas();
  Q2_R->Draw("AP");
  Q2_R->Fit("pol0","W","",1,4);
  qRwien[2]->Fit("pol0","W","",5,9);
  TF1 *fr1 = Q2_R->GetFunction("pol0");
  fr1->SetLineColor(1); fr1->SetLineWidth(1);
  TPaveText *qr1 = new TPaveText(1,0.029,4,0.0295);
  qr1->SetTextColor(kRed);
  qr1->AddText(Form("p0 = %0.5e#pm%0.1e",fr1->GetParameter(0),fr1->GetParError(0)));
  TF1 *fr2 = qRwien[2]->GetFunction("pol0");
  fr2->SetLineColor(1); fr2->SetLineWidth(1);
  TPaveText *qr2 = new TPaveText(5,0.029,9,0.0295);
  qr2->SetTextColor(kBlue);
  qr2->AddText(Form("p0 = %0.5e#pm%0.1e",fr2->GetParameter(0),fr2->GetParError(0)));
  qr2->SetTextAlign(12);
  qr1->Draw("sames");
  qr2->Draw("sames");
  lav->Draw();
  tq->Draw();
  cqR->SaveAs("Q2usRfit.png");

  TCanvas *cthL = new TCanvas();
  Th_L->Draw("AP");
  Th_L->Fit("pol0","W","",1,4);
  thLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fthl1 = Th_L->GetFunction("pol0");
  fthl1->SetLineColor(1); fthl1->SetLineWidth(1);
  TPaveText *thl1 = new TPaveText(1,-0.01,4,-0.008);
  thl1->SetTextColor(kRed);
  thl1->AddText(Form("p0 = %0.5e#pm%0.1e",fthl1->GetParameter(0),fthl1->GetParError(0)));
  TF1 *fthl2 = thLwien[2]->GetFunction("pol0");
  fthl2->SetLineColor(1); fthl2->SetLineWidth(1);
  TPaveText *thl2 = new TPaveText(5,-0.01,9,-0.008);
  thl2->SetTextColor(kBlue);
  thl2->AddText(Form("p0 = %0.5e#pm%0.1e",fthl2->GetParameter(0),fthl2->GetParError(0)));
  thl2->SetTextAlign(12);
  thl1->Draw("sames");
  thl2->Draw("sames");
  lav->Draw();
  tth->Draw(); 
  cthL->SaveAs("ThusLfit.png"); 

  TCanvas *cthR = new TCanvas();
  Th_R->Draw("AP");
  Th_R->Fit("pol0","W","",1,4);
  thRwien[2]->Fit("pol0","W","",5,9);
  TF1 *fthr1 = Th_R->GetFunction("pol0");
  fthr1->SetLineColor(1); fthr1->SetLineWidth(1);
  TPaveText *thr1 = new TPaveText(1,-0.01,4,-0.008);
  thr1->SetTextColor(kRed);
  thr1->AddText(Form("p0 = %0.5e#pm%0.1e",fthr1->GetParameter(0),fthr1->GetParError(0)));
  TF1 *fthr2 = thRwien[2]->GetFunction("pol0");
  fthr2->SetLineColor(1); fthr2->SetLineWidth(1);
  TPaveText *thr2 = new TPaveText(5,-0.01,9,-0.008);
  thr2->SetTextColor(kBlue);
  thr2->AddText(Form("p0 = %0.5e#pm%0.1e",fthr2->GetParameter(0),fthr2->GetParError(0)));
  thr2->SetTextAlign(12);
  thr1->Draw("sames");
  thr2->Draw("sames");
  lav->Draw(); 
  tth->Draw();
  cthR->SaveAs("ThusRfit.png");


  TCanvas *cphL = new TCanvas();
  Ph_L->Draw("AP");
  Ph_L->Fit("pol0","W","",1,4);
  phLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fphl1 = Ph_L->GetFunction("pol0");
  fphl1->SetLineColor(1); fphl1->SetLineWidth(1);
  TPaveText *phl1 = new TPaveText(1,-0.01,4,-0.008);
  phl1->SetTextColor(kRed);
  phl1->AddText(Form("p0 = %0.5e#pm%0.1e",fphl1->GetParameter(0),fphl1->GetParError(0)));
  TF1 *fphl2 = thLwien[2]->GetFunction("pol0");
  fphl2->SetLineColor(1); fphl2->SetLineWidth(1);
  TPaveText *phl2 = new TPaveText(5,-0.01,9,-0.008);
  phl2->SetTextColor(kBlue);
  phl2->AddText(Form("p0 = %0.5e#pm%0.1e",fphl2->GetParameter(0),fphl2->GetParError(0)));
  phl2->SetTextAlign(12);
  phl1->Draw("sames"); 
  phl2->Draw("sames"); 
  lav->Draw();
  tphL->Draw();
  cphL->SaveAs("PhusLfit.png");

  TCanvas *cphR = new TCanvas();
  Ph_R->Draw("AP");
  Ph_R->Fit("pol0","W","",1,4);
  phLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fphr1 = Ph_R->GetFunction("pol0");
  fphr1->SetLineColor(1); fphr1->SetLineWidth(1);
  TPaveText *phr1 = new TPaveText(1,0.,4,0.002);
  phr1->SetTextColor(kRed);
  phr1->AddText(Form("p0 = %0.5e#pm%0.1e",fphr1->GetParameter(0),fphr1->GetParError(0)));
  TF1 *fphr2 = thLwien[2]->GetFunction("pol0");
  fphr2->SetLineColor(1); fphr2->SetLineWidth(1);
  TPaveText *phr2 = new TPaveText(5,0.,9,0.002);
  phr2->SetTextColor(kBlue);
  phr2->AddText(Form("p0 = %0.5e#pm%0.1e",fphr2->GetParameter(0),fphr2->GetParError(0)));
  phr2->SetTextAlign(12);
  phr1->Draw("sames"); 
  phr2->Draw("sames");
  lav->Draw();
  tphR->Draw();
  cphR->SaveAs("PhusRfit.png");


  TCanvas *clL = new TCanvas();
  Lab_L->Draw("AP");
  clL->SaveAs("LabangleL.png");

  //Lab_L->Fit("pol0","W","",1,4);
  
   // lbLwien[2]->Fit("pol0","W","",5,9);

  TF1 *flb1 = Lab_L->GetFunction("pol0");
  flb1->SetLineColor(1); flb1->SetLineWidth(1);
  TPaveText *lbl1 = new TPaveText(1,4.5,4,4.51);
  lbl1->SetTextColor(kRed);
  lbl1->AddText(Form("p0 = %0.5e#pm%0.1e",flb1->GetParameter(0),flb1->GetParError(0)));
  TF1 *flb2 = lbLwien[2]->GetFunction("pol0");
  flb2->SetLineColor(1); flb2->SetLineWidth(1);
  TPaveText *lbl2 = new TPaveText(5,4.5,9,4.51);
  lbl2->SetTextColor(kBlue);
  lbl2->AddText(Form("p0 = %0.5e#pm%0.1e",flb2->GetParameter(0),flb2->GetParError(0)));
  lbl2->SetTextAlign(12);
  lbl1->Draw("sames");
  lbl2->Draw("sames");
  lav->Draw();
  tlab->Draw();

 // clL->SaveAs("LabusLfit.png"); 

  TCanvas *clR = new TCanvas();
  Lab_R->Draw("AP");
  clR->SaveAs("LabangleR.png");


  Lab_R->Fit("pol0","W","",1,4);
  lbRwien[2]->Fit("pol0","W","",5,9);
  TF1 *flbR1 = Lab_R->GetFunction("pol0");
  flbR1->SetLineColor(1); flbR1->SetLineWidth(1);
  TPaveText *lbr1 = new TPaveText(1,4.54,4,4.55);
  lbr1->SetTextColor(kRed);
  lbr1->AddText(Form("p0 = %0.5e#pm%0.1e",flbR1->GetParameter(0),flbR1->GetParError(0)));
  TF1 *flbR2 = lbRwien[2]->GetFunction("pol0");
  flbR2->SetLineColor(1); flbR2->SetLineWidth(1);
  TPaveText *lbr2 = new TPaveText(5,4.54,9,4.55);//hopefully this can be changed
  lbr2->SetTextColor(kBlue);
  lbr2->AddText(Form("p0 = %0.5e#pm%0.1e",flbR2->GetParameter(0),flbR2->GetParError(0)));
  lbr2->SetTextAlign(12);
  lbr1->Draw("sames");
  lbr2->Draw("sames");
  lav->Draw();
  tlab->Draw();

//  clR->SaveAs("LabusRfit.png");

  TCanvas *caL = new TCanvas();
  APV_L->Draw("AP");
  APV_L->Fit("pol0","W","",1,4);
  asLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fasL1 = APV_L->GetFunction("pol0");
  fasL1->SetLineColor(1); fasL1->SetLineWidth(1);
  TPaveText *asl1 = new TPaveText(1,2.4,4,2.405);
  asl1->SetTextColor(kRed);
  asl1->AddText(Form("p0 = %0.5e#pm%0.1e",fasL1->GetParameter(0),fasL1->GetParError(0)));
  TF1 *fasL2 = asLwien[2]->GetFunction("pol0");
  fasL2->SetLineColor(1); fasL2->SetLineWidth(1);
  TPaveText *asl2 = new TPaveText(5,2.4,9,2.405);//hopefully this can be changed
  asl2->SetTextColor(kBlue);
  asl2->AddText(Form("p0 = %0.5e#pm%0.1e",fasL2->GetParameter(0),fasL2->GetParError(0)));
  asl2->SetTextAlign(12);
  asl1->Draw("sames"); 
  asl2->Draw("sames"); 
  lav->Draw();
  taL->Draw();
  caL->SaveAs("APVusLfit.png");

  TCanvas *caR = new TCanvas();
  APV_R->Draw("AP");
  APV_R->Fit("pol0","W","",1,4);
  asRwien[2]->Fit("pol0","W","",5,9);
  TF1 *fasR1 = APV_R->GetFunction("pol0");
  fasR1->SetLineColor(1); fasR1->SetLineWidth(1);
  TPaveText *asr1 = new TPaveText(1,2.4,4,2.405);
  asr1->SetTextColor(kRed);
  asr1->AddText(Form("p0 = %0.5e#pm%0.1e",fasR1->GetParameter(0),fasR1->GetParError(0)));
  TF1 *fasR2 = asRwien[2]->GetFunction("pol0");
  fasR2->SetLineColor(1); fasR2->SetLineWidth(1);
  TPaveText *asr2 = new TPaveText(5,2.4,9,2.405);//hopefully this can be changed
  asr2->SetTextColor(kBlue);
  asr2->AddText(Form("p0 = %0.5e#pm%0.1e",fasR2->GetParameter(0),fasR2->GetParError(0)));
  asr2->SetTextAlign(12);
  asr1->Draw("sames"); 
  asr2->Draw("sames"); 
  lav->Draw();
  ta->Draw();
  caR->SaveAs("APVusRfit.png");

  TCanvas *cqA = new TCanvas();
  Q2Av->Draw("AP");
  Q2Av->Fit("pol0","W","",1,4);
  qavwien[2]->Fit("pol0","W","",5,9);
  TF1 *fav1 = Q2Av->GetFunction("pol0");
  fav1->SetLineColor(1); fav1->SetLineWidth(1);
  TPaveText *qav1 = new TPaveText(1,0.029,4,0.0295);
  qav1->SetTextColor(kRed);
  qav1->AddText(Form("p0 = %0.5e#pm%0.1e",fav1->GetParameter(0),fav1->GetParError(0)));
  TF1 *fav2 = qavwien[2]->GetFunction("pol0");
  fav2->SetLineColor(1); fav2->SetLineWidth(1);
  TPaveText *qav2 = new TPaveText(5,0.029,9,0.0295);
  qav2->SetTextColor(kBlue);
  qav2->AddText(Form("p0 = %0.5e#pm%0.1e",fav2->GetParameter(0),fav2->GetParError(0)));
  qav2->SetTextAlign(12);
  qav1->Draw("sames");
  qav2->Draw("sames");
  lav->Draw(); 
  tq->Draw();
  cqA->SaveAs("Q2Avusfit.png");
 
  TCanvas *cqD = new TCanvas();
  Q2DD->Draw("AP");
  Q2DD->Fit("pol0","W","",1,4);
  qdwien[2]->Fit("pol0","W","",5,9);
  TF1 *fdd1 = Q2DD->GetFunction("pol0");
  fdd1->SetLineColor(1); fdd1->SetLineWidth(1);
  TPaveText *qdd1 = new TPaveText(1,-0.002,4,-0.001);
  qdd1->SetTextColor(kRed);
  qdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fdd1->GetParameter(0),fdd1->GetParError(0)));
  TF1 *fdd2 = qdwien[2]->GetFunction("pol0");
  fdd2->SetLineColor(1); fdd2->SetLineWidth(1);
  TPaveText *qdd2 = new TPaveText(5,-0.002,9,-0.001);
  qdd2->SetTextColor(kBlue);
  qdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fdd2->GetParameter(0),fdd2->GetParError(0)));
  qdd2->SetTextAlign(12);
  qdd1->Draw("sames");
  qdd2->Draw("sames");
  lav->Draw();
  tqd->Draw();
  cqD->SaveAs("Q2DDusfit.png");
  


  TCanvas *caA = new TCanvas();
  APVAvg->Draw("AP");
  APVAvg->Fit("pol0","W","",1,4);
  asavwien[2]->Fit("pol0","W","",5,9);
  TF1 *fasav1 = APVAvg->GetFunction("pol0");
  fasav1->SetLineColor(1); fasav1->SetLineWidth(1);
  TPaveText *asav1 = new TPaveText(1,2.4,4,2.405);
  asav1->SetTextColor(kRed);
  asav1->AddText(Form("p0 = %0.5e#pm%0.1e",fasav1->GetParameter(0),fasav1->GetParError(0)));
  TF1 *fasav2 = asavwien[2]->GetFunction("pol0");
  fasav2->SetLineColor(1); fasav2->SetLineWidth(1);
  TPaveText *asav2 = new TPaveText(5,2.4,9,2.405);
  asav2->SetTextColor(kBlue);
  asav2->AddText(Form("p0 = %0.5e#pm%0.1e",fasav2->GetParameter(0),fasav2->GetParError(0)));
  asav2->SetTextAlign(12);
  asav1->Draw("sames");
  asav2->Draw("sames");
  lav->Draw();
  ta->Draw();
  caA->SaveAs("APVusavfit.png");


  TCanvas *caD = new TCanvas();
  APVDD->Draw("AP");
  APVDD->Fit("pol0","W","",1,4);
  asdwien[2]->Fit("pol0","W","",5,9);
  TF1 *fasad1 = APVDD->GetFunction("pol0");
  fasad1->SetLineColor(1); fasad1->SetLineWidth(1);
  TPaveText *asad1 = new TPaveText(1,-0.01,4,-0.008);
  asad1->SetTextColor(kRed);
  asad1->AddText(Form("p0 = %0.5e#pm%0.1e",fasad1->GetParameter(0),fasad1->GetParError(0)));
  TF1 *fasad2 = asdwien[2]->GetFunction("pol0");
  fasad2->SetLineColor(1); fasad2->SetLineWidth(1);
  TPaveText *asad2 = new TPaveText(5,-0.01,9,-0.008);
  asad2->SetTextColor(kBlue);
  asad2->AddText(Form("p0 = %0.5e#pm%0.1e",fasad2->GetParameter(0),fasad2->GetParError(0)));
  asad2->SetTextAlign(12);
  asad1->Draw("sames");
  asad2->Draw("sames");
  lav->Draw();
  tad->Draw();
  caD->SaveAs("APVusddfit.png");
 


  TCanvas *ctA = new TCanvas();
  ThAv->Draw("AP");
  ThAv->Fit("pol0","W","",1,4);
  thavwien[2]->Fit("pol0","W","",5,9);
  TF1 *fthav1 = ThAv->GetFunction("pol0");
  fthav1->SetLineColor(1); fthav1->SetLineWidth(1);
  TPaveText *thav1 = new TPaveText(1,-0.01,4,-0.008);
  thav1->SetTextColor(kRed);
  thav1->AddText(Form("p0 = %0.5e#pm%0.1e",fthav1->GetParameter(0),fthav1->GetParError(0)));
  TF1 *fthav2 = thavwien[2]->GetFunction("pol0");
  fthav2->SetLineColor(1); fthav2->SetLineWidth(1);
  TPaveText *thav2 = new TPaveText(5,-0.01,9,-0.008);
  thav2->SetTextColor(kBlue);
  thav2->AddText(Form("p0 = %0.5e#pm%0.1e",fthav2->GetParameter(0),fthav2->GetParError(0)));
  thav2->SetTextAlign(12);
  thav1->Draw("sames");
  thav2->Draw("sames");
  lav->Draw();
  tth->Draw();
  ctA->SaveAs("Thusavfit.png");

  TCanvas *ctD = new TCanvas();
  ThDD->Draw("AP");
  ThDD->Fit("pol0","W","",1,4);
  thdwien[2]->Fit("pol0","W","",5,9);
  TF1 *fthdd1 = ThDD->GetFunction("pol0");
  fthdd1->SetLineColor(1); fthdd1->SetLineWidth(1);
  TPaveText *thdd1 = new TPaveText(1,-0.005,4,-0.003);
  thdd1->SetTextColor(kRed);
  thdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fthdd1->GetParameter(0),fthdd1->GetParError(0)));
  TF1 *fthdd2 = thdwien[2]->GetFunction("pol0");
  fthdd2->SetLineColor(1); fthdd2->SetLineWidth(1);
  TPaveText *thdd2 = new TPaveText(5,-0.005,9,-0.003);
  thdd2->SetTextColor(kBlue);
  thdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fthdd2->GetParameter(0),fthdd2->GetParError(0)));
  thdd2->SetTextAlign(12);
  thdd1->Draw("sames");
  thdd2->Draw("sames");
  lav->Draw();
  tthd->Draw();
  ctD->SaveAs("Thusddfit.png");


  TCanvas *cpA = new TCanvas();
  PhAv->Draw("AP");
  PhAv->Fit("pol0","W","",1,4);
  phavwien[2]->Fit("pol0","W","",5,9);
  TF1 *fphav1 = PhAv->GetFunction("pol0");
  fphav1->SetLineColor(1); fphav1->SetLineWidth(1);
  TPaveText *phav1 = new TPaveText(1,-5e-4,4,-3e-4);
  phav1->SetTextColor(kRed);
  phav1->AddText(Form("p0 = %0.5e#pm%0.1e",fphav1->GetParameter(0),fphav1->GetParError(0)));
  TF1 *fphav2 = phavwien[2]->GetFunction("pol0");
  fphav2->SetLineColor(1); fphav2->SetLineWidth(1);
  TPaveText *phav2 = new TPaveText(5,-5e-4,9,-3e-4);
  phav2->SetTextColor(kBlue);
  phav2->AddText(Form("p0 = %0.5e#pm%0.1e",fphav2->GetParameter(0),fphav2->GetParError(0)));
  phav2->SetTextAlign(12);
  phav1->Draw("sames");
  phav2->Draw("sames");
  lav->Draw();
  tphA->Draw();
 cpA->SaveAs("Phusavfit.png");

  TCanvas *cpD = new TCanvas();
  PhDD->Draw("AP");
  PhDD->Fit("pol0","W","",1,4);
  phdwien[2]->Fit("pol0","W","",5,9);
  TF1 *fphd1 = PhDD->GetFunction("pol0");
  fphd1->SetLineColor(1); fphd1->SetLineWidth(1);
  TPaveText *phd1 = new TPaveText(1,-0.01,4,-0.008);
  phd1->SetTextColor(kRed);
  phd1->AddText(Form("p0 = %0.5e#pm%0.1e",fphd1->GetParameter(0),fphd1->GetParError(0)));
  TF1 *fphd2 = phdwien[2]->GetFunction("pol0");
  fphd2->SetLineColor(1); fphd2->SetLineWidth(1);
  TPaveText *phd2 = new TPaveText(5,-0.01,9,-0.008);
  phd2->SetTextColor(kBlue);
  phd2->AddText(Form("p0 = %0.5e#pm%0.1e",fphd2->GetParameter(0),fphd2->GetParError(0)));
  phd2->SetTextAlign(12);
  phd1->Draw("sames");
  phd2->Draw("sames");
  lav->Draw();  
  tphL->Draw();
  cpD->SaveAs("Phusddfit.png");


  TCanvas *cxi = new TCanvas();
  Xifac->Draw("AP");
  trans()->Draw("F");
  dit->Draw();
  qdit->Draw();
  cxi->SaveAs("xiComp_thphmoment.png");



  TCanvas *cxi = new TCanvas();
  Xifac->Draw("AP");
  lav->Draw();
  cxi->SaveAs("xiMeasurements.png");


  TCanvas *cxi1 = new TCanvas();
  Xifac->Draw("AP");
  Xifac->Fit("pol0","W","",1,4);
  xiwien[2]->Fit("pol0","W","",5,9);
  TF1 *fx1 = Xifac->GetFunction("pol0");
  fx1->SetLineColor(1); fx1->SetLineWidth(1);
  TPaveText *qx1 = new TPaveText(1,-0.5,4,-0.4);
  qx1->SetTextColor(kRed);
  qx1->AddText(Form("p0 = %0.5e#pm%0.1e",fx1->GetParameter(0),fx1->GetParError(0)));
  TF1 *fx2 = xiwien[2]->GetFunction("pol0");
  fx2->SetLineColor(1); fx2->SetLineWidth(1);
  TPaveText *qx2 = new TPaveText(5,-0.5,9,-0.4);
  qx2->SetTextColor(kBlue);
  qx2->AddText(Form("p0 = %0.5e#pm%0.1e",fx2->GetParameter(0),fx2->GetParError(0)));
  qx2->SetTextAlign(12);
  lav->Draw();
  qx1->Draw("sames");
  qx2->Draw("sames");
  tx->Draw(); 
  cxi1->SaveAs("xiFits_phmoment.png");
  */


  TCanvas *c_verL = new TCanvas();
  Ver_L->Draw("AP");
/*
  Ver_L->Fit("pol0","W","",1,4);
  verLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fvL1 = Ver_L->GetFunction("pol0");
  fvL1->SetLineColor(1); fvL1->SetLineWidth(1);
  TPaveText *qvL1 = new TPaveText(1,0.05,4,0.055);
  qvL1->SetTextColor(kRed);
  qvL1->AddText(Form("p0 = %0.5e#pm%0.1e",fvL1->GetParameter(0),fvL1->GetParError(0)));
  TF1 *fvL2 = verLwien[2]->GetFunction("pol0");
  fvL2->SetLineColor(1); fvL2->SetLineWidth(1);
  TPaveText *qvL2 = new TPaveText(5,0.05,9,0.055);
  qvL2->SetTextColor(kBlue);
  qvL2->AddText(Form("p0 = %0.5e#pm%0.1e",fvL2->GetParameter(0),fvL2->GetParError(0)));
  lav->Draw();
  tv->Draw();
  qvL1->Draw("sames");
  qvL2->Draw("sames");
*/
  c_verL->SaveAs("VerLHRSnofit.png");
 
  TCanvas *c_verR = new TCanvas();
  Ver_R->Draw("AP");
/*
  Ver_R->Fit("pol0","W","",1,4);
  verRwien[2]->Fit("pol0","W","",5,9);
  TF1 *fvR1 = Ver_R->GetFunction("pol0");
  fvR1->SetLineColor(1); fvR1->SetLineWidth(1);
  TPaveText *qvR1 = new TPaveText(1,0.5,4,0.55);
  qvR1->SetTextColor(kRed);
  qvR1->AddText(Form("p0 = %0.5e#pm%0.1e",fvR1->GetParameter(0),fvR1->GetParError(0)));
  TF1 *fvR2 = verRwien[2]->GetFunction("pol0");
  fvR2->SetLineColor(1); fvR2->SetLineWidth(1);
  TPaveText *qvR2 = new TPaveText(5,0.05,9,0.055);
  qvR2->SetTextColor(kBlue);
  qvR2->AddText(Form("p0 = %0.5e#pm%0.1e",fvR2->GetParameter(0),fvR2->GetParError(0)));
  lav->Draw();
  tv->Draw();
  qvR1->Draw("sames");
  qvR2->Draw("sames");
*/
  c_verR->SaveAs("VerRHRSnofit.png");
 
  TCanvas *c_vav = new TCanvas();
  VerAv->Draw("AP");
/*
  VerAv->Fit("pol0","W","",1,4);
  veravwien[2]->Fit("pol0","W","",5,9);
  TF1 *fvav1 = VerAv->GetFunction("pol0");
  fvav1->SetLineColor(1); fvav1->SetLineWidth(1);
  TPaveText *qvav1 = new TPaveText(1,0.05,4,0.055);
  qvav1->SetTextColor(kRed);
  qvav1->AddText(Form("p0 = %0.5e#pm%0.1e",fvav1->GetParameter(0),fvav1->GetParError(0)));
  TF1 *fvav2 = veravwien[2]->GetFunction("pol0");
  fvav2->SetLineColor(1); fvav2->SetLineWidth(1);
  TPaveText *qvav2 = new TPaveText(5,0.05,9,0.055);
  qvav2->SetTextColor(kBlue);
  qvav2->AddText(Form("p0 = %0.5e#pm%0.1e",fvav2->GetParameter(0),fvav2->GetParError(0)));
  lav->Draw();
  tv->Draw();
  qvav1->Draw("sames");
  qvav2->Draw("sames");
*/
  c_vav->SaveAs("VerAvgnofit.png");


  TCanvas *c_vdd = new TCanvas();
  VerDD->Draw("AP");
/*
  VerDD->Fit("pol0","W","",1,4);
  verdwien[2]->Fit("pol0","W","",5,9);
  TF1 *fvdd1 = VerDD->GetFunction("pol0");
  fvdd1->SetLineColor(1); fvdd1->SetLineWidth(1);
  TPaveText *qvdd1 = new TPaveText(1,-7e-4,4,-6e-4);
  qvdd1->SetTextColor(kRed);
  qvdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fvdd1->GetParameter(0),fvdd1->GetParError(0)));
  TF1 *fvdd2 = verdwien[2]->GetFunction("pol0");
  fvdd2->SetLineColor(1); fvdd2->SetLineWidth(1);
  TPaveText *qvdd2 = new TPaveText(5,-7e-4,9,-6e-4);
  qvdd2->SetTextColor(kBlue);
  qvdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fvdd2->GetParameter(0),fvdd2->GetParError(0)));
  lav->Draw();
  tvv->Draw();
  qvdd1->Draw("sames");
  qvdd2->Draw("sames");
*/ 
  c_vdd->SaveAs("VerDDnofit.png");
/*

  TCanvas *c_horL = new TCanvas();
  Hor_L->Draw("AP");
  Hor_L->Fit("pol0","W","",1,4);
  horLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fhL1 = Hor_L->GetFunction("pol0");
  fhL1->SetLineColor(1); fhL1->SetLineWidth(1);
  TPaveText *qhL1 = new TPaveText(1,-0.005,4,-0.0045);
  qhL1->SetTextColor(kRed);
  qhL1->AddText(Form("p0 = %0.5e#pm%0.1e",fhL1->GetParameter(0),fhL1->GetParError(0)));
  TF1 *fhL2 = horLwien[2]->GetFunction("pol0");
  fhL2->SetLineColor(1); fhL2->SetLineWidth(1);
  TPaveText *qhL2 = new TPaveText(5,-0.005,9,-0.0045);
  qhL2->SetTextColor(kBlue);
  qhL2->AddText(Form("p0 = %0.5e#pm%0.1e",fhL2->GetParameter(0),fhL2->GetParError(0)));
  lav->Draw();
  thL->Draw();
  qhL1->Draw("sames");
  qhL2->Draw("sames");
  c_horL->SaveAs("HorLHRSfit.png");

  TCanvas *c_horR = new TCanvas();
  Hor_R->Draw("AP");
  Hor_R->Fit("pol0","W","",1,4);
  horRwien[2]->Fit("pol0","W","",5,9);
  TF1 *fhR1 = Hor_R->GetFunction("pol0");
  fhR1->SetLineColor(1); fhR1->SetLineWidth(1);
  TPaveText *qhR1 = new TPaveText(1,0,4,0.0005);
  qhR1->SetTextColor(kRed);
  qhR1->AddText(Form("p0 = %0.5e#pm%0.1e",fhR1->GetParameter(0),fhR1->GetParError(0)));
  TF1 *fhR2 = horRwien[2]->GetFunction("pol0");
  fhR2->SetLineColor(1); fhR2->SetLineWidth(1);
  TPaveText *qhR2 = new TPaveText(5,0.0,9,0.0005);
  qhR2->SetTextColor(kBlue);
  qhR2->AddText(Form("p0 = %0.5e#pm%0.1e",fhR2->GetParameter(0),fhR2->GetParError(0)));
  lav->Draw();
  thR->Draw();
  qhR1->Draw("sames");
  qhR2->Draw("sames");
  c_horR->SaveAs("HorRHRSfit.png");


 TCanvas *c_hav = new TCanvas();
 HorAv->Draw("AP");
 HorAv->Fit("pol0","W","",1,4);
 horavwien[2]->Fit("pol0","W","",5,9);
 TF1 *fhav1 = HorAv->GetFunction("pol0");
 fhav1->SetLineColor(1); fhav1->SetLineWidth(1);
 TPaveText *qhav1 = new TPaveText(1,-5e-4,4,-4.5e-4);
 qhav1->SetTextColor(kRed);
 qhav1->AddText(Form("p0 = %0.5e#pm%0.1e",fhav1->GetParameter(0),fhav1->GetParError(0)));
 TF1 *fhav2 = horavwien[2]->GetFunction("pol0");
 fhav2->SetLineColor(1); fhav2->SetLineWidth(1);
 TPaveText *qhav2 = new TPaveText(5,-5e-4,9,-4.5e-4);
 qhav2->SetTextColor(kBlue);
 qhav2->AddText(Form("p0 = %0.5e#pm%0.1e",fhav2->GetParameter(0),fhav2->GetParError(0)));
 lav->Draw();
 thh->Draw();
 qhav1->Draw("sames");
 qhav2->Draw("sames");
 c_hav->SaveAs("HorAvgfit.png");


  TCanvas *c_hdd = new TCanvas();
  HorDD->Draw("AP");
  HorDD->Fit("pol0","W","",1,4);
  hordwien[2]->Fit("pol0","W","",5,9);
  TF1 *fhdd1 = HorDD->GetFunction("pol0");
  fhdd1->SetLineColor(1); fhdd1->SetLineWidth(1);
  TPaveText *qhdd1 = new TPaveText(1,-5e-3,4,-4.5e-3);
  qhdd1->SetTextColor(kRed);
  qhdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fhdd1->GetParameter(0),fhdd1->GetParError(0)));
  TF1 *fhdd2 = hordwien[2]->GetFunction("pol0");
  fhdd2->SetLineColor(1); fhdd2->SetLineWidth(1);
  TPaveText *qhdd2 = new TPaveText(5,-5e-3,9,-4.5e-3);
  qhdd2->SetTextColor(kBlue);
  qhdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fhdd2->GetParameter(0),fhdd2->GetParError(0)));
  lav->Draw();
  tv->Draw();
  thL->Draw();
  qhdd1->Draw("sames");
  qhdd2->Draw("sames");
 c_hdd->SaveAs("HorDDfit.png");


 TCanvas *c_cosL = new TCanvas();
 Cos_L->Draw("AP");
 
 Cos_L->Fit("pol0","W","",1,4);
  cosLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fcL1 = Cos_L->GetFunction("pol0");
  fcL1->SetLineColor(1); fcL1->SetLineWidth(1);
  TPaveText *qcL1 = new TPaveText(1,-0.06,4,-0.055);
  qcL1->SetTextColor(kRed);
  qcL1->AddText(Form("p0 = %0.5e#pm%0.1e",fcL1->GetParameter(0),fcL1->GetParError(0)));
  TF1 *fcL2 = cosLwien[2]->GetFunction("pol0");
  fcL2->SetLineColor(1); fcL2->SetLineWidth(1);
  TPaveText *qcL2 = new TPaveText(5,-0.06,9,-0.055);
  qcL2->SetTextColor(kBlue);
  qcL2->AddText(Form("p0 = %0.5e#pm%0.1e",fcL2->GetParameter(0),fcL2->GetParError(0)));
  lav->Draw();
  qcL1->Draw("sames");
  qcL2->Draw("sames");
  c_cosL->SaveAs("cosmainLfit.png");

  c_cosL->SaveAs("LHRScosmain.png"); 

 TCanvas *c_cosR = new TCanvas();
 Cos_R->Draw("AP");

  Cos_R->Fit("pol0","W","",1,4);
  cosRwien[2]->Fit("pol0","W","",5,9);
  TF1 *fcR1 = Cos_R->GetFunction("pol0");
  fcR1->SetLineColor(1); fcR1->SetLineWidth(1);
  TPaveText *qcR1 = new TPaveText(1,0.0,4,0.005);
  qcR1->SetTextColor(kRed);
  qcR1->AddText(Form("p0 = %0.5e#pm%0.1e",fcR1->GetParameter(0),fcR1->GetParError(0)));
  TF1 *fcR2 = cosRwien[2]->GetFunction("pol0");
  fcR2->SetLineColor(1); fcR2->SetLineWidth(1);
  TPaveText *qcR2 = new TPaveText(5,0.,9,0.005);
  qcR2->SetTextColor(kBlue);
  qcR2->AddText(Form("p0 = %0.5e#pm%0.1e",fcR2->GetParameter(0),fcR2->GetParError(0)));
  lav->Draw();
  qcR1->Draw("sames");
  qcR2->Draw("sames");
  c_cosR->SaveAs("cosmainRfit.png");

  c_cosR->SaveAs("RHRScosmain.png");
*/

  TCanvas *c_sinL = new TCanvas();
  Sin_L->Draw("AP");
/*
  Sin_L->Fit("pol0","W","",1,4);
  sinLwien[2]->Fit("pol0","W","",5,9);
  TF1 *fsL1 = Sin_L->GetFunction("pol0");
  fsL1->SetLineColor(1); fsL1->SetLineWidth(1);
  TPaveText *qssL1 = new TPaveText(1,0.9,4,0.91);
  qssL1->SetTextColor(kRed);
  qssL1->AddText(Form("p0 = %0.5e#pm%0.1e",fsL1->GetParameter(0),fsL1->GetParError(0)));
  TF1 *fsL2 = sinLwien[2]->GetFunction("pol0");
  fsL2->SetLineColor(1); fsL2->SetLineWidth(1);
  TPaveText *qssL2 = new TPaveText(5,0.9,9,0.91);
  qssL2->SetTextColor(kBlue);
  qssL2->AddText(Form("p0 = %0.5e#pm%0.1e",fsL2->GetParameter(0),fsL2->GetParError(0)));
  lav->Draw();
  qssL1->Draw("sames");
  qssL2->Draw("sames");
  c_sinL->SaveAs("sinmainLfit.png");
*/
  c_sinL->SaveAs("sinmainL_nofit.png"); 

  TCanvas *c_sinR = new TCanvas();
  Sin_R->Draw("AP");
/*  Sin_R->Fit("pol0","W","",1,4);
  sinRwien[2]->Fit("pol0","W","",5,9);
  TF1 *fsR1 = Sin_R->GetFunction("pol0");
  fsR1->SetLineColor(1); fsR1->SetLineWidth(1);
  TPaveText *qssR1 = new TPaveText(1,0.9,4,0.91);
  qssR1->SetTextColor(kRed);
  qssR1->AddText(Form("p0 = %0.5e#pm%0.1e",fsR1->GetParameter(0),fsR1->GetParError(0)));
  TF1 *fsR2 = sinRwien[2]->GetFunction("pol0");
  fsR2->SetLineColor(1); fsR2->SetLineWidth(1);
  TPaveText *qssR2 = new TPaveText(5,0.9,9,0.91);
  qssR2->SetTextColor(kBlue);
  qssR2->AddText(Form("p0 = %0.5e#pm%0.1e",fsR2->GetParameter(0),fsR2->GetParError(0)));
  lav->Draw();
  qssR1->Draw("sames");
  qssR2->Draw("sames");
  c_sinR->SaveAs("sinmainRfit.png");
*/
  c_sinR->SaveAs("sinmainR_nofit.png");
/*
 TCanvas *c_cosav = new TCanvas();
 CosAv->Draw("AP");
 CosAv->Fit("pol0","W","",1,4);
  cosavwien[2]->Fit("pol0","W","",5,9);
  TF1 *fcavg1 = CosAv->GetFunction("pol0");
  fcavg1->SetLineColor(1); fcavg1->SetLineWidth(1);
  TPaveText *qcavg1 = new TPaveText(1,-0.006,4,-0.005);
  qcavg1->SetTextColor(kRed);
  qcavg1->AddText(Form("p0 = %0.5e#pm%0.1e",fcavg1->GetParameter(0),fcavg1->GetParError(0)));
  TF1 *fcavg2 = cosavwien[2]->GetFunction("pol0");
  fcavg2->SetLineColor(1); fcavg2->SetLineWidth(1);
  TPaveText *qcavg2 = new TPaveText(5,-0.006,9,-0.005);
  qcavg2->SetTextColor(kBlue);
  qcavg2->AddText(Form("p0 = %0.5e#pm%0.1e",fcavg2->GetParameter(0),fcavg2->GetParError(0)));
  lav->Draw();
  qcavg1->Draw("sames");
  qcavg2->Draw("sames");
 c_cosav->SaveAs("cosavmainfit.png");
*/

/*
 TCanvas *c_cosd = new TCanvas();
 CosDD->Draw("AP");
 CosDD->Fit("pol0","W","",1,4);
  cosdwien[2]->Fit("pol0","W","",5,9);
  TF1 *fcdd1 = CosDD->GetFunction("pol0");
  fcdd1->SetLineColor(1); fcdd1->SetLineWidth(1);
  TPaveText *qcdd1 = new TPaveText(1,-0.08,4,-0.07);
  qcdd1->SetTextColor(kRed);
  qcdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fcdd1->GetParameter(0),fcdd1->GetParError(0)));
  TF1 *fcdd2 = cosdwien[2]->GetFunction("pol0");
  fcdd2->SetLineColor(1); fcdd2->SetLineWidth(1);
  TPaveText *qcdd2 = new TPaveText(5,-0.08,9,-0.07);
  qcdd2->SetTextColor(kBlue);
  qcdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fcdd2->GetParameter(0),fcdd2->GetParError(0)));
  lav->Draw();
  qcdd1->Draw("sames");
  qcdd2->Draw("sames");
 c_cosd->SaveAs("cosddmainfit.png");
*/

  TCanvas *c_sinav = new TCanvas();
  SinAv->Draw("AP");
/*
  SinAv->Fit("pol0","W","",1,4);
  sinavwien[2]->Fit("pol0","W","",5,9);
  TF1 *fsavg1 = SinAv->GetFunction("pol0");
  fsavg1->SetLineColor(1); fsavg1->SetLineWidth(1);
  TPaveText *qssavg1 = new TPaveText(1,0.9,4,0.91);
  qssavg1->SetTextColor(kRed);
  qssavg1->AddText(Form("p0 = %0.5e#pm%0.1e",fsavg1->GetParameter(0),fsavg1->GetParError(0)));
  TF1 *fsavg2 = sinavwien[2]->GetFunction("pol0");
  fsavg2->SetLineColor(1); fsavg2->SetLineWidth(1);
  TPaveText *qssavg2 = new TPaveText(5,0.9,9,0.91);
  qssavg2->SetTextColor(kBlue);
  qssavg2->AddText(Form("p0 = %0.5e#pm%0.1e",fsavg2->GetParameter(0),fsavg2->GetParError(0)));
  lav->Draw();
  qssavg1->Draw("sames");
  qssavg2->Draw("sames");
*/
// c_sinav->SaveAs("sinavmainfit.png");
  c_sinav->SaveAs("sinavmain_nofit.png");


 TCanvas *c_sind = new TCanvas();
 SinDD->Draw("AP");
/*
  SinDD->Fit("pol0","W","",1,4);
  sindwien[2]->Fit("pol0","W","",5,9);
  TF1 *fsdd1 = SinDD->GetFunction("pol0");
  fsdd1->SetLineColor(1); fsdd1->SetLineWidth(1);
  TPaveText *qssdd1 = new TPaveText(1,-0.006,4,-0.005);
  qssdd1->SetTextColor(kRed);
  qssdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fsdd1->GetParameter(0),fsdd1->GetParError(0)));
  TF1 *fsdd2 = sindwien[2]->GetFunction("pol0");
  fsdd2->SetLineColor(1); fsdd2->SetLineWidth(1);
  TPaveText *qssdd2 = new TPaveText(5,-0.006,9,-0.005);
  qssdd2->SetTextColor(kBlue);
  qssdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fsdd2->GetParameter(0),fsdd2->GetParError(0)));
  lav->Draw();
  qssdd1->Draw("sames");
  qssdd2->Draw("sames");
*/
// c_sind->SaveAs("sinddmainfit.png");    
 c_sind->SaveAs("sinddmain_nofit.png");

/*
 TCanvas *c_xixi = new TCanvas();
 Xifac->Draw("AP");
 c_xixi->SaveAs("Newxi_calc.png");
*/



}

void ATCorr::PlotATDetPlots(arm Arm){
 
 int nn, mm;// This allows calling the same function for both arms despite the difference in the data set
 double var1, var2;
 string hrs;
 string atin = "ATin";
 string atout = "ATout";

 double qin1, qin2, qout1, qout2;
 double apvin1, apvin2, apvout1, apvout2;
 double labin1, labin2, labout1, labout2;
 double apvd1, apvd2,qd1,qd2, qa1, qa2;

 int color[2] = {2,4};

 if(Arm == kLHRS){
 ReadFiles("ThtgLHRSmain_newDB.txt");
 nn = 3; mm = 6;
 var1 = -0.015; var2 = 0;
 qin1 = 0.0285; qin2 = 0.03; qout1 = 0.029; qout2 = 0.0305;
 apvin1 = 2.4; apvin2 = 2.44; apvout1 = 2.4; apvout2 = 2.44;
 labin1 = 4.45; labin2 = 4.49; labout1 = 4.50; labout2 = 4.55;
 apvd1 = -0.01; apvd2 = 0.005; qd1 = -0.002;qd2 = 0.002; 
 qa1 = qout1; qa2 = qin2;
 hrs = "LHRS";
 } if(Arm == kRHRS){
 ReadFiles("ThtgRHRSmain_newDB.txt");
 nn = 1; mm = 5;
 var1 = 0; var2 = 0.015;
 qin1 = 0.03; qin2 = 0.0315; qout1 = 0.295; qout2 = 0.031;
 apvin1 = 2.4; apvin2 = 2.44; apvout1 = 2.4; apvout2 = 2.44;
 labin1 = 4.45; labin2 = 4.49; labout1 = 4.50; labout2 = 4.55;
 apvd1 = -0.01; apvd2 = 0.005; qd1 = -0.002;qd2 = 0.002;
 qa1 = qin1; qa2 = qout2;
 hrs = "RHRS";
 }


 ComputeChi();
 
 double wien2[nn], wien2err[nn], wien3[mm], wien3err[mm];
 
 //A_T Detector combos
  double qsqavg2[nn], qsqavg2err[nn], qsqavg3[mm], qsqavg3err[mm];
  double qsqdd2[nn], qsqdd2err[nn], qsqdd3[mm], qsqdd3err[mm];
  double thtgdd2[nn], thtgdd2err[nn], thtgdd3[mm], thtgdd3err[mm];
  double thtgavg2[nn], thtgavg2err[nn], thtgavg3[mm], thtgavg3err[mm];   
  double phtgdd2[nn], phtgdd2err[nn], phtgdd3[mm], phtgdd3err[mm];
  double phtgavg2[nn], phtgavg2err[nn], phtgavg3[mm], phtgavg3err[mm];
  double asymavg2[nn], asymavg2err[nn], asymavg3[mm], asymavg3err[mm];     
  double asymdd2[nn], asymdd2err[nn], asymdd3[mm], asymdd3err[mm];
  double horavg2[nn], horavg2err[nn], horavg3[mm], horavg3err[mm];
  double hordd2[nn], hordd2err[nn], hordd3[mm], hordd3err[mm];
  double veravg2[nn], veravg2err[nn], veravg3[mm], veravg3err[mm];
  double verdd2[nn], verdd2err[nn], verdd3[mm], verdd3err[mm];
  double cosavg2[nn], cosavg2err[nn], cosavg3[mm], cosavg3err[mm];
  double cosdd2[nn], cosdd2err[nn], cosdd3[mm], cosdd3err[mm];
  double sinavg2[nn], sinavg2err[nn], sinavg3[mm], sinavg3err[mm];
  double sindd2[nn], sindd2err[nn], sindd3[mm], sindd3err[mm];

 //Chi
 double chi2[nn], chi2err[nn], chi3[mm], chi3err[mm];

 //Individual A_Ts
  double qsL2[nn], qsL2er[nn], qsL3[mm], qsL3er[mm];
  double qsR2[nn], qsR2er[nn], qsR3[mm], qsR3er[mm];
  double thsL2[nn], thsL2er[nn], thsL3[mm], thsL3er[mm];
  double thsR2[nn], thsR2er[nn], thsR3[mm], thsR3er[mm];
  double phsL2[nn], phsL2er[nn], phsL3[mm], phsL3er[mm];
  double phsR2[nn], phsR2er[nn], phsR3[mm], phsR3er[mm];
  double asyL2[nn], asyL2er[nn], asyL3[mm], asyL3er[mm];
  double asyR2[nn], asyR2er[nn], asyR3[mm], asyR3er[mm];
  double labL2[nn], labL2er[nn], labL3[mm], labL3er[mm];
  double labR2[nn], labR2er[nn], labR3[mm], labR3er[mm];
  double verL2[nn], verL3[mm], horL2[nn], horL3[mm];
  double verR2[nn], verR3[mm], horR2[nn], horR3[mm];
  double verL2er[nn], verL3er[mm], horL2er[nn], horL3er[mm];
  double verR2er[nn], verR3er[mm], horR2er[nn], horR3er[mm];
  double cosL2[nn], cosL3[mm], sinL2[nn], sinL3[mm];
  double cosR2[nn], cosR3[mm],  sinR2[nn], sinR3[mm];
  double cosL2er[nn], cosL3er[mm], sinL2er[nn], sinL3er[mm];
  double cosR2er[nn], cosR3er[mm], sinR2er[nn], sinR3er[mm];
 

  int ii = 0; 
  //Separating at the Wien level
  for(int i = 0; i < nn+mm; i++){
   if(i < nn) { 
                qsqavg2[i] = QsqAvg[i]; qsqavg2err[i] = QsqAvgerr[i]; qsqdd2[i] = QsqDiff[i]; qsqdd2err[i] = QsqDifferr[i];
                thtgavg2[i] = ThtgAvg[i]; thtgavg2err[i] = ThtgAvgerr[i]; thtgdd2[i] = ThtgDiff[i]; thtgdd2err[i] = ThtgDifferr[i];
                phtgavg2[i] = PhtgAvg[i]; phtgavg2err[i] = PhtgAvgerr[i]; phtgdd2[i] = PhtgDiff[i]; phtgdd2err[i] = PhtgDifferr[i];
                asymavg2[i] = AsymAvg[i]; asymavg2err[i] = AsymAvgerr[i]; asymdd2[i] = AsymDiff[i]; asymdd2err[i] = AsymDifferr[i];
                qsL2[i] = QsqL[i]; qsL2er[i] = QsqerL[i]; qsR2[i] = QsqR[i]; qsR2er[i] = QsqerR[i];
                thsL2[i] = ThtgL[i]; thsL2er[i] = ThtgerL[i]; thsR2[i] = ThtgR[i]; thsR2er[i] = ThtgerR[i];
                phsL2[i] = PhtgL[i]; phsL2er[i] = PhtgerL[i]; phsR2[i] = PhtgR[i]; phsR2er[i] = PhtgerR[i];
                asyL2[i] = asymL[i]; asyL2er[i] = asymerL[i]; asyR2[i] = asymR[i]; asyR2er[i] = asymerR[i];
                labL2[i] = AngleL[i]; labL2er[i] = AngleerL[i]; labR2[i] = AngleR[i]; labR2er[i] = AngleerR[i];
                verL2[i] = VermomL[i]; verL2er[i] = VermomerL[i]; verR2[i] = VermomR[i]; verR2er[i] = VermomerR[i];
                horL2[i] = HormomL[i]; horL2er[i] = HormomerL[i]; horR2[i] = HormomR[i]; horR2er[i] = HormomerR[i];
                cosL2[i] = cosphL[i]; cosL2er[i] = cospherL[i]; cosR2[i] = cosphR[i]; cosR2er[i] = cospherR[i];
                sinL2[i] = sinphL[i]; sinL2er[i] = sinpherL[i]; sinR2[i] = sinphR[i]; sinR2er[i] = sinpherR[i];
                cosavg2[i] = CosAvg[i]; cosavg2err[i] = CosAvgerr[i]; cosdd2[i] = CosDiff[i]; cosdd2err[i] = CosDifferr[i];
                sinavg2[i] = SinAvg[i]; sinavg2err[i] = SinAvgerr[i]; sindd2[i] = SinDiff[i]; sindd2err[i] = SinDifferr[i];
                veravg2[i] = VerMomAvg[i]; veravg2err[i] = VerMomAvgerr[i]; verdd2[i] = VerMomDiff[i]; verdd2err[i] = VerMomDifferr[i];
                horavg2[i] = HorMomAvg[i]; horavg2err[i] = HorMomAvgerr[i]; hordd2[i] = HorMomDiff[i]; hordd2err[i] = HorMomDifferr[i];
                chi2[i] = Chi[i]; chi2err[i] = Chierr[i];                                                                                                        
                wien2[i] = i+1; wien2err[i] = 0;

    } else { 
                qsqavg3[ii] = QsqAvg[i]; qsqavg3err[ii] = QsqAvgerr[i]; qsqdd3[ii] = QsqDiff[i]; qsqdd3err[ii] = QsqDifferr[i];
                thtgavg3[ii] = ThtgAvg[i]; thtgavg3err[ii] = ThtgAvgerr[i]; thtgdd3[ii] = ThtgDiff[i]; thtgdd3err[ii] = ThtgDifferr[i];
                phtgavg3[ii] = PhtgAvg[i]; phtgavg3err[ii] = PhtgAvgerr[i]; phtgdd3[ii] = PhtgDiff[i]; phtgdd3err[ii] = PhtgDifferr[i];
                asymavg3[ii] = AsymAvg[i]; asymavg3err[ii] = AsymAvgerr[i]; asymdd3[ii] = AsymDiff[i]; asymdd3err[ii] = AsymDifferr[i];
                qsL3[ii] = QsqL[i]; qsL3er[ii] = QsqerL[i]; qsR3[ii] = QsqR[i]; qsR3er[ii] = QsqerR[i];
                thsL3[ii] = ThtgL[i]; thsL3er[ii] = ThtgerL[i]; thsR3[ii] = ThtgR[i]; thsR3er[ii] = ThtgerR[i];
                phsL3[ii] = PhtgL[i]; phsL3er[ii] = PhtgerL[i]; phsR3[ii] = PhtgR[i]; phsR3er[ii] = PhtgerR[i];
                asyL3[ii] = asymL[i]; asyL3er[ii] = asymerL[i]; asyR3[ii] = asymR[i]; asyR3er[ii] = asymerR[i];
                labL3[ii] = AngleL[i]; labL3er[ii] = AngleerL[i]; labR3[ii] = AngleR[i]; labR3er[ii] = AngleerR[i];
                verL3[ii] = VermomL[i]; verL3er[ii] = VermomerL[i]; verR3[ii] = VermomR[i]; verR3er[ii] = VermomerR[i];
                horL3[ii] = HormomL[i]; horL3er[ii] = HormomerL[i]; horR3[ii] = HormomR[i]; horR3er[ii] = HormomerR[i];
                cosL3[ii] = cosphL[i]; cosL3er[ii] = cospherL[i]; cosR3[ii] = cosphR[i]; cosR3er[ii] = cospherR[i];
                sinL3[ii] = sinphL[i]; sinL3er[ii] = sinpherL[i]; sinR3[ii] = sinphR[i]; sinR3er[ii] = sinpherR[i];
                cosavg3[ii] = CosAvg[i]; cosavg3err[ii] = CosAvgerr[i]; cosdd3[ii] = CosDiff[i]; cosdd3err[ii] = CosDifferr[i];
                sinavg3[ii] = SinAvg[i]; sinavg3err[ii] = SinAvgerr[i]; sindd3[ii] = SinDiff[i]; sindd3err[ii] = SinDifferr[i];
                veravg3[ii] = VerMomAvg[i]; veravg3err[ii] = VerMomAvgerr[i]; verdd3[ii] = VerMomDiff[i]; verdd3err[ii] = VerMomDifferr[i];
                horavg3[ii] = HorMomAvg[i]; horavg3err[ii] = HorMomAvgerr[i]; hordd3[ii] = HorMomDiff[i]; hordd3err[ii] = HorMomDifferr[i];
                chi3[ii] = Chi[i]; chi3err[ii] = Chierr[i];
                wien3[ii] = 1+i; wien3err[ii] = 0;
           //    cout << chi3[ii] << "  " << chi3err[ii] << endl;;
               ii++;
   }

  }


   Chifac = new TMultiGraph(); 
   Chifac->SetTitle(Form("#chi_{%s} (Up/Down Acceptance Factor); Run Number; #chi_{%s}",hrs.c_str(),hrs.c_str()));
   Chifac->SetMinimum(-10); Chifac->SetMaximum(0);   

   Q2_L = new TMultiGraph();
   Q2_L->SetTitle(Form("%s %s Q^{2} (GeV/c)^{2}; Run Number; Q^{2} (GeV/c)^{2}",hrs.c_str(),atin.c_str()));
   Q2_L->SetMinimum(qin1); Q2_L->SetMaximum(qin2);

   Q2_R = new TMultiGraph();
   Q2_R->SetTitle(Form("%s %s Q^{2} (GeV/c)^{2}; Run Number; Q^{2} (GeV/c)^{2}",hrs.c_str(),atout.c_str()));
   Q2_R->SetMinimum(qout1); Q2_R->SetMaximum(qout2);

   Q2Av = new TMultiGraph();
   Q2Av->SetTitle(Form("%s A_{T} Average Q^{2} (GeV/c)^{2}; Run Number; Q^{2} (GeV/c)^{2}",hrs.c_str()));
   Q2Av->SetMinimum(qa1); Q2Av->SetMaximum(qa2);

   Q2DD = new TMultiGraph();
   Q2DD->SetTitle(Form("%s A_{T} Q^{2} Double Difference (GeV/c)^{2}; Run Number; Q^{2} (GeV/c)^{2}",hrs.c_str()));
   Q2DD->SetMinimum(qd1); Q2DD->SetMaximum(qd2);

   Ph_L = new TMultiGraph();
   Ph_L->SetTitle(Form("%s %s  <#phi_{tg}>; Run Number; <#phi_{tg}>",hrs.c_str(),atin.c_str()));
   Ph_L->SetMinimum(var1); Ph_L->SetMaximum(var2);

   Ph_R = new TMultiGraph();
   Ph_R->SetTitle(Form("%s %s <#phi_{tg}>; Run Number; <#phi_{tg}>",hrs.c_str(),atout.c_str()));
   Ph_R->SetMinimum(var1); Ph_R->SetMaximum(var2);

   PhAv = new TMultiGraph();
   PhAv->SetTitle(Form("%s A_{T} Average <#phi_{tg}>; Run Number; <#phi_{tg}>",hrs.c_str()));
   PhAv->SetMinimum(var1); PhAv->SetMaximum(var2);

   PhDD = new TMultiGraph();
   PhDD->SetTitle(Form("%s <#phi_{tg}> Double Difference; Run Number; <#phi_{tg}>",hrs.c_str()));
   PhDD->SetMinimum(-0.01); PhDD->SetMaximum(0.01);

   Th_L = new TMultiGraph();
   Th_L->SetTitle(Form("%s %s  <#theta_{tg}>; Run Number; <#theta_{tg}>",hrs.c_str(),atin.c_str()));
   Th_L->SetMinimum(0.0); Th_L->SetMaximum(0.03);

   Th_R = new TMultiGraph();
   Th_R->SetTitle(Form("%s %s  <#theta_{tg}>; Run Number; <#theta_{tg}>",hrs.c_str(),atout.c_str()));
   Th_R->SetMinimum(-0.03); Th_R->SetMaximum(0.0);

   ThAv = new TMultiGraph();
   ThAv->SetTitle(Form("%s A_{T} Average <#theta_{tg}>; Run Number; <#theta_{tg}>",hrs.c_str()));
   ThAv->SetMinimum(-0.01); ThAv->SetMaximum(0.01);

   ThDD = new TMultiGraph();
   ThDD->SetTitle(Form("%s A_{T} <#theta_{tg}> Double Difference; Run Number; <#theta_{tg}>",hrs.c_str()));
   ThDD->SetMinimum(0); ThDD->SetMaximum(0.03);


   APV_L = new TMultiGraph();
   APV_L->SetTitle(Form("%s %s <A_{PV}> (ppm); Run Number; <A_{PV}> (ppm)",hrs.c_str(),atin.c_str()));
   APV_L->SetMinimum(apvin1); APV_L->SetMaximum(apvin2);

   APV_R = new TMultiGraph();
   APV_R->SetTitle(Form("%s %s <A_{PV}> (ppm); Run Number; <A_{PV}> (ppm)",hrs.c_str(),atout.c_str()));
   APV_R->SetMinimum(apvout1); APV_R->SetMaximum(apvout2);


   APVAvg = new TMultiGraph();
   APVAvg->SetTitle(Form("%s A_{T} Average <A_{PV}> (ppm); Run Number; <A_{PV}> (ppm)",hrs.c_str()));
   APVAvg->SetMinimum(apvin1); APVAvg->SetMaximum(apvin2);

   APVDD = new TMultiGraph();
   APVDD->SetTitle(Form("%s A_{T} <A_{PV}> Double Difference (ppm); Run Number; <A_{PV}> (ppm)",hrs.c_str()));
   APVDD->SetMinimum(apvd1); APVDD->SetMaximum(apvd2);

   Lab_L = new TMultiGraph();
   Lab_L->SetTitle(Form("%s %s #theta_{lab} (deg); Run Number; #theta_{lab} (deg)",hrs.c_str(),atin.c_str()));
   Lab_L->SetMinimum(labin1); Lab_L->SetMaximum(labin2);

   Lab_R = new TMultiGraph();
   Lab_R->SetTitle(Form("%s %s #theta_{lab} (deg); Run Number; #theta_{lab} (deg)",hrs.c_str(),atout.c_str()));
   Lab_R->SetMinimum(labout1); Lab_R->SetMaximum(labout2);

   Ver_L = new TMultiGraph();
   Ver_L->SetTitle(Form("%s %s sin(#theta)sin(#phi); Run Number; sin(#theta)sin(#phi)",hrs.c_str(),atin.c_str()));
//   Ver_L->SetMinimum(0.05); Ver_L->SetMaximum(0.1);

   Ver_R = new TMultiGraph();
   Ver_R->SetTitle(Form("%s %s sin(#theta)sin(#phi); Run Number; sin(#theta)sin(#phi)",hrs.c_str(),atout.c_str()));
 
//  Ver_R->SetMinimum(0.05); Ver_R->SetMaximum(0.1);

   VerAv = new TMultiGraph();
   VerAv->SetTitle(Form("%s A_{T} Average sin(#theta)sin(#phi); Run Number; sin(#theta)sin(#phi)",hrs.c_str()));
//   VerAv->SetMinimum(0.05); VerAv->SetMaximum(0.1);

   VerDD = new TMultiGraph();
   VerDD->SetTitle(Form("%s A_{T} sin(#theta)sin(#phi) Double Difference; Run Number; sin(#theta)sin(#phi)",hrs.c_str()));
//   VerDD->SetMinimum(-7e-4); VerDD->SetMaximum(0);

   Hor_L = new TMultiGraph();
   Hor_L->SetTitle(Form("%s %s sin(#theta)cos(#phi); Run Number; sin(#theta)cos(#phi)",hrs.c_str(),atin.c_str()));
//   Hor_L->SetMinimum(-0.005); Hor_L->SetMaximum(0);

   Hor_R = new TMultiGraph();
   Hor_R->SetTitle(Form("%s %s sin(#theta)cos(#phi); Run Number; sin(#theta)cos(#phi)",hrs.c_str(),atout.c_str()));
//   Hor_R->SetMinimum(0); Hor_R->SetMinimum(0.005);

   HorAv = new TMultiGraph();
   HorAv->SetTitle(Form("%s A_{T} Average sin(#theta)cos(#phi); Run Number; sin(#theta)cos(#phi)",hrs.c_str()));
//   HorAv->SetMinimum(-5e-4); HorAv->SetMaximum(5e-4);

   HorDD = new TMultiGraph();
   HorDD->SetTitle(Form("%s A_{T} sin(#theta)cos(#phi) Double Difference; Run Number; sin(#theta)cos(#phi)", hrs.c_str()));
//   HorDD->SetMinimum(-0.005); HorDD->SetMaximum(0);

   Cos_L = new TMultiGraph();
   Cos_L->SetTitle(Form("%s %s cos(#phi); Run Number; cos(#phi)",hrs.c_str(),atin.c_str()));
   Cos_L->SetMinimum(0); Cos_L->SetMaximum(0.4);

   Cos_R = new TMultiGraph();
   Cos_R->SetTitle(Form("%s %s cos(#phi); Run Number; cos(#phi)",hrs.c_str(),atout.c_str()));
   Cos_R->SetMinimum(-0.4); Cos_R->SetMaximum(0);

   CosAv = new TMultiGraph();
   CosAv->SetTitle(Form("%s A_{T} Average cos(#phi); Run Number; cos(#phi)",hrs.c_str()));
   CosAv->SetMinimum(-0.006); CosAv->SetMaximum(0.006);

   CosDD = new TMultiGraph();
   CosDD->SetTitle(Form("%s A_{T} cos(#phi) Double Difference; Run Number; cos(#phi)",hrs.c_str()));
//   CosDD->SetMinimum(0); CosDD->SetMaximum(0.4);


   Sin_L = new TMultiGraph();
   Sin_L->SetTitle(Form("%s %s sin(#phi); Run Number; sin(#phi)",hrs.c_str(),atin.c_str()));
   Sin_L->SetMinimum(0.9);Sin_L->SetMaximum(1.05);

   Sin_R = new TMultiGraph();
   Sin_R->SetTitle(Form("%s %s sin(#phi); Run Number; sin(#phi)",hrs.c_str(),atout.c_str()));
   Sin_R->SetMinimum(0.9); Sin_R->SetMaximum(1.05);

   SinAv = new TMultiGraph();
   SinAv->SetTitle(Form("%s A_{T} Average sin(#phi); Run Number; sin(#phi)",hrs.c_str()));
   SinAv->SetMinimum(0.9); SinAv->SetMaximum(1.05);

   SinDD = new TMultiGraph();
   SinDD->SetTitle(Form("%s A_{T} sin(#phi) Double Difference; Run Number; sin(#phi)", hrs.c_str()));
   SinDD->SetMinimum(-0.006); SinDD->SetMaximum(0.006);



   //Combos
   TGraphErrors *chiwien[2], *qavwien[2], *qdwien[2], *thavwien[2], *thdwien[2], *phavwien[2], *phdwien[2], *asavwien[2], *asdwien[2];   
   TGraphErrors *veravwien[2], *horavwien[2], *cosavwien[2], *sinavwien[2], *verdwien[2], *hordwien[2], *cosdwien[2], *sindwien[2];

  //Individual
   TGraphErrors *qLwien[2], *qRwien[2], *thLwien[2], *thRwien[2], *phLwien[2], *phRwien[2], *asLwien[2], *asRwien[2], *lbLwien[2], *lbRwien[2];
   TGraphErrors *cosLwien[2], *cosRwien[2], *sinLwien[2], *sinRwien[2], *verLwien[2], *verRwien[2], *horLwien[2], *horRwien[2];

   chiwien[0] = new TGraphErrors(nn,wien2,chi2,wien2err,chi2err); chiwien[1] = new TGraphErrors(mm,wien3,chi3,wien3err,chi3err); 

   qavwien[0] = new TGraphErrors(nn,wien2,qsqavg2,wien2err,qsqavg2err); qdwien[0] = new TGraphErrors(nn,wien2,qsqdd2,wien2err,qsqdd2err);
   qavwien[1] = new TGraphErrors(mm,wien3,qsqavg3,wien3err,qsqavg3err); qdwien[1] = new TGraphErrors(mm,wien3,qsqdd3,wien3err,qsqdd3err);

   thavwien[0] = new TGraphErrors(nn,wien2,thtgavg2,wien2err,thtgavg2err); thdwien[0] = new TGraphErrors(nn,wien2,thtgdd2,wien2err,thtgdd2err);
   thavwien[1] = new TGraphErrors(mm,wien3,thtgavg3,wien3err,thtgavg3err); thdwien[1] = new TGraphErrors(mm,wien3,thtgdd3,wien3err,thtgdd3err);

   phavwien[0] = new TGraphErrors(nn,wien2,phtgavg2,wien2err,phtgavg2err); phdwien[0] = new TGraphErrors(nn,wien2,phtgdd2,wien2err,phtgdd2err);
   phavwien[1] = new TGraphErrors(mm,wien3,phtgavg3,wien3err,phtgavg3err); phdwien[1] = new TGraphErrors(mm,wien3,phtgdd3,wien3err,phtgdd3err);

   asavwien[0] = new TGraphErrors(nn,wien2,asymavg2,wien2err,asymavg2err); asdwien[0] = new TGraphErrors(nn,wien2,asymdd2,wien2err,asymdd2err);
   asavwien[1] = new TGraphErrors(mm,wien3,asymavg3,wien3err,asymavg3err); asdwien[1] = new TGraphErrors(mm,wien3,asymdd3,wien3err,asymdd3err);

   veravwien[0] = new TGraphErrors(nn,wien2,veravg2,wien2err,veravg2err); verdwien[0] = new TGraphErrors(nn,wien2,verdd2,wien2err,verdd2err);
   veravwien[1] = new TGraphErrors(mm,wien3,veravg3,wien3err,veravg3err); verdwien[1] = new TGraphErrors(mm,wien3,verdd3,wien3err,verdd3err);

   horavwien[0] = new TGraphErrors(nn,wien2,horavg2,wien2err,horavg2err); hordwien[0] = new TGraphErrors(nn,wien2,hordd2,wien2err,hordd2err);
   horavwien[1] = new TGraphErrors(mm,wien3,horavg3,wien3err,horavg3err); hordwien[1] = new TGraphErrors(mm,wien3,hordd3,wien3err,hordd3err);
  
   cosavwien[0] = new TGraphErrors(nn,wien2,cosavg2,wien2err,cosavg2err); cosdwien[0] = new TGraphErrors(nn,wien2,cosdd2,wien2err,cosdd2err);
   cosavwien[1] = new TGraphErrors(mm,wien3,cosavg3,wien3err,cosavg3err); cosdwien[1] = new TGraphErrors(mm,wien3,cosdd3,wien3err,cosdd3err);

   sinavwien[0] = new TGraphErrors(nn,wien2,sinavg2,wien2err,sinavg2err); sindwien[0] = new TGraphErrors(nn,wien2,sindd2,wien2err,sindd2err);
   sinavwien[1] = new TGraphErrors(mm,wien3,sinavg3,wien3err,sinavg3err); sindwien[1] = new TGraphErrors(mm,wien3,sindd3,wien3err,sindd3err);

   qLwien[0] = new TGraphErrors(nn,wien2,qsL2,wien2err,qsL2er); qRwien[0] = new TGraphErrors(nn,wien2,qsR2,wien2err,qsR2er);
   qLwien[1] = new TGraphErrors(mm,wien3,qsL3,wien3err,qsL3er); qRwien[1] = new TGraphErrors(mm,wien3,qsR3,wien3err,qsR3er);

   thLwien[0] = new TGraphErrors(nn,wien2,thsL2,wien2err,thsL2er); thRwien[0] = new TGraphErrors(nn,wien2,thsR2,wien2err,thsR2er);
   thLwien[1] = new TGraphErrors(mm,wien3,thsL3,wien3err,thsL3er); thRwien[1] = new TGraphErrors(mm,wien3,thsR3,wien3err,thsR3er);

   phLwien[0] = new TGraphErrors(nn,wien2,phsL2,wien2err,phsL3er); phRwien[0] = new TGraphErrors(nn,wien2,phsR2,wien2err,phsR2er);
   phLwien[1] = new TGraphErrors(mm,wien3,phsL3,wien3err,phsL3er); phRwien[1] = new TGraphErrors(mm,wien3,phsR3,wien3err,phsR3er);

   asLwien[0] = new TGraphErrors(nn,wien2,asyL2,wien2err,asyL3er); asRwien[0] = new TGraphErrors(nn,wien2,asyR2,wien2err,asyR2er);
   asLwien[1] = new TGraphErrors(mm,wien3,asyL3,wien3err,asyL3er); asRwien[1] = new TGraphErrors(mm,wien3,asyR3,wien3err,asyR3er);

   lbLwien[0] = new TGraphErrors(nn,wien2,labL2,wien2err,labL3er); lbRwien[0] = new TGraphErrors(nn,wien2,labR2,wien2err,labR2er);
   lbLwien[1] = new TGraphErrors(mm,wien3,labL3,wien3err,labL3er); lbRwien[1] = new TGraphErrors(mm,wien3,labR3,wien3err,labR3er);

   cosLwien[0] = new TGraphErrors(nn,wien2,cosL2,wien2err,cosL3er); cosRwien[0] = new TGraphErrors(nn,wien2,cosR2,wien2err,cosR2er);
   cosLwien[1] = new TGraphErrors(mm,wien3,cosL3,wien3err,cosL3er); cosRwien[1] = new TGraphErrors(mm,wien3,cosR3,wien3err,cosR3er);

   sinLwien[0] = new TGraphErrors(nn,wien2,sinL2,wien2err,sinL3er); sinRwien[0] = new TGraphErrors(nn,wien2,sinR2,wien2err,sinR2er);
   sinLwien[1] = new TGraphErrors(mm,wien3,sinL3,wien3err,sinL3er); sinRwien[1] = new TGraphErrors(mm,wien3,sinR3,wien3err,sinR3er);

   verLwien[0] = new TGraphErrors(nn,wien2,verL2,wien2err,verL3er); verRwien[0] = new TGraphErrors(nn,wien2,verR2,wien2err,verR2er);
   verLwien[1] = new TGraphErrors(mm,wien3,verL3,wien3err,verL3er); verRwien[1] = new TGraphErrors(mm,wien3,verR3,wien3err,verR3er);

   horLwien[0] = new TGraphErrors(nn,wien2,horL2,wien2err,horL2er); horRwien[0] = new TGraphErrors(nn,wien2,horR2,wien2err,horR2er);
   horLwien[1] = new TGraphErrors(mm,wien3,horL3,wien3err,horL3er); horRwien[1] = new TGraphErrors(mm,wien3,horR3,wien3err,horR3er);

   TLegend *lav = new TLegend(0.55,0.7,0.9,0.9);
  

  for(int i = 0; i < 2; i++){
      chiwien[i]->SetMarkerStyle(21); chiwien[i]->SetMarkerColor(color[i]);
      qavwien[i]->SetMarkerStyle(21); qavwien[i]->SetMarkerColor(color[i]);
      qdwien[i]->SetMarkerStyle(21); qdwien[i]->SetMarkerColor(color[i]);
      thavwien[i]->SetMarkerStyle(21); thavwien[i]->SetMarkerColor(color[i]);
      thdwien[i]->SetMarkerStyle(21); thdwien[i]->SetMarkerColor(color[i]);
      phavwien[i]->SetMarkerStyle(21); phavwien[i]->SetMarkerColor(color[i]);
      phdwien[i]->SetMarkerStyle(21); phdwien[i]->SetMarkerColor(color[i]);
      asavwien[i]->SetMarkerStyle(21); asavwien[i]->SetMarkerColor(color[i]);
      veravwien[i]->SetMarkerStyle(21); veravwien[i]->SetMarkerColor(color[i]);
      verdwien[i]->SetMarkerStyle(21); verdwien[i]->SetMarkerColor(color[i]);
      horavwien[i]->SetMarkerStyle(21); horavwien[i]->SetMarkerColor(color[i]);
      hordwien[i]->SetMarkerStyle(21); hordwien[i]->SetMarkerColor(color[i]);
      cosavwien[i]->SetMarkerStyle(21); cosavwien[i]->SetMarkerColor(color[i]);
      cosdwien[i]->SetMarkerStyle(21); cosdwien[i]->SetMarkerColor(color[i]);
      sinavwien[i]->SetMarkerStyle(21); sinavwien[i]->SetMarkerColor(color[i]);
      sindwien[i]->SetMarkerStyle(21); sindwien[i]->SetMarkerColor(color[i]);
      asdwien[i]->SetMarkerStyle(21); asdwien[i]->SetMarkerColor(color[i]);
      qLwien[i]->SetMarkerStyle(21); qLwien[i]->SetMarkerColor(color[i]);
      qRwien[i]->SetMarkerStyle(21); qRwien[i]->SetMarkerColor(color[i]);
      thLwien[i]->SetMarkerStyle(21); thLwien[i]->SetMarkerColor(color[i]);
      thRwien[i]->SetMarkerStyle(21); thRwien[i]->SetMarkerColor(color[i]);
      phLwien[i]->SetMarkerStyle(21); phLwien[i]->SetMarkerColor(color[i]);
      phRwien[i]->SetMarkerStyle(21); phRwien[i]->SetMarkerColor(color[i]);
      asLwien[i]->SetMarkerStyle(21); asLwien[i]->SetMarkerColor(color[i]);
      asRwien[i]->SetMarkerStyle(21); asRwien[i]->SetMarkerColor(color[i]); 
      lbLwien[i]->SetMarkerStyle(21); lbLwien[i]->SetMarkerColor(color[i]);
      lbRwien[i]->SetMarkerStyle(21); lbRwien[i]->SetMarkerColor(color[i]);
      cosLwien[i]->SetMarkerStyle(21); cosLwien[i]->SetMarkerColor(color[i]);
      cosRwien[i]->SetMarkerStyle(21); cosRwien[i]->SetMarkerColor(color[i]);
      sinLwien[i]->SetMarkerStyle(21); sinLwien[i]->SetMarkerColor(color[i]);
      sinRwien[i]->SetMarkerStyle(21); sinRwien[i]->SetMarkerColor(color[i]);
      verLwien[i]->SetMarkerStyle(21); verLwien[i]->SetMarkerColor(color[i]);
      verRwien[i]->SetMarkerStyle(21); verRwien[i]->SetMarkerColor(color[i]);
      horLwien[i]->SetMarkerStyle(21); horLwien[i]->SetMarkerColor(color[i]);
      horRwien[i]->SetMarkerStyle(21); horRwien[i]->SetMarkerColor(color[i]);
      



      Chifac->Add(chiwien[i]);
      Q2_L->Add(qLwien[i]); Q2_R->Add(qRwien[i]); 
      Lab_L->Add(lbLwien[i]); Lab_R->Add(lbRwien[i]); 
      Th_L->Add(thLwien[i]); Th_R->Add(thRwien[i]);
      Ph_L->Add(phLwien[i]); Ph_R->Add(phRwien[i]);
      APV_L->Add(asLwien[i]); APV_R->Add(asRwien[i]);
      Cos_L->Add(cosLwien[i]); Cos_R->Add(cosRwien[i]);
      Sin_L->Add(sinLwien[i]); Sin_R->Add(sinRwien[i]);
      Ver_L->Add(verLwien[i]); Ver_R->Add(verRwien[i]);
      Hor_L->Add(horLwien[i]); Hor_R->Add(horRwien[i]);
      Q2Av->Add(qavwien[i]); Q2DD->Add(qdwien[i]);
      ThAv->Add(thavwien[i]); ThDD->Add(thdwien[i]);
      PhAv->Add(phavwien[i]); PhDD->Add(phdwien[i]);
      APVAvg->Add(asavwien[i]); APVDD->Add(asdwien[i]);
      VerAv->Add(veravwien[i]); VerDD->Add(verdwien[i]);
      HorAv->Add(horavwien[i]); HorDD->Add(hordwien[i]);
      CosAv->Add(cosavwien[i]); CosDD->Add(cosdwien[i]);
      SinAv->Add(sinavwien[i]); SinDD->Add(sindwien[i]);

    if(i == 0) { lav->AddEntry(chiwien[i],"Wien Left","p"); } else { lav->AddEntry(chiwien[i],"Wien Right","p"); }

     }

/*
  TCanvas *cqL = new TCanvas();
  Q2_L->Draw("AP");
//  Q2_L->Fit("pol0","W","",1,nn);
  qLwien[1]->Fit("pol0","W","",nn+1,nn+mm);
//  TF1 *fl1 = Q2_L->GetFunction("pol0");
//  fl1->SetLineColor(1); fl1->SetLineWidth(1);
//  TPaveText *ql1 = new TPaveText(1,qin1,nn,qin1+0.002);
//  ql1->SetTextColor(kRed);
//  ql1->AddText(Form("p0 = %0.5e#pm%0.1e",fl1->GetParameter(0),fl1->GetParError(0)));
  TF1 *fl2 = qLwien[1]->GetFunction("pol0");
  fl2->SetLineColor(1); fl2->SetLineWidth(1);
  TPaveText *ql2 = new TPaveText(nn+1,qin1,nn+mm,qin1+0.0002);
  ql2->SetTextColor(kBlue);
  ql2->AddText(Form("p0 = %0.5e#pm%0.1e",fl2->GetParameter(0),fl2->GetParError(0)));
  ql2->SetTextAlign(12);
//  ql1->Draw("sames");
  ql2->Draw("sames");
  lav->Draw();
  cqL->SaveAs(Form("Q2%s%sfit.png",atin.c_str(),hrs.c_str()));
 
 
  TCanvas *cqR = new TCanvas();
  Q2_R->Draw("AP");
//  Q2_R->Fit("pol0","W","",1,nn);
  qRwien[1]->Fit("pol0","W","",nn+1,nn+mm);
//  TF1 *fr1 = Q2_R->GetFunction("pol0");
//  fr1->SetLineColor(1); fr1->SetLineWidth(1);
//  TPaveText *qr1 = new TPaveText(1,qout1,4,qout1+0.002);
//  qr1->SetTextColor(kRed);
//  qr1->AddText(Form("p0 = %0.5e#pm%0.1e",fr1->GetParameter(0),fr1->GetParError(0)));
  TF1 *fr2 = qRwien[1]->GetFunction("pol0");
  fr2->SetLineColor(1); fr2->SetLineWidth(1);
  TPaveText *qr2 = new TPaveText(nn+1,qout1,nn+mm,qout1+0.0002);
  qr2->SetTextColor(kBlue);
  qr2->AddText(Form("p0 = %0.5e#pm%0.1e",fr2->GetParameter(0),fr2->GetParError(0)));
  qr2->SetTextAlign(12);
//  qr1->Draw("sames");
  qr2->Draw("sames");
  lav->Draw();
  cqR->SaveAs(Form("Q2%s%sfit.png",atout.c_str(),hrs.c_str()));


  TCanvas *cthL = new TCanvas();
  Th_L->Draw("AP");
//  Th_L->Fit("pol0","W","",1,nn);
  thLwien[1]->Fit("pol0","W","",nn+1,nn+mm);
//  TF1 *fthl1 = Th_L->GetFunction("pol0");
//  fthl1->SetLineColor(1); fthl1->SetLineWidth(1);
//  TPaveText *thl1 = new TPaveText(1,0,nn,0.002);
//  thl1->SetTextColor(kRed);
//  thl1->AddText(Form("p0 = %0.5e#pm%0.1e",fthl1->GetParameter(0),fthl1->GetParError(0)));
  TF1 *fthl2 = thLwien[1]->GetFunction("pol0");
  fthl2->SetLineColor(1); fthl2->SetLineWidth(1);
  TPaveText *thl2 = new TPaveText(nn+1,0.0,nn+mm,0.002);
  thl2->SetTextColor(kBlue);
  thl2->AddText(Form("p0 = %0.5e#pm%0.1e",fthl2->GetParameter(0),fthl2->GetParError(0)));
  thl2->SetTextAlign(12);
//  thl1->Draw("sames");
  thl2->Draw("sames");
  lav->Draw();
 // cthL->SaveAs(Form("Th%s%sfit.png",atin.c_str(),hrs.c_str()));

  TCanvas *cthR = new TCanvas();
  Th_R->Draw("AP");
//  Th_R->Fit("pol0","W","",1,nn);
  thRwien[1]->Fit("pol0","W","",nn,nn+mm);
//  TF1 *fthr1 = Th_R->GetFunction("pol0");
//  fthr1->SetLineColor(1); fthr1->SetLineWidth(1);
//  TPaveText *thr1 = new TPaveText(1,-0.03,nn,-0.028);
//  thr1->SetTextColor(kRed);
//  thr1->AddText(Form("p0 = %0.5e#pm%0.1e",fthr1->GetParameter(0),fthr1->GetParError(0)));
  TF1 *fthr2 = thRwien[1]->GetFunction("pol0");
  fthr2->SetLineColor(1); fthr2->SetLineWidth(1);
  TPaveText *thr2 = new TPaveText(nn+1,-0.03,nn+mm,-0.028);
  thr2->SetTextColor(kBlue);
  thr2->AddText(Form("p0 = %0.5e#pm%0.1e",fthr2->GetParameter(0),fthr2->GetParError(0)));
  thr2->SetTextAlign(12);
//  thr1->Draw("sames");
  thr2->Draw("sames");
  lav->Draw();
  cthR->SaveAs(Form("Th%s%sfit.png",atout.c_str(),hrs.c_str()));

  TCanvas *cphL = new TCanvas();
  Ph_L->Draw("AP");
  Ph_L->Fit("pol0","W","",1,nn);
  phLwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fphl1 = Ph_L->GetFunction("pol0");
  fphl1->SetLineColor(1); fphl1->SetLineWidth(1);
  TPaveText *phl1 = new TPaveText(1,var1,nn,var1+0.002);
  phl1->SetTextColor(kRed);
  phl1->AddText(Form("p0 = %0.5e#pm%0.1e",fphl1->GetParameter(0),fphl1->GetParError(0)));
  TF1 *fphl2 = phLwien[1]->GetFunction("pol0");
  fphl2->SetLineColor(1); fphl2->SetLineWidth(1);
  TPaveText *phl2 = new TPaveText(nn+1,var1,nn+mm,var1+0.002);
  phl2->SetTextColor(kBlue);
  phl2->AddText(Form("p0 = %0.5e#pm%0.1e",fphl2->GetParameter(0),fphl2->GetParError(0)));
  phl2->SetTextAlign(12);
//  phl1->Draw("sames");
  phl2->Draw("sames");
  lav->Draw();
 // cphL->SaveAs(Form("Ph%s%sfit.png",atin.c_str(),hrs.c_str()));

  TCanvas *cphR = new TCanvas();
  Ph_R->Draw("AP");
  Ph_R->Fit("pol0","W","",1,nn);
  phRwien[1]->Fit("pol0","W","",nn,nn+mm);
  TF1 *fphr1 = Ph_R->GetFunction("pol0");
  fphr1->SetLineColor(1); fphr1->SetLineWidth(1);
  TPaveText *phr1 = new TPaveText(1,var1,nn,var1+0.002);
  phr1->SetTextColor(kRed);
  phr1->AddText(Form("p0 = %0.5e#pm%0.1e",fphr1->GetParameter(0),fphr1->GetParError(0)));
  TF1 *fphr2 = phRwien[1]->GetFunction("pol0");
  fphr2->SetLineColor(1); fphr2->SetLineWidth(1);
  TPaveText *phr2 = new TPaveText(nn+1,var1,nn+mm,var1+0.002);
  phr2->SetTextColor(kBlue);
  phr2->AddText(Form("p0 = %0.5e#pm%0.1e",fphr2->GetParameter(0),fphr2->GetParError(0)));
  phr2->SetTextAlign(12);
//  phr1->Draw("sames");
  phr2->Draw("sames");
  lav->Draw();
  cphR->SaveAs(Form("Ph%s%sfit.png",atout.c_str(),hrs.c_str()));


  TCanvas *caL = new TCanvas();
  APV_L->Draw("AP");
  APV_L->Fit("pol0","W","",1,nn);
  asLwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fasl1 = APV_L->GetFunction("pol0");
  fasl1->SetLineColor(1); fasl1->SetLineWidth(1);
  TPaveText *asl1 = new TPaveText(1,apvin1,nn,apvin1+0.002);
  asl1->SetTextColor(kRed);
  asl1->AddText(Form("p0 = %0.5e#pm%0.1e",fasl1->GetParameter(0),fasl1->GetParError(0)));
  TF1 *fasl2 = asLwien[1]->GetFunction("pol0");
  fasl2->SetLineColor(1); fasl2->SetLineWidth(1);
  TPaveText *asl2 = new TPaveText(nn+1,apvin1,nn+mm,apvin1+0.002);
  asl2->SetTextColor(kBlue);
  asl2->AddText(Form("p0 = %0.5e#pm%0.1e",fasl2->GetParameter(0),fasl2->GetParError(0)));
  asl2->SetTextAlign(12);
//  asl1->Draw("sames");
  asl2->Draw("sames");
  lav->Draw();
  caL->SaveAs(Form("apv%s%sfit.png",atin.c_str(),hrs.c_str()));

  TCanvas *caR = new TCanvas();
  APV_R->Draw("AP");
  APV_R->Fit("pol0","W","",1,nn);
  asRwien[1]->Fit("pol0","W","",nn,nn+mm);
  TF1 *fasr1 = APV_R->GetFunction("pol0");
  fasr1->SetLineColor(1); fasr1->SetLineWidth(1);
  TPaveText *asr1 = new TPaveText(1,apvout1,nn,apvout1+0.002);
  asr1->SetTextColor(kRed);
  asr1->AddText(Form("p0 = %0.5e#pm%0.1e",fasr1->GetParameter(0),fasr1->GetParError(0)));
  TF1 *fasr2 = asRwien[1]->GetFunction("pol0");
  fasr2->SetLineColor(1); fasr2->SetLineWidth(1);
  TPaveText *asr2 = new TPaveText(nn+1,apvout1,nn+mm,apvout1+0.002);
  asr2->SetTextColor(kBlue);
  asr2->AddText(Form("p0 = %0.5e#pm%0.1e",fasr2->GetParameter(0),fasr2->GetParError(0)));
  asr2->SetTextAlign(12);
//  asr1->Draw("sames");
  asr2->Draw("sames");
  lav->Draw();
  caR->SaveAs(Form("apv%s%sfit.png",atout.c_str(),hrs.c_str()));


  TCanvas *clL = new TCanvas();
  Lab_L->Draw("AP");
  Lab_L->Fit("pol0","W","",1,nn);
  lbLwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *flbl1 = Lab_L->GetFunction("pol0");
  flbl1->SetLineColor(1); flbl1->SetLineWidth(1);
  TPaveText *lbl1 = new TPaveText(1,labin1,nn,labin1+0.002);
  lbl1->SetTextColor(kRed);
  lbl1->AddText(Form("p0 = %0.5e#pm%0.1e",flbl1->GetParameter(0),flbl1->GetParError(0)));
  TF1 *flbl2 = lbLwien[1]->GetFunction("pol0");
  flbl2->SetLineColor(1); flbl2->SetLineWidth(1);
  TPaveText *lbl2 = new TPaveText(nn+1,labin1,nn+mm,labin1+0.002);
  lbl2->SetTextColor(kBlue);
  lbl2->AddText(Form("p0 = %0.5e#pm%0.1e",flbl2->GetParameter(0),flbl2->GetParError(0)));
  lbl2->SetTextAlign(12);
//  lbl1->Draw("sames");
  lbl2->Draw("sames");
  lav->Draw();
  clL->SaveAs(Form("lab%s%sfit.png",atin.c_str(),hrs.c_str()));


  TCanvas *clR = new TCanvas();
  Lab_R->Draw("AP");
  Lab_R->Fit("pol0","W","",1,nn);
  lbRwien[1]->Fit("pol0","W","",nn,nn+mm);
  TF1 *flbr1 = Lab_R->GetFunction("pol0");
  flbr1->SetLineColor(1); flbr1->SetLineWidth(1);
  TPaveText *lbr1 = new TPaveText(1,labout1,nn,labout1+0.002);
  lbr1->SetTextColor(kRed);
  lbr1->AddText(Form("p0 = %0.5e#pm%0.1e",flbr1->GetParameter(0),flbr1->GetParError(0)));
  TF1 *flbr2 = lbRwien[1]->GetFunction("pol0");
  flbr2->SetLineColor(1); flbr2->SetLineWidth(1);
  TPaveText *lbr2 = new TPaveText(nn+1,labout1,nn+mm,labout1+0.002);
  lbr2->SetTextColor(kBlue);
  lbr2->AddText(Form("p0 = %0.5e#pm%0.1e",flbr2->GetParameter(0),flbr2->GetParError(0)));
  lbr2->SetTextAlign(12);
//  lbr1->Draw("sames");
  lbr2->Draw("sames");
  lav->Draw();
  clR->SaveAs(Form("lab%s%sfit.png",atout.c_str(),hrs.c_str()));
 

  TCanvas *ctA = new TCanvas();
  ThAv->Draw("AP");
  ThAv->Fit("pol0","W","",1,nn);
  thavwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fthav1 = ThAv->GetFunction("pol0");
  fthav1->SetLineColor(1); fthav1->SetLineWidth(1);
  TPaveText *thav1 = new TPaveText(1,-0.01,nn,-0.008);
  thav1->SetTextColor(kRed);
  thav1->AddText(Form("p0 = %0.5e#pm%0.1e",fthav1->GetParameter(0),fthav1->GetParError(0)));
  TF1 *fthav2 = thavwien[1]->GetFunction("pol0");
  fthav2->SetLineColor(1); fthav2->SetLineWidth(1);
  TPaveText *thav2 = new TPaveText(nn+1,-0.01,nn+mm,-0.008);
  thav2->SetTextColor(kBlue);
  thav2->AddText(Form("p0 = %0.5e#pm%0.1e",fthav2->GetParameter(0),fthav2->GetParError(0)));
  thav2->SetTextAlign(12);
//  thav1->Draw("sames");
  thav2->Draw("sames");
  lav->Draw();
  ctA->SaveAs(Form("Thatav%sfit.png",hrs.c_str()));

  TCanvas *ctD = new TCanvas();
  ThDD->Draw("AP");
  ThDD->Fit("pol0","W","",1,nn);
  thdwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fthdd1 = ThDD->GetFunction("pol0");
  fthdd1->SetLineColor(1); fthdd1->SetLineWidth(1);
  TPaveText *thdd1 = new TPaveText(1,0.0,nn,0.003);
  thdd1->SetTextColor(kRed);
  thdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fthdd1->GetParameter(0),fthdd1->GetParError(0)));
  TF1 *fthdd2 = thdwien[1]->GetFunction("pol0");
  fthdd2->SetLineColor(1); fthdd2->SetLineWidth(1);
  TPaveText *thdd2 = new TPaveText(nn+1,0.0,nn+mm,0.003);
  thdd2->SetTextColor(kBlue);
  thdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fthdd2->GetParameter(0),fthdd2->GetParError(0)));
  thdd2->SetTextAlign(12);
//  thdd1->Draw("sames");
  thdd2->Draw("sames");
  lav->Draw();
  ctD->SaveAs(Form("Thatdd%sfit.png",hrs.c_str()));


  TCanvas *cpA = new TCanvas();
  PhAv->Draw("AP");
  PhAv->Fit("pol0","W","",1,nn);
  phavwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fphav1 = PhAv->GetFunction("pol0");
  fphav1->SetLineColor(1); fphav1->SetLineWidth(1);
  TPaveText *phav1 = new TPaveText(1,var1,4,var1+0.002);
  phav1->SetTextColor(kRed);
  phav1->AddText(Form("p0 = %0.5e#pm%0.1e",fphav1->GetParameter(0),fphav1->GetParError(0)));
  TF1 *fphav2 = phavwien[1]->GetFunction("pol0");
  fphav2->SetLineColor(1); fphav2->SetLineWidth(1);
  TPaveText *phav2 = new TPaveText(5,var1,9,var1+0.002);
  phav2->SetTextColor(kBlue);
  phav2->AddText(Form("p0 = %0.5e#pm%0.1e",fphav2->GetParameter(0),fphav2->GetParError(0)));
  phav2->SetTextAlign(12);
//  phav1->Draw("sames");
  phav2->Draw("sames");
  lav->Draw();
  cpA->SaveAs(Form("Phatav%sfit.png",hrs.c_str()));

  TCanvas *cpD = new TCanvas();
  PhDD->Draw("AP");
  PhDD->Fit("pol0","W","",1,nn);
  phdwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fphd1 = PhDD->GetFunction("pol0");
  fphd1->SetLineColor(1); fphd1->SetLineWidth(1);
  TPaveText *phd1 = new TPaveText(1,-0.01,nn,-0.008);
  phd1->SetTextColor(kRed);
  phd1->AddText(Form("p0 = %0.5e#pm%0.1e",fphd1->GetParameter(0),fphd1->GetParError(0)));
  TF1 *fphd2 = phdwien[1]->GetFunction("pol0");
  fphd2->SetLineColor(1); fphd2->SetLineWidth(1);
  TPaveText *phd2 = new TPaveText(nn+1,-0.01,nn+mm,-0.008);
  phd2->SetTextColor(kBlue);
  phd2->AddText(Form("p0 = %0.5e#pm%0.1e",fphd2->GetParameter(0),fphd2->GetParError(0)));
  phd2->SetTextAlign(12);
//  phd1->Draw("sames");
  phd2->Draw("sames");
  lav->Draw();
  cpD->SaveAs(Form("Phatdd%sfit.png",hrs.c_str()));


  TCanvas *caA = new TCanvas();
  APVAvg->Draw("AP");
  APVAvg->Fit("pol0","W","",1,nn);
  asavwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fasav1 = APVAvg->GetFunction("pol0");
  fasav1->SetLineColor(1); fasav1->SetLineWidth(1);
  TPaveText *asav1 = new TPaveText(1,apvin1,nn,apvin1+0.005);
  asav1->SetTextColor(kRed);
  asav1->AddText(Form("p0 = %0.5e#pm%0.1e",fasav1->GetParameter(0),fasav1->GetParError(0)));
  TF1 *fasav2 = asavwien[1]->GetFunction("pol0");
  fasav2->SetLineColor(1); fasav2->SetLineWidth(1);
  TPaveText *asav2 = new TPaveText(nn,apvin1,nn+mm,apvin1+0.005);
  asav2->SetTextColor(kBlue);
  asav2->AddText(Form("p0 = %0.5e#pm%0.1e",fasav2->GetParameter(0),fasav2->GetParError(0)));
  asav2->SetTextAlign(12);
//  asav1->Draw("sames");
  asav2->Draw("sames");
  lav->Draw();
  caA->SaveAs(Form("APVatav%sfit.png",hrs.c_str()));
*/

  TCanvas *caD = new TCanvas();
  APVDD->Draw("AP");
//  APVDD->Fit("pol0","W","",1,nn);
  asdwien[1]->Fit("pol0","W","",nn+1,nn+mm);
//  TF1 *fasad1 = APVDD->GetFunction("pol0");
 // fasad1->SetLineColor(1); fasad1->SetLineWidth(1);
//  TPaveText *asad1 = new TPaveText(1,apvd1,nn,apvd1+0.003);
//  asad1->SetTextColor(kRed);
//  asad1->AddText(Form("p0 = %0.5e#pm%0.1e",fasad1->GetParameter(0),fasad1->GetParError(0)));
  TF1 *fasad2 = asdwien[1]->GetFunction("pol0");
  fasad2->SetLineColor(1); fasad2->SetLineWidth(1);
  TPaveText *asad2 = new TPaveText(nn+1,apvd1,nn+mm,apvd1+0.003);
  asad2->SetTextColor(kBlue);
  asad2->AddText(Form("p0 = %0.5e#pm%0.1e",fasad2->GetParameter(0),fasad2->GetParError(0)));
  asad2->SetTextAlign(12);
//  asad1->Draw("sames");
  asad2->Draw("sames");
  lav->Draw();
  caD->SaveAs(Form("APVatdd%sfit.png",hrs.c_str()));
/*
  TCanvas *cqA = new TCanvas();
  Q2Av->Draw("AP");
  Q2Av->Fit("pol0","W","",1,nn);
  qavwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fav1 = Q2Av->GetFunction("pol0");
  fav1->SetLineColor(1); fav1->SetLineWidth(1);
  TPaveText *qav1 = new TPaveText(1,qa1,nn,qa1+0.0001);
  qav1->SetTextColor(kRed);
  qav1->AddText(Form("p0 = %0.5e#pm%0.1e",fav1->GetParameter(0),fav1->GetParError(0)));
  TF1 *fav2 = qavwien[1]->GetFunction("pol0");
  fav2->SetLineColor(1); fav2->SetLineWidth(1);
  TPaveText *qav2 = new TPaveText(nn+1,qa1,nn+mm,qa1+0.0001);
  qav2->SetTextColor(kBlue);
  qav2->AddText(Form("p0 = %0.5e#pm%0.1e",fav2->GetParameter(0),fav2->GetParError(0)));
  qav2->SetTextAlign(12);
//  qav1->Draw("sames");
  qav2->Draw("sames");
  lav->Draw();
  cqA->SaveAs(Form("Q2atav%sfit.png",hrs.c_str()));

  TCanvas *cqD = new TCanvas();
  Q2DD->Draw("AP");
  Q2DD->Fit("pol0","W","",1,nn);
  qdwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fdd1 = Q2DD->GetFunction("pol0");
  fdd1->SetLineColor(1); fdd1->SetLineWidth(1);
  TPaveText *qdd1 = new TPaveText(1,qd1,nn,qd1+0.001);
  qdd1->SetTextColor(kRed);
  qdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fdd1->GetParameter(0),fdd1->GetParError(0)));
  TF1 *fdd2 = qdwien[1]->GetFunction("pol0");
  fdd2->SetLineColor(1); fdd2->SetLineWidth(1);
  TPaveText *qdd2 = new TPaveText(nn+1,qd1,nn+mm,qd1+0.001);
  qdd2->SetTextColor(kBlue);
  qdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fdd2->GetParameter(0),fdd2->GetParError(0)));
  qdd2->SetTextAlign(12);
//  qdd1->Draw("sames");
  qdd2->Draw("sames");
  lav->Draw();
  cqD->SaveAs(Form("Q2atdd%sfit.png",hrs.c_str()));
*/

/*  
  TCanvas *cchi = new TCanvas();
  Chifac->Draw("AP");
//  Chifac->Fit("pol0","W","",1,nn); 
  chiwien[1]->Fit("pol0","W","",nn+1,nn+mm);
//  TF1 *fchi1 = Chifac->GetFunction("pol0");
//  fchi1->SetLineColor(1); fchi1->SetLineWidth(1);
//  TPaveText *qchi1 = new TPaveText(1,-10,nn,-9);  
//  qchi1->SetTextColor(kRed);
//  qchi1->AddText(Form("p0 =%0.5e#pm%0.1e",fchi1->GetParameter(0),fchi1->GetParError(0)));  
  TF1 *fchi2 = chiwien[1]->GetFunction("pol0");
  fchi2->SetLineColor(1); fchi2->SetLineWidth(1);
  TPaveText *qchi2 = new TPaveText(nn+1,-10,mm+nn,-9);
  qchi2->SetTextColor(kBlue);
  qchi2->AddText(Form("p0 =%0.5e#pm%0.1e",fchi2->GetParameter(0),fchi2->GetParError(0)));
  qchi2->SetTextAlign(12);
//  qchi1->Draw("sames");
  qchi2->Draw("sames");
  lav->Draw();
  cchi->SaveAs(Form("Chi%sfit.png",hrs.c_str()));


 TCanvas *c_cosL = new TCanvas();
 Cos_L->Draw("AP");
 Cos_L->Fit("pol0","W","",1,nn);
  cosLwien[1]->Fit("pol0","W","",nn+1,nn+mm);
  TF1 *fcL1 = Cos_L->GetFunction("pol0");
  fcL1->SetLineColor(1); fcL1->SetLineWidth(1);
  TPaveText *qcL1 = new TPaveText(1,-0.06,nn,-0.055);
  qcL1->SetTextColor(kRed);
  qcL1->AddText(Form("p0 = %0.5e#pm%0.1e",fcL1->GetParameter(0),fcL1->GetParError(0)));
  TF1 *fcL2 = cosLwien[1]->GetFunction("pol0");
  fcL2->SetLineColor(1); fcL2->SetLineWidth(1);
  TPaveText *qcL2 = new TPaveText(5,-0.06,nn+mm,-0.055);
  qcL2->SetTextColor(kBlue);
  qcL2->AddText(Form("p0 = %0.5e#pm%0.1e",fcL2->GetParameter(0),fcL2->GetParError(0)));
  lav->Draw();
  qcL1->Draw("sames");
  qcL2->Draw("sames");

 c_cosL->SaveAs(Form("%s%scosfit.png",hrs.c_str(),atin.c_str()));

 TCanvas *c_cosR = new TCanvas();
 Cos_R->Draw("AP");

  Cos_R->Fit("pol0","W","",1,4);
  cosRwien[1]->Fit("pol0","W","",5,9);
  TF1 *fcR1 = Cos_R->GetFunction("pol0");
  fcR1->SetLineColor(1); fcR1->SetLineWidth(1);
  TPaveText *qcR1 = new TPaveText(1,0.0,4,0.005);
  qcR1->SetTextColor(kRed);
  qcR1->AddText(Form("p0 = %0.5e#pm%0.1e",fcR1->GetParameter(0),fcR1->GetParError(0)));
  TF1 *fcR2 = cosRwien[1]->GetFunction("pol0");
  fcR2->SetLineColor(1); fcR2->SetLineWidth(1);
  TPaveText *qcR2 = new TPaveText(5,0.,9,0.005);
  qcR2->SetTextColor(kBlue);
  qcR2->AddText(Form("p0 = %0.5e#pm%0.1e",fcR2->GetParameter(0),fcR2->GetParError(0)));
  lav->Draw();
  qcR1->Draw("sames");
  qcR2->Draw("sames");
*/
// c_cosR->SaveAs(Form("%s%scos.png",hrs.c_str(),atout.c_str()));

/*
  TCanvas *c_sinL = new TCanvas();
  Sin_L->Draw("AP");
  Sin_L->Fit("pol0","W","",1,4);
  sinLwien[1]->Fit("pol0","W","",5,9);
  TF1 *fsL1 = Sin_L->GetFunction("pol0");
  fsL1->SetLineColor(1); fsL1->SetLineWidth(1);
  TPaveText *qssL1 = new TPaveText(1,0.9,4,0.91);
  qssL1->SetTextColor(kRed);
  qssL1->AddText(Form("p0 = %0.5e#pm%0.1e",fsL1->GetParameter(0),fsL1->GetParError(0)));
  TF1 *fsL2 = sinLwien[1]->GetFunction("pol0");
  fsL2->SetLineColor(1); fsL2->SetLineWidth(1);
  TPaveText *qssL2 = new TPaveText(5,0.9,9,0.91);
  qssL2->SetTextColor(kBlue);
  qssL2->AddText(Form("p0 = %0.5e#pm%0.1e",fsL2->GetParameter(0),fsL2->GetParError(0)));
  lav->Draw();
  qssL1->Draw("sames");
  qssL2->Draw("sames");

  TCanvas *c_sinR = new TCanvas();
  Sin_R->Draw("AP");
  Sin_R->Fit("pol0","W","",1,4);
  sinRwien[1]->Fit("pol0","W","",5,9);
  TF1 *fsR1 = Sin_R->GetFunction("pol0");
  fsR1->SetLineColor(1); fsR1->SetLineWidth(1);
  TPaveText *qssR1 = new TPaveText(1,0.9,4,0.91);
  qssR1->SetTextColor(kRed);
  qssR1->AddText(Form("p0 = %0.5e#pm%0.1e",fsR1->GetParameter(0),fsR1->GetParError(0)));
  TF1 *fsR2 = sinRwien[1]->GetFunction("pol0");
  fsR2->SetLineColor(1); fsR2->SetLineWidth(1);
  TPaveText *qssR2 = new TPaveText(5,0.9,9,0.91);
  qssR2->SetTextColor(kBlue);
  qssR2->AddText(Form("p0 = %0.5e#pm%0.1e",fsR2->GetParameter(0),fsR2->GetParError(0)));
  lav->Draw();
  qssR1->Draw("sames");
  qssR2->Draw("sames");
*/
/*
 TCanvas *c_cosav = new TCanvas();
 CosAv->Draw("AP");
 CosAv->Fit("pol0","W","",1,4);
  cosavwien[1]->Fit("pol0","W","",5,9);
  TF1 *fcavg1 = CosAv->GetFunction("pol0");
  fcavg1->SetLineColor(1); fcavg1->SetLineWidth(1);
  TPaveText *qcavg1 = new TPaveText(1,-0.006,4,-0.005);
  qcavg1->SetTextColor(kRed);
  qcavg1->AddText(Form("p0 = %0.5e#pm%0.1e",fcavg1->GetParameter(0),fcavg1->GetParError(0)));
  TF1 *fcavg2 = cosavwien[1]->GetFunction("pol0");
  fcavg2->SetLineColor(1); fcavg2->SetLineWidth(1);
  TPaveText *qcavg2 = new TPaveText(5,-0.006,9,-0.005);
  qcavg2->SetTextColor(kBlue);
  qcavg2->AddText(Form("p0 = %0.5e#pm%0.1e",fcavg2->GetParameter(0),fcavg2->GetParError(0)));
  lav->Draw();
  qcavg1->Draw("sames");
  qcavg2->Draw("sames");
*/

 TCanvas *c_cosd = new TCanvas();
 CosDD->Draw("AP");



/*
 CosDD->Fit("pol0","W","",1,4);
  cosdwien[1]->Fit("pol0","W","",5,9);
  TF1 *fcdd1 = CosDD->GetFunction("pol0");
  fcdd1->SetLineColor(1); fcdd1->SetLineWidth(1);
  TPaveText *qcdd1 = new TPaveText(1,-0.06,4,-0.05);
  qcdd1->SetTextColor(kRed);
  qcdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fcdd1->GetParameter(0),fcdd1->GetParError(0)));
  TF1 *fcdd2 = cosdwien[1]->GetFunction("pol0");
  fcdd2->SetLineColor(1); fcdd2->SetLineWidth(1);
  TPaveText *qcdd2 = new TPaveText(5,-0.06,9,-0.05);
  qcdd2->SetTextColor(kBlue);
  qcdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fcdd2->GetParameter(0),fcdd2->GetParError(0)));
  lav->Draw();
  qcdd1->Draw("sames");
  qcdd2->Draw("sames");
*/
//  c_cosd->SaveAs(Form("%scosdd.png",hrs.c_str()));

/*
 TCanvas *c_sinav = new TCanvas();
  SinAv->Draw("AP");
  SinAv->Fit("pol0","W","",1,4);
  sinavwien[1]->Fit("pol0","W","",5,9);
  TF1 *fsavg1 = SinAv->GetFunction("pol0");
  fsavg1->SetLineColor(1); fsavg1->SetLineWidth(1);
  TPaveText *qssavg1 = new TPaveText(1,0.9,4,0.91);
  qssavg1->SetTextColor(kRed);
  qssavg1->AddText(Form("p0 = %0.5e#pm%0.1e",fsavg1->GetParameter(0),fsavg1->GetParError(0)));
  TF1 *fsavg2 = sinavwien[1]->GetFunction("pol0");
  fsavg2->SetLineColor(1); fsavg2->SetLineWidth(1);
  TPaveText *qssavg2 = new TPaveText(5,0.9,9,0.91);
  qssavg2->SetTextColor(kBlue);
  qssavg2->AddText(Form("p0 = %0.5e#pm%0.1e",fsavg2->GetParameter(0),fsavg2->GetParError(0)));
  lav->Draw();
  qssavg1->Draw("sames");
  qssavg2->Draw("sames");

 TCanvas *c_sind = new TCanvas();
 SinDD->Draw("AP");
  SinDD->Fit("pol0","W","",1,4);
  sindwien[1]->Fit("pol0","W","",5,9);
  TF1 *fsdd1 = SinDD->GetFunction("pol0");
  fsdd1->SetLineColor(1); fsdd1->SetLineWidth(1);
  TPaveText *qssdd1 = new TPaveText(1,-0.006,4,-0.005);
  qssdd1->SetTextColor(kRed);
  qssdd1->AddText(Form("p0 = %0.5e#pm%0.1e",fsdd1->GetParameter(0),fsdd1->GetParError(0)));
  TF1 *fsdd2 = sindwien[1]->GetFunction("pol0");
  fsdd2->SetLineColor(1); fsdd2->SetLineWidth(1);
  TPaveText *qssdd2 = new TPaveText(5,-0.006,9,-0.005);
  qssdd2->SetTextColor(kBlue);
  qssdd2->AddText(Form("p0 = %0.5e#pm%0.1e",fsdd2->GetParameter(0),fsdd2->GetParError(0)));
  lav->Draw();
  qssdd1->Draw("sames");
  qssdd2->Draw("sames");
*/
}

void ATCorr::PlotATDetPlots(int i){

 //These files have different number of entries 
 ReadFiles("ThtgLHRSmain_newDB_combo.txt","ThtgRHRSmain_newDB.txt");

 //Erasing the first three elements
// thtgmain1.erase(thtgmain1.begin(),thtgmain1.begin()+3);
// thtgmain1err.erase(thtgmain1err.begin(),thtgmain1err.begin()+3);
 
 
 //Now we can compute 
 ComputeChi(4,3);
 //  ComputeChi(4);

 //I am only going to draw three plots
 Chifac = new TMultiGraph();
 Chifac->SetTitle("#chi Factor; Run Number; #chi");
 Chifac->SetMinimum(-10); Chifac->SetMaximum(0);

 APVAvg = new TMultiGraph();
 APVAvg->SetTitle("Average of <A_{PV}>_{atldd} and <A_{PV}>_{atrdd}; Run Number; Asymmetry (ppm)");
 APVAvg->SetMinimum(-0.015); APVAvg->SetMaximum(0.005);

 ThAv = new TMultiGraph();
 ThAv->SetTitle("Average of <#theta_{tg}>_{atldd} and <#theta_{tg}>_{atrdd}; Run Number; <#theta>");
 ThAv->SetMinimum(0); ThAv->SetMaximum(0.3); 

 int color[2] = {2,4};
 //Preparing Data 
 double apvwien2[1],apvwien2err[1], apvwien3[5], apvwien3err[5], thwien2[1], thwien2err[1], thwien3[5], thwien3err[5];
 double wien2[1] = {1}; double wien2err[1] = {0}; 
 double wien3[5] = { 2,3,4,5,6};  
 double wien3err[5] = {0};
 double chi2[1], chi2err[1], chi3[5], chi3err[5];

 int ii = 0;
 for(int i = 0; i < thtgmain1.size(); i++){
  if(i == 0){ 
     apvwien2[i] = AsymAll[i]; apvwien2err[i] = AsymAllerr[i];
     chi2[i] = Chi[i]; chi2err[i] = Chierr[i];
     thwien2[i] = CosAll[i]; thwien2err[i] = CosAllerr[i];
     
    cout << thwien2[i] << "   " << thwien2err[i] << endl;
 //    cout << apvwien2[i] << "   " << apvwien2err[i] << endl;
     cout << chi2[i] << "  " << chi2err[i] << endl;
  } else {         
     apvwien3[ii] = AsymAll[i]; apvwien3err[ii] = AsymAllerr[i];
     chi3[ii] = Chi[i]; chi3err[i] = Chierr[i];
     thwien3[ii] = CosAll[i]; thwien3err[ii] = CosAllerr[i];
  //   cout << chi3[ii] <<  "  " << chi3err[ii] << endl;
   //  cout << apvwien3[ii] << "    " << apvwien3err[ii] << endl;
     ii++; 
  }
 
 } 

  TLegend *lav = new TLegend(0.55,0.7,0.9,0.9);

  TGraphErrors *chiwien[2], *apvwien[2], *thwien[2];
  chiwien[0] = new TGraphErrors(1,wien2,chi2,wien2err,chi2err); chiwien[1] = new TGraphErrors(5,wien3,chi3,wien3err,chi3err);
  apvwien[0] = new TGraphErrors(1,wien2,apvwien2,wien2err,apvwien2err); apvwien[1] = new TGraphErrors(5,wien3,apvwien3,wien3err,apvwien3err);
  thwien[0] = new TGraphErrors(1,wien2,thwien2,wien2err,thwien2err); thwien[1] = new TGraphErrors(5,wien3,thwien3,wien3err,thwien3err);

  for(int i = 0; i < 2; i++){ 
   chiwien[i]->SetMarkerStyle(21); apvwien[i]->SetMarkerStyle(21); thwien[i]->SetMarkerStyle(21);
   chiwien[i]->SetMarkerColor(color[i]); apvwien[i]->SetMarkerColor(color[i]); thwien[i]->SetMarkerColor(color[i]);

   Chifac->Add(chiwien[i]); APVAvg->Add(apvwien[i]); 
   ThAv->Add(thwien[i]);

  if(i == 0) { lav->AddEntry(chiwien[i],"Wien Left","p"); } else { lav->AddEntry(chiwien[i],"Wien Right","p"); }
  }
/*
  
  TCanvas *cchi = new TCanvas();
  Chifac->Draw("AP");
  Chifac->Fit("pol0","W","",1,nn); 
  chiwien[1]->Fit("pol0","W","",2,6);
  TF1 *fchi1 = Chifac->GetFunction("pol0");
  fchi1->SetLineColor(1); fchi1->SetLineWidth(1);
  TPaveText *qchi1 = new TPaveText(1,-10,nn,-9);  
  qchi1->SetTextColor(kRed);
  qchi1->AddText(Form("p0 =%0.5e#pm%0.1e",fchi1->GetParameter(0),fchi1->GetParError(0)));  
  TF1 *fchi2 = chiwien[1]->GetFunction("pol0");
  fchi2->SetLineColor(1); fchi2->SetLineWidth(1);
  TPaveText *qchi2 = new TPaveText(2,-10,6,-9);
  qchi2->SetTextColor(kBlue);
  qchi2->AddText(Form("p0 =%0.5e#pm%0.1e",fchi2->GetParameter(0),fchi2->GetParError(0)));
  qchi2->SetTextAlign(12);
  qchi1->Draw("sames");
   qchi2->Draw("sames");
  lav->Draw();
  cchi->SaveAs("Chifit_thetadb1.png");

   
                   
 TCanvas *cA = new TCanvas();
 APVAvg->Draw("AP");
 apvwien[1]->Fit("pol0","W","",2,6);
 TF1 *fas2 = apvwien[1]->GetFunction("pol0");
 fas2->SetLineColor(1); fas2->SetLineWidth(1);
 TPaveText *qas2 = new TPaveText(2,-0.015,6,-0.013);
 qas2->SetTextColor(kBlue);
 qas2->AddText(Form("p0 =%0.5e#pm%0.1e",fas2->GetParameter(0),fas2->GetParError(0)));
 qas2->SetTextAlign(12);
 lav->Draw();
 qas2->Draw("sames");
 cA->SaveAs("APV_AVg.png");
*/



 TCanvas *cT = new TCanvas();
 ThAv->Draw("AP");
 thwien[1]->Fit("pol0","W","",2,6);
 TF1 *fth2 = thwien[1]->GetFunction("pol0");
 fth2->SetLineColor(1); fth2->SetLineWidth(1);
 TPaveText *qth2 = new TPaveText(2,0.0,6,0.05);
 qth2->SetTextColor(kBlue);
 qth2->AddText(Form("p0 =%0.5e#pm%0.1e",fth2->GetParameter(0),fth2->GetParError(0)));
 qth2->SetTextAlign(12);
 lav->Draw();
 qth2->Draw("sames");
 cT->SaveAs("Cos_Avgdb1b.png");

 
 



}
void ATCorr::ComputeXi(){
  

  for(int i = 0; i < QsqL.size(); i++){
     double num = PhtgL[i]+PhtgR[i]; double den = PhtgL[i]-PhtgR[i];
     double err = sqrt(PhtgerL[i]*PhtgerL[i]+PhtgerR[i]*PhtgerR[i]);//errors are the same
  
     Xi[i] = (PhtgL[i]+PhtgR[i])/(PhtgL[i]-PhtgR[i]);
     Xierr[i] = fabs(Xi[i])*sqrt( pow( (err/num),2 ) + pow( (err/den),2 ) );   
  }

}

void ATCorr::ComputeXi(int i){

   for(int i = 0; i < QsqL.size(); i++){
      double num = sinphL[i]-sinphR[i];
      double den = sinphL[i]+sinphR[i];
      double err = sqrt(sinpherL[i]*sinpherL[i]+sinpherR[i]*sinpherR[i]); 
   
      Xi[i] = num/den;
      Xierr[i] = fabs(Xi[i])*sqrt( pow( (err/num),2 ) + pow( (err/den),2 ) );
   

  }

}

void ATCorr::ComputeXi(int i, int j){
  for(int i = 0; i < QsqL.size(); i++){
   double num = VermomL[i]-VermomR[i];
   double den = VermomR[i]+VermomR[i];
   double err = sqrt(VermomerL[i]*VermomerL[i]+VermomerR[i]*VermomerR[i]);

   Xi[i] = num/den;
   Xierr[i] = fabs(Xi[i])*sqrt( pow( (err/num),2 ) + pow( (err/den),2 ) );
   

   }
}

void ATCorr::ComputeQsqAvg(){

  for(int i = 0; i < QsqL.size(); i++){
    QsqAvg[i] = 0.5*(QsqL[i]+QsqR[i]); 
    QsqAvgerr[i] = 0.5*sqrt( pow( QsqerL[i],2 ) + pow( QsqerR[i],2) );

  }

}

void ATCorr::ComputeQsqDiff(){

   for(int i = 0; i < QsqL.size(); i++){
    QsqDiff[i] = 0.5*(QsqL[i]-QsqR[i]);
    QsqDifferr[i] = 0.5*sqrt( pow( QsqerL[i],2 ) + pow( QsqerR[i],2 ) );

  }


}

void ATCorr::ComputeAsymAvg(){

  for(int i = 0; i < asymL.size(); i++){
    AsymAvg[i] = 0.5*(asymL[i]+asymR[i]);
    AsymAvgerr[i] = 0.5*sqrt( pow( (asymerL[i]),2 ) + pow( asymerR[i],2) );

  }

}

void ATCorr::ComputeAsymDiff(){

   for(int i = 0; i < asymL.size(); i++){
    AsymDiff[i] = 0.5*(asymL[i]-asymR[i]);
    AsymDifferr[i] = 0.5*sqrt( pow( asymerL[i],2 ) + pow( asymerR[i],2) );

  }


}

void ATCorr::ComputeThtgAvg(){

  for(int i = 0; i < ThtgL.size(); i++){
    ThtgAvg[i] = 0.5*(ThtgL[i]+ThtgR[i]);
    ThtgAvgerr[i] = 0.5*ThtgAvg[i]*sqrt( pow( ThtgerL[i],2 ) + pow( ThtgerR[i],2) );

  }

}

void ATCorr::ComputeThtgDiff(){

   for(int i = 0; i < ThtgL.size(); i++){
    ThtgDiff[i] = 0.5*(ThtgL[i]-ThtgR[i]);
    ThtgDifferr[i] = 0.5*sqrt( pow( ThtgerL[i],2 ) + pow( ThtgerR[i],2) );

  }


}

void ATCorr::ComputePhtgAvg(){

  for(int i = 0; i < PhtgL.size(); i++){
    PhtgAvg[i] = 0.5*(PhtgL[i]+PhtgR[i]);
    PhtgAvgerr[i] = 0.5*sqrt( pow( PhtgerL[i],2 ) + pow( PhtgerR[i],2) );

  }

}

void ATCorr::ComputePhtgDiff(){

   for(int i = 0; i < PhtgL.size(); i++){
    PhtgDiff[i] = 0.5*(PhtgL[i]-PhtgR[i]);
    PhtgDifferr[i] = 0.5*sqrt( pow( PhtgerL[i],2 ) + pow( PhtgerR[i],2) );

  }


}


void ATCorr::ReadFiles(string f1){
  
  ifstream F;
  F.open(f1.c_str());
  
  double run, tg, tgerr, cos, coserr;

  while(F >> run >> tg >> tgerr >> cos >> coserr){ thtgmain1.push_back(tg); thtgmain1err.push_back(tgerr); cosmain1.push_back(cos); cosmain1err.push_back(coserr);}

}

void ATCorr::ReadFiles(string f1, string f2){

  ifstream F1, F2;
  F1.open(f1.c_str());
  F2.open(f2.c_str());

  double run1, run2, tg1, tg2, tgerr1, tgerr2, cos1, cos2, coserr1, coserr2;

 while(F1 >> run1 >> tg1 >> tgerr1 >> cos1 >> coserr1){ thtgmain1.push_back(tg1); thtgmain1err.push_back(tgerr1); cosmain1.push_back(cos1); cosmain1err.push_back(coserr1); }
 while(F2 >> run2 >> tg2 >> tgerr2 >> cos2 >> coserr2){ thtgmain2.push_back(tg2); thtgmain2err.push_back(tgerr2); cosmain2.push_back(cos2); cosmain2err.push_back(coserr2); }

}

void ATCorr::ComputeChi(){

   for(int i = 0; i < ThtgL.size(); i++){
//   Chi[i] = 0.5*(ThtgL[i]-ThtgR[i])/thtgmain1[i];
//   double num = 0.5*(ThtgL[i]-ThtgR[i]);
//   double numerr = 0.5*sqrt(ThtgerL[i]*ThtgerL[i]+ThtgerR[i]*ThtgerR[i]);  
     double num = 0.5*(cosphL[i]-cosphR[i]);
     double numerr = 0.5*sqrt(cospherL[i]*cospherL[i] + cospherR[i]*cospherR[i]);
     Chi[i] = 0.5*(cosphL[i] - cosphR[i])/cosmain1[i];
     Chierr[i] = fabs(Chi[i])*sqrt(pow((cosmain1err[i]/cosmain1[i]),2) + pow( (numerr/num),2));
  }


}


void ATCorr::ComputeChi(int j){

  for(int i = 0; i < ThtgL.size(); i++){ 
    double num1 = 0.5*(ThtgL[i]-ThtgR[i]); double num2 = 0.5*(ThtgL1[i]-ThtgR1[i]);
    double num = (num1+num2);
    double numerr = 0.5*sqrt(ThtgerL[i]*ThtgerL[i]+ThtgerR[i]*ThtgerR[i]+ThtgerL1[i]*ThtgerL1[i]+ThtgerR1[i]*ThtgerR1[i]); 
    double den = thtgmain1[i]+thtgmain2[i];
    double denerr = sqrt(thtgmain1err[i]*thtgmain1err[i]+thtgmain2err[i]*thtgmain2err[i]);

   Chi[i] = (0.5*(ThtgL[i]-ThtgR[i])+0.5*(ThtgL1[i]-ThtgR1[i]))/(thtgmain1[i]-thtgmain2[i]);
   Chierr[i] = fabs(Chi[i])*sqrt( pow((numerr/num),2)+pow((denerr/den),2) );

  }


}

void ATCorr::ComputeChi(int j, int k){
 
    for(int i = 0; i < cosphL.size(); i++){
     double num1 = 0.5*(cosphL[i]-cosphR[i]); double num2 = 0.5*(cosphL1[i]-cosphR1[i]);   
     double numerr = 0.5*sqrt(cospherL[i]*cospherL[i]+cospherR[i]*cospherR[i]+cospherL1[i]*cospherL1[i]+cospherR1[i]*cospherR1[i]);
     double num = num1-num2;
     double den = cosmain1[i]-cosmain2[i];
     double denerr = sqrt(cosmain1err[i]*cosmain1err[i]+cosmain2err[i]*cosmain2err[i]);
     Chi[i] =  num/den;
     Chierr[i] = fabs(Chi[i])*sqrt(pow((numerr/num),2)+pow((denerr/den),2));

  
   }
 

}

void ATCorr::ComputeTh(){
  
  for(int i = 0; i < ThtgL.size(); i++){
 //  ThAll[i] = 0.25*(ThtgL[i]-ThtgR[i]+(ThtgL1[i]-ThtgR1[i]));
 //  ThAllerr[i] = 0.25*sqrt(ThtgerL[i]*ThtgerL[i]+ThtgerR[i]*ThtgerR[i]+ThtgerL1[i]*ThtgerL1[i]+ThtgerR1[i]*ThtgerR1[i]);
     CosAll[i] = 0.25*(cosphL[i]-cosphR[i]-(cosphL1[i]-cosphR1[i]));
    CosAllerr[i] = 0.25*sqrt(cospherL[i]*cospherL[i]+cospherR[i]*cospherR[i]+cospherL1[i]*cospherL1[i]+cospherR1[i]*cospherR1[i]);

  
  }
}

void ATCorr::ComputeAsym(){
  for(int i = 0; i < ThtgL.size(); i++){
    AsymAll[i] = 0.25*((asymL[i]-asymR[i])+(asymL1[i]-asymR1[i]));
    AsymAllerr[i] = 0.25*sqrt(asymerL[i]*asymerL[i]+asymerR[i]*asymerR[i]+asymerL1[i]*asymerL1[i]+asymerR1[i]*asymerR1[i]);
  }



}


void ATCorr::ComputeVerAvg(){
 
   for(int i = 0; i < VermomL.size(); i++){
    VerMomAvg[i] = 0.5*(VermomL[i]+VermomR[i]);
    VerMomAvgerr[i] = 0.5*sqrt(VermomerL[i]*VermomerL[i] + VermomerR[i]*VermomerR[i]);
  }

}

void ATCorr::ComputeVerDiff(){
  
  for(int i = 0; i < VermomL.size(); i++){
  VerMomDiff[i] = 0.5*(VermomL[i]-VermomR[i]); 
  VerMomDifferr[i] = 0.5*sqrt(VermomerL[i]*VermomerL[i] + VermomerR[i]*VermomerR[i]);
 
  }

}

void ATCorr::ComputeHorAvg(){

   for(int i = 0; i < HormomL.size(); i++){
    HorMomAvg[i] = 0.5*(HormomL[i]+HormomR[i]);
    HorMomAvgerr[i] = 0.5*sqrt(HormomerL[i]*HormomerL[i] + HormomerR[i]*HormomerR[i]);
  }

}

void ATCorr::ComputeHorDiff(){

  for(int i = 0; i < HormomL.size(); i++){
  HorMomDiff[i] = 0.5*(HormomL[i]-HormomR[i]);
  HorMomDifferr[i] = 0.5*sqrt(HormomerL[i]*HormomerL[i] + HormomerR[i]*HormomerR[i]);

  }

}

void ATCorr::ComputeCosDiff(){

 for(int i = 0; i < cosphL.size(); i++){
  CosDiff[i] = 0.5*(cosphL[i]-cosphR[i]);
  CosDifferr[i] = 0.5*sqrt(cospherL[i]*cospherL[i] + cospherR[i]*cospherR[i]);


  }
}

void ATCorr::ComputeCosAvg(){

 for(int i = 0; i < cosphL.size(); i++){
  CosAvg[i] = 0.5*(cosphL[i]+cosphR[i]);
  CosAvgerr[i] = 0.5*sqrt(cospherL[i]*cospherL[i] + cospherR[i]*cospherR[i]);

  }
}

void ATCorr::ComputeSinDiff(){

 for(int i = 0; i < sinphL.size(); i++){
  SinDiff[i] = 0.5*(sinphL[i]-sinphR[i]);
  SinDifferr[i] = 0.5*sqrt(sinpherL[i]*sinpherL[i] + sinpherR[i]*sinpherR[i]);

  }
}

void ATCorr::ComputeSinAvg(){

 for(int i = 0; i < sinphL.size(); i++){
  SinAvg[i] = 0.5*(sinphL[i]+sinphR[i]);
  SinAvgerr[i] = 0.5*sqrt(sinpherL[i]*sinpherL[i] + sinpherR[i]*sinpherR[i]);

  }
}


TGraph* ATCorr::trans(){

  const int nn = 9;

  TGraph *t = new TGraph(2*nn);

  for(int i = 0; i < nn; i++){
    t->SetPoint(i,i+1,0.071);
    t->SetPoint(nn+i,nn-i,-0.145); 
    
   }
  


  t->SetFillStyle(3013); t->SetFillColor(16);
  

  return t;
 

}



