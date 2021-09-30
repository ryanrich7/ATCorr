#ifndef ATCorr_H
#define ATCorr_H

#include "TGraph.h"
#include <vector>
#include <string>

 class TGraph;

 class ATCorr  {
    public:

    enum arm { kLHRS, kRHRS };

    ATCorr(string, string); //reads text files for A_T correction. Form Vertical, LHRS and RHRS main detector files, for A_T correction use files on same arm
//  ATCorr(string, string, string, string); //for 4 detector combo correction
    ~ATCorr();
    
    void PlotMainDetPlots();
    void PlotATDetPlots(arm Arm = kLHRS);
    void PlotATDetPlots(int );//4 detector combo
    int Init();
 

    private:   

    //Computes kinematics
    void ComputeXi(); //left right apparatus asymmetry -- vertical A_T correction
    void ComputeQsqAvg();
    void ComputeQsqDiff();
    void ComputeAsymAvg();
    void ComputeAsymDiff();
    void ComputeThtgAvg();
    void ComputeThtgDiff();
    void ComputePhtgAvg();
    void ComputePhtgDiff();
    void ComputeChi(); //up down acceptance factor -- horizontal A_T correction per arm
    //void ComputeChi(int);//up down acceptance factor -- correct the average
    
    TMultiGraph *Xifac, *Chifac; 

    //Other Kinematic Quantities
    TMultiGraph *Q2_L, *Q2_R, *Q2Av, *Q2DD, *Th_L, *Th_R, *ThAv, *ThDD, *Ph_L, *Ph_R, *PhAv, *PhDD, *APV_L, *APV_R, *APVAvg, *APVDD, *Lab_L, *Lab_R;
   

    //detectors
    vector<double> QsqL, AngleL, asymL, ThtgL, PhtgL;
    vector<double> QsqR, AngleR, asymR, ThtgR, PhtgR;
    vector<double> QsqerL, AngleerL, asymerL, ThtgerL, PhtgerL;
    vector<double> QsqerR, AngleerR, asymerR, ThtgerR, PhtgerR;
    vector<int> RunL,RunR;
    //Need for chi calculations
    vector<double> thtgmain1, thtgmain2, thtgmain1err, thtgmain2err;
   
 
    double Xi[20],Xierr[20],QsqAvg[20],QsqAvgerr[20],QsqDiff[20],QsqDifferr[20];
    double AsymAvg[20],AsymAvgerr[20],AsymDiff[20],AsymDifferr[20];
    double ThtgAvg[20],ThtgAvgerr[20],ThtgDiff[20],ThtgDifferr[20];
    double PhtgAvg[20],PhtgAvgerr[20],PhtgDiff[20],PhtgDifferr[20];

    //For Chi calculations -- average for four A_T det combo
    double ThtgmainL[20],ThtgmainR[20], ThtgmainAv[20];
  
    double Chi[20], Chierr[20];

    TGraph *trans(); //this for the intergrating transverse data

    //only needed for A_T detectors
    void ReadFiles(string);
    void ReadFiles(string,string);


};

#endif
