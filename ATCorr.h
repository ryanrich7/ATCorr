#ifndef ATCorr_H
#define ATCorr_H

#include <vector>
#include <string>
#include "TGraph.h"

 class TGraph;

 class ATCorr  {
    public:

    ATCorr(string, string); //reads text files for A_T correction. Form Vertical, LHRS and RHRS main detector files, for A_T correction use files on same arm
    ~ATCorr();

    
    void PlotMainDetPlots();
    void PlotATDetPlots();
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
    void ComputeChi(); //up down acceptance factor -- horizontal A_T correction

    
   //Additional function needed for the Chi computation
 //   void Load(); 

    //Correction factors
    TMultiGraph *Xifac, *Chifac; 

    //Other Kinematic Quantities
    TMultiGraph *Q2_L, *Q2_R, *Q2Av, *Q2DD, *Th_L, *Th_R, *ThAv, *ThDD, *Ph_L, *Ph_R, *PhAv, *PhDD, *APV_L, *APV_R, *APVAvg, *APVDD, *Lab_L, *Lab_R;
   

    //detectors
    vector<double> QsqL, AngleL, asymL, ThtgL, PhtgL;
    vector<double> QsqR, AngleR, asymR, ThtgR, PhtgR;
    vector<double> QsqerL, AngleerL, asymerL, ThtgerL, PhtgerL;
    vector<double> QsqerR, AngleerR, asymerR, ThtgerR, PhtgerR;
    vector<int> RunL,RunR;

    double Xi[20],Xierr[20],QsqAvg[20],QsqAvgerr[20],QsqDiff[20],QsqDifferr[20];
    double AsymAvg[20],AsymAvgerr[20],AsymDiff[20],AsymDifferr[20];
    double ThtgAvg[20],ThtgAvgerr[20],ThtgDiff[20],ThtgDifferr[20];
    double PhtgAvg[20],PhtgAvgerr[20],PhtgDiff[20],PhtgDifferr[20];

    double Chi[20], Chierr[20];

    TGraph *trans(); //this for the intergrating transverse data




};

#endif
