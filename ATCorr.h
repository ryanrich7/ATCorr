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
    ATCorr(string, string, string, string); //for 4 detector combo correction
    ~ATCorr();
    
    void PlotMainDetPlots();
    void PlotATDetPlots(arm Arm = kLHRS);
    void PlotATDetPlots(int );//4 detector combo
    int Init();
    int Init(int );// 4 detector combo 
    //Making public right now

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
    void ComputeVerAvg();
    void ComputeVerDiff();
    void ComputeHorAvg();
    void ComputeHorDiff();
    void ComputeChi(); //up down acceptance factor -- horizontal A_T correction per arm
    void ComputeChi(int);//up down acceptance factor -- correct the average
    void ComputeAsym();//for all 4 detectors   
    void ComputeTh();//for all 4 detectors
    void ComputeXi(int);//calculation of left/right apparatus asymmetry using sin(phi)
    void ComputeCosAvg();
    void ComputeSinAvg();
    void ComputeCosDiff();
    void ComputeSinDiff();
    void ComputeXi(int,int);//calculation of left/right apparatus asymmetry using sin(phi)*sin(theta)
    void ComputeChi(int,int);
 
 
    TMultiGraph *Xifac, *Chifac; 

    //Other Kinematic Quantities
    TMultiGraph *Q2_L, *Q2_R, *Q2Av, *Q2DD, *Th_L, *Th_R, *ThAv, *ThDD, *Ph_L, *Ph_R, *PhAv, *PhDD, *APV_L, *APV_R, *APVAvg, *APVDD, *Lab_L, *Lab_R;
    TMultiGraph *Ver_L, *Ver_R, *VerAv, *VerDD, *Hor_L, *Hor_R, *HorAv, *HorDD; 
    TMultiGraph *Cos_L, *Cos_R, *Sin_L, *Sin_R, *CosAv, *CosDD, *SinAv, *SinDD;

    //detectors
    vector<double> QsqL, AngleL, asymL, ThtgL, PhtgL;
    vector<double> QsqR, AngleR, asymR, ThtgR, PhtgR;
    vector<double> QsqerL, AngleerL, asymerL, ThtgerL, PhtgerL;
    vector<double> QsqerR, AngleerR, asymerR, ThtgerR, PhtgerR;
    vector<double> cosphL,cosphR, cospherL, cospherR;
    vector<double> sinphL, sinphR, sinpherL, sinpherR;
    vector<double> VermomL, VermomR,HormomL,HormomR;
    vector<double> VermomerL, VermomerR, HormomerL, HormomerR;
   
    vector<int> RunL,RunR;
    //Need for chi calculations 
    vector<double> thtgmain1, thtgmain2, thtgmain1err, thtgmain2err, cosmain1, cosmain2, cosmain1err, cosmain2err;
    //Need for 4 detector combo calculations
    vector<double> QsqL1,AngleL1,asymL1,ThtgL1,PhtgL1;
    vector<double> QsqR1, AngleR1,asymR1,ThtgR1,PhtgR1;  
    vector<double> QsqerL1,AngleerL1,asymerL1,ThtgerL1,PhtgerL1;
    vector<double> QsqerR1, AngleerR1,asymerR1,ThtgerR1,PhtgerR1;
    vector<double> cosphL1,cosphR1, cospherL1, cospherR1;
    vector<double> sinphL1, sinphR1, sinpherL1, sinpherR1;
    vector<double> VermomL1, VermomR1,HormomL1,HormomR1;
    vector<double> VermomerL1, VermomerR1, HormomerL1, HormomerR1;
    vector<int> RunL1, RunR1;


 
    double Xi[20],Xierr[20],QsqAvg[20],QsqAvgerr[20],QsqDiff[20],QsqDifferr[20];
    double AsymAvg[20],AsymAvgerr[20],AsymDiff[20],AsymDifferr[20];
    double ThtgAvg[20],ThtgAvgerr[20],ThtgDiff[20],ThtgDifferr[20];
    double PhtgAvg[20],PhtgAvgerr[20],PhtgDiff[20],PhtgDifferr[20];
    double VerMomAvg[20], VerMomAvgerr[20], VerMomDiff[20], VerMomDifferr[20];
    double HorMomAvg[20], HorMomAvgerr[20], HorMomDiff[20], HorMomDifferr[20];  
    double CosAvg[20], CosAvgerr[20], CosDiff[20],CosDifferr[20];
    double SinAvg[20], SinAvgerr[20], SinDiff[20], SinDifferr[20];

    //For Chi calculations -- average for four A_T det combo
    double ThtgmainL[20],ThtgmainR[20], ThtgmainAv[20];
    double AsymAll[20], AsymAllerr[20], ThAll[20], ThAllerr[20], CosAll[20], CosAllerr[20];
    double Chi[20], Chierr[20];

    TGraph *trans(); //this for the intergrating transverse data

    //only needed for A_T detectors
    void ReadFiles(string);
    void ReadFiles(string,string);


};

#endif
