#include        <iostream>
#include        <fstream>
#include        <iomanip>
#include        <cmath>
#include        <math.h>
#include        "Library/TNT/tnt.h"
#include        "Library/MyFun/myFun.h"
using namespace std;
using namespace TNT;

extern "C" {

class arrayCGH2d {

// NANCY: move obs out as public.

  public:
    double      p, a, b, c, base1, base2, oldtarget, emtarget;
    Vector<int>     mixState;
    Vector<double>  paraMu, estBS;
    Matrix<double>  obs, estPara, paraV, paraSigAA, paraSigAB, 
            paraSigBA, paraSigBB;
    arrayCGH2d(Matrix<double>, double, double, double, double, 
    double, Vector<double>, Matrix<double>, Matrix<double>, 
    Matrix<double>, Matrix<double>, Matrix<double>, int, int);
    ~arrayCGH2d(){};
    void        BcmixSmoothEst();
    void        FindState();
    void        SelectAllHyper(int paraSet);
    void        PrintAllHyper();

// NANCY: add mystream.

    ofstream mystream;  // file stream for debugging output.


  private:
    int         N, K, M;
    double      detV, iv11, iv12, iv22, ivmu1, ivmu2, 
            tmpLog00, aAA, bAA, dAA, aAB, bAB, dAB,
            aBA, bBA, dBA, aBB, bBB, dBB, detSigAA,
            detSigAB, detSigBA, detSigBB,
            y1, y2, pstar, sumvec;
    Vector<double>  fP, bP;
    Matrix<int>     fFilter, bFilter;
    Matrix<double>  cenObs, fQ, fLog, bQ, bLog, probseq;
    Array3D<double> fIV, fB, bIV, bB;

    void        GetCenteredIntensity();
    void        AddBackBaseline();      // NANCY: add this function.
    void        CompInternalPara();
    void        BcmixForward();
    void        BcmixBackward();
    void        BcmixSmooth();
    void        EM_allhyper_BcmixSM(Matrix<double>&,Matrix<double>&,
    Matrix<double>&, Matrix<double>&, Matrix<double>&);
    void                EM_allhyper_Update(Matrix<double>,Matrix<double>,
        Matrix<double>, Matrix<double>, Matrix<double>, int);
    void    CompSmoothWeight(double, double, double, double, double, 
    double, double, double, double, double, double&, double&, 
    double&, double&, double&, double&);
    void        UpdateFun(double&, double&, double&, double&, double&,
    double&, double, double, double, double, double, double, double,
    int);

    // sub functions 
    double      GetTarget(double, double, double, double, double, double);
    void        EM_Update_paraSig(Matrix<double>&, double, double,
    double, double);
    void    SaveEMOldPara(double&, double&, double&, Vector<double>&,
    Matrix<double>&, Matrix<double>&, Matrix<double>&, 
    Matrix<double>&, Matrix<double>&, double, double, double, 
    Vector<double>, Matrix<double>, Matrix<double>, Matrix<double>, 
    Matrix<double>, Matrix<double>);
    void    SubFunOneInEM(double&, double &, double&, double&, double&,
    double, double, double, double, double, double);
    void        Get2dSymMatInv(Matrix<double>, double&, double&,double&,
    double&);
    void        CompIVinvMultiB(double, double, double, double, double,
    double&, double&, double&, double&, double&);

}; /* end of -- class arrayCGH */

arrayCGH2d::arrayCGH2d(Matrix<double> obs1, double p1, double a1,
    double b1, double base11, double base21, Vector<double> paraMu1, 
    Matrix<double> paraV1, Matrix<double> paraSigAA1, 
    Matrix<double> paraSigAB1, Matrix<double> paraSigBA1,
    Matrix<double> paraSigBB1, int K1, int M1)
{
  K=K1;  M=M1;  p=p1;   a=a1;   b=b1;   c=1-a-b;    N=obs1.num_rows()-1;
  paraMu=paraMu1;       paraV=paraV1;
  paraSigAA=paraSigAA1;     paraSigAB=paraSigAB1;
  paraSigBA=paraSigBA1;     paraSigBB=paraSigBB1;
  fP.newsize(N+1);      bP.newsize(N+1);
  fFilter.newsize(N+1,K+1); bFilter.newsize(N+1,K+1);
  fQ.newsize(N+1,K+1);      fLog.newsize(N+1,K+1);
  bQ.newsize(N+1,K+1);      bLog.newsize(N+1,K+1);
  estPara.newsize(N+1,2);   mixState.newsize(N+1);
  estBS.newsize(N+1);       obs.newsize(N+1,2);
  cenObs.newsize(N+1,2);
  
  
//NANCY:  
//  for(int i=1;i<=N;i++) for(int j=0;j<2;j++)    obs[i][j]=obs1[i-1][j];
  for(int i=1;i<=N;i++) for(int j=0;j<2;j++)    obs[i][j]=obs1[i][j];
  
  Array3D<double>   nfIV(N+1,K+1,3), nfB(N+1,K+1,3), 
            nbIV(N+1,K+1,3), nbB(N+1,K+1,3);
  fIV=nfIV; fB=nfB;     bIV=nbIV;   bB=nbB;   fP=0.0;
  fLog=0.0; bLog=0.0;   fQ=0.0;     bQ=0.0;   bP=0.0;
  base1=base11;     base2=base21;
  emtarget = - exp(100);
}


/********************************************************************/

void   arrayCGH2d::BcmixSmoothEst()
{
  GetCenteredIntensity();
  CompInternalPara();
  BcmixForward();
  BcmixBackward();
  BcmixSmooth();
  // NANCY: Need to add baseline to estPara[0], so that subsequent call to FindState() would work.
  // Thus, estPara is the estimated parameters for obs, mixstate, and not for cenObs, mixstate.
  // "cenObs" exists only internally to BcmixSmoothEst, and outside of it (e.g. FindState()) estPara
  // should have baseline added in.
  AddBackBaseline();
}

void  arrayCGH2d::FindState()
{
  int   tmpmin;         double   alpha, beta, z1, z2;
  Vector<double>        target(4);
  for (int t=1; t<N+1; t++) {   y1=obs[t][0];   y2=obs[t][1];
    alpha=estPara[t][0];    beta=estPara[t][1];    tmpmin=0;
    z1=y1-2.0*alpha;    z2=y2;
    target[0]=GetTarget(aAA, bAA, dAA, z1, z2, detSigAA);
    z1=y1-(alpha+beta); z2=y2-(alpha-beta);
    target[1]=GetTarget(aAB, bAB, dAB, z1, z2, detSigAB);
    z1=y1-(alpha-beta); z2=y2-(alpha+beta);
    target[2]=GetTarget(aBA, bBA, dBA, z1, z2, detSigBA);
    z1=y1;      z2=y2-2.0*alpha;
    target[3]=GetTarget(aBB, bBB, dBB, z1, z2, detSigBB);
    for (int i=0; i<4; i++)  if (target[i]<target[tmpmin])  tmpmin=i;
    mixState[t] = tmpmin+1;
  }
}


void    arrayCGH2d::SelectAllHyper(int paraSet)
{
  int           iter=1, emIterMax=1;
  double        emOldp, emOlda, emOldb, oldtarget;
  Vector<double>    emOldMu, emBS(N+1);
  Matrix<double>    emOldV,emOldSigAA,emOldSigAB,emOldSigBA,emOldSigBB;
  Matrix<double>        emSig(N+1,2), emSig2(N+1,3), emPSig(N+1,2),
            emPSig2(N+1,3), emSigAA(N+1,3), emSigAB(N+1,3),
            emSigBA(N+1,3), emSigBB(N+1,3), probseq(N+1,6);

  GetCenteredIntensity();
  SaveEMOldPara(emOldp, emOlda, emOldb, emOldMu, emOldV, emOldSigAA, 
    emOldSigAB, emOldSigBA, emOldSigBB, p, a, b, paraMu, paraV, 
    paraSigAA, paraSigAB, paraSigBA, paraSigBB);
  EM_allhyper_BcmixSM(emSig,emSig2,emPSig,emPSig2,probseq);
  EM_allhyper_Update(emSig,emSig2,emPSig,emPSig2,probseq, paraSet);
  oldtarget = emtarget-1;
  //  while ((oldtarget<emtarget)&(iter<emIterMax)) {
  while (iter<emIterMax) {
    SaveEMOldPara(emOldp, emOlda, emOldb, emOldMu, emOldV, emOldSigAA,
      emOldSigAB, emOldSigBA, emOldSigBB, p, a, b, paraMu, paraV,
      paraSigAA, paraSigAB, paraSigBA, paraSigBB);
    EM_allhyper_BcmixSM(emSig,emSig2,emPSig,emPSig2,probseq);
    EM_allhyper_Update(emSig,emSig2,emPSig,emPSig2,probseq, paraSet);
    iter++; oldtarget = emtarget;
  }
}

void    arrayCGH2d::PrintAllHyper()
{
  cout << "p=" << p << "\t" << "a="<< a << "\t" << "b=" << b <<endl;
  cout << "paraMu = " << paraMu[0] << " " << paraMu[1] << endl;
  cout << "paraV = " << paraV << endl;
  cout << paraSigAA << "\n" << paraSigAB << "\n" << paraSigBA <<
      "\n" << paraSigBB << endl;
}


/******************* private functions *******************************/

// estPara = E(theta_t|Data), estBS = Prob(theta_t=B|Data)
// probseq[t][1] = Prob( theta_t = theta_{t+1} = B |Data )
// probseq[t][2] = Prob( theta_t = B, theta_{t+1} != B |Data)
// probseq[t][3] = Prob( theta_t !=B, theta_{t+1} = B |Data)
// probseq[t][4] = Prob( B != theta_t != theta_{t+1} != B |Data)
// probseq[t][5] = Prob( theta_t = theta_{t+1} != B |Data)
// emSig = E(theta_t|Data),   emSig2=E(theta_t^2|Data)
// emPSig = E(theta_t * 1_{B != theta_t != theta_{t-1}|Data)
// emPSig2 = E( (theta_t-mu)^2 * 1_{B != theta_t != theta_{t-1}|Data)

/* The function can only be called by "SelectAllHyper" */
void   arrayCGH2d::EM_allhyper_BcmixSM(Matrix<double> &emSig,
        Matrix<double> &emSig2, Matrix<double> &emPSig,
        Matrix<double> &emPSig2, Matrix<double> &probseq)
{
  CompInternalPara();     BcmixForward();     BcmixBackward();
  double    alpha, tmpa, total,
        tmp, vb1, vb2, tmpv11, tmpv12, tmpv22, logij, incre;

  emSig=0.0;   emSig2=0.0;   emPSig=0.0;  emPSig2=0.0;  probseq=0.0;
  for (int t=1; t<N; t++) {
    probseq[t][1]=fP[t]*bP[t+1]*(1-p)/c;  probseq[t][2]=fP[t]*(1-bP[t+1]);
    probseq[t][3]=(1-fP[t])*bP[t+1];
    probseq[t][4]=(1-fP[t])*(1-bP[t+1])*b/p;
    alpha=probseq[t][1]+probseq[t][2];   tmpa=(b+(p-b)*bP[t+1])/p;
    for (int i=1;i<=min(t,K);i++) {
      CompIVinvMultiB(fIV[t][i][0], fIV[t][i][1], fIV[t][i][2],
    fB[t][i][0], fB[t][i][1], tmpv11, tmpv12, tmpv22, vb1, vb2);
      SubFunOneInEM(emSig[t][0],emSig[t][1], emSig2[t][0],emSig2[t][1],
    emSig2[t][2], fQ[t][i]*tmpa, vb1,vb2, tmpv11,tmpv12,tmpv22);
    }
    for (int i=1;i<=min(t,K);i++)  for (int j=1;j<=min(N-t,K);j++) {
      CompSmoothWeight(fIV[t][i][0],fIV[t][i][1],fIV[t][i][2],
    bIV[t+1][j][0],bIV[t+1][j][1],bIV[t+1][j][2], fB[t][i][0],
    fB[t][i][1], bB[t+1][j][0],bB[t+1][j][1], tmpv11, tmpv12,
    tmpv22, vb1, vb2, logij);
      probseq[t][5] +=incre =a/p*fQ[t][i]*bQ[t+1][j]*exp(fLog[t][i]+
    bLog[t+1][j]-logij-tmpLog00);
      SubFunOneInEM(emSig[t][0],emSig[t][1], emSig2[t][0],emSig2[t][1],
    emSig2[t][2], incre, vb1, vb2, tmpv11, tmpv12, tmpv22);
    }
    total=alpha+probseq[t][3]+probseq[t][4]+probseq[t][5];
    for (int i=1; i<6; i++) probseq[t][i] /=total;
    tmpa = (1+(b/p-1)*(1-fP[t]))/total;
    for (int j=1;j<=min(N-t,K);j++) {   tmp=tmpa*bQ[t+1][j];
      CompIVinvMultiB(bIV[t+1][j][0],bIV[t+1][j][1],bIV[t+1][j][2],
    bB[t+1][j][0],bB[t+1][j][1], tmpv11,tmpv12,tmpv22, vb1,vb2);
      emPSig[t][0] += tmp*vb1;          emPSig[t][1] += tmp*vb2;
      emPSig2[t][0] +=tmp*((vb1-paraMu[0])*(vb1-paraMu[0])+tmpv11);
      emPSig2[t][1] +=tmp*((vb1-paraMu[0])*(vb2-paraMu[1])+tmpv12);
      emPSig2[t][2] +=tmp*((vb2-paraMu[1])*(vb2-paraMu[1])+tmpv22);
    }
    estBS[t]=alpha/total;
    for (int k=0;k<2;k++) { emSig[t][k]/=total; emPSig[t][k]/=total; }
    for (int k=0;k<3;k++) { emSig2[t][k]/=total; emPSig2[t][k]/=total; }
  }
//cout << probseq << endl;
}

/* The function can only be called by "SelectAllHyper" */
void    arrayCGH2d::EM_allhyper_Update(Matrix<double> emSig,
        Matrix<double> emSig2, Matrix<double> emPSig,
        Matrix<double> emPSig2, Matrix<double> probseq,
        int paraSet)
{
  double    inAAnu11, inAAnu12, inAAnu22, inABnu11, inABnu12,
        inABnu22, inBAnu11, inBAnu12, inBAnu22, inBBnu11,
        inBBnu12, inBBnu22, AAde, ABde, BAde, BBde, pnu, pde,
        vde, bnu, bde, cnu, mu1nu, mu2nu, v11nu, v12nu, v22nu;
  inAAnu11=inAAnu12=inAAnu22=inABnu11=inABnu12=inABnu22=inBAnu11=
    inBAnu12=inBAnu22=inBBnu11=inBBnu12=inBBnu22=AAde=ABde=BAde=
    BBde=pnu=pde=vde=bnu=bde=cnu=mu1nu=mu2nu=v11nu=v12nu=v22nu=0.0;
  for (int t=1; t<N; t++) {
    pnu +=probseq[t][2];  bnu +=probseq[t][4];  cnu +=probseq[t][3];
    pde +=probseq[t][1]+probseq[t][2];
    bde +=probseq[t][3]+probseq[t][4]+probseq[t][5];
    vde +=probseq[t][2]+probseq[t][4];
    mu1nu +=emPSig[t][0];      mu2nu +=emPSig[t][1];
    v11nu +=emPSig2[t][0];  v12nu +=emPSig2[t][1];  v22nu +=emPSig2[t][2];
    y1=cenObs[t][0];        y2=cenObs[t][1];
    switch ( mixState[t] ) {
      case 1:   AAde += 1.0;        inAAnu22 += y2*y2;
    inAAnu11 += y1*y1-4.0*y1*emSig[t][0] + 4.0*emSig2[t][0];
    inAAnu12 += y1*y2-2.0*y2*emSig[t][0];           break;
      case 2:   ABde += 1.0;
    inABnu11 += y1*y1-2.0*y1*(emSig[t][0]+emSig[t][1]) +
      emSig2[t][0] + 2.0*emSig2[t][1] + emSig2[t][2];
    inABnu12 += y1*y2-(emSig[t][0]-emSig[t][1])*y1
      -(emSig[t][0]+emSig[t][1])*y2 + emSig2[t][0] - emSig2[t][2];
    inABnu22 += y2*y2-2.0*y2*(emSig[t][0]-emSig[t][1]) +
      emSig2[t][0]-2.0*emSig2[t][1]+emSig2[t][2];   break;
      case 3:   BAde += 1.0;
    inBAnu11 += y1*y1-2.0*y1*(emSig[t][0]-emSig[t][1]) +
      emSig2[t][0] - 2.0*emSig2[t][1] + emSig2[t][2];
    inBAnu12 += y1*y2-(emSig[t][0]-emSig[t][1])*y2
      -(emSig[t][0]+emSig[t][1])*y1 + emSig2[t][0] - emSig2[t][2];
    inBAnu22 += y2*y2-2.0*y2*(emSig[t][0]+emSig[t][1]) +
      emSig2[t][0]+2.0*emSig2[t][1]+emSig2[t][2];   break;
      case 4:   BBde += 1.0;        inBBnu11 += y1*y1;
        inBBnu22 += y2*y2-4.0*y2*emSig[t][0] + 4.0*emSig2[t][0];
        inBBnu12 += y1*y2-2.0*y1*emSig[t][0];           break;
    }
  }
  switch( paraSet ) {
    case 1: // update the transition matrix
      p=pnu/pde;    b=bnu/bde;   c=cnu/bde;   a=1-b-c;   break;
    case 2: // update other parameters
      paraMu[0]=mu1nu/vde;      paraMu[1]=mu2nu/vde;
      EM_Update_paraSig(paraV, v11nu, v12nu, v22nu, vde);
      EM_Update_paraSig(paraSigAA,inAAnu11,inAAnu12,inAAnu22,AAde);
      EM_Update_paraSig(paraSigAB,inABnu11,inABnu12,inABnu22,ABde);
      EM_Update_paraSig(paraSigBA,inBAnu11,inBAnu12,inBAnu22,BAde);
      EM_Update_paraSig(paraSigBB,inBBnu11,inBBnu12,inBBnu22,BBde);
      break;
    case 3: // update all hypers
      p=pnu/pde;    b=bnu/bde;   c=cnu/bde;    a=1-b-c;
      paraMu[0]=mu1nu/vde;  paraMu[1]=mu2nu/vde;
      EM_Update_paraSig(paraV, v11nu, v12nu, v22nu, vde);
      EM_Update_paraSig(paraSigAA, inAAnu11, inAAnu12, inAAnu22, AAde);
      EM_Update_paraSig(paraSigAB, inABnu11, inABnu12, inABnu22, ABde);
      EM_Update_paraSig(paraSigBA, inBAnu11, inBAnu12, inBAnu22, BAde);
      EM_Update_paraSig(paraSigBB, inBBnu11, inBBnu12, inBBnu22, BBde);
      break;
    case 4: // NANCY: clamp everything except paraSig**.
      EM_Update_paraSig(paraSigAA, inAAnu11, inAAnu12, inAAnu22, AAde);
      EM_Update_paraSig(paraSigAB, inABnu11, inABnu12, inABnu22, ABde);
      EM_Update_paraSig(paraSigBA, inBAnu11, inBAnu12, inBAnu22, BAde);
      EM_Update_paraSig(paraSigBB, inBBnu11, inBBnu12, inBBnu22, BBde);
      break;
      
   }
   CompInternalPara();
}




//  The function computes the smooth estimate given the forward and 
// backward results.
void   arrayCGH2d::BcmixSmooth()
{
  double     alpha, tmpa,
    tmp, tmpv11, tmpv12, tmpv22, vb1, vb2, logij, incre, total;
  estPara=0.0;
  for (int t=1; t<N+1; t++) {
    total = alpha = fP[t]*(c+(1-p-c)*bP[t+1])/c;
    // Note that we did not add the baseline to the estimate.
    tmpa = (b+(p-b)*bP[t+1])/p;     total+=(1-fP[t])*tmpa;
    for (int i=1;i<=min(t,K);i++) { tmp = tmpa*fQ[t][i];
      CompIVinvMultiB(fIV[t][i][0], fIV[t][i][1], fIV[t][i][2],
    fB[t][i][0], fB[t][i][1], tmpv11, tmpv12, tmpv22, vb1, vb2);
      estPara[t][0]+= tmp*vb1;  estPara[t][1]+= tmp*vb2;
    }
    for (int i=1;i<=min(t,K);i++)  for (int j=1;j<=min(N-t,K);j++) {
      CompSmoothWeight(fIV[t][i][0], fIV[t][i][1], fIV[t][i][2],
    bIV[t+1][j][0], bIV[t+1][j][1], bIV[t+1][j][2], fB[t][i][0],
    fB[t][i][1], bB[t+1][j][0], bB[t+1][j][1], tmpv11, tmpv12,
    tmpv22, vb1, vb2, logij);
      total += incre = a/p*fQ[t][i]*bQ[t+1][j] * exp( fLog[t][i]
    + bLog[t+1][j]-logij-tmpLog00);
      estPara[t][0] +=incre*vb1;     estPara[t][1] +=incre*vb2;
    }
    estPara[t][0]/=total;   estPara[t][1]/=total;
    estBS[t]=alpha/total;
  }
  for (int i=1;i<=K;i++) {  tmp = fQ[N][i];
    CompIVinvMultiB(fIV[N][i][0], fIV[N][i][1], fIV[N][i][2], 
      fB[N][i][0], fB[N][i][1], tmpv11, tmpv12, tmpv22, vb1, vb2);
    estPara[N][0] += tmp*vb1;   estPara[N][1] += tmp*vb2;
  }
  estBS[N]=fP[N];
}

void   arrayCGH2d::BcmixBackward()
{
  Vector<int>       tmpBackFil(K+2);    int s;
  Vector<double>    bQstar(K+2), tmpLog(K+2), numTmp(K+2);
  Matrix<double>    tmpIV(K+2,3), tmpB(K+2,2);

  bFilter[N][1]=1;  bQ[N][1]=p/(p+c);   bP[N]=1-bQ[N][1];
  UpdateFun(bIV[N][1][0],bIV[N][1][1],bIV[N][1][2], bB[N][1][0],
    bB[N][1][1], bLog[N][1], iv11, iv12, iv22, ivmu1, ivmu2, 
    cenObs[N][0],cenObs[N][1], mixState[N]);
  for (int t=N-1; t>=N-K+1; t--) { y1=cenObs[t][0]; y2=cenObs[t][1]; s=N+1-t;
    UpdateFun(bIV[t][s][0],bIV[t][s][1],bIV[t][s][2], bB[t][s][0],
      bB[t][s][1], bLog[t][s], iv11, iv12, iv22, ivmu1, ivmu2, y1, y2, 
      mixState[t]);
    numTmp[s] = tmpLog00-bLog[t][s];        bFilter[t][s]=t;
    for (int i=1; i<s; i++) {
      UpdateFun(bIV[t][i][0],bIV[t][i][1],bIV[t][i][2], bB[t][i][0],
    bB[t][i][1], bLog[t][i],bIV[t+1][i][0],bIV[t+1][i][1],bIV[t+1][i][2],
    bB[t+1][i][0], bB[t+1][i][1], y1, y2, mixState[t]);
      numTmp[i] = bLog[t+1][i]-bLog[t][i];  bFilter[t][i]=N+1-i;
    }
    double      numTmpMax=numTmp[1];
    for (int i=2;i<=s;i++) if (numTmpMax<numTmp[i]) numTmpMax=numTmp[i];
    for (int i=1;i<=s;i++) numTmp[i] -= numTmpMax;
    sumvec=pstar=((1-p-c)*bP[t+1]+c)*exp(-numTmpMax);
    sumvec+=bQstar[s]=(b+(p-b)*bP[t+1])*exp(numTmp[s]);
    for (int i=1;i<s;i++)  sumvec+=bQstar[i]=a*bQ[t+1][i]*exp(numTmp[i]);
    bP[t] = pstar/sumvec;
    for (int i=1; i<=s; i++)    bQ[t][i]=bQstar[i]/sumvec;
  }

  for (int t=N-K; t>=1; t--) { y1=cenObs[t][0]; y2=cenObs[t][1];
    UpdateFun(tmpIV[K+1][0],tmpIV[K+1][1],tmpIV[K+1][2], tmpB[K+1][0],
      tmpB[K+1][1], tmpLog[K+1], iv11, iv12, iv22, ivmu1, ivmu2, 
      y1, y2, mixState[t]);
    numTmp[K+1]=tmpLog00-tmpLog[K+1];
    for (int i=1; i<K+1; i++) {
      UpdateFun(tmpIV[i][0],tmpIV[i][1],tmpIV[i][2],tmpB[i][0],tmpB[i][1],
    tmpLog[i], bIV[t+1][i][0],bIV[t+1][i][1],bIV[t+1][i][2],
    bB[t+1][i][0],bB[t+1][i][1], y1, y2, mixState[t]);
      numTmp[i] = bLog[t+1][i]-tmpLog[i];  tmpBackFil[i]=bFilter[t+1][i];
    }
    double      numTmpMax=numTmp[1];
    for (int i=1; i<=K+1; i++) if (numTmpMax<numTmp[i]) numTmpMax=numTmp[i];
    for (int i=1; i<=K+1; i++) numTmp[i] -= numTmpMax;
    sumvec = pstar = ((1-p-c)*bP[t+1]+c)*exp(-numTmpMax);
    sumvec += bQstar[K+1]=(b+(p-b)*bP[t+1])*exp(numTmp[K+1]);
    for (int i=1;i<K+1;i++)   sumvec+=bQstar[i]=a*bQ[t+1][i]*exp(numTmp[i]);
    bP[t] = pstar/sumvec;   int  minind = 1;
    for (int i=1;i<=K+1-M;i++) if (bQstar[minind]>bQstar[i]) minind=i;
    double delta=(1-bP[t])/(sumvec-pstar-bQstar[minind]);
    for (int i=1; i<minind; i++) {
      for (int j=0; j<3; j++)   bIV[t][i][j] = tmpIV[i][j];
      for (int j=0; j<2; j++)   bB[t][i][j] = tmpB[i][j];
      bLog[t][i]=tmpLog[i];     bFilter[t][i]=tmpBackFil[i]; 
      bQ[t][i]=bQstar[i]*delta;
    }
    for (int i=minind+1; i<=K+1; i++) {
      for (int j=0; j<3; j++)   bIV[t][i-1][j] = tmpIV[i][j];
      for (int j=0; j<2; j++)   bB[t][i-1][j] = tmpB[i][j];
      bLog[t][i-1]=tmpLog[i];   bFilter[t][i-1]=tmpBackFil[i]; 
      bQ[t][i-1]=bQstar[i]*delta;
    }
  }
}


void   arrayCGH2d::BcmixForward()
{
  Vector<int>       tmpForFil(K+2);
  Vector<double>    fQstar(K+2), tmpLog(K+2), numTmp(K+2);
  Matrix<double>    tmpIV(K+2,3), tmpB(K+2,2);
  fFilter[1][1]=1;  fQ[1][1]=p/(p+c);   fP[1]=1-fQ[1][1];
  UpdateFun(fIV[1][1][0],fIV[1][1][1],fIV[1][1][2], fB[1][1][0],fB[1][1][1],
    fLog[1][1],iv11,iv12,iv22,ivmu1,ivmu2, cenObs[1][0],cenObs[1][1],
    mixState[1]);
  for (int t=2; t<=K; t++) {    y1=cenObs[t][0];     y2=cenObs[t][1];
    UpdateFun(fIV[t][t][0],fIV[t][t][1],fIV[t][t][2], fB[t][t][0],fB[t][t][1],
      fLog[t][t], iv11, iv12, iv22, ivmu1, ivmu2, y1, y2, mixState[t]);
    numTmp[t]=tmpLog00-fLog[t][t];      fFilter[t][t]=t;
    for (int i=1; i<t; i++) {
      UpdateFun(fIV[t][i][0],fIV[t][i][1],fIV[t][i][2], fB[t][i][0],
    fB[t][i][1], fLog[t][i], fIV[t-1][i][0], fIV[t-1][i][1],
    fIV[t-1][i][2], fB[t-1][i][0], fB[t-1][i][1], y1, y2, mixState[t]);
      numTmp[i] = fLog[t-1][i]-fLog[t][i];  fFilter[t][i]=i;
    } // Note that leaving fQ[t-1][i] in the end can void NaNs.
    double      numTmpMax=numTmp[1];
    for (int i=2;i<=t;i++) if (numTmpMax<numTmp[i]) numTmpMax=numTmp[i];
    for (int i=1;i<=t;i++) numTmp[i] -= numTmpMax;
    sumvec = pstar = ((1-p-c)*fP[t-1]+c)*exp(-numTmpMax);
    sumvec+=fQstar[t]=(b+(p-b)*fP[t-1])*exp( numTmp[t] );
    for (int i=1;i<t;i++)  sumvec+=fQstar[i]=a*fQ[t-1][i]*exp(numTmp[i]);
    fP[t] = pstar/sumvec;
    for (int i=1; i<=t; i++)    fQ[t][i]=fQstar[i]/sumvec;
  }

  for (int t=K+1; t<N+1; t++) { y1=cenObs[t][0];    y2=cenObs[t][1];
    UpdateFun(tmpIV[K+1][0],tmpIV[K+1][1],tmpIV[K+1][2], tmpB[K+1][0],
      tmpB[K+1][1], tmpLog[K+1], iv11,iv12,iv22,ivmu1,ivmu2, y1, y2, 
      mixState[t]);
    numTmp[K+1]=tmpLog00-tmpLog[K+1];
    for (int i=1; i<K+1; i++) {
      UpdateFun(tmpIV[i][0],tmpIV[i][1],tmpIV[i][2], tmpB[i][0],tmpB[i][1],
    tmpLog[i], fIV[t-1][i][0],fIV[t-1][i][1],fIV[t-1][i][2],
    fB[t-1][i][0],fB[t-1][i][1],y1, y2, mixState[t]);
      numTmp[i]=fLog[t-1][i]-tmpLog[i];  tmpForFil[i]=fFilter[t-1][i];
    }
    double  numTmpMax=numTmp[1];
    for (int i=2;i<=K+1;i++) if (numTmpMax<numTmp[i]) numTmpMax=numTmp[i];
    for (int i=1; i<=K+1; i++)  numTmp[i] -= numTmpMax;
    sumvec = pstar = ((1-p-c)*fP[t-1]+c)*exp(-numTmpMax);
    sumvec += fQstar[K+1] = (b+(p-b)*fP[t-1])*exp( numTmp[K+1] );
    for (int i=1;i<K+1;i++)  sumvec+=fQstar[i]=a*fQ[t-1][i]*exp(numTmp[i]);
    fP[t] = pstar/sumvec;   int  minind = 1;
    for (int i=1;i<=K+1-M;i++) if (fQstar[minind]>fQstar[i]) minind=i;
    double  delta=(1-fP[t])/(sumvec-pstar-fQstar[minind]);
    for (int i=1; i<minind; i++) {
      for (int j=0; j<3; j++)   fIV[t][i][j] = tmpIV[i][j];
      for (int j=0; j<2; j++)   fB[t][i][j] = tmpB[i][j];
      fLog[t][i]=tmpLog[i]; fFilter[t][i]=tmpForFil[i];  
      fQ[t][i]=fQstar[i]*delta;
    }
    for (int i=minind+1; i<=K+1; i++) {
      for (int j=0; j<3; j++)   fIV[t][i-1][j] = tmpIV[i][j];
      for (int j=0; j<2; j++)   fB[t][i-1][j] = tmpB[i][j];
      fLog[t][i-1]=tmpLog[i];   fFilter[t][i-1]=tmpForFil[i];
      fQ[t][i-1]=fQstar[i]*delta;
    }
  }
}


void    arrayCGH2d::CompInternalPara()
{
  Get2dSymMatInv(paraV, detV, iv11, iv12, iv22);
  Get2dSymMatInv(paraSigAA, detSigAA, aAA, bAA, dAA);
  Get2dSymMatInv(paraSigAB, detSigAB, aAB, bAB, dAB);
  Get2dSymMatInv(paraSigBA, detSigBA, aBA, bBA, dBA);
  Get2dSymMatInv(paraSigBB, detSigBB, aBB, bBB, dBB);
  ivmu1 = iv22*paraMu[0]-iv12*paraMu[1];
  ivmu2 = iv11*paraMu[1]-iv12*paraMu[0];
  tmpLog00 = -0.5*( paraMu[0]*ivmu1+paraMu[1]*ivmu2+log(detV) );
}

void    arrayCGH2d::GetCenteredIntensity()
{
  for (int t=1; t<=N; t++)  switch(mixState[t]) {
    case 1: cenObs[t][0]=obs[t][0]-2.0*base1;
    cenObs[t][1]=obs[t][1];         break;
    case 2: cenObs[t][0]=obs[t][0]-(base1+base2);
    cenObs[t][1]=obs[t][1]-(base1-base2);   break;
    case 3: cenObs[t][0]=obs[t][0]-(base1-base2);
    cenObs[t][1]=obs[t][1]-(base1+base2);   break;
    case 4: cenObs[t][0]=obs[t][0];
    cenObs[t][1]=obs[t][1]-2.0*base1;   break;
  }
}



//NANCY:
void   arrayCGH2d::AddBackBaseline()
{
  for (int t=1; t<N+1; t++) {
    estPara[t][0] = estPara[t][0]+base1;
    estPara[t][1] = estPara[t][1]+base2;
  }
}

void    arrayCGH2d::CompSmoothWeight(double fiv11, double fiv12,
        double fiv22, double biv11, double biv12, double biv22,
        double fb1, double fb2, double bb1, double bb2,
        double &tmpv11, double &tmpv12, double &tmpv22,
        double &vb1, double &vb2, double &logij)
{
  double    tmpiv11, tmpiv12, tmpiv22, bij1, bij2, tmp;
  tmpiv11=fiv11+biv11-iv11;     tmpiv12=fiv12+biv12-iv12;
  tmpiv22=fiv22+biv22-iv22;
  bij1=fb1+bb1-ivmu1;       bij2=fb2+bb2-ivmu2;
  tmp = 1.0/(tmpiv11*tmpiv22-tmpiv12*tmpiv12);
  tmpv11 = tmp*tmpiv22;     tmpv12 = -tmp*tmpiv12;
  tmpv22 = tmp*tmpiv11;
  vb1 = tmpv11*bij1+tmpv12*bij2;    vb2 = tmpv12*bij1+tmpv22*bij2;
  logij = -0.5*(bij1*vb1+bij2*vb2+log(tmp));
}


/*  the function implement the following updating schemes: given V, 
*   mu, compute     newV, newMu, newB, newLogDen
    iv = V^{-1} = ( iv11, iv12;  iv12, iv22);
    invSigma = (a, b; b, d);
    ivmu = V^{-1}Mu = (ivmu1, ivmu2)
    M = newV^{-1} = iv + At'*invSigma*At
    B = ivmu + At'invSigma*yt
    logden = 0.5*(log|M| - B*M^{-1}B)
    increlogden = Yt'*Sig^{-1}*At*Base-Base'*At'*Sig^{-1}*At*Base/2
*/
void   arrayCGH2d::UpdateFun(double &M11, double &M12, double &M22,
        double &B1, double &B2, double &logden, double iv11,
        double iv12, double iv22, double ivmu1, double ivmu2,
        double y1, double y2, int mState)
{
  switch (mState) {
    case 1:      M11=iv11+4.0*aAA;      M12=iv12;       M22=iv22;
      B1=ivmu1+2.0*(aAA*y1+bAA*y2);     B2=ivmu2;
    break;    case 2:
      M11=iv11+aAB+dAB+bAB+bAB;     M12=iv12+aAB-dAB;
      M22=iv22+aAB+dAB-bAB-bAB;
      B1=ivmu1+(aAB+bAB)*y1+(bAB+dAB)*y2;
      B2=ivmu2+(aAB-bAB)*y1+(bAB-dAB)*y2;
    break;    case 3:
      M11=iv11+aBA+dBA+bBA+bBA; M12=iv12+dBA-aBA;
      M22=iv22+aBA+dBA-bBA-bBA;
      B1=ivmu1+(aBA+bBA)*y1+(bBA+dBA)*y2;
      B2=ivmu2+(bBA-aBA)*y1+(dBA-bBA)*y2;
    break;    case 4:
      M11=iv11+4.0*dBB;  M12=iv12;      M22=iv22;
      B1=ivmu1+2.0*(bBB*y1+dBB*y2);     B2=ivmu2;
    break;
  }
  double    detM=M11*M22-M12*M12;
  logden = 0.5*(log(detM)-(M22*B1*B1-2*M12*B1*B2+M11*B2*B2)/detM);
}

// detabd = det(paraSig)
double arrayCGH2d::GetTarget(double a, double b, double d, double z1,
        double z2, double detabd)
{  return       a*z1*z1+2*b*z1*z2+d*z2*z2-log(detabd);  }

void    arrayCGH2d::EM_Update_paraSig(Matrix<double>  &paraSigMat,
        double increnu11, double increnu12, double increnu22,
        double de)
{
  paraSigMat[0][0]=increnu11/de;
  paraSigMat[0][1]=paraSigMat[1][0]=increnu12/de;
  paraSigMat[1][1]=increnu22/de;
}

void    arrayCGH2d::SaveEMOldPara(double &oldp, double &olda, 
        double &oldb, Vector<double> &oldmu, Matrix<double> &oldV,
        Matrix<double> &oldSigAA, Matrix<double> &oldSigAB, 
        Matrix<double> &oldSigBA, Matrix<double> &oldSigBB,
        double newp, double newa, double newb, 
        Vector<double> newmu, Matrix<double> newV,
        Matrix<double> newSigAA, Matrix<double> newSigAB,
        Matrix<double> newSigBA, Matrix<double> newSigBB)
{
  oldp=newp;    olda=newa;      oldb=newb;  oldmu=newmu;    
  oldV=newV;    oldSigAA=newSigAA;
  oldSigAB=newSigAB;    oldSigBA=newSigBA;      oldSigBB=newSigBB;
}

/**************************** sub functions **************************/
void arrayCGH2d::SubFunOneInEM(double &emsig1, double &emsig2,
    double &emsig21, double &emsig22, double &emsig23, double tmp,
    double vb1,double vb2, double tmpv11,double tmpv12,double tmpv22)
{
   emsig1 += tmp*vb1;           emsig2 += tmp*vb2;
   emsig21+=tmp*(vb1*vb1+tmpv11);
   emsig22+=tmp*(vb1*vb2+tmpv12);
   emsig23+=tmp*(vb2*vb2+tmpv22);
}


/* The function calculates the inverse of a 2X2 symmetric matrix 'mat'
    (m11   m12) = mat^{-1}
    (m12   m22)
*/
void arrayCGH2d::Get2dSymMatInv(Matrix<double> mat, double &detmat,
        double &im11, double &im12, double &im22)
{
  detmat=mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
  im11=mat[1][1]/detmat; im12=-mat[0][1]/detmat; im22=mat[0][0]/detmat;
}

// This function compute        iv^{-1} * b.
void    arrayCGH2d::CompIVinvMultiB(double iv11, double iv12,
        double iv22, double b1, double b2, double &v11,
        double &v12, double &v22, double &vb1, double &vb2)
{
  double        tmp=1.0/(iv11*iv22-iv12*iv12);
  vb1 = tmp*(iv22*b1-iv12*b2);  vb2 = tmp*(iv11*b2-iv12*b1);
  v11 = tmp*iv22;       v12=-tmp*iv12;          v22=tmp*iv11;
}


void    cppEstHyper( double *obs, int *mixState, int *n, 
                        double *p, double *a, double *b, 
                        double *base1, double *base2, 
                        double *mu, double *v, double *sigAA, 
                        double *sigAB, double *sigBA, double *sigBB, 
                        int *K, int *M, int *niters, int *nBCMIXiters,
                        double *estSig, double *estBS, int *selectHyper) 
{
  
    ofstream progstream;
    progstream.open("cppEstHyperprogress.txt"); // Outputs debugging comments into a file.
    
    progstream << "In cppEstHyper 2.\n";
    
    
    // -------------- STEP 1: Unpack arguments passed in from R. -------------------
    
    int                 N = *n;
    Matrix<double>      obsvec(N+1,2);
    
    obsvec[0][0]=0; obsvec[0][1]=0;
    for (int i=1; i<N+1; i++){ obsvec[i][0]=obs[2*(i-1)]; obsvec[i][1] = obs[2*(i-1)+1];}
    
    progstream << "Copied in obsvec.\n";
    
    Vector<double>   paraMu(2); 
    paraMu[0] = mu[0]; paraMu[1] = mu[1]; 
    
    progstream << "Copied in paraMu:" << mu[0]<<" "<<mu[1]<<endl;
    
    int              paraK=*K, paraM=*M;
    double           parap = *p;
    double           paraa = *a;    
    double           parab = *b;
    double           base11 = *base1;
    double           base21 = *base2;

    progstream << "Copied in K="<<paraK<<", M="<<paraM<<", p="<<parap<<".\n";
    
    Matrix<double>   paraSigAA(2,2), paraSigAB(2,2), paraSigBA(2,2), 
                    paraSigBB(2,2), paraV(2,2);
    paraV[0][0] =  v[0]; paraV[0][1] = v[1]; paraV[1][0] = v[2]; paraV[1][1] = v[3];
    paraSigAA[0][0] =  sigAA[0]; paraSigAA[0][1] = sigAA[1]; paraSigAA[1][0] = sigAA[2]; paraSigAA[1][1] = sigAA[3];
    paraSigAB[0][0] =  sigAB[0]; paraSigAB[0][1] = sigAB[1]; paraSigAB[1][0] = sigAB[2]; paraSigAB[1][1] = sigAB[3];
    paraSigBA[0][0] =  sigBA[0]; paraSigBA[0][1] = sigBA[1]; paraSigBA[1][0] = sigBA[2]; paraSigBA[1][1] = sigBA[3];
    paraSigBB[0][0] =  sigBB[0]; paraSigBB[0][1] = sigBB[1]; paraSigBB[1][0] = sigBB[2]; paraSigBB[1][1] = sigBB[3];
    
    progstream << "Copied in paraV: "<<paraV[0][0]<<" "<<paraV[0][1]<<" "<<paraV[1][0]<<" "<<paraV[1][1]<<endl;
    progstream << "Copied in paraSigAA: "<<paraSigAA[0][0]<<" "<<paraSigAA[0][1]<<" "<<paraSigAA[1][0]<<" "<<paraSigAA[1][1]<<endl;
    progstream << "Copied in paraSigAB: "<<paraSigAB[0][0]<<" "<<paraSigAB[0][1]<<" "<<paraSigAB[1][0]<<" "<<paraSigAB[1][1]<<endl;
    progstream << "Copied in paraSigBA: "<<paraSigBA[0][0]<<" "<<paraSigBA[0][1]<<" "<<paraSigBA[1][0]<<" "<<paraSigBA[1][1]<<endl;
    progstream << "Copied in paraSigBB: "<<paraSigBB[0][0]<<" "<<paraSigBB[0][1]<<" "<<paraSigBB[1][0]<<" "<<paraSigBB[1][1]<<endl;
    
    arrayCGH2d       arr( obsvec, parap, paraa, parab, base11, base21, paraMu, paraV, paraSigAA, paraSigAB,
                            paraSigBA, paraSigBB, paraK, paraM);
    
    
    progstream << "Created arrayCGH2d object.\n";  
    progstream << "Reading in mixState:\n";
    
    arr.mixState[0] = 0;
    for (int i=1; i<N+1; i++){ 
        arr.mixState[i] = mixState[i-1];
    }
  
    
    arr.mystream.open("R_internal.txt");
    //arr.mystream << "Started arr.mystream...any debugging comments internal to arrayCGH2d object are piped here.\n";


   // -------------- STEP 2: Run BCMIX / FindState / SelectAllHyper Loop. -------------------


    progstream << "\n\nStarting BcmixSmooth...\n";
    progstream << "**** niters= "<<*niters << ", nBCMIXiters= "<<*nBCMIXiters<<endl;

    for (int k=1; k<= *niters; k++) {
        progstream << "\n\n\nk="<<k<<endl;

        for( int kk=1; kk<=*nBCMIXiters; kk++){
             progstream << "\n  kk="<<kk<<endl;
             arr.BcmixSmoothEst();    
             arr.FindState();            
        }
        
        for (int t=1; t<20; t++)  progstream << "t=" << t << " " <<arr.estPara[t][0]<<" "<<arr.estPara[t][1]<<" "<<arr.mixState[t]<<endl;    
        
        progstream << "selectHyper = "<< *selectHyper <<endl;
        
        if(*selectHyper > 0){
            progstream << "Computing hyperparameters ..."<<endl;
            arr.SelectAllHyper(*selectHyper);    
            
            progstream << "p=" << arr.p << "\t" << "a="<< arr.a << "\t" << "b=" << arr.b <<endl;
            progstream << "paraMu = " << arr.paraMu[0] << " " << arr.paraMu[1] << endl;
            progstream<< "paraV = " << arr.paraV << endl;
            progstream << arr.paraSigAA << "\n" << arr.paraSigAB << "\n" << arr.paraSigBA <<"\n" << arr.paraSigBB << endl;
            for (int t=1; t<20; t++)  progstream << "t=" << t << " " <<arr.estPara[t][0]<<" "<<arr.estPara[t][1]<<" "<<arr.mixState[t]<<endl;    
        }
        
    }
    progstream << "... BcmixSmooth Done.\n";

  // -------------- STEP 3: Package variables for returning to R.    -------------------

    for(int i=0; i<N; i++){
        estSig[2*i] = arr.estPara[i+1][0];
        estSig[2*i+1] = arr.estPara[i+1][1];
    }
    
    progstream << "Stored estSig.\n";
    
    for (int i=0; i<N; i++) estBS[i]=arr.estBS[i+1];
    
    progstream << "Stored estCP.\n";

    for (int i=0; i<N; i++) mixState[i] = arr.mixState[i+1];
    
    progstream << "Stored arr.mixState.\n";

    *p = arr.p;
    *a = arr.a;
    *b = arr.b;
    *base1 = arr.base1;
    *base2 = arr.base2;
    mu[0] = arr.paraMu[0];
    mu[1] = arr.paraMu[1];
    v[0] = arr.paraV[0][0]; v[1] = arr.paraV[0][1]; v[2] = arr.paraV[1][0]; v[3] = arr.paraV[1][1];
    sigAA[0] = arr.paraSigAA[0][0]; sigAA[1] = arr.paraSigAA[0][1]; sigAA[2] = arr.paraSigAA[1][0]; sigAA[3] = arr.paraSigAA[1][1];
    sigAB[0] = arr.paraSigAB[0][0]; sigAB[1] = arr.paraSigAB[0][1]; sigAB[2] = arr.paraSigAB[1][0]; sigAB[3] = arr.paraSigAB[1][1];
    sigBA[0] = arr.paraSigBA[0][0]; sigBA[1] = arr.paraSigBA[0][1]; sigBA[2] = arr.paraSigBA[1][0]; sigBA[3] = arr.paraSigBA[1][1];
    sigBB[0] = arr.paraSigBB[0][0]; sigBB[1] = arr.paraSigBB[0][1]; sigBB[2] = arr.paraSigBB[1][0]; sigBB[3] = arr.paraSigBB[1][1];
   
    progstream << "Stored estimated hyperparameters.\n";

    arr.mystream.close();  
    progstream.close();  
  
}




}
