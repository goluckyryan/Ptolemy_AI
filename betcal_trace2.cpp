// Use math_utils for CG
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <complex>
#include <vector>
#include "include/math_utils.h"

int32_t lcg_state=12345;
float lcg_next(){lcg_state=(int32_t)(1664525LL*(int64_t)lcg_state+1013904223LL);return(float)((double)lcg_state*4.656612873e-10);}

int main(){
    const int LMN=0,LMX=20,LXMN=2,LXMX=2,LDELMX=2;
    const int JPBASE=1,JPMX=3,NLO=LMX-LMN+1;
    const int NMX=LXMX+1,NMLX=LXMX-LXMN+1;
    const double UK=1.23280,FACTOR=0.5/UK,PI=3.14159265358979320;
    const double ETA_A=7.416,ETA_B=0.612;
    int nact=0;
    struct Slot{int ldel,lx,jp;};
    std::vector<Slot> jtocs;
    for(int jp=JPBASE;jp<=JPMX;jp+=2)
        for(int ldel=-LDELMX;ldel<=LDELMX;ldel++)
            {jtocs.push_back({ldel,LXMX,jp});nact++;}
    
    // SMAG/SPHASE
    std::vector<std::vector<float>> SMAG(nact,std::vector<float>(NLO));
    std::vector<std::vector<float>> SPHASE(nact,std::vector<float>(NLO));
    for(int li=LMN;li<=LMX;li++){
        int ili=li-LMN;
        for(int k=0;k<nact;k++){
            SMAG[k][ili]=fabsf(lcg_next()*0.5f);
            SPHASE[k][ili]=lcg_next()*(float)PI;
        }
    }
    
    // Coulomb phases
    std::vector<double> sigin(LMX+3,0),sigot(LMX+3,0);
    for(int l=1;l<=LMX+1;l++){sigin[l]=sigin[l-1]+atan2(ETA_A,(double)l);sigot[l]=sigot[l-1]+atan2(ETA_B,(double)l);}
    
    // BETCAL — only for LO=0 (to trace)
    std::vector<std::vector<std::complex<double>>> BETAS(nact,std::vector<std::complex<double>>(NLO,0.0));
    
    int lo=0;
    int ilo=lo-LMN;
    int LIMN=lo-LDELMX;
    int MXMX=std::min(LXMX,lo);
    int NTEMPS=500;
    std::vector<double> TEMPS(NTEMPS,0.0);
    
    // Build TEMPS for LO=0
    for(int li=LIMN;li<=lo+LDELMX;li+=2){
        if(li<LMN||li>LMX) continue;
        int LX1=std::max(std::abs(li-lo),LXMN);
        for(int lx=LX1;lx<=LXMX;lx++){
            int mxz=(lx+li-lo)%2; if(mxz<0)mxz+=2;
            int base=NMX*(lx-LXMN+NMLX*(li-LIMN)/2);
            for(int mx=mxz;mx<=lx;mx++){
                double cg=ClebschGordan((double)li,0.0,(double)lx,(double)mx,(double)lo,(double)mx);
                TEMPS[base+mx]=FACTOR*(2*li+1)*cg;
                printf("TEMPS[%d]=%.8f (li=%d,lx=%d,mx=%d,lo=%d,CG=%.8f)\n",base+mx,TEMPS[base+mx],li,lx,mx,lo,cg);
            }
        }
    }
    
    // KOFFS loop for LO=0
    int LIPREV=-100000,MXZ_cur=0,KOFFZ_cur=0;
    for(int k=0;k<nact;k++){
        int lx=jtocs[k].lx,ldel=jtocs[k].ldel;
        int li=lo-ldel;
        if(li>=LIPREV){MXZ_cur=lx+ldel;KOFFZ_cur=k-std::max(0,MXZ_cur);}
        LIPREV=li;
        if(li<LMN||li>LMX) continue;
        int ili_li=li-LMN;
        float amag=SMAG[k][ili_li];
        if(amag==0.0f) continue;
        float phase_f=SPHASE[k][ili_li];
        double phase=(double)phase_f+sigin[li]+sigot[lo];
        double smatr=(double)amag*sin(phase);
        double smati=-(double)amag*cos(phase);
        int temps_base=NMX*(lx-LXMN+NMLX*(li-LIMN)/2);
        printf("k=%d ldel=%d li=%d: amag=%.6f smatr=%.6f smati=%.6f\n",k,ldel,li,amag,smatr,smati);
        for(int mx=std::max(0,MXZ_cur);mx<=lx;mx++){
            double t=TEMPS[temps_base+mx];
            int kidx=KOFFZ_cur+mx;
            if(t!=0.0) printf("  MX=%d TEMPS[%d]=%.8f kidx=%d contrib=(%.8f,%.8f)\n",mx,temps_base+mx,t,kidx,t*smatr,t*smati);
            if(kidx>=0&&kidx<nact) BETAS[kidx][ilo]+=std::complex<double>(t*smatr,t*smati);
        }
    }
    
    // sqfac second pass
    std::vector<double> sqfac(LXMX+2,1.0);
    for(int mx=MXMX+1;mx<=LXMX+1;mx++) sqfac[mx]=0.0;
    printf("MXMX=%d sqfac[0..3]=%f %f %f %f\n",MXMX,sqfac[0],sqfac[1],sqfac[2],sqfac[3]);
    
    printf("BETAS before sqfac:\n");
    for(int k=0;k<nact;k++) if(BETAS[k][ilo]!=std::complex<double>(0,0))
        printf("  k=%d: (%.8f,%.8f)\n",k,BETAS[k][ilo].real(),BETAS[k][ilo].imag());
    
    for(int k=0;k<nact;k++){
        int lx=jtocs[k].lx,ldel=jtocs[k].ldel;
        int mx=(ldel+lx+1)/2;
        BETAS[k][ilo]*=sqfac[mx];
    }
    
    printf("BETAS after sqfac (LO=0):\n");
    for(int k=0;k<nact;k++) if(BETAS[k][ilo]!=std::complex<double>(0,0))
        printf("  k=%d: (%.8f,%.8f)\n",k,BETAS[k][ilo].real(),BETAS[k][ilo].imag());
    printf("Expected: k=0: (-0.00421, -0.04060), k=5: (0.24275, -0.17778)\n");
}
