// Same as betcal_faithful_test.cpp but with debug prints for ILO=0 accumulation
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
    const int JPBASE=1,JPMX=3,NLO=LMX-LMN+1,NMX=LXMX+1,NMLX=LXMX-LXMN+1;
    const double UK=1.23280,FACTOR=0.5/UK,PI=3.14159265358979320;
    const double ETA_A=7.416,ETA_B=0.612;
    struct Slot{int ldel,lx,jp;bool active;};
    std::vector<Slot> jtocs;
    int nact=0;
    for(int jp=JPBASE;jp<=JPMX;jp+=2)
        for(int ldel=-LDELMX;ldel<=LDELMX;ldel++)
            {jtocs.push_back({ldel,LXMX,jp,true});nact++;}
    
    std::vector<std::vector<float>> SMAG(nact,std::vector<float>(NLO));
    std::vector<std::vector<float>> SPHASE(nact,std::vector<float>(NLO));
    for(int li=LMN;li<=LMX;li++){
        int ili=li-LMN;
        for(int k=0;k<nact;k++){
            SMAG[k][ili]=fabsf(lcg_next()*0.5f);
            SPHASE[k][ili]=lcg_next()*(float)PI;
        }
    }
    
    std::vector<double> sigin(LMX+3,0),sigot(LMX+3,0);
    for(int l=1;l<=LMX+1;l++){sigin[l]=sigin[l-1]+atan2(ETA_A,(double)l);sigot[l]=sigot[l-1]+atan2(ETA_B,(double)l);}
    
    std::vector<std::vector<std::complex<double>>> BETAS(nact,std::vector<std::complex<double>>(NLO,0.0));
    
    for(int lo=LMN;lo<=LMX;lo++){
        int ilo=lo-LMN;
        int LIMN=lo-LDELMX;
        int MXMX=std::min(LXMX,lo);
        std::vector<double> TEMPS(500,0.0);
        
        for(int li=LIMN;li<=lo+LDELMX;li+=2){
            if(li<LMN||li>LMX) continue;
            int LX1=std::max(std::abs(li-lo),LXMN);
            for(int lx=LX1;lx<=LXMX;lx++){
                int mxz=(lx+li-lo)%2; if(mxz<0)mxz+=2;
                int base=NMX*(lx-LXMN+NMLX*(li-LIMN)/2);
                for(int mx=mxz;mx<=lx;mx++){
                    double cg=ClebschGordan((double)li,0.0,(double)lx,(double)mx,(double)lo,(double)mx);
                    TEMPS[base+mx]=FACTOR*(2*li+1)*cg;
                }
            }
        }
        
        int LIPREV=-100000,MXZ_cur=0,KOFFZ_cur=0;
        for(int k=0;k<nact;k++){
            if(!jtocs[k].active) continue;
            int lx=jtocs[k].lx,ldel=jtocs[k].ldel,li=lo-ldel;
            if(li>=LIPREV){MXZ_cur=lx+ldel;KOFFZ_cur=k-std::max(0,MXZ_cur);}
            LIPREV=li;
            if(li<LMN||li>LMX) continue;
            float amag=SMAG[k][li-LMN];
            if(amag==0.0f) continue;
            double phase=(double)SPHASE[k][li-LMN]+sigin[li]+sigot[lo];
            double smatr=(double)amag*sin(phase),smati=-(double)amag*cos(phase);
            int tbase=NMX*(lx-LXMN+NMLX*(li-LIMN)/2);
            for(int mx=std::max(0,MXZ_cur);mx<=lx;mx++){
                double t=TEMPS[tbase+mx];
                int kidx=KOFFZ_cur+mx;
                if(t!=0.0&&kidx>=0&&kidx<nact){
                    auto old=BETAS[kidx][ilo];
                    BETAS[kidx][ilo]+=std::complex<double>(t*smatr,t*smati);
                    if(ilo==0&&(kidx==0||kidx==2)) // trace k=0,2 for ILO=0
                        printf("LO=%d k=%d->kidx=%d MX=%d t=%.6f smatr=%.6f: BETAS[%d][0] %.6f->%.6f\n",
                               lo,k,kidx,mx,t,smatr,kidx,old.real(),BETAS[kidx][0].real());
                }
            }
        }
        
        std::vector<double> sqfac(LXMX+2,1.0);
        for(int mx=MXMX+1;mx<=LXMX+1;mx++) sqfac[mx]=0.0;
        for(int k=0;k<nact;k++){
            if(!jtocs[k].active) continue;
            int lx=jtocs[k].lx,ldel=jtocs[k].ldel,mx=(ldel+lx+1)/2;
            if(mx<0||mx>LXMX+1) continue;
            if(ilo==0&&(k==0||k==2)&&BETAS[k][0]!=std::complex<double>(0,0))
                printf("LO=%d sqfac pass: k=%d MX=%d sqfac=%.6f BETAS[k][0] %.6f->%.6f\n",
                       lo,k,mx,sqfac[mx],BETAS[k][0].real(),BETAS[k][0].real()*sqfac[mx]);
            BETAS[k][ilo]*=sqfac[mx];
        }
    }
    
    printf("\nFinal BETAS[k][ILO=0]:\n");
    for(int k=0;k<nact;k++) printf("k=%d: (%.8f,%.8f)\n",k,BETAS[k][0].real(),BETAS[k][0].imag());
}
