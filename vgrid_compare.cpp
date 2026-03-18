// vgrid_compare.cpp: compare C++ VMIN/VMAX/VMID per U vs Fortran vgrid_test output
// Build: g++ -O2 -std=c++17 vgrid_compare.cpp -o vgrid_compare
// Run:   ../fortran_testing/vgrid_test > /tmp/fort_vgrid.txt && ./vgrid_compare < /tmp/fort_vgrid.txt
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

static void GaussLegendre(int n, std::vector<double>& x, std::vector<double>& w) {
    x.resize(n); w.resize(n);
    if(n==1){x[0]=0;w[0]=2;return;}
    int m=(n+1)/2;
    double e1=n*(n+1),e2=M_PI/(4*n+2),e3=1.-(1.-1./n)/(8.*n*n);
    for(int i=1;i<=m;i++){
        double t=(4*i-1)*e2,x0=e3*cos(t),pkm1=1,pk=x0,fk=1;
        for(int k=2;k<=n;k++){fk++;double t1=x0*pk,pkp1=t1-pkm1-(t1-pkm1)/fk+t1;pkm1=pk;pk=pkp1;}
        double den=1-x0*x0,d1=n*(pkm1-x0*pk),dpn=d1/den;
        double d2=( 2*x0*dpn-e1*pk)/den,d3=(4*x0*d2+(2-e1)*dpn)/den,d4=(6*x0*d3+(6-e1)*d2)/den;
        double u=pk/dpn,v=d2/dpn,h1=1+.5*u*(v+u*(v*v-u*d3/(3*dpn))),h=-u*h1;
        double p=pk+h*(dpn+.5*h*(d2+h/3*(d3+.25*h*d4)));
        double dp=dpn+h*(d2+.5*h*(d3+h*d4/3));h=h-p/dp;
        x[n-i]=x0+h;x[i-1]=-(x0+h);
        double fx=d1-h*e1*(pk+.5*h*(dpn+h/3*(d2+.25*h*(d3+.2*h*d4))));
        w[n-i]=2*(1-x[i-1]*x[i-1])/(fx*fx);w[i-1]=w[n-i];
    }
    if((m+m)>n)x[m-1]=0;
}
static void CubMap(int mt,double xlo,double xmid,double xhi,double gamma,
                   std::vector<double>& a,std::vector<double>& w){
    int n=a.size();
    double tau=(gamma>1e-6)?log(gamma+sqrt(gamma*gamma+1)):gamma*(1-gamma*gamma/6);
    double xlen=xhi-xlo,xadd=xlo+xhi;
    if(mt==0){for(int i=0;i<n;i++){a[i]=xlo+.5*xlen*(a[i]+1);w[i]*=.5*xlen;}}
    else if(mt==1){xmid=std::max(xmid,xlo+xlen/7);xmid=std::min(xmid,.5*xadd);double A=.5*xadd-xmid,B=.5*xlen,C=.5*xadd;for(int i=0;i<n;i++){double tu=tau*a[i],xi=sinh(tu)/gamma;a[i]=A*(xi*xi-1)*(xi+1)+B*xi+C;w[i]*=(tau/gamma)*cosh(tu)*((3*xi-1)*(xi+1)*A+B);}}
    else if(mt==2){double A=-xmid*xlen,B=xlen,C=xmid*xadd-2*xlo*xhi,D=xadd-2*xmid;for(int i=0;i<n;i++){double tu=tau*a[i],sh=sinh(tu),denom=B-(D/gamma)*sh;a[i]=(-A+(C/gamma)*sh)/denom;w[i]*=(tau/gamma)*cosh(tu)*((B*C-A*D)/(denom*denom));}}
    else if(mt==3){for(int i=0;i<n;i++){double tu=tau*a[i];a[i]=xlo+.5*xlen*(sinh(tu)/gamma+1);w[i]*=.5*xlen*(tau/gamma)*cosh(tu);}}
}
static double Aitlag(double x,double si,const std::vector<double>& t,int nait){
    int n=t.size(); if(x<0)return t[0]; double xh=x*si; int ii=int(xh);
    if(ii>=n)return t[n-1];
    static const double inv[]={0,1,.5,1./3,.25,.2,1./6,1./7,.125,1./9,.1};
    int ns=ii-nait/2; if(ns<0)ns=0; int ne=ns+nait; if(ne>=n)ne=n-1; ns=ne-nait;
    std::vector<double> fs(nait+2),de(nait+2);
    fs[1]=t[ns]; double d=xh-ns; de[1]=d; double f=0;
    for(int i=1;i<=nait;i++){f=t[ns+i];d-=1;for(int j=1;j<=i;j++)f=(f*de[j]-fs[j]*d)*inv[i+1-j];de[i+1]=d;fs[i+1]=f;}
    return f;
}
static double bsprod(const std::vector<double>& tT,const std::vector<double>& tP,
                     double hT,double hP,double bmxT,double bmxP,
                     double S1,double T1,double S2,double T2,
                     double ri,double ro,double X){
    double rt2=S1*S1*ri*ri+T1*T1*ro*ro+2*S1*T1*ri*ro*X;
    double rp2=S2*S2*ri*ri+T2*T2*ro*ro+2*S2*T2*ri*ro*X;
    if(rt2<0)rt2=0; if(rp2<0)rp2=0;
    double rt=sqrt(rt2),rp=sqrt(rp2);
    if(rt>bmxT||rp>bmxP)return 0;
    return Aitlag(rt,1/hT,tT,4)*Aitlag(rp,1/hP,tP,4);
}

int main(){
    // Same setup as Fortran
    const int NPSUM=40,NTAB=641,NPDIF=40,LOOKST=100;
    const double hT=0.05,hP=0.05,bmxT=(NTAB-1)*hT,bmxP=(NTAB-1)*hP;
    std::vector<double> tabT(NTAB),tabP(NTAB);
    for(int i=0;i<NTAB;i++){double r=i*hT;tabT[i]=(i==0)?0:r*exp(-0.433*r);}
    for(int i=0;i<NTAB;i++){double r=i*hP;tabP[i]=-0.8*exp(-0.3*r*r);}
    // Kinematics
    double tmp=1.0/(1.0+(1.0/16.0)*(1.0+1.0));
    double S1=(1+1)*(1+1.0/16)*tmp,T1=-(1+1.0/16)*tmp,S2=(1+1)*tmp,T2=-S1;
    double ST1=2*S1*T1,ST2=2*S2*T2;
    double DXV=2.0/(LOOKST*LOOKST);
    double XS[2]={1.0-DXV,-1.0+DXV};
    double WVWMAX=75,DWCUT=2e-6,RVRLIM=WVWMAX*DWCUT;
    // U grid
    std::vector<double> xi_s(NPSUM),wi_s(NPSUM);
    GaussLegendre(NPSUM,xi_s,wi_s);
    double SUMMIN=0,SUMMAX=30.4,SUMMID=SUMMAX/2,GAMSUM=1;
    CubMap(2,SUMMIN,SUMMID,SUMMAX,GAMSUM,xi_s,wi_s);

    // Stage 1: VMAX/VMIN scan
    std::vector<double> vmax_u(NPSUM),vmin_u(NPSUM);
    double DV=1.0/LOOKST;
    double VMAX=1,VMIN=1;
    for(int IU=0;IU<NPSUM;IU++){
        double U=xi_s[IU],VLEN=2*U;
        if(U<1){vmax_u[IU]=VLEN;vmin_u[IU]=-VLEN;continue;}
        // Positive scan
        double vval=std::min(1.0,VMAX+3*DV);
        while(vval>0.5*DV){
            double ri=U+vval*0.5*VLEN,ro=U-vval*0.5*VLEN;
            if(ri<1e-6||ro<1e-6){vval=0;break;}
            double ulim=RVRLIM/std::max(1e-2,ri*ro);
            bool ok=true;
            for(int k=0;k<2;k++){double f=std::abs(bsprod(tabT,tabP,hT,hP,bmxT,bmxP,S1,T1,S2,T2,ri,ro,XS[k]));if(f>ulim){ok=false;break;}}
            if(!ok){vval=std::min(1.0,vval+DV);break;}
            vval-=DV;
        }
        vval=std::min(1.0,vval);
        {double ri=U+vval*0.5*VLEN,ro=U-vval*0.5*VLEN;
         double rt=sqrt(S1*S1*ri*ri+T1*T1*ro*ro+ST1*ri*ro*(1-DXV));
         double rp=sqrt(S2*S2*ri*ri+T2*T2*ro*ro+ST2*ri*ro*(1-DXV));
         if(rt>bmxT||rp>bmxP)vval-=DV;}
        VMAX=std::max(0.0,vval);
        vmax_u[IU]=VMAX*VLEN;
        // Negative scan
        vval=std::min(1.0,VMIN+3*DV);
        while(vval>0.5*DV){
            double ri=U-vval*0.5*VLEN,ro=U+vval*0.5*VLEN;
            if(ri<0||ro<0){vval=0;break;}
            double ulim=RVRLIM/std::max(1e-2,std::abs(ri*ro));
            bool ok=true;
            for(int k=0;k<2;k++){double f=std::abs(bsprod(tabT,tabP,hT,hP,bmxT,bmxP,S1,T1,S2,T2,ri,ro,XS[k]));if(f>ulim){ok=false;break;}}
            if(!ok){vval=std::min(1.0,vval+DV);break;}
            vval-=DV;
        }
        vval=std::min(1.0,vval);
        {double ri=U-vval*0.5*VLEN,ro=U+vval*0.5*VLEN;
         double rt=sqrt(S1*S1*ri*ri+T1*T1*ro*ro+ST1*ri*ro*(1-DXV));
         double rp=sqrt(S2*S2*ri*ri+T2*T2*ro*ro+ST2*ri*ro*(1-DXV));
         if(rt>bmxT||rp>bmxP)vval-=DV;}
        VMIN=std::max(0.0,vval);
        vmin_u[IU]=-VMIN*VLEN;
    }

    // Stage 2: VMID (V of max |bsprod|)
    std::vector<double> vmid_u(NPSUM);
    for(int IU=0;IU<NPSUM;IU++){
        double U=xi_s[IU],vmn=vmin_u[IU],vmx=vmax_u[IU];
        int IMAX=2*NPDIF;
        double dv2=(vmx-vmn)/(IMAX+1),vval=vmn,pvpmax=0,vofmax=0.5*(vmn+vmx);
        for(int IV=1;IV<=IMAX;IV++){
            vval+=dv2;
            double ri=U+0.5*vval,ro=U-0.5*vval,temp=0;
            for(int k=0;k<2;k++) temp+=std::abs(bsprod(tabT,tabP,hT,hP,bmxT,bmxP,S1,T1,S2,T2,ri,ro,XS[k]));
            if(temp>pvpmax){pvpmax=temp;vofmax=vval;}
        }
        double vm=vofmax,tmp2=0.3*(vmx-vmn);
        vm=std::min(std::max(vm,vmn+tmp2),vmx-tmp2);
        vmid_u[IU]=vm;
    }

    // Read Fortran output
    struct FRow{double U,vmin,vmid,vmax;};
    std::vector<FRow> frows;
    std::string line;
    while(std::getline(std::cin,line)){
        if(line.find("IU")==std::string::npos&&line.find("=")==std::string::npos){
            std::istringstream ss(line);int iu;double u,vn,vm,vx;
            if(ss>>iu>>u>>vn>>vm>>vx)frows.push_back({u,vn,vm,vx});
        }
    }
    // Compare
    std::cout<<std::setw(4)<<"IU"<<std::setw(8)<<"U"
             <<std::setw(12)<<"F_VMIN"<<std::setw(12)<<"C_VMIN"<<std::setw(8)<<"Δ%"
             <<std::setw(12)<<"F_VMID"<<std::setw(12)<<"C_VMID"<<std::setw(8)<<"Δ%"
             <<std::setw(12)<<"F_VMAX"<<std::setw(12)<<"C_VMAX"<<std::setw(8)<<"Δ%"<<"\n";
    double max_vmax=0,max_vmin=0,max_vmid=0;
    for(int i=0;i<(int)frows.size()&&i<NPSUM;i++){
        auto&f=frows[i];
        double cvn=vmin_u[i],cvm=vmid_u[i],cvx=vmax_u[i];
        auto rel=[](double a,double b){return std::abs(b)>1e-6?100*(a-b)/b:0;};
        double rvn=rel(cvn,f.vmin),rvm=rel(cvm,f.vmid),rvx=rel(cvx,f.vmax);
        max_vmin=std::max(max_vmin,std::abs(rvn));
        max_vmid=std::max(max_vmid,std::abs(rvm));
        max_vmax=std::max(max_vmax,std::abs(rvx));
        std::string flag=(std::abs(rvx)>5||std::abs(rvm)>5)?" <<<":"";
        std::cout<<std::setw(4)<<i+1<<std::fixed<<std::setw(8)<<std::setprecision(2)<<f.U
                 <<std::setw(12)<<std::setprecision(3)<<f.vmin<<std::setw(12)<<cvn<<std::setw(8)<<std::setprecision(1)<<rvn
                 <<std::setw(12)<<f.vmid<<std::setw(12)<<cvm<<std::setw(8)<<rvm
                 <<std::setw(12)<<f.vmax<<std::setw(12)<<cvx<<std::setw(8)<<rvx<<flag<<"\n";
    }
    std::cout<<"\nMax VMIN err: "<<max_vmin<<"%, VMID err: "<<max_vmid<<"%, VMAX err: "<<max_vmax<<"%\n";
}
