#ifndef ISOTOPE_H
#define ISOTOPE_H

#include <string>
#include <vector>
#include <cmath>

class Isotope {
public:
  int A, Z;
  double Mass, MassError, BEA;
  std::string Name, Symbol;
  std::string dataSource;
  
  Isotope(){findHeliosPath();};
  Isotope(int a, int z){ findHeliosPath();SetIso(a,z);  };
  Isotope(std::string name){ findHeliosPath(); SetIsoByName(name); };

  void SetIso(int a, int z);
  void SetIsoByName(std::string name);

  double CalSp(int Np, int Nn);
  double CalSp2(int a, int z);

  double CalBeta(double T){
    double Etot = Mass + T;
    double gamma = 1 + T/Mass;
    double beta = sqrt(1 - 1 / gamma / gamma ) ;
    return beta;
  }
   
private:
  void FindMassByAZ(int a, int z);
  void FindMassByName(std::string name);

  int fileStartLine;
  int fileEndLine;
  int lineMass050_099;
  int lineMass100_149;
  int lineMass150_199;
  int lineMass200;
  
  void setFileLines(){
    fileStartLine = 37;
    fileEndLine = 3594;
    
    lineMass050_099 = 466;
    lineMass100_149 = 1160;
    lineMass150_199 = 1994;
    lineMass200     = 2774;
  }
  
  void findHeliosPath(){
      dataSource = "data/mass20.txt";
  }
};

#endif
