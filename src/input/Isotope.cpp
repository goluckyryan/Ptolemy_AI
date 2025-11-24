#include "Isotope.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>

// Constants
const double amu = 931.494028;
const double mp = 938.272013;
const double mn = 939.565346;

void Isotope::SetIso(int a, int z){
    this->A = a;
    this->Z = z;
    FindMassByAZ(a,z); 
}

void Isotope::SetIsoByName(std::string name){
    FindMassByName(name); 
}

void Isotope::FindMassByAZ(int A, int Z){
  std::string line;
  int    lineNum=0;
  int    list_A, list_Z;

  std::ifstream myfile;
  int    flag=0;

  setFileLines();

  int numLineStart = fileStartLine;
  int numLineEnd  = fileEndLine;

  if ( A >= 50 && A < 100) numLineStart = lineMass050_099;
  if ( A >=100 && A < 150) numLineStart = lineMass100_149;
  if ( A >=150 && A < 200) numLineStart = lineMass150_199;
  if ( A >=200 ) numLineStart           = lineMass200;

  myfile.open(dataSource.c_str());

  if (myfile.is_open()) {
    while (flag == 0 && lineNum <numLineEnd){
      lineNum ++ ;
      getline (myfile,line);

      if (lineNum >= numLineStart ){
        list_Z     = atoi((line.substr(10,5)).c_str());
      	list_A     = atoi((line.substr(15,5)).c_str());

      	if ( A == list_A && Z == list_Z) {
            this->BEA       = atof((line.substr(54,11)).c_str());
      		this->Mass      = list_Z*mp + (list_A-list_Z)*mn - this->BEA/1000*list_A;
            this->MassError = atof((line.substr(65,7)).c_str());
            std::string str = line.substr(20,3);
            str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
            this->Symbol    = str;
            
            std::ostringstream ss;
            ss << A << this->Symbol;
            this->Name      = ss.str();
     		   flag = 1;
      	}else if ( list_A > A) {
          this->BEA  = -404;
          this->Mass = -404;
          this->MassError = -404;
          this->Symbol = "non";
          this->Name   = "non";
          break;
        }

      }
    }
    
    if( this->Name == "1H" ) this->Name = "p";
    if( this->Name == "2H" ) this->Name = "d";
    if( this->Name == "3H" ) this->Name = "t";
    if( this->Name == "4He" ) this->Name = "a";
    
    myfile.close();
  }else {
    printf("Unable to open %s\n", dataSource.c_str());
  }
}

void Isotope::FindMassByName(std::string name){

    // done seperate the Mass number and the name 
  if( name == "n" ) {
    this->Name = "1n";
    this->BEA       = 0;
    this->Mass      = mn;
    this->MassError = 0;
    this->Name      = "n";
    this->A         = 1;
    this->Z         = 0;
    return;
  }
    if( name == "p" ) name = "1H";
    if( name == "d" ) name = "2H";
    if( name == "t" ) name = "3H";
    if( name == "a" ) name = "4He";
    
    std::string temp = name;
    int lastDigit = 0;

    for(int i=0; temp[i]; i++){
      if(temp[i] == '0') lastDigit =  i; 
      if(temp[i] == '1') lastDigit =  i; 
      if(temp[i] == '2') lastDigit =  i; 
      if(temp[i] == '3') lastDigit =  i; 
      if(temp[i] == '4') lastDigit =  i; 
      if(temp[i] == '5') lastDigit =  i; 
      if(temp[i] == '6') lastDigit =  i; 
      if(temp[i] == '7') lastDigit =  i; 
      if(temp[i] == '8') lastDigit =  i; 
      if(temp[i] == '9') lastDigit =  i; 
    }

    this->Symbol = temp.substr(lastDigit + 1);
    //check is Symbol is 2 charaters, if not, add " " at the end
    if( this->Symbol.length() == 1 ){
      this->Symbol = this->Symbol + " ";
    }


    temp = name;
    int len = temp.length();
    temp = temp.substr(0, lastDigit + 1);
    
    this->A = atoi(temp.c_str());    

    // find the nucleus in the data
    std::string line;
    int    lineNum=0;
    int    list_A;
    std::string list_symbol;

    std::ifstream myfile;
    int    flag=0;

    setFileLines();

    int numLineStart = fileStartLine;
    int numLineEnd  = fileEndLine;

    if ( A >= 50 && A < 100) numLineStart = lineMass050_099;
    if ( A >=100 && A < 150) numLineStart = lineMass100_149;
    if ( A >=150 && A < 200) numLineStart = lineMass150_199;
    if ( A >=200 ) numLineStart           = lineMass200;

    myfile.open(dataSource.c_str());

    if (myfile.is_open()) {
      while (flag == 0 && lineNum <numLineEnd){
        lineNum ++ ;
        getline (myfile,line);

        if (lineNum >= numLineStart ){
          list_symbol  = line.substr(20,3);
          // Trim spaces from list_symbol
          list_symbol.erase(std::remove(list_symbol.begin(), list_symbol.end(), ' '), list_symbol.end());
          
          list_A       = atoi((line.substr(15,5)).c_str());

          // Trim spaces from this->Symbol for comparison
          std::string currentSymbol = this->Symbol;
          currentSymbol.erase(std::remove(currentSymbol.begin(), currentSymbol.end(), ' '), currentSymbol.end());

          if ( this->A == list_A &&  currentSymbol == list_symbol) {
            this->Z         = atoi((line.substr(10,5)).c_str());
            this->BEA       = atof((line.substr(54,11)).c_str());
       		this->Mass      = this->Z*mp + (list_A-this->Z)*mn - this->BEA/1000*list_A;
            this->MassError = atof((line.substr(65,7)).c_str());
            
            // Update Symbol to match file (or keep trimmed)
            this->Symbol    = list_symbol; 
            
            std::ostringstream ss;
            ss << this->A << this->Symbol;
            this->Name      = ss.str();
            flag = 1;
          }else if ( list_A > this->A) {
            this->BEA  = -404;
            this->Mass = -404;
            this->MassError = -404;
            this->Symbol = "non";
            this->Name   = "non";
            break;
          }

        }
      }
      myfile.close();
    }else {
      printf("Unable to open %s\n", dataSource.c_str());
    }
}

double Isotope::CalSp(int Np, int Nn){
  Isotope nucleusD(A - Np - Nn, Z - Np);  

  if( nucleusD.Mass != -404){
    return nucleusD.Mass + Nn*mn + Np*mp - this->Mass;
  }else{
    return -404;
  }
}

double Isotope::CalSp2(int a, int z){
  Isotope nucleusD(A - a , Z - z);
  Isotope nucleusS(a,z);  

  if( nucleusD.Mass != -404 && nucleusS.Mass != -404){
    return nucleusD.Mass + nucleusS.Mass - this->Mass;
  }else{
    return -404;
  }
}
