#ifndef PTOLEMY_PARSER_H
#define PTOLEMY_PARSER_H

#include <string>
#include "dwba.h"

class PtolemyParser {
public:
    PtolemyParser();
    ~PtolemyParser();
    
    void ParseFile(const std::string& filename, DWBA& dwba);
    
private:
    void ParseLine(const std::string& line, DWBA& dwba);
};

#endif
