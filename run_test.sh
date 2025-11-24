#!/bin/bash
#run the original ptolemy using Cleopatra and the file "DWBA" in working directory
if [ ! -f "cleopatra.txt" ]; then
    echo "Running Cleopatra (Ptolemy) and copying results..."
    cd ../digios/analysis/working
    ../Cleopatra/Cleopatra.sh DWBA 1 1 1 0 0 0 180 5
    #cp the results to the Cpp_AI directory
    cd ../../../Cpp_AI
    cp ../digios/analysis/working/DWBA cleopatra.txt
    cp ../digios/analysis/working/DWBA.in ptolemy_input.txt
    cp ../digios/analysis/working/DWBA.out ptolemy_output.txt
    cp ../digios/analysis/working/DWBA.Xsec.txt ptolemy_Xsec.txt
else
    echo "Cleopatra results (cleopatra.txt) already exist. Skipping Ptolemy run and file copy."
fi

# Compile C++ code
echo "Compiling C++ code..."
g++ -std=c++17 -Iinclude src/main.cpp src/dwba/*.cpp src/input/*.cpp -o dwba
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

# Run C++ DWBA
echo "Running C++ DWBA..."
./dwba test_input.txt 0 180 5 > cpp_output.txt 2>&1