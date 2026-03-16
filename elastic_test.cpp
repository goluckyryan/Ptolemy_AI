// elastic_test.cpp — Test the unified ElasticSolver for two cases
#include <iostream>
#include "elastic.h"

int main() {
    // ===== Test 1: 148Sm(α,α) at 50 MeV =====
    std::cout << "===== 148Sm(α,α) at 50 MeV =====\n";
    {
        ElasticSolver s;
        s.SetSystem(4, 2, 148, 62, 50.0);
        s.AddVolumeWS({-65.5, -29.8}, 1.427, 0.671);  // V<0 = attractive
        s.AddCoulomb(1.4);
        s.SetLmax(90);
        s.CalcKinematics();
        s.PrintKinematics();
        s.CalcScatteringMatrix();
        s.PrintSMatrix(12);
        s.SaveSMatrix("sm148_aa_smat_cpp.txt");
        s.SaveDCS("sm148_aa_xsec_cpp.txt");
        std::cout << "\nDCS at key angles:\n";
        for (double th : {10.0, 20.0, 40.0, 60.0, 90.0, 120.0, 150.0})
            std::cout << "  theta=" << th << ": " << s.DCSUnpolarized(th) << " mb/sr\n";
    }

    // ===== Test 2: 60Ni(p,p) at 30 MeV =====
    std::cout << "\n===== 60Ni(p,p) at 30 MeV =====\n";
    {
        ElasticSolver s;
        s.SetSystem(1, 1, 60, 28, 30.0);
        s.AddVolumeWS ({-47.937, -2.853}, 1.200, 0.669);
        s.AddSurfaceWS({    0.0, -6.878}, 1.280, 0.550);
        s.AddSpinOrbit({ -5.250,  0.162}, 1.020, 0.590);
        s.AddCoulomb(1.258);
        s.SetLmax(25);
        s.CalcKinematics();
        s.PrintKinematics();
        s.CalcScatteringMatrix();
        s.PrintSMatrix(12);
        s.SaveSMatrix("ni60pp_smat_cpp.txt");
        s.SaveDCS("ni60pp_xsec_cpp.txt");
        std::cout << "\nDCS at key angles:\n";
        for (double th : {10.0, 20.0, 30.0, 40.0, 60.0, 90.0, 120.0})
            std::cout << "  theta=" << th << ": " << s.DCSUnpolarized(th) << " mb/sr\n";
    }

    return 0;
}
