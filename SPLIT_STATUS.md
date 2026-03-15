# SPLIT_STATUS.md — dwba.cpp extraction progress

xsectn.cpp: DONE — XSectn (line 1698-2267) extracted
ineldc.cpp: DONE — GrdSet (line 568-575) + InelDc (line 576-1697) extracted
  - File: src/dwba/ineldc.cpp (1213 lines)
  - Includes CubMap static helper (reproduced from dwba.cpp lines 19-77)
  - Compile test: PASSED (no errors, /tmp/ineldc.o 53K)
  - Symbols verified: DWBA::GrdSet, DWBA::InelDc both present
