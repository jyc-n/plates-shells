#ifndef PLATES_SHELLS_PRE_PROCESSOR_H
#define PLATES_SHELLS_PRE_PROCESSOR_H

class Parameters;
class Geometry;
class Boundary;

class PreProcessorImpl {

public:
    PreProcessorImpl(Parameters* SimPar, Geometry* SimGeo);

    // main pre-processing function
    void PreProcess(
                   const bool AR_FLAG, const bool K_FLAG,
                   int num_len, int num_wid,
                   double k_s, double k_sh, double k_b
                   );

private:

    // subroutines
    void readInput();
    void readGeoFile();
    void buildNodes(bool AR_FLAG);
    void buildMesh();
    void buildNodeElementList();
    void buildEdgeHingeList();

    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
};

#endif //PLATES_SHELLS_PRE_PROCESSOR_H
