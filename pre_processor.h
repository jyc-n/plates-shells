#ifndef PLATES_SHELLS_PRE_PROCESSOR_H
#define PLATES_SHELLS_PRE_PROCESSOR_H

class Parameters;
class Geometry;
class Boundary;

class PreProcessorImpl {

public:
    PreProcessorImpl(Parameters* SimPar, Geometry* SimGeo);

    // main pre-processing function
    void PreProcess();

private:

    // subroutines
    void readInput();
    void readGeoFile();
    void buildNodes();
    void buildMesh();
    void buildNodeElementList();
    void buildEdgeHingeList();

    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
};

#endif //PLATES_SHELLS_PRE_PROCESSOR_H
