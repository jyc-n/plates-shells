#ifndef PLATES_SHELLS_PRE_PROCESSOR_H
#define PLATES_SHELLS_PRE_PROCESSOR_H

class Arguments;
class Parameters;
class Geometry;
class Boundary;

class PreProcessorImpl {

public:
    PreProcessorImpl(Parameters* SimPar, Geometry* SimGeo);

    // main pre-processing function
    void PreProcess(Arguments t_args);

private:

    // subroutines
    void readInput();
    void readGeoFile();
    void buildNodes(int opt);
    void buildMesh();
    void buildNodeElementList();
    void buildEdgeHingeList();
    
    // other functions
    void print_parameters();

    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
};

#endif //PLATES_SHELLS_PRE_PROCESSOR_H
