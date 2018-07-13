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

    // subroutines
    void readInput();
    void initCoord();
    void initConn();
    void initGeoList();

private:
    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
};

#endif //PLATES_SHELLS_PRE_PROCESSOR_H
