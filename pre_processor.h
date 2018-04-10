#ifndef PLATES_SHELLS_PRE_PROCESSOR_H
#define PLATES_SHELLS_PRE_PROCESSOR_H

class Parameters;
class Geometry;
class Element;

// pre-processing functions
void read_input(Parameters& Params);
std::vector<Element> pre_processor(Geometry& Geo, Parameters& Params);

// other helper function for parameters
std::string& trim(std::string& str);

#endif //PLATES_SHELLS_PRE_PROCESSOR_H
