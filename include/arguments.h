#ifndef PLATES_SHELLS_ARGUMENTS_H
#define PLATES_SHELLS_ARGUMENTS_H

#include <string>

enum Opts {
    DEFAULT = 1,
    NUMS = 3,
    NUMS_DIMS = 5
};

class Arguments {

public:
    Arguments(int t_argc);
    void getArgs(int t_nl, int t_nw, double t_l=0.0, double t_w=0.0);
    // TODO: implement -p for E,t, -r for refinement
    int argc;
    int num_len;
    int num_wid;
    double len;
    double wid;
    std::string inputPath;
    std::string outputPath;
};

inline Arguments::Arguments(int t_argc)
    : argc(t_argc), num_len(0), num_wid(0), len(0.0), wid(0.0),
    //   inputPath("/Users/chenjingyu/Dropbox/Research/Codes/plates-shells/"),
    //   outputPath("/Users/chenjingyu/Dropbox/Research/Codes/plates-shells/results/")
      inputPath("/home/sci04/Codes/plates-shells/"),
      outputPath("/home/sci04/Codes/plates-shells/results/")
{}

inline void Arguments::getArgs(int t_nl, int t_nw, double t_l, double t_w) {
    num_len = t_nl;
    num_wid = t_nw;
    len = t_l;
    wid = t_w;
}

#endif	// PLATES_SHELLS_ARGUMENTS_H