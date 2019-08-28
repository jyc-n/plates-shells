#ifndef PLATES_SHELLS_ARGUMENTS_H
#define PLATES_SHELLS_ARGUMENTS_H

#include <string>
// #include <filesystem>

enum Opts {
    DEFAULT = 1,
    OPTIONS = 4,
};

class Arguments {

public:
    Arguments(int t_argc, char* argv[]);
    
    int type;                       // 0-default, 1-paramter, 2-refinement
    std::string flag;               // -p for parameter, -r for refinement
    std::string data1;              // -p for E, -r for #nodes on length
    std::string data2;              // -p for t, -r for #nodes on width
    std::string inputPath;
    std::string outputPath;
};

inline Arguments::Arguments(int t_argc, char* argv[])
{
    if (t_argc == DEFAULT) {
        type = 0;
    }
    else if (t_argc == OPTIONS) {
        flag  = argv[1];
        if (flag == "-p") {
            type = 1;
        }
        else if(flag == "-r") {
            type = 2;
        }
        else {
            throw "wrong options!";
        }
        data1 = argv[2];
        data2 = argv[3];
    }
    else {
        throw "wrong number of arguments!";
    }

    // std::string tmp = std::filesystem::current_path().string();
    // std::string tmp = "/Users/rla/Dropbox/Research/Codes/plates-shells";
    std::string tmp = "/Users/chenjingyu/Dropbox/Research/Codes/plates-shells";
    // std::string tmp = "/home/jingyuchen/code/plates-shells";
    inputPath = tmp+"/";
    outputPath = tmp+"/results/";
    std::cout << "cwd is " << inputPath << std::endl;
}

#endif	// PLATES_SHELLS_ARGUMENTS_H