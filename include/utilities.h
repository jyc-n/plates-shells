#ifndef PLATES_SHELLS_UTILITIES_H
#define PLATES_SHELLS_UTILITIES_H

#include <chrono>

class Timer
{
    public:
        Timer(bool LONGTIME = false)
        {
            start(LONGTIME);
        }
        void start(bool LONGTIME = false)
        {
            if (LONGTIME)
                m_ltime = std::chrono::system_clock::now();
            else
                m_time = std::chrono::high_resolution_clock::now();
        }
        double elapsed(bool LONGTIME = false) const
        {
            if (LONGTIME)
            {
                std::chrono::duration<double> diff =
                        std::chrono::system_clock::now() - m_ltime;
                return diff.count();
            }
            else
            {
                std::chrono::duration<double,std::milli> diff =
                        std::chrono::high_resolution_clock::now() - m_time;
                return diff.count();
            }
        }
    private:
        std::chrono::high_resolution_clock::time_point m_time;
        std::chrono::system_clock::time_point m_ltime;
};

#endif //PLATES_SHELLS_UTILITIES_H
