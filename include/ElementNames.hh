#ifndef ElementNames_HH
#define ElementNames_HH

#include <string>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

class ElementNames
{
    public:
        ElementNames();
        virtual ~ElementNames();
        static void ClearStore();
        static void SetElementNames();
        static string GetName(int Z)
        {
            return elementName[Z];
        }
        static bool CheckName(string name);
        static bool CheckName(string name, int Z);
        static double StringPercentMatch(string base, string secondary);
        static string ExtractString(std::stringstream &stream, char delim, int outType=7);
        static string *elementName;
    protected:
    private:

};

#endif // ElementNames_HH
