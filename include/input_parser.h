/* input_parser.h
 * Author : iain
 * Input parser code found here:
 * https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
 */
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include <vector>

class InputParser
{
    public:
        InputParser (int &argc, char **argv);
        const std::string& getCmdOption(const std::string &option) const;
        bool cmdOptionExists(const std::string &option) const;
    private:
        std::vector <std::string> tokens;
};

/* Usage :
int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        // Do stuff
    }
    const std::string &filename = input.getCmdOption("-f");
    if (!filename.empty()){
        // Do interesting things ...
    }
    return 0;
}
*/
#endif
