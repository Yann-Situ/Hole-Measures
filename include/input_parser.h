/* input_parser.h
 * Author : iain
 * Input parser code found here:
 * https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
 */
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include <vector>
#include <algorithm>

/// @author iain
class InputParser
{
private:
    std::vector <std::string> tokens;

public:
    InputParser (int &argc, char **argv) {
        for (int i=1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }

    const std::string& getCmdOption(const std::string &option) const {
    std::vector<std::string>::const_iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }
    static const std::string empty_string("");
    return empty_string;
}

bool cmdOptionExists(const std::string &option) const {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
    != this->tokens.end();
}

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
