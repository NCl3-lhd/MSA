#ifndef STARALIGN_H
#define STARALIGN_H
#include <vector>
#include <array>
#include <string>
int findCenterSeq(const std::vector<std::string>&);
std::vector<std::array<std::string, 2>> getAlignedSeqs(int, const std::vector<std::string>&);
void getMarkInsertion(const std::vector<std::array<std::string, 2>>&, std::vector<int>&);
std::string insertGap(const std::vector<int>&, const std::string&);
std::vector<std::string> insertGapToAlignedSeqs(const std::vector<std::string>&, int, const std::vector<std::array<std::string, 2>>&, const std::vector<int>&);
void starAlign(const std::vector<std::string>&, std::vector<std::string>&);
#endif