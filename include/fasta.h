#ifndef FASTA_H
#define FASTA_H

#include<vector>
#include<string>
void readFasta(const std::string&, std::vector<std::string>&, std::vector<std::string>&);
void writeFasta(const std::string&, std::vector<std::string>&, std::vector<std::string>&);
#endif
