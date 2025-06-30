#ifndef  FASTA_H
#define FASTA_H

#include<vector>
#include<string>
void readFasta(const std::string& path, std::vector<std::string>& seqs, std::vector<std::string>& labels);
void writeFasta(const std::string& path, std::vector<std::string>& seqs, std::vector<std::string>& labels);
#endif
