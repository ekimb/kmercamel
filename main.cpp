#include "greedy.cpp"
#include "generalized_simplitigs.cpp"
#include "generalized_simplitigs_ac.cpp"
#include "kmers.cpp"

#include <iostream>
#include <string>

#include "bioio.hpp"

void WriteSuperstring(KMerSet result, std::string name) {
    std::cout << name << std::endl;
    std::string superstring = "";
    for (int i = 0; i < result.superstring.length(); ++i) {
        superstring += result.mask[i] ? result.superstring[i] : std::tolower(result.superstring[i]);
    }
    std::cout << superstring << std::endl;
}

void WriteStats(KMerSet result, std::vector<KMer> kMers, std::string data, std::string name) {
    std::cout << "name:                       " << name << std::endl;
    std::cout << "superstring length:         " << result.superstring.length() << std::endl;
    std::cout << "k-mers count:               " << kMers.size() << std::endl;
    std::cout << "length of scanned sequence: " << data.length() << std::endl;
    std::cout << "coefficient:                " << result.superstring.length() / (double)kMers.size() << std::endl;
    std::cout << "========================================="  << std::endl;
}

bool HasFlag(char **begin, char **end, std::string param) {
    return std::find(begin, end, param) != end;
}

std::string GetFlag(char **begin, char **end, std::string param, std::string def) {
    auto flag = std::find(begin, end, param);
    if (flag != end && ++flag != end) return *flag;
    return def;
}

int GetFlagAsInt(char **begin, char **end, std::string param, int def) {
    auto ret = GetFlag(begin, end, param, "");
    if (ret == "") return def;
    return std::stoi(ret);
}

void Help() {
    std::cerr << "Example usage:       ./kmers path_to_fasta -k 13 -d 5 -a greedy" << std::endl;
    std::cerr << "Possible algorithms: greedy pseudosimplitigs pseudosimplitigsAC" << std::endl;
}

int main(int argc, char **argv) {
    std::string path;
    int k;
    int d_max;
    std::string algorithm;
    bool printStats = false;
    try {
        path = argv[1];
        k = GetFlagAsInt(argv, argv + argc, "-k", 13);
        d_max = GetFlagAsInt(argv, argv + argc, "-d", 5);
        algorithm = GetFlag(argv, argv + argc, "-a", "greedy");
        printStats = HasFlag(argv, argv + argc, "-s");
    } catch (std::exception) {
        Help();
        return 1;
    };
    auto data = bioio::read_fasta(path);
    for (auto record : data) {
        auto kMers = ConstructKMers(record.sequence, k);
        KMerSet result;
        if (algorithm == "greedy")
            result = Greedy(kMers);
        else if (algorithm == "pseudosimplitigs")
            result = GreedyGeneralizedSimplitigs(kMers, k, d_max);
        else if (algorithm == "pseudosimplitigsAC")
            result = GreedyGeneralizedSimplitigsAC(kMers, k, d_max);
        else {
            Help();
            return 1;
        }
        if (printStats)
            WriteStats(result, kMers, record.sequence, record.name);
        else
            WriteSuperstring(result, record.name);
        return 0;
    }
}
