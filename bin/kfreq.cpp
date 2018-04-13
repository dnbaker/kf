#include <getopt.h>
#include <thread>
#include <omp.h>
#include "kfreq.h"

#ifndef FLOAT_TYPE
#define FLOAT_TYPE double
#endif

using namespace kf;

void usage(char **argv) {
    std::fprintf(stderr, "Usage: %s [flags] [genome1] [genome2] ...\n"
                         "Flags:\n"
                         "-k\tSet kmer size [4]\n"
                         "-b\tEmit binary [false]\n"
                         "-p\tSet number of threads [1]. Using -1 will result in all available cores being used\n"
                         "-o\tSet output file for distance table, if produced.\n"
                         "-c\tSketch only, don't calculate distances.\n"
                         "-R\tDo not reverse complement. [Default: always reverse complement.]\n"
                 , *argv);
    std::fflush(stderr);
    std::exit(EXIT_FAILURE);
}

std::string canonicalize(const char *str) {
    if(const char *slash = std::strrchr(str, '/'))
        return slash + 1;
    return str;
}

std::vector<std::vector<FLOAT_TYPE>> pairwise_pearson(const std::vector<std::vector<FLOAT_TYPE>> &profiles) {
    std::vector<std::vector<FLOAT_TYPE>> ret;
    ret.reserve(profiles.size());
    while(ret.size() < profiles.size()) ret.emplace_back(profiles.size());
    for(size_t i(0); i < profiles.size(); ++i)
        for(size_t j(i); j < profiles.size(); ++j)
            ret[i][j] = ret[j][i] = freq::pearsonr_naive(profiles[i], profiles[j]);
    return ret;
}

void print_distmat(std::FILE *ofp, const std::vector<std::vector<FLOAT_TYPE>> &profiles, const std::vector<std::string> &paths) {
    std::fprintf(ofp, "#Path");
    for(const auto &path: paths) std::fprintf(ofp, "\t%s", path.data());
    std::fputc('\n', ofp);
    for(size_t i(0); i < profiles.size(); ++i) {
        std::fputs(paths[i].data(), ofp);
        for(const auto &val: profiles[i]) std::fprintf(ofp, "\t%f", val);
        std::fputc('\n', ofp);
    }
}

void emit_zscores(const std::string &path, const std::vector<FLOAT_TYPE> &data) {
    std::FILE *ofp = std::fopen(path.data(), "wb");
    if(ofp == nullptr) throw std::runtime_error(std::string("Could not open file at ") + path);
    for(const auto &val: data) std::fprintf(ofp, "%f\n", val);
    std::fclose(ofp);
}

int main(int argc, char *argv[]) {
    if(argc == 1) usage(argv);

    using KFType = freq::KFC;
    std::vector<std::string> paths;
    bool rc = true, calculate_distances = true;
    unsigned ks = 4;
    int c, nthreads = 1;
    std::FILE *ofp = stdout;
    while((c = getopt(argc, argv, "Rcbo:k:h?")) >= 0) {
        switch(c) {
            case 'o': ofp = std::fopen(optarg, "wb"); break;
            case 'k': ks = std::atoi(optarg); break;
            case 'p': nthreads = std::atoi(optarg); break;
            case 'c': calculate_distances = false; break;
            case 'R': rc = false; break;
            case 'h': case '?': usage(argv);
        }
    }
    if(nthreads < 0) nthreads = std::thread::hardware_concurrency();
    omp_set_num_threads(nthreads);
    if(ks > 16 || ks < 2) {
        std::fprintf(stderr, "ks: %u. Max supported: 16. Min: 2. Setting to 4.\n", ks);
        usage(argv);
    }
    for(char **p(argv + optind); *p; paths.emplace_back(*p++));
    std::vector<std::vector<FLOAT_TYPE>> profiles;
    if(calculate_distances) profiles.resize(paths.size());
    std::vector<KFType> kfcs; kfcs.reserve(nthreads);
    std::vector<kseq_t> kseqs; kseqs.reserve(nthreads);
    while(kseqs.size() < (unsigned)nthreads) kseqs.emplace_back(kseq_init_stack());
    while(kfcs.size() < (unsigned)nthreads) kfcs.emplace_back(ks);
    #pragma omp parallel for
    for(unsigned i = 0; i < paths.size(); ++i) {
        const auto tid =  omp_get_thread_num();
        assert(tid <= kfcs.size());
        auto &kfc = kfcs[tid];
        kfc.add(paths[i].data(), kseqs.data() + tid);
        if(rc) rc_collapse(kfc);
        auto zs = calc_zscores(kfc);
        emit_zscores(canonicalize(paths[i].data()) + ".k" + std::to_string(ks) + ".txt", zs);
        if(calculate_distances) profiles[i] = zs;
        kfc.clear();
    }
    for(auto &ks: kseqs) kseq_destroy_stack(ks);
    if(calculate_distances) {
        std::fprintf(stderr, "calculating distances\n");
        print_distmat(ofp, pairwise_pearson(profiles), paths);
    }
    if(ofp != stdout) std::fclose(ofp);
}
