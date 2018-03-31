#include <getopt.h>
#include <thread>
#include <omp.h>
#include "kfreq.h"

using namespace kf;

void usage(char **argv) {
    std::fprintf(stderr, "Usage: %s [max_kmer_size] [genome1] [genome2] ...\n"
                         "Flags:\n"
                         "-k\tSet kmer size [4]\n"
                         "-b\tEmit binary [false]\n"
                         "-p\tSet number of threads [1]. Using -1 will result in all available cores being used\n"
                 , *argv);
    std::exit(EXIT_FAILURE);
}

std::string canonicalize(const char *str) {
    if(const char *slash = std::strrchr(str, '/'))
        return slash + 1;
    return str;
}

int main(int argc, char *argv[]) {
    if(argc == 1) usage(argv);
    std::vector<std::string> paths;
    bool emit_binary = false;
    unsigned ks = 4;
    int c, nthreads = 1;
    while((c = getopt(argc, argv, "bk:")) >= 0) {
        switch(c) {
            case 'k': ks = std::atoi(optarg); break;
            case 'b': emit_binary = true; break;
            case 'p': nthreads = std::atoi(optarg); break;
        }
    }
    if(nthreads < 0)
        nthreads = std::thread::hardware_concurrency();
    if(ks > 16 || ks < 2) {
        std::fprintf(stderr, "ks: %u. Max supported: 16. Min: 2\n", ks);
        usage(argv);
    }
    for(char **p(argv + optind); *p; paths.emplace_back(*p++));
    if(nthreads == 1) {
        kseq_t seq = kseq_init_stack();
        freq::KFC kc(ks);
        for(const auto &path: paths) {
            std::string outpath = canonicalize(path.data()) + ".k" + std::to_string(ks) + ".bin";
            kc.add(path.data(), &seq);
            kc.write(outpath.data(), emit_binary);
            kc.clear();
        }
        kseq_destroy_stack(seq);
    } else {
        std::vector<kseq_t> kseqs; kseqs.reserve(nthreads); while(kseqs.size() < (unsigned)nthreads) kseqs.emplace_back(kseq_init_stack());
        std::vector<freq::KFC> kfcs; kfcs.reserve(nthreads); while(kfcs.size() < (unsigned)nthreads) kfcs.emplace_back(ks);
        std::vector<std::string> pathnames; for(const auto &path: paths) pathnames.emplace_back(canonicalize(path.data()) + ".k" + std::to_string(ks) + ".bin");
        #pragma omp parallel for
        for(unsigned i = 0; i < paths.size(); ++i) {
            auto &kfc = kfcs[omp_get_thread_num()];
            kfc.add(paths[i].data(), kseqs.data() + omp_get_thread_num());
            kfc.write(pathnames[i].data(), emit_binary);
            kfc.clear();
        }
        for(auto &ks: kseqs) kseq_destroy_stack(ks);
    }
}
