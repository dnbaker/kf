#include <getopt.h>
#include "kfreq.h"

using namespace kf;

void usage(char **argv) {
    std::fprintf(stderr, "Usage: %s [max_kmer_size] [genome1] [genome2] ...\n", *argv);
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
    int c;
    while((c = getopt(argc, argv, "bk:")) >= 0) {
        switch(c) {
            case 'k': ks = std::atoi(optarg); break;
            case 'b': emit_binary = true; break;
        }
    }
    if(ks > 16 || ks < 2) {
        std::fprintf(stderr, "ks: %u. Max supported: 16. Min: 2\n", ks);
        usage(argv);
    }
    for(char **p(argv + optind); *p; paths.emplace_back(*p++));
    kseq_t seq = kseq_init_stack();
    freq::KFC kc(ks);
    for(const auto &path: paths) {
        std::string outpath = canonicalize(path.data()) + ".k" + std::to_string(ks) + ".bin";
        kc.add(path.data(), &seq);
        kc.write(outpath.data(), emit_binary);
        kc.clear();
    }
    kseq_destroy_stack(seq);
}
