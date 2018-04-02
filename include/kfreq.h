#pragma once
#include "kmerutil.h"
#include <vector>
#include <cstring>
#include <stdexcept>
#include <limits>

namespace kf {

static inline void kseq_assign(kseq_t *ks, gzFile fp) {
    if(!ks->f) {
        ks->f = (kstream_t*)calloc(1, sizeof(kstream_t));
        ks->f->buf = (unsigned char*)malloc(KSTREAM_SIZE);
    } else ks->f->is_eof = ks->f->begin = ks->f->end = 0;
    ks->f->f = fp;
}

static inline kseq_t kseq_init_stack() {
    kseq_t ret;
    std::memset(&ret, 0, sizeof(ret));
    return ret;
}
static inline void kseq_destroy_stack(kseq_t &ks) {
    free(ks.name.s); free(ks.comment.s); free(ks.seq.s); free(ks.qual.s);
    ks_destroy(ks.f);
}



namespace freq {

enum Mode {
    TEXT   = 0,
    BINARY = 1,
    DETECT = 2
};

static const char KF_BIN [] {'#', 'k', 'f', 'b', 'i', 'n', '\n'};
static const char KF_TEXT [] {'#', 'k', 'f', 't', 'x', 't', '\n'};

template<typename SizeType=std::size_t>
inline SizeType str2kmer(const std::string &str) {
    const unsigned k = str.size();
    auto it(str.cbegin());
    SizeType v = cstr_lut[*it++];
    do {
        v <<= 2;
        if((v |= cstr_lut[*it++]) == SizeType(-1))
            throw std::runtime_error("Illegal character!");
    } while(it < str.cend());
    return v;
}

// Counts short kmer occurrences using arrays. (Supported: up to 16)
template<typename SizeType, typename=typename std::enable_if<std::is_integral<SizeType>::value && std::is_unsigned<SizeType>::value>::type>
class KFreqArray {
    struct SubKFreq {
        const unsigned k_; // Kmer size
        u32            v_; // Value
        uint8_t        f_; // How full?
        std::vector<SizeType> data_;
        SubKFreq(unsigned k): k_(k), v_(0), f_(0), data_(1ull << (k << 1)) {
        }
        void clear_kmer() {
            f_ = v_ = 0;
        }
        void clear() {
            clear_kmer();
            std::fill(std::begin(data_), std::end(data_), 0);
        }
        void write(gzFile fp) const {
            gzwrite(fp, (void *)data_.data(), data_.size() * sizeof(SizeType));
#if !NDEBUG
            std::fprintf(stderr, "For k = %u:", k_);
            for(const auto &el: data_) {
                std::fprintf(stderr, "%u|", (unsigned)el);
            }
            std::fputc('\n', stderr);
#endif
        }
    };
    unsigned maxk_;
    std::vector<SubKFreq> freqs_;
public:
    KFreqArray(unsigned k): maxk_(k) {
        if(std::numeric_limits<SizeType>::max() < (1ull << (k << 1))) {
            throw std::runtime_error(std::string("SizeType with width ") + std::to_string(sizeof(SizeType) * CHAR_BIT) + " is not long enough for k = " + std::to_string(maxk_));
        }
        while(freqs_.size() < maxk_) {
            freqs_.emplace_back(freqs_.size() + 1);
        }
        //if(maxk_ != 4) throw std::runtime_error("I'm making it for only k == 4 for now because I'm lazy.");
    }
    KFreqArray(const char *path) {
        gzFile fp = gzopen(path, "rb");
        if(!fp) throw std::runtime_error("Could not open file.");
        char buf[sizeof(KF_BIN)];
        gzread(fp, buf, sizeof(buf));
        bool read_binary;
        if(std::string(buf) == KF_BIN) read_binary == true;
        else if(std::string(buf) == KF_TEXT) read_binary = false;
        else {
            buf[sizeof(buf) - 1] = '\0';
            throw std::runtime_error(std::string("Unexpected magic string: ") + buf);
        }
        if(read_binary) {
            gzread(fp, &maxk_, sizeof(maxk_));
            while(freqs_.size() < maxk_) freqs_.emplace_back(freqs_.size() + 1);
            for(auto &sf: freqs_) gzread(fp, sf.data_.data(), sf.data_.size() * sizeof(SizeType));
        } else {
            char *line, *p;
            std::vector<char> linebuf(256);
            if((line = gzgets(fp, linebuf.data(), linebuf.size())) == nullptr) throw std::runtime_error("Could not read from file.");
            for(p = line; *p; ++p);
            while(!std::isdigit(*p)) --p; // In case the newline is attached
            while(std::isdigit(*p)) --p;
            maxk_ = std::atoi(p + 1);
            while(freqs_.size() < maxk_) freqs_.emplace_back(freqs_.size() + 1);
            for(unsigned i(0); i < maxk_;++i) {
                auto &freq = freqs_[i];
                if((line = gzgets(fp, linebuf.data(), linebuf.size())) == nullptr) throw std::runtime_error("Could not read from file.");
                if(!(p = std::strchr(line, '['))) goto fail;
                unsigned j(0);
                do {
                    freq.data_[j++] = static_cast<SizeType>(std::strtoull(++p, nullptr, 10));
                    while(std::isdigit(*p)) ++p;
                } while(j < freq.data_.size());
            }
            if(false) {
                fail:
                throw std::runtime_error("Error in parsing.");
            }
        }
        gzclose(fp);
    }
    void clear_kmers() {
        for(auto &freq: freqs_) freq.clear_kmer();
    }
#define __kmask32(k) (UINT32_C(-1) >> (32 - (k << 1)))
    void process_seq(const char *s, size_t l) {
        clear_kmers();
        u32 cc;
        size_t i = 0;
        while(i < l) {
            if((cc = cstr_lut[s[i++]]) == UINT32_C(-1)) {
                clear_kmers();
                //std::fprintf(stderr, "i is now %zu. Kmers are cleared. Continue.\n", i);
                continue;
            }
            //std::fprintf(stderr, "Now incrementing counts\n");
            ++freqs_[0].data_[cc];
            for(auto it(freqs_.begin() + 1); it < freqs_.end(); ++it) {
                //std::fprintf(stderr, "Now doing stuff at pos %zu of %zu\n", std::distance(freqs_.begin(), it), freqs_.size());
                auto &sf(*it);
                sf.v_ <<= 2;
                sf.v_ |= cc;
                if(sf.f_ == sf.k_ - 1) {
                    sf.v_ &= __kmask32(sf.k_);
                    ++sf.data_[sf.v_];
                } else ++sf.f_;
            }
            //std::fprintf(stderr, "Done incrementing counts\n");
        }
    }
#undef __kmask32
    void add(const char *path, kseq_t *ks=nullptr) {
        const bool destroy = (ks == nullptr);
        gzFile fp(gzopen(path, "rb"));
        if(destroy) ks = kseq_init(fp);
        else       kseq_assign(ks, fp);

        while(kseq_read(ks) >= 0) process_seq(ks->seq.s, ks->seq.l);

        if(destroy) kseq_destroy(ks);
        gzclose(fp);
    }
    void clear() {
        for(auto &freq: freqs_)
            freq.clear();
    }
    void write(const char *path, bool emit_binary=false) {
        gzFile fp = gzopen(path, "wb");
        if(fp == nullptr) throw std::runtime_error("Could not open file for output.");
        if(emit_binary) {
            gzwrite(fp, (void *)KF_BIN, sizeof(KF_BIN));
            gzwrite(fp, (void *)&maxk_, sizeof(maxk_));
            for(const auto &freq: freqs_) {
                freq.write(fp);
            }
        } else {
            gzwrite(fp, (void *)KF_TEXT, sizeof(KF_TEXT));
            gzprintf(fp, "#Max k: %u\n", maxk_);
            for(const auto &sf: freqs_) {
                gzprintf(fp, "%u: [", sf.k_);
                for(size_t i(0); i < sf.data_.size() - 1; gzprintf(fp, "%zu|", size_t(sf.data_[i++])));
                gzprintf(fp, "%zu]\n", size_t(sf.data_.back()));
            }
        }
        gzclose(fp);
    }
    SizeType count(const std::string &str) {
        return count(str.size(), str2kmer<SizeType>(str));
    }
    SizeType count(unsigned k, SizeType value) {
        return freqs_[k - 1].data_[value];
    }
};


using KFC = KFreqArray<u32>;

}

}
