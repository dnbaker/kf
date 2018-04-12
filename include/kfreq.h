#pragma once
#include "kmerutil.h"
#include <numeric>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

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

template<typename FloatType, typename=typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
double pearsonr_naive(const std::vector<FloatType> &v1, const std::vector<FloatType> &v2) {
    // This is far from the fastest way to perform this task. See https://github.com/dnbaker/vec/blob/master/stats.h for an accelerated pearsonr implementation.
    const auto m1(std::accumulate(v1.cbegin(), v1.cend(), static_cast<FloatType>(0))),
               m2(std::accumulate(v2.cbegin(), v2.cend(), static_cast<FloatType>(0)));
    FloatType sd = 0., s1s = 0., s2s = 0., val1, val2;
    for(size_t i(0); i < v1.size(); ++i) {
        val1 = v1[i] - m1, val2 = v2[i] - m2;
        s1s += val1 * val1, s2s += val2 * val2, sd  += val1 * val2;
    }
    const FloatType rden = std::sqrt(s1s) * std::sqrt(s2s); // two square roots for better floating-point accuracy.
    return std::min(std::max(sd / rden, FloatType(-1.)), FloatType(1.));
}

enum Mode {
    TEXT   = 0,
    BINARY = 1,
    DETECT = 2
};

static const char KF_BIN [] {'#', 'k', 'f', 'b', 'i', 'n', '\n'};
static const char KF_TEXT [] {'#', 'k', 'f', 't', 'x', 't', '\n'};
static const char KFL_BIN [] {'#', 'k', 'f', 'l', 'b', 'i', 'n', '\n'};
static const char KFL_TEXT [] {'#', 'k', 'f', 'l', 't', 'x', 't', '\n'};

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
// From Jellyfish, by Guillaume Marcais
static INLINE u32 reverse_complement(u32 kmer, uint8_t n) {
    kmer = ((kmer >> 2)  & 0x33333333U) | ((kmer & 0x33333333U) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0FU) | ((kmer & 0x0F0F0F0FU) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FFU) | ((kmer & 0x00FF00FFU) << 8);
    kmer = ( kmer >> 16               ) | ( kmer               << 16);
    return (((u32)-1) - kmer) >> (CHAR_BIT * sizeof(kmer) - (n << 1));
}

template<typename SizeType, typename=typename std::enable_if<std::is_integral<SizeType>::value && std::is_unsigned<SizeType>::value>::type>
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
    void rc_collapse() {
        u32 k, rc;
        for(k = 0; k < k_; ++k) {
            if((rc = reverse_complement(k, k_)) < k)
                data_[rc] += data_[k], data_[k] = 0;
            else
                data_[k] += data_[rc], data_[rc] = 0;
        }
    }
};

template<typename KFType>
void rc_collapse(KFType &kf) {
    for(auto &fs: kf.freqs())
        fs.rc_collapse();
}

// Counts short kmer occurrences using arrays. (Supported: up to 16)
template<typename SizeType, typename=typename std::enable_if<std::is_integral<SizeType>::value && std::is_unsigned<SizeType>::value>::type>
class KFreqArray {
    unsigned maxk_;
    std::vector<SubKFreq<SizeType>> freqs_;
    using FreqType = std::vector<SubKFreq<SizeType>>;
public:
    FreqType       &freqs()       {return freqs_;}
    const FreqType &freqs() const {return freqs_;}
    using size_type = SizeType;
    KFreqArray(unsigned k): maxk_(k) {
        if(std::numeric_limits<SizeType>::max() < (1ull << (k << 1)))
            throw std::runtime_error(std::string("SizeType with width ") + std::to_string(sizeof(SizeType) * CHAR_BIT) + " is not long enough for k = " + std::to_string(maxk_));
        while(freqs_.size() < maxk_) freqs_.emplace_back(freqs_.size() + 1);
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
#define __kmask32(k) (UINT32_C(-1) >> (32 - ((k) << 1)))
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
    SizeType count(const std::string &str) const {
        return count(str.size(), str2kmer<SizeType>(str));
    }
    SizeType count(unsigned k, SizeType value) const {
        return freqs_[k - 1].data_[value];
    }
    unsigned maxk() const {return maxk_;}
};

// Counts short kmer occurrences using arrays. (Supported: up to 16)
template<typename SizeType, typename=typename std::enable_if<std::is_integral<SizeType>::value && std::is_unsigned<SizeType>::value>::type>
class KFreqList {
    const uint16_t maxk_;
    const uint16_t   nk_;
    std::vector<SubKFreq<SizeType>> freqs_;
    using FreqType = std::vector<SubKFreq<SizeType>>;
public:
    FreqType       &freqs()       {return freqs_;}
    const FreqType &freqs() const {return freqs_;}
    using size_type = SizeType;
    KFreqList(unsigned k, unsigned num_kmers=3): maxk_(k), nk_(num_kmers) {
        if(std::numeric_limits<SizeType>::max() > (1ull << (k << 1)))
            throw std::runtime_error(std::string("SizeType with width ") + std::to_string(sizeof(SizeType) * CHAR_BIT) + " is not long enough for k = " + std::to_string(maxk_));
        for(k = maxk_ - nk_; k < maxk_;freqs_.emplace_back(k+++1));
    }
    KFreqList(const char *path) {
        gzFile fp = gzopen(path, "rb");
        if(!fp) throw std::runtime_error("Could not open file.");
        char buf[sizeof(KF_BIN)];
        gzread(fp, buf, sizeof(buf));
        bool read_binary;
        if(std::string(buf) == KFL_BIN) read_binary == true;
        else if(std::string(buf) == KFL_TEXT) read_binary = false;
        else {
            buf[sizeof(buf) - 1] = '\0';
            throw std::runtime_error(std::string("Unexpected magic string: ") + buf);
        }
        if(read_binary) {
            gzread(fp, &maxk_, sizeof(maxk_));
            gzread(fp, &nk_, sizeof(nk_));
            for(unsigned k = maxk_ - nk_; k < maxk_;freqs_.emplace_back(k+++1));
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
            auto it = freqs_.begin();
            for(unsigned i(maxk_ - nk_); i < maxk_;++i) {
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
            for(auto it(freqs_.begin()); it < freqs_.end(); ++it) {
                //std::fprintf(stderr, "Now doing stuff at pos %zu of %zu\n", std::distance(freqs_.begin(), it), freqs_.size());
                auto &sf(*it);
                sf.v_ <<= 2;
                sf.v_ |= cc;
                sf.v_ &= __kmask32(sf.k_);
                sf.f_ += (sf.f_ != sf.k_ - 1);
                sf.data_[sf.v_] += (sf.f_ == sf.k_ - 1);
            }
            //std::fprintf(stderr, "Done incrementing counts\n");
        }
    }
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
            gzwrite(fp, (void *)KFL_BIN, sizeof(KFL_BIN));
            gzwrite(fp, (void *)&maxk_, sizeof(maxk_));
            gzwrite(fp, (void *)&nk_, sizeof(nk_));
            for(const auto &freq: freqs_) freq.write(fp);
        } else {
            gzwrite(fp, (void *)KFL_TEXT, sizeof(KF_TEXT));
            gzprintf(fp, "#Max k: %u\n", maxk_);
            gzprintf(fp, "#nk: %u\n", nk_);
            for(const auto &sf: freqs_) {
                gzprintf(fp, "%u: [", sf.k_);
                for(size_t i(0); i < sf.data_.size() - 1; gzprintf(fp, "%zu|", size_t(sf.data_[i++])));
                gzprintf(fp, "%zu]\n", size_t(sf.data_.back()));
            }
        }
        gzclose(fp);
    }
    SizeType count(const std::string &str) const {
        return count(str.size(), str2kmer<SizeType>(str));
    }
    SizeType count(unsigned k, SizeType value) const {
        return freqs_[k - (maxk_ - nk_ + 1)].data_[value];
    }
    unsigned maxk() const {return maxk_;}
};
using KFC = KFreqArray<u32>;
using KFL = KFreqList<u32>;

template<typename KFType, typename FloatType=double,
         typename=typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
std::vector<FloatType> calc_zscores(const KFType &kf) {
#if 0
    Simple port of the below code:
def calc_zscores(c, k):
    # for each kmer, calc expected frequencies with a
    # maximal order Markov chain model (k, k-1, k-2)
    mer_z = []
    for mer in gen_kmers(k): # even if ambig nucs in Counters, we skip them
        mid = c[k-2][mer[1:k-1]]
        try:
            exp = 1. * c[k-1][mer[:k-1]] * c[k-1][mer[1:]] / mid
            std = sqrt(exp * (mid-c[k-1][mer[:k-1]]) * (mid-c[k-1][mer[1:]]) / (mid**2))
            mer_z.append((c[k][mer]-exp) / std)
        except ZeroDivisionError:
            if mid > 0: # if a (k-1)mer is absent
                mer_z.append(1 / (mid**2))
            elif mid == 0: # if a (k-2)"mid"mer is absent
                mer_z.append(0)
    return mer_z
#endif
    const unsigned k = kf.maxk();
    std::vector<FloatType> ret;
    ret.reserve(1u << (k << 1));
    using sz_t = typename KFType::size_type;
    sz_t mid;
    u32 km1l, km1r;
    FloatType xp, std, fmid;
    for(unsigned i(0), max_kmer(1u << (k << 1)); i < max_kmer; ++i) {
        if((mid = kf.count(k - 2, (i >> 2) & __kmask32(k - 2))) == 0) {
            ret.push_back(0.); continue;
        }
        km1l = kf.count(k - 1, i & __kmask32(k - 1));
        km1r = kf.count(k - 1, (i>>2) & __kmask32(k - 1));
        fmid = 1. / mid;
        xp = static_cast<FloatType>(km1l * km1r) * fmid;
        if((std = xp * (mid - km1l) * (mid - km1r) * fmid * fmid) == 0.)
            ret.push_back(1./(static_cast<FloatType>(mid) * static_cast<FloatType>(mid)));
        else ret.push_back((kf.count(k, i) - xp) / std);
    }
    assert(ret.size() == (1u << (k << 1)));
    return ret;
}

} // namespace freq

#undef __kmask32
} // namespace kf
