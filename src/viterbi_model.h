#pragma once

#include <algorithm>
#include <cctype>
#include <vector>
#include <string>
#include <array>
#include <cmath>

namespace uniann {

// Coarse states: N, E0, E1, E2, I0, I1, I2.
constexpr int NUM_STATES = 7;

inline const std::array<std::string, NUM_STATES> state_name = {
    "N", "E0", "E1", "E2", "I0", "I1", "I2"
};

constexpr double NEG_INF = -1e9;
constexpr int MIN_INTRON = 40;
constexpr int MIN_EXON   = 3;
constexpr int MIN_INTER  = 30;
constexpr int MIN_SINGLE = 100;

inline bool is_exon(int s) {
    return (s == 1 || s == 2 || s == 3);
}

inline bool is_intron(int s) {
    return (s == 4 || s == 5 || s == 6);
}

struct PathMetadata {
    int intron_len = 0;
    int exon_len = 0;
    int inter_len = 0;
    int exon_from = 0;
    // The start-codon rule requires the N predecessor to be N itself.
    int predecessor_state = -1;
};

// Emission columns: N, E0, E1, E2 and one intron column shared by I0-I2.
constexpr int EMISSION_COLUMNS = 5;
using EmissionRow = std::array<float, EMISSION_COLUMNS>;

inline int emission_column(int state) {
    return is_intron(state) ? EMISSION_COLUMNS - 1 : state;
}

// Motif starting at a 0-based sequence position. The four kinds differ in
// their first two bases, so a position carries at most one of them.
enum class SiteMotif { None, Donor, Acceptor, Start, Stop };

inline SiteMotif motif_starting_at(const std::vector<char> &seq, int start) {
    const int length = seq.size();
    if (start < 0 || start + 1 >= length) return SiteMotif::None;
    const char first = std::toupper(seq[start]);
    const char second = std::toupper(seq[start + 1]);
    if (first == 'G' && second == 'T') return SiteMotif::Donor;
    if (first == 'A' && second == 'G') return SiteMotif::Acceptor;
    if (start + 2 >= length) return SiteMotif::None;
    const char third = std::toupper(seq[start + 2]);
    if (first == 'A' && second == 'T' && third == 'G') return SiteMotif::Start;
    if (first == 'T' && ((second == 'A' && (third == 'A' || third == 'G')) ||
                         (second == 'G' && third == 'A')))
        return SiteMotif::Stop;
    return SiteMotif::None;
}

struct ModelInputs {
    const std::vector<EmissionRow> &emit;
    // One score per position, keyed by the motif the sequence carries there:
    // GT or AG at the dinucleotide start, ATG or a stop codon at the codon
    // start. The transition rules only read a score where the sequence has the
    // corresponding motif, so four dense arrays collapse into one without
    // changing any decision.
    const std::vector<double> &site_score;
    const std::vector<char> &seq;
    const std::vector<std::vector<double>> &trans;
};

struct TransitionResult {
    bool allowed = false;
    double transition_score = NEG_INF;
    double emission_score = NEG_INF;
    PathMetadata next_metadata;
};

inline int splice_motif_start_0based(int dp_position) {
    return dp_position - 1;
}

inline int codon_start_0based(int dp_position) {
    return std::max(0, dp_position - 2);
}

inline bool is_stop_codon(const std::string &codon) {
    return codon == "TAA" || codon == "TAG" || codon == "TGA";
}

// Sequence motifs at one DP position. Every transition out of a position
// inspects the same bases, so callers compute this once per position instead
// of once per (from, to) pair.
struct SiteContext {
    bool donor = false;     // GT starting at position - 1
    bool acceptor = false;  // AG starting at position - 1
    bool start = false;     // ATG ending at position
    bool stop = false;      // TAA, TAG or TGA ending at position
};

inline SiteContext site_context(const std::vector<char> &seq, int position) {
    SiteContext site;
    if (position <= 0 || position >= static_cast<int>(seq.size())) return site;
    const char current = std::toupper(seq[position]);
    const char previous = std::toupper(seq[position - 1]);
    site.donor = previous == 'G' && current == 'T';
    site.acceptor = previous == 'A' && current == 'G';
    if (position >= 2) {
        const char first = std::toupper(seq[position - 2]);
        site.start = first == 'A' && previous == 'T' && current == 'G';
        site.stop = first == 'T' &&
            ((previous == 'A' && (current == 'A' || current == 'G')) ||
             (previous == 'G' && current == 'A'));
    }
    return site;
}

inline double destination_emission(int position, int to, const ModelInputs &inputs) {
    const int column = emission_column(to);
    double emission = inputs.emit[position][column];
    // Preserve the legacy TAG stop / AG acceptor overlap convention.
    if (emission <= -1e6 && position > 0 && position + 1 < static_cast<int>(inputs.seq.size()) &&
        std::toupper(inputs.seq[position - 1]) == 'A' &&
        std::toupper(inputs.seq[position]) == 'G') {
        emission = inputs.emit[position + 1][column];
    }
    return emission;
}

inline TransitionResult evaluate_transition(
    int position,
    int from,
    int to,
    const PathMetadata &previous,
    const ModelInputs &inputs,
    const SiteContext &site
) {
    TransitionResult result;
    if (position <= 0 || position >= static_cast<int>(inputs.seq.size()) ||
        from < 0 || from >= NUM_STATES || to < 0 || to >= NUM_STATES) {
        return result;
    }

    double log_t = inputs.trans[from][to];
    const bool stop = site.stop;

    if (is_exon(to) && is_exon(from) && from == to && stop &&
        codon_start_0based(position) % 3 == to - 1) {
        log_t = NEG_INF;
    }

    if (is_exon(from) && is_intron(to) && site.donor) {
        const int length = previous.exon_len - 2;
        log_t = NEG_INF;
        if ((to - 4) == (from - 1) && (length >= MIN_EXON || previous.exon_from == 0)) {
            log_t = inputs.site_score[splice_motif_start_0based(position)];
        }
    }

    if (is_intron(from) && is_exon(to) && site.acceptor) {
        const int length = previous.intron_len + 2;
        log_t = NEG_INF;
        if (length >= MIN_INTRON) {
            const int intron_frame = from - 4;
            const int exon_frame = to - 1;
            if ((exon_frame - intron_frame + 3) % 3 == length % 3) {
                log_t = inputs.site_score[splice_motif_start_0based(position)];
            }
        }
    }

    if (is_exon(from) && to == 0 && stop) {
        const int length = previous.exon_len - 2;
        const int frame = from - 1;
        if (codon_start_0based(position) % 3 == frame) {
            if (is_intron(previous.exon_from) ||
                (previous.exon_from == 0 && length > MIN_SINGLE)) {
                log_t = inputs.site_score[codon_start_0based(position)];
            } else {
                // A finite model penalty, not a forbidden-transition sentinel.
                log_t = -1e3;
            }
        }
    }

    if (is_exon(to) && from == 0 && position >= 2) {
        const int length = previous.inter_len + 2;
        if (site.start &&
            (length >= MIN_INTER || position < MIN_INTER) &&
            previous.predecessor_state == 0) {
            const int frame = to - 1;
            if (codon_start_0based(position) % 3 == frame) {
                log_t = position < 25 ? 1.0 : inputs.site_score[codon_start_0based(position)];
            }
        }
    }

    result.transition_score = log_t;
    result.emission_score = destination_emission(position, to, inputs);
    result.allowed = log_t > NEG_INF && result.emission_score > NEG_INF;

    PathMetadata next;
    next.predecessor_state = from;
    if (is_intron(to)) {
        next.intron_len = is_intron(from) ? previous.intron_len + 1 : 1;
    }
    if (is_exon(to)) {
        if (is_exon(from)) {
            next.exon_len = previous.exon_len + 1;
            next.exon_from = previous.exon_from;
        } else {
            next.exon_len = from == 0 ? 3 : 1;
            next.exon_from = from;
        }
    }
    if (to == 0) {
        next.inter_len = from == 0 ? previous.inter_len + 1 : 1;
    }
    result.next_metadata = next;
    return result;
}

inline TransitionResult evaluate_transition(
    int position,
    int from,
    int to,
    const PathMetadata &previous,
    const ModelInputs &inputs
) {
    return evaluate_transition(position, from, to, previous, inputs,
                               site_context(inputs.seq, position));
}

} // namespace uniann
