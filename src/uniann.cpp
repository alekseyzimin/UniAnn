#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <cerrno>
#include <cstdlib>

#include "viterbi_model.h"
#include <limits>

using namespace std;
using namespace uniann;

//------------------------------------------------------------
// Helpers
//------------------------------------------------------------

inline bool is_label_exon(string s) {
    return (s == "E0" || s == "E1" || s == "E2");
}

inline bool is_label_intron(string s) {
    return (s == "I0" || s == "I1" || s == "I2");
}


//------------------------------------------------------------
// Read FASTA (single sequence)
//------------------------------------------------------------
string read_fasta(const string &file) {
    ifstream in(file);
    if (!in) {
        cerr << "Cannot open FASTA " << file << "\n";
        exit(1);
    }

    string seq, line;
    while (getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line) {
            if (!isspace(c)) seq.push_back(c);
        }
    }
    return seq;
}

//------------------------------------------------------------
// Whitespace-separated numeric fields parsed in place. A stringstream per
// line costs more than the decoder itself on a chromosome-sized emission
// file; strtol/strtof/strtod accept the same decimal text.
//------------------------------------------------------------
static bool parse_int_field(const char *&cursor, int &value) {
    char *end = nullptr;
    errno = 0;
    const long parsed = strtol(cursor, &end, 10);
    if (end == cursor || errno == ERANGE ||
        parsed < numeric_limits<int>::min() || parsed > numeric_limits<int>::max())
        return false;
    cursor = end;
    value = static_cast<int>(parsed);
    return true;
}

static bool parse_float_field(const char *&cursor, float &value) {
    char *end = nullptr;
    value = strtof(cursor, &end);
    if (end == cursor) return false;
    cursor = end;
    return true;
}

static bool parse_double_field(const char *&cursor, double &value) {
    char *end = nullptr;
    value = strtod(cursor, &end);
    if (end == cursor) return false;
    cursor = end;
    return true;
}

//------------------------------------------------------------
// Load emissions: pos \t 5 values (states 5 and 6 share the intron column)
//------------------------------------------------------------
vector<EmissionRow> load_emissions(const string &file, int L) {
    vector<EmissionRow> emit(L);
    for (int i = 0; i < L; i++)
        emit[i].fill(NEG_INF);

    ifstream in(file);
    if (!in) {
        cerr << "Cannot open emissions " << file << "\n";
        exit(1);
    }

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const char *cursor = line.c_str();
        int pos;
        if (!parse_int_field(cursor, pos) || pos < 0 || pos >= L)
            throw invalid_argument("Invalid emission position in " + file);
        for (int s = 0; s < EMISSION_COLUMNS; s++) {
            if (!parse_float_field(cursor, emit[pos][s]) || !isfinite(emit[pos][s]))
                throw invalid_argument("Invalid emission value in " + file);
        }
    }
    return emit;
}

//------------------------------------------------------------
// Load sparse ATG/GT/AG/stop scores into one per-position array
//------------------------------------------------------------
// A score is kept only where the sequence has the motif the file describes;
// the transition rules never read a score anywhere else. Positions whose
// sequence differs are reported, as the former separate safety check did.
static const char *motif_label(SiteMotif motif) {
    switch (motif) {
        case SiteMotif::Donor: return "GT";
        case SiteMotif::Acceptor: return "AG";
        case SiteMotif::Start: return "ATG";
        case SiteMotif::Stop: return "STOP";
        default: return "";
    }
}

void load_site_scores(const string &file, SiteMotif motif, const vector<char> &seq,
                      vector<double> &site_score) {
    const int L = seq.size();
    ifstream in(file);
    if (!in) {
        cerr << "Cannot open scores " << file << "\n";
        exit(1);
    }

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const char *cursor = line.c_str();
        int pos;
        double score;
        if (!parse_int_field(cursor, pos) || pos < 0 || pos >= L ||
            !parse_double_field(cursor, score) || !isfinite(score))
            throw invalid_argument("Invalid site score in " + file);
        if (motif_starting_at(seq, pos) == motif) {
            site_score[pos] = score;
        } else if (score > NEG_INF) {
            const int width = motif == SiteMotif::Donor || motif == SiteMotif::Acceptor ? 2 : 3;
            string found;
            for (int i = pos; i < pos + width && i < L; i++) found.push_back(toupper(seq[i]));
            cerr << "WARNING: " << motif_label(motif) << " score " << score << " at position "
                 << pos << " but sequence has " << found << "\n";
        }
    }
}

//------------------------------------------------------------
// Extract FASTA header
//------------------------------------------------------------
string get_fasta_header(const string &file) {
    ifstream in(file);
    if (!in) return "sequence";

    string line;
    while (getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            stringstream ss(line.substr(1));
            string id;
            ss >> id;
            return id;
        }
    }
    return "sequence";
}

//------------------------------------------------------------
// Baseline transition matrix (log-space), excluding GT/AG overrides
//------------------------------------------------------------
vector<vector<double>> init_transitions() {
    vector<vector<double>> trans(NUM_STATES, vector<double>(NUM_STATES, NEG_INF));

    double self_prob = 0.0;

    // Noncoding self
    trans[0][0] = self_prob;

    // Exon frame cycling: none (self only)
    trans[1][1] = self_prob;
    trans[2][2] = self_prob;
    trans[3][3] = self_prob;

    // Intron states: self only
    trans[4][4] = self_prob;
    trans[5][5] = self_prob;
    trans[6][6] = self_prob;

    return trans;
}

//------------------------------------------------------------
// DP structures
//------------------------------------------------------------
// Field order keeps the cell at 24 bytes; the previous order padded it to 32.
struct DPCell {
    double dp;
    int intron_len;
    int exon_len;
    int inter_len;
    short int bt;
    short int exon_from;
};

// One contiguous block instead of one heap allocation per position: the
// per-row vector header and allocator overhead cost about 40 bytes per base.
class DPMatrix {
    vector<DPCell> cells;
public:
    explicit DPMatrix(int L) : cells(static_cast<size_t>(L) * NUM_STATES) {}
    DPCell *operator[](int position) { return cells.data() + static_cast<size_t>(position) * NUM_STATES; }
    const DPCell *operator[](int position) const { return cells.data() + static_cast<size_t>(position) * NUM_STATES; }
};

DPMatrix init_dp(int L, const vector<EmissionRow> &emit)
{
    DPMatrix dp(L);

    // Initialization at position 0 — force start in N
    for (int s = 0; s < NUM_STATES; s++) {
        double start_prob = (s == 0 ? 0.0 : NEG_INF);
        double e = emit[0][emission_column(s)];

        dp[0][s].dp = start_prob + e;
        dp[0][s].bt = -1;
        dp[0][s].intron_len = is_intron(s) ? 1 : 0;
        dp[0][s].exon_len   = is_exon(s)   ? 1 : 0;
        dp[0][s].inter_len  = (s == 0)   ? 1 : 0;
        dp[0][s].exon_from   = is_exon(s)   ? 0 : 0;
    }

    return dp;
}

//------------------------------------------------------------
// Full Viterbi DP
//------------------------------------------------------------
void run_viterbi(
    DPMatrix &dp,
    const ModelInputs &inputs
) {
    const vector<char> &seq = inputs.seq;
    const int L = seq.size();

    for (int i = 1; i < L; i++) {
        const SiteContext site = site_context(seq, i);

        for (int to = 0; to < NUM_STATES; to++) {
            double best = -1e18;
            int best_from = -1;

            for (int from = 0; from < NUM_STATES; from++) {
                PathMetadata previous_metadata;
                previous_metadata.intron_len = dp[i - 1][from].intron_len;
                previous_metadata.exon_len = dp[i - 1][from].exon_len;
                previous_metadata.inter_len = dp[i - 1][from].inter_len;
                previous_metadata.exon_from = dp[i - 1][from].exon_from;
                previous_metadata.predecessor_state = dp[i - 1][from].bt;

                const auto edge = evaluate_transition(i, from, to, previous_metadata, inputs, site);
                // Preserve upstream propagation of the finite NEG_INF sentinel:
                // a forbidden edge still competes with a very low score.
                double cand = dp[i - 1][from].dp + edge.transition_score + edge.emission_score;

                if (cand > best) {
                    best = cand;
                    best_from = from;
                }
            }

            //--------------------------------------------------------
            // Store best transition
            //--------------------------------------------------------
            dp[i][to].dp = best;
            dp[i][to].bt = best_from;

            // Track intron length
            if (is_intron(to)) {
                if (best_from >= 0 && is_intron(best_from))
                    dp[i][to].intron_len = dp[i - 1][best_from].intron_len + 1;
                else
                    dp[i][to].intron_len = 1;
            } else {
                dp[i][to].intron_len = 0;
            }

            // Track exon length
            if (is_exon(to)) {
                if (best_from >= 0 && is_exon(best_from))
                    dp[i][to].exon_len = dp[i - 1][best_from].exon_len + 1;
                else {
                  if (best_from == 0)
                    dp[i][to].exon_len = 3; //count the ATG
                  else 
                    dp[i][to].exon_len = 1;   
                }
            } else {
                dp[i][to].exon_len = 0;
            }

            // Track if the first exon
            if (is_exon(to)) {
              if (best_from >= 0) {
                if(is_exon(best_from)){
                  dp[i][to].exon_from=dp[i-1][to].exon_from;
                } else {
                  dp[i][to].exon_from = best_from;
                }
              }else{
                dp[i][to].exon_from = -1;
              }
            }
  
            // Track inter length
            if (to == 0) {
                if (best_from == 0)
                    dp[i][to].inter_len = dp[i - 1][best_from].inter_len + 1;
                else
                    dp[i][to].inter_len = 1;
            } else {
                dp[i][to].inter_len = 0;
            }
        }
    }
}

//------------------------------------------------------------
// Termination: find best final state
//------------------------------------------------------------
pair<double, int> viterbi_termination(
    const DPMatrix &dp,
    int L
) {
    double best_final = -1e18;
    int best_state = -1;

    for (int s = 0; s < NUM_STATES; s++) {
        if (dp[L - 1][s].dp > best_final) {
            best_final = dp[L - 1][s].dp;
            best_state = s;
        }
    }
    return {best_final, best_state};
}

//------------------------------------------------------------
// Backtrace: reconstruct optimal path
//------------------------------------------------------------
vector<int> viterbi_backtrace(
    const DPMatrix &dp,
    int L,
    int best_state
) {
    vector<int> path_states(L);
    int cur = best_state;

    for (int i = L - 1; i >= 0; i--) {
        path_states[i] = cur;
        cur = dp[i][cur].bt;
        if (cur < 0) break;
    }
    return path_states;
}

//------------------------------------------------------------
// Backtrace: get scores for optimal path
//------------------------------------------------------------
vector<double> viterbi_backtrace_scores(
    const DPMatrix &dp,
    int L,
    int best_state
) {
    vector<double> path_scores(L);
    int cur = best_state;

    for (int i = L - 1; i >= 0; i--) {
        path_scores[i] = dp[i][cur].dp;
        cur = dp[i][cur].bt;
        if (cur < 0) break;
    }
    return path_scores;
}

//------------------------------------------------------------
// Global index counter (mirrors Perl's $index)
//------------------------------------------------------------
int gff_index = 0;

//------------------------------------------------------------
// Write a single GFF3 feature
//------------------------------------------------------------
void write_gff_feature(
    const string &seqid,
    const string &state,
    int start0,
    int end0,
    const string &f_fasta,
    double best_final
) {
    // Extract filename from path
    size_t slash = f_fasta.find_last_of("/\\");
    string tname = (slash == string::npos)
        ? f_fasta
        : f_fasta.substr(slash + 1);

    int start = start0;   // 0-based
    int end   = end0;

    string type;

    if (state[0] == 'E') {
        type = "CDS";
        start = start0;
    }
    else if (state[0] == 'I') {
        type = "intron";
        end = end0;
    }
    else {
        type = "region";   // noncoding
        end = end0;
        gff_index++;
    }

    // Perl rounds best_final to 2 decimals
    double bf = floor(best_final * 100.0) / 100.0;
    
    if(type != "region")
      cout << seqid << "\t"
         << "UniAnn" << "\t"
         << type << "\t"
         << start << "\t"
         << end << "\t"
         << bf << "\t"
         << "+" << "\t"
         << "." << "\t"
         << "Parent=" << tname << "." << gff_index
         << ";state=" << state
         << "\n";
}

//------------------------------------------------------------
// Emit all GFF3 features from the state path
//------------------------------------------------------------
// Consume state changes so both a per-base path and compressed runs use the
// same upstream GFF boundaries and interval-score convention.
class GffPathWriter {
    const string &seqid;
    const string &fasta;
    string current_state;
    string previous_state;
    int start = 0;
    double score_offset;

public:
    GffPathWriter(const string &sequence_id, const string &fasta_filename,
                  const string &initial_state, double initial_score)
        : seqid(sequence_id), fasta(fasta_filename), current_state(initial_state),
          previous_state(initial_state), score_offset(initial_score) {
        cout << "##gff-version 3\n";
    }

    void transition(const string &next_state, int position, double score) {
        int end = position - 1;
        if (current_state == "N") end = position - 2;
        if (next_state == "N") end = position + 1;
        if (end > start) {
            const double interval_score = (score - score_offset) / (end - start + 1);
            if (is_label_exon(current_state) && is_label_intron(previous_state)) {
                write_gff_feature(seqid, current_state, start + 2, end, fasta, interval_score);
            } else if (is_label_intron(current_state) && is_label_exon(previous_state)) {
                write_gff_feature(seqid, current_state, start, end + 2, fasta, interval_score);
            } else {
                write_gff_feature(seqid, current_state, start, end, fasta, interval_score);
            }
        }
        if (current_state == "N") {
            start = position - 1;
            score_offset = score;
        } else if (next_state == "N") {
            start = position + 2;
        } else {
            start = position;
        }
        previous_state = current_state;
        current_state = next_state;
    }

    void finish(int length, double best_final) {
        // Preserve upstream's special final-segment score and boundaries.
        write_gff_feature(seqid, current_state, start, length - 1, fasta, best_final);
    }
};

void write_gff_from_path(
    const vector<int> &path_states, const vector<double> &scores,
    const string &seqid, const string &f_fasta, double best_final
) {
    GffPathWriter writer(seqid, f_fasta, state_name[path_states[0]], scores[0]);
    for (size_t i = 1; i < path_states.size(); ++i) {
        if (path_states[i] != path_states[i - 1])
            writer.transition(state_name[path_states[i]], i, scores[i]);
    }
    writer.finish(path_states.size(), best_final);
}

//------------------------------------------------------------
// MAIN
//------------------------------------------------------------
int run_uniann(int argc, char** argv) {

    if (argc < 7) {
        cerr << "Usage: " << argv[0]
             << " seq.fasta emissions.txt gt.txt ag.txt atg.txt stop.txt [--no-dp-dump]\n";
        return 1;
    }
    bool no_dp_dump = false;
    for (int argument = 7; argument < argc; ++argument) {
        const string option = argv[argument];
        if (option == "--no-dp-dump") no_dp_dump = true;
        else throw invalid_argument("Unknown option: " + option);
    }

    string f_fasta = argv[1];
    string f_emit  = argv[2];
    string f_gt    = argv[3];
    string f_ag    = argv[4];
    string f_atg   = argv[5];
    string f_stop   = argv[6];

    //--------------------------------------------------------
    // Load FASTA
    //--------------------------------------------------------
    string seq_str = read_fasta(f_fasta);
    if (seq_str.empty() || seq_str.size() > static_cast<size_t>(numeric_limits<int>::max()))
        throw invalid_argument("FASTA must contain 1..INT_MAX sequence bases");
    const int L = seq_str.size();

    vector<char> seq(seq_str.begin(), seq_str.end());
    string().swap(seq_str);

    //--------------------------------------------------------
    // Load emissions and site scores; the score loader also reports
    // coordinates whose sequence lacks the expected motif.
    //--------------------------------------------------------
    auto emit = load_emissions(f_emit, L);
    vector<double> site_score(L, NEG_INF);
    load_site_scores(f_gt, SiteMotif::Donor, seq, site_score);
    load_site_scores(f_ag, SiteMotif::Acceptor, seq, site_score);
    load_site_scores(f_atg, SiteMotif::Start, seq, site_score);
    load_site_scores(f_stop, SiteMotif::Stop, seq, site_score);

    //--------------------------------------------------------
    // Initialize transitions
    //--------------------------------------------------------
    auto trans = init_transitions();
    const ModelInputs inputs{emit, site_score, seq, trans};

    //--------------------------------------------------------
    // Initialize DP
    //--------------------------------------------------------
    auto dp = init_dp(L, emit);

    //--------------------------------------------------------
    // Run full Viterbi
    //--------------------------------------------------------
    run_viterbi(dp, inputs);

    //print DP and BT matrices unless the caller asked to skip the dump
    if (!no_dp_dump) {
      for (int i = 0; i < L; i++) {
      fprintf(stderr,"%d\tdp\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",i,int(dp[i][0].dp),int(dp[i][1].dp),int(dp[i][2].dp),int(dp[i][3].dp),int(dp[i][4].dp),int(dp[i][5].dp),int(dp[i][6].dp));
      fprintf(stderr,"%d\tbt\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",i,int(dp[i][0].bt),int(dp[i][1].bt),int(dp[i][2].bt),int(dp[i][3].bt),int(dp[i][4].bt),int(dp[i][5].bt),int(dp[i][6].bt));
      }
    }


    //--------------------------------------------------------
    // Termination
    //--------------------------------------------------------
    auto [best_final, best_state] = viterbi_termination(dp, L);

    //--------------------------------------------------------
    // Backtrace
    //--------------------------------------------------------
    auto path_states = viterbi_backtrace(dp, L, best_state);
    auto path_scores = viterbi_backtrace_scores(dp, L, best_state);

    //--------------------------------------------------------
    // Write GFF3 output
    //--------------------------------------------------------
    string seqid = get_fasta_header(f_fasta);
    write_gff_from_path(path_states, path_scores, seqid, f_fasta, best_final);

    return 0;
}


int main(int argc, char **argv) {
    try {
        return run_uniann(argc, argv);
    } catch (const exception &error) {
        cerr << error.what() << '\n';
        return 1;
    }
}
