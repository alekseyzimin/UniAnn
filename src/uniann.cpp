#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>

using namespace std;

//------------------------------------------------------------
// States:
// 0: N
// 1: E0  2: E1  3: E2
// 4: I0  5: I1  6: I2
//------------------------------------------------------------

static const int NUM_STATES = 7;

static const array<string, NUM_STATES> state_name = {
    "N", "E0", "E1", "E2", "I0", "I1", "I2"
};

static const double NEG_INF = -1e9;
static const int MIN_INTRON = 30;
static const int MIN_EXON   = 40;
static const int MIN_INTER  = 200;

//------------------------------------------------------------
// Helpers
//------------------------------------------------------------

inline bool is_exon(int s) {
    return (s == 1 || s == 2 || s == 3);
}

inline bool is_intron(int s) {
    return (s == 4 || s == 5 || s == 6);
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
            if (!isspace(c)) seq.push_back(toupper(c));
        }
    }
    return seq;
}

//------------------------------------------------------------
// Load emissions: pos \t 7 values
//------------------------------------------------------------
vector<array<double, NUM_STATES>> load_emissions(const string &file, int L) {
    vector<array<double, NUM_STATES>> emit(L);
    for (int i = 0; i < L; i++)
        for (int s = 0; s < NUM_STATES; s++)
            emit[i][s] = NEG_INF;

    ifstream in(file);
    if (!in) {
        cerr << "Cannot open emissions " << file << "\n";
        exit(1);
    }

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        int pos;
        ss >> pos;
        for (int s = 0; s < NUM_STATES; s++) {
            ss >> emit[pos][s];
        }
    }
    return emit;
}

//------------------------------------------------------------
// Load sparse GT/AG scores
//------------------------------------------------------------
vector<double> load_sparse_scores(const string &file, int L) {
    vector<double> scores(L, NEG_INF);

    ifstream in(file);
    if (!in) {
        cerr << "Cannot open scores " << file << "\n";
        exit(1);
    }

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        int pos;
        double score;
        ss >> pos >> score;
        scores[pos] = score;
    }
    return scores;
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
// Safety check: verify GT/AG coordinates match the sequence
//------------------------------------------------------------
void safety_check_gt_ag(
    const vector<char> &seq,
    const vector<double> &gt_score,
    const vector<double> &ag_score
) {
    int L = seq.size();
    for (int pos = 0; pos < L - 1; pos++) {

        // Check GT
        if (gt_score[pos] > -1e8) {
            string dinuc;
            dinuc.push_back(seq[pos]);
            dinuc.push_back(seq[pos + 1]);
            if (dinuc != "GT") {
                cerr << "WARNING: GT score at position "
                     << pos << " but sequence has " << dinuc << "\n";
            }
        }

        // Check AG
        if (ag_score[pos] > -1e8) {
            string dinuc;
            dinuc.push_back(seq[pos]);
            dinuc.push_back(seq[pos + 1]);
            if (dinuc != "AG") {
                cerr << "WARNING: AG score at position "
                     << pos << " but sequence has " << dinuc << "\n";
            }
        }
    }
}

//------------------------------------------------------------
// Baseline transition matrix (log-space), excluding GT/AG overrides
//------------------------------------------------------------
vector<vector<double>> init_transitions() {
    vector<vector<double>> trans(NUM_STATES, vector<double>(NUM_STATES, -1e9));

    double max_log_prob = 1.0;

    // Noncoding self
    trans[0][0] = log(1.0);

    // Exon frame cycling: none (self only)
    trans[1][1] = max_log_prob;
    trans[2][2] = max_log_prob;
    trans[3][3] = max_log_prob;

    // Intron states: self only
    trans[4][4] = max_log_prob;
    trans[5][5] = max_log_prob;
    trans[6][6] = max_log_prob;

    return trans;
}

//------------------------------------------------------------
// DP structures
//------------------------------------------------------------
struct DPCell {
    double dp;
    int bt;
    int intron_len;
    int exon_len;
    int inter_len;
};

vector<vector<DPCell>> init_dp(int L,
                               const vector<array<double, NUM_STATES>> &emit)
{
    vector<vector<DPCell>> dp(L, vector<DPCell>(NUM_STATES));

    // Initialization at position 0 — force start in N
    for (int s = 0; s < NUM_STATES; s++) {
        double start_prob = (s == 0 ? 0.0 : -1e9);
        double e = emit[0][s];

        dp[0][s].dp = start_prob + e;
        dp[0][s].bt = -1;
        dp[0][s].intron_len = is_intron(s) ? 1 : 0;
        dp[0][s].exon_len   = is_exon(s)   ? 1 : 0;
        dp[0][s].inter_len  = (s == 0)   ? 1 : 0;
    }

    return dp;
}

//------------------------------------------------------------
// Full Viterbi DP
//------------------------------------------------------------
void run_viterbi(
    vector<vector<DPCell>> &dp,
    const vector<array<double, NUM_STATES>> &emit,
    const vector<double> &gt_score,
    const vector<double> &ag_score,
    const vector<char> &seq,
    const vector<vector<double>> &trans
) {
    int L = seq.size();

    for (int i = 1; i < L; i++) {
        char b_prev = seq[i - 1];
        char b      = seq[i];

        // on stop do not allow to continue in the same exon
        

        for (int to = 0; to < NUM_STATES; to++) {

            double emit_log = emit[i][to];

            // Fix for TAG stop where AG is acceptor
            if (emit_log <= -1e6 && b_prev == 'A' && b == 'G' && i + 1 < L) {
                emit_log = emit[i + 1][to];
            }

            double best = -1e18;
            int best_from = -1;

            for (int from = 0; from < NUM_STATES; from++) {

                double log_t = trans[from][to];

                //--------------------------------------------------------
                // Exon → Intron only at GT AND only if exon length ≥ MIN_EXON
                //--------------------------------------------------------
                if (is_exon(from) && is_intron(to) &&
                    b_prev == 'G' && b == 'T')
                {
                    int len = dp[i - 1][from].exon_len - 2;

                    cerr << "DEBUG at " << i
                         << " trying transition " << state_name[from]
                         << " " << state_name[to]
                         << " score " << dp[i - 1][from].dp
                         << " emit " << emit_log
                         << " length " << len << "\n";

                    log_t = -1e9;

                    if ((to - 4) == (from - 1) && len >= MIN_EXON) {
                        log_t = gt_score[i - 1];
                    }

                    cerr << "DEBUG probability " << log_t << "\n";
                }

                //--------------------------------------------------------
                // Intron → Exon only at AG AND only if intron length ≥ MIN_INTRON
                //--------------------------------------------------------
                if (is_intron(from) && is_exon(to) &&
                    b_prev == 'A' && b == 'G')
                {
                    int len = dp[i - 1][from].intron_len + 2;

                    cerr << "DEBUG at " << i
                         << " trying transition " << state_name[from]
                         << " " << state_name[to]
                         << " score " << dp[i - 1][from].dp
                         << " emit " << emit_log
                         << " length " << len << "\n";

                    log_t = -1e9;

                    if (len >= MIN_INTRON) {
                        int f = from - 4; // intron frame 0,1,2
                        int e = to - 1;   // exon frame 0,1,2
                        int mod = len % 3;

                        // Frame‑compatible transitions
                        if ((f == 0 && e == 0) || (f == 1 && e == 1) || (f == 2 && e == 2)) {
                            if (mod == 0) log_t = ag_score[i - 1];
                        }
                        else if ((f == 0 && e == 1) || (f == 1 && e == 2) || (f == 2 && e == 0)) {
                            if (mod == 1) log_t = ag_score[i - 1];
                        }
                        else if ((f == 0 && e == 2) || (f == 1 && e == 0) || (f == 2 && e == 1)) {
                            if (mod == 2) log_t = ag_score[i - 1];
                        }
                    }

                    cerr << "DEBUG probability " << log_t << "\n";
                }

                //--------------------------------------------------------
                // Exon → Noncoding after STOP codon (TAA, TAG, TGA)
                //--------------------------------------------------------
                if (is_exon(from) && to == 0 && i >= 2) {
                    string codon;
                    codon.push_back(seq[i - 2]);
                    codon.push_back(seq[i - 1]);
                    codon.push_back(seq[i]);

                    if (codon == "TAA" || codon == "TAG" || codon == "TGA") {

                        cerr << "DEBUG at " << i
                             << " trying transition " << state_name[from]
                             << " " << state_name[to]
                             << " score " << dp[i - 1][from].dp
                             << " emission " << emit_log << "\n";

                        int frame = from - 1;
                        if (((i - 2) % 3) == frame) {
                            log_t = 1.0; // max_log_prob
                        }

                        cerr << "DEBUG probability " << log_t << "\n";
                    }
                }

                //--------------------------------------------------------
                // Noncoding → Exon after START codon (ATG)
                //--------------------------------------------------------
                if (is_exon(to) && from == 0 && i >= 2) {
                    string codon;
                    codon.push_back(seq[i - 2]);
                    codon.push_back(seq[i - 1]);
                    codon.push_back(seq[i]);
                    int len = dp[i - 1][from].inter_len + 2;
                    if (codon == "ATG" && (len >= MIN_INTER || i < MIN_INTER)) {

                        cerr << "DEBUG at " << i
                             << " trying transition " << state_name[from]
                             << " " << state_name[to]
                             << " score " << dp[i - 1][from].dp
                             << " emission " << emit_log << "\n";

                        int frame = to - 1;
                        if (((i - 2) % 3) == frame) {
                            log_t = 1.0;
                        }

                        cerr << "DEBUG probability " << log_t << "\n";
                    }
                }

                //--------------------------------------------------------
                // Candidate score
                //--------------------------------------------------------
                double cand = dp[i - 1][from].dp + log_t + emit_log;

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
                else
                    dp[i][to].exon_len = 1;
            } else {
                dp[i][to].exon_len = 0;
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
    const vector<vector<DPCell>> &dp,
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
    const vector<vector<DPCell>> &dp,
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
// Convert state numbers to labels
//------------------------------------------------------------
vector<string> states_to_labels(const vector<int> &path_states) {
    vector<string> labels(path_states.size());
    for (size_t i = 0; i < path_states.size(); i++) {
        labels[i] = state_name[path_states[i]];
    }
    return labels;
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

    // Perl logic: if start0 == 0, decrement start0
    if (start0 == 0) start0--;

    string type;

    if (state[0] == 'E') {
        type = "CDS";
        start = start0 + 2;
    }
    else if (state[0] == 'I') {
        type = "intron";
        end = end0 + 2;
    }
    else {
        type = "region";   // noncoding
        end = end0 + 2;
        gff_index++;
    }

    // Perl rounds best_final to 2 decimals
    double bf = floor(best_final * 100.0) / 100.0;

    cout << seqid << "\t"
         << "HMM" << "\t"
         << type << "\t"
         << start << "\t"
         << end << "\t"
         << bf << "\t"
         << "." << "\t"
         << "." << "\t"
         << "Parent=" << tname << "." << gff_index
         << ";state=" << state
         << "\n";
}

//------------------------------------------------------------
// Emit all GFF3 features from the state path
//------------------------------------------------------------
void write_gff_from_path(
    const vector<string> &labels,
    const string &seqid,
    const string &f_fasta,
    double best_final
) {
    cout << "##gff-version 3\n";

    string current_state = labels[0];
    int start = 0;
    int end = 0;

    for (size_t i = 1; i < labels.size(); i++) {
        if (labels[i] != current_state) {
          end = i - 1;
          if (current_state == "N")
            end = i - 2;

          if (end > start) 
              write_gff_feature(seqid, current_state, start, end, f_fasta, best_final);
          if (current_state == "N"){
            start = i - 2;
          } else {
            start = i;
          }
          current_state = labels[i];
        }
    }

    // Final segment
    write_gff_feature(seqid, current_state, start, labels.size() - 1, f_fasta, best_final);
}

//------------------------------------------------------------
// MAIN
//------------------------------------------------------------
int main(int argc, char** argv) {

    if (argc < 5) {
        cerr << "Usage: " << argv[0]
             << " seq.fasta emissions.txt gt.txt ag.txt\n";
        return 1;
    }

    string f_fasta = argv[1];
    string f_emit  = argv[2];
    string f_gt    = argv[3];
    string f_ag    = argv[4];

    //--------------------------------------------------------
    // Load FASTA
    //--------------------------------------------------------
    string seq_str = read_fasta(f_fasta);
    int L = seq_str.size();

    vector<char> seq(L);
    for (int i = 0; i < L; i++)
        seq[i] = seq_str[i];

    //--------------------------------------------------------
    // Load emissions and splice scores
    //--------------------------------------------------------
    auto emit     = load_emissions(f_emit, L);
    auto gt_score = load_sparse_scores(f_gt, L);
    auto ag_score = load_sparse_scores(f_ag, L);

    //--------------------------------------------------------
    // Safety check: GT/AG coordinates match sequence
    //--------------------------------------------------------
    safety_check_gt_ag(seq, gt_score, ag_score);

    //--------------------------------------------------------
    // Initialize transitions
    //--------------------------------------------------------
    auto trans = init_transitions();

    //--------------------------------------------------------
    // Initialize DP
    //--------------------------------------------------------
    auto dp = init_dp(L, emit);

    //--------------------------------------------------------
    // Run full Viterbi
    //--------------------------------------------------------
    run_viterbi(dp, emit, gt_score, ag_score, seq, trans);

    //--------------------------------------------------------
    // Termination
    //--------------------------------------------------------
    auto [best_final, best_state] = viterbi_termination(dp, L);

    //--------------------------------------------------------
    // Backtrace
    //--------------------------------------------------------
    auto path_states = viterbi_backtrace(dp, L, best_state);
    auto path_labels = states_to_labels(path_states);

    //--------------------------------------------------------
    // Write GFF3 output
    //--------------------------------------------------------
    string seqid = get_fasta_header(f_fasta);
    write_gff_from_path(path_labels, seqid, f_fasta, best_final);

    return 0;
}

