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
#include <limits>
#include <stdexcept>

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
static const int MIN_INTRON = 40;
static const int MIN_EXON   = 3;
static const int MIN_INTER  = 30;
static const int MIN_SINGLE = 100;

//------------------------------------------------------------
// Helpers
//------------------------------------------------------------

inline bool is_exon(int s) {
    return (s == 1 || s == 2 || s == 3);
}

inline bool is_intron(int s) {
    return (s == 4 || s == 5 || s == 6);
}

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
vector<array<float, NUM_STATES>> load_emissions(const string &file, int L) {
    vector<array<float, NUM_STATES>> emit(L);
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
        const char *cursor = line.c_str();
        int pos;
        if (!parse_int_field(cursor, pos) || pos < 0 || pos >= L)
            throw invalid_argument("Invalid emission position in " + file);
        for (int s = 0; s < NUM_STATES-2; s++) {
            if (!parse_float_field(cursor, emit[pos][s]) || !isfinite(emit[pos][s]))
                throw invalid_argument("Invalid emission value in " + file);
        }
        emit[pos][5]=emit[pos][4];
        emit[pos][6]=emit[pos][4];
    }
    return emit;
}

//------------------------------------------------------------
// Load sparse ATG/GT/AG scores
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
        const char *cursor = line.c_str();
        int pos;
        double score;
        if (!parse_int_field(cursor, pos) || pos < 0 || pos >= L ||
            !parse_double_field(cursor, score) || !isfinite(score))
            throw invalid_argument("Invalid site score in " + file);
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
void safety_check_gt_ag_atg(
    const vector<char> &seq,
    const vector<double> &gt_score,
    const vector<double> &ag_score,
    const vector<double> &atg_score,
    const vector<double> &stop_score
) {
    int L = seq.size();
    for (int pos = 0; pos < L - 1; pos++) {

        // Check GT
        if (gt_score[pos] > NEG_INF) {
            string dinuc;
            dinuc.push_back(toupper(seq[pos]));
            dinuc.push_back(toupper(seq[pos + 1]));
            if (dinuc != "GT") {
                cerr << "WARNING: GT score "<< gt_score[pos] << " at position "
                     << pos << " but sequence has " << dinuc << "\n";
            }
        }

        // Check AG
        if (ag_score[pos] > NEG_INF) {
            string dinuc;
            dinuc.push_back(toupper(seq[pos]));
            dinuc.push_back(toupper(seq[pos + 1]));
            if (dinuc != "AG") {
                cerr << "WARNING: AG score "<< ag_score[pos] << " at position "
                     << pos << " but sequence has " << dinuc << "\n";
            }
        }

        // Check ATG
        if (atg_score[pos] > NEG_INF) {
            string trinuc;
            trinuc.push_back(toupper(seq[pos]));
            trinuc.push_back(toupper(seq[pos + 1]));
            trinuc.push_back(toupper(seq[pos + 2]));
            if (trinuc != "ATG") {
                cerr << "WARNING: ATG score "<< atg_score[pos] << " at position "
                     << pos << " but sequence has " << trinuc << "\n";
            }
        }

        // Check STOP
        if (stop_score[pos] > NEG_INF) {
            string trinuc;
            trinuc.push_back(toupper(seq[pos]));
            trinuc.push_back(toupper(seq[pos + 1]));
            trinuc.push_back(toupper(seq[pos + 2]));
            if (trinuc != "TAG"  && trinuc != "TGA" && trinuc != "TAA") {
                cerr << "WARNING: STOP score "<< stop_score[pos] << " at position "
                     << pos << " but sequence has " << trinuc << "\n";
            }
        }
    }
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

DPMatrix init_dp(int L, const vector<array<float, NUM_STATES>> &emit)
{
    DPMatrix dp(L);

    // Initialization at position 0 — force start in N
    for (int s = 0; s < NUM_STATES; s++) {
        double start_prob = (s == 0 ? 0.0 : NEG_INF);
        double e = emit[0][s];

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
    const vector<array<float, NUM_STATES>> &emit,
    const vector<double> &gt_score,
    const vector<double> &ag_score,
    const vector<double> &atg_score,
    const vector<double> &stop_score,
    const vector<char> &seq,
    const vector<vector<double>> &trans
) {
    int L = seq.size();

    for (int i = 1; i < L; i++) {
        char b_prev = seq[i - 1];
        char b      = seq[i];
        bool is_stop = false;

        // on stop do not allow to continue in the same exon
        if ( i >= 2 ) {
          string codon;
          codon.push_back(toupper(seq[i - 2]));
          codon.push_back(toupper(seq[i - 1]));
          codon.push_back(toupper(seq[i]));

          if (codon == "TAA" || codon == "TAG" || codon == "TGA") 
            is_stop = true;
        }

        for (int to = 0; to < NUM_STATES; to++) {
            double emit_log = emit[i][to];

            // Fix for TAG stop where AG is acceptor
            if (emit_log <= -1e6 && toupper(b_prev) == 'A' && toupper(b) == 'G' && i + 1 < L) {
                emit_log = emit[i + 1][to];
            }

            double best = -1e18;
            int best_from = -1;

            for (int from = 0; from < NUM_STATES; from++) {

                double log_t = trans[from][to];

                //prohibit staying in the exon if a stop is found
                if (is_exon(to) && is_exon(from) && from==to && is_stop) {
                  if (((i - 2) % 3) == to - 1) {
                    log_t = NEG_INF;
                  }
                }

                //--------------------------------------------------------
                // Exon → Intron only at GT AND only if exon length ≥ MIN_EXON or it is the first exon
                //--------------------------------------------------------
                if (is_exon(from) && is_intron(to) &&
                    toupper(b_prev) == 'G' && toupper(b) == 'T')
                {
                    int len = dp[i - 1][from].exon_len - 2;
                    log_t = NEG_INF;

                    if ((to - 4) == (from - 1) && (len >= MIN_EXON || dp[i - 1][from].exon_from == 0)) {
                    // go into intron in the same frame
                    //cerr << "DEBUG at " << i
                    //     << " trying transition " << state_name[from]
                    //     << " " << state_name[to]
                    //     << " score " << dp[i - 1][from].dp
                    //     << " emit " << emit_log
                    //     << " length " << len << "\n";
                        log_t = gt_score[i - 1];
                    }

                    //cerr << "DEBUG GT probability " << log_t << "\n";
                }

                //--------------------------------------------------------
                // Intron → Exon only at AG AND only if intron length ≥ MIN_INTRON
                //--------------------------------------------------------
                if (is_intron(from) && is_exon(to) &&
                    toupper(b_prev) == 'A' && toupper(b) == 'G')
                {
                    int len = dp[i - 1][from].intron_len + 2;
                    log_t = NEG_INF;

                    if (len >= MIN_INTRON) {
                        int f = from - 4; // intron frame 0,1,2
                        int e = to - 1;   // exon frame 0,1,2
                        int mod = len % 3;
                        //cerr << "DEBUG at " << i
                        // << " trying transition " << state_name[from]
                        // << " " << state_name[to]
                        // << " score " << dp[i - 1][from].dp
                        // << " emit " << emit_log
                        // << " length " << len << "\n";

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

                    //cerr << "DEBUG AG probability " << log_t << "\n";
                }

                //--------------------------------------------------------
                // Exon → Noncoding after STOP codon (TAA, TAG, TGA)
                //--------------------------------------------------------
                if (is_exon(from) && to == 0 && is_stop) {
                    int len = dp[i - 1][from].exon_len - 2;
                    int frame = from - 1;
                    if (((i - 2) % 3) == frame) {
                        //cerr << "DEBUG at " << i
                        //     << " trying transition " << state_name[from]
                        //     << " " << state_name[to]
                        //     << " origin "<< dp[i-1][from].exon_from
                        //     << " score " << dp[i - 1][from].dp
                        //     << " emission " << emit_log << "\n";
                        if(is_intron(dp[i-1][from].exon_from) || (dp[i-1][from].exon_from == 0 && len > MIN_SINGLE)){ // if came from intron or came from noncoding and min length satisfied
                          log_t = stop_score[i - 2]; 
                        }else{
                          log_t = -1e3;
                        }
                    }
                    //cerr << "DEBUG STOP probability " << log_t << "\n";
                }

                //--------------------------------------------------------
                // Noncoding → Exon after START codon (ATG)
                //--------------------------------------------------------
                if (is_exon(to) && from == 0 && i >= 2) {
                    string codon;
                    codon.push_back(toupper(seq[i - 2]));
                    codon.push_back(toupper(seq[i - 1]));
                    codon.push_back(toupper(seq[i]));
                    int len = dp[i - 1][from].inter_len + 2;
                    if (codon == "ATG" && (len >= MIN_INTER || i < MIN_INTER) && dp[i-1][from].bt == 0) { //must come from non-coding

                        int frame = to - 1;
                        if (((i - 2) % 3) == frame) { //only in the right frame and if psauron score is positive in this frame
                        //cerr << "DEBUG at " << i
                        //     << " trying transition " << state_name[from]
                        //     << " " << state_name[to]
                        //     << " score " << dp[i - 1][from].dp
                        //     << " emission " << emit_log << "\n";

                            if (i < 25) {
                              log_t = 1.0;
                            } else {
                              log_t = atg_score[i - 2];
                              //log_t = 1.0;
                            }
                        }

                        //cerr << "DEBUG ATG probability " << log_t << "\n";
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
void write_gff_from_path(
    const vector<int> &path_states,
    const vector<double> &scores,
    const string &seqid,
    const string &f_fasta,
    double best_final
) { 
  // Look labels up per position instead of materializing one string per base.
  const auto labels = [&](size_t i) -> const string & { return state_name[path_states[i]]; };
  cout << "##gff-version 3\n";
  string current_state = labels(0);
  string previous_state = labels(0);
  double score_offset = scores[0];

  int start = 0;
  int end = 0;

  for (size_t i = 1; i < path_states.size(); i++) {
    if (labels(i) != current_state) {
      end = i - 1;
      if (current_state == "N")  { // transitioned to exon
        end = i - 2;
      }
      if (labels(i) == "N") { //transitioned to non-coding
        end = i + 1;
      }
      if(end > start) {
        if(is_label_exon(current_state) && is_label_intron(previous_state)) {
          write_gff_feature(seqid, current_state, start+2, end, f_fasta, (scores[i]-score_offset)/(end-start+1));
        }else if(is_label_intron(current_state) && is_label_exon(previous_state)) {
          write_gff_feature(seqid, current_state, start, end+2, f_fasta, (scores[i]-score_offset)/(end-start+1));
        }else{
          write_gff_feature(seqid, current_state, start, end, f_fasta, (scores[i]-score_offset)/(end-start+1)); 
        }
      }
      if (current_state == "N"){
        start = i - 1;
        score_offset = scores[i];
      } else if (labels(i) == "N" ){
        start = i + 2;
      } else {
        start = i;
      }
      previous_state = current_state;
      current_state = labels(i);
    }
  }

  // Final segment
  write_gff_feature(seqid, current_state, start, path_states.size() - 1, f_fasta, best_final);
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
    // Load emissions and splice scores
    //--------------------------------------------------------
    auto emit     =  load_emissions(f_emit, L);
    auto gt_score =  load_sparse_scores(f_gt, L);
    auto ag_score =  load_sparse_scores(f_ag, L);
    auto atg_score = load_sparse_scores(f_atg, L);
    auto stop_score = load_sparse_scores(f_stop, L);

    //--------------------------------------------------------
    // Safety check: GT/AG coordinates match sequence
    //--------------------------------------------------------
    safety_check_gt_ag_atg(seq, gt_score, ag_score, atg_score, stop_score);

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
    run_viterbi(dp, emit, gt_score, ag_score, atg_score, stop_score, seq, trans);

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
