// compute_markov_scores.cpp
#include <bits/stdc++.h>
using namespace std;

// Trim helpers
static inline std::string ltrim(const std::string &s) {
    size_t i = 0;
    while (i < s.size() && isspace((unsigned char)s[i])) ++i;
    return s.substr(i);
}
static inline std::string rtrim(const std::string &s) {
    if (s.empty()) return s;
    size_t i = s.size();
    while (i > 0 && isspace((unsigned char)s[i-1])) --i;
    return s.substr(0, i);
}
static inline std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}
static inline std::string toupper_str(const std::string &s) {
    std::string r = s;
    for (char &c : r) c = toupper((unsigned char)c);
    return r;
}

// split by whitespace
static inline vector<string> split_ws(const string &s) {
    vector<string> out;
    istringstream iss(s);
    string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage:\ncompute_markov_scores genome.fa positive_model.pwm negative_model.pwm\n";
        return 1;
    }

    string fasta_path = argv[1];
    string pos_model_path = argv[2];
    string neg_model_path = argv[3];

    // --- nucleotide code maps ---
    vector<char> narray = {'A','C','G','T'};
    unordered_map<string,int> code2;
    unordered_map<string,int> code3;
    unordered_map<char,int> code1;
    int n = 0, n2 = 0, n3 = 0;
    for (int i = 0; i < 4; ++i) {
        code1[narray[i]] = n++;
        for (int j = 0; j < 4; ++j) {
            string s2; s2.push_back(narray[i]); s2.push_back(narray[j]);
            code2[s2] = n2++;
            for (int k = 0; k < 4; ++k) {
                string s3; s3.push_back(narray[i]); s3.push_back(narray[j]); s3.push_back(narray[k]);
                code3[s3] = n3++;
            }
        }
    }

    // --- variables to hold HMM tables ---
    vector<vector<double>> donor_freq;        // donor 0HMM (4 columns)
    vector<vector<double>> donor_hmm_freq;    // donor 1HMM (16 columns)
    vector<vector<double>> donor_hmm2_freq;   // donor 2HMM (64 columns)

    vector<vector<double>> acceptor_freq;     // acceptor 0HMM
    vector<vector<double>> acceptor_hmm_freq; // acceptor 1HMM
    vector<vector<double>> acceptor_hmm2_freq;// acceptor 2HMM

    vector<vector<double>> donor_nfreq;
    vector<vector<double>> donor_hmm_nfreq;
    vector<vector<double>> donor_hmm2_nfreq;

    vector<vector<double>> acceptor_nfreq;
    vector<vector<double>> acceptor_hmm_nfreq;
    vector<vector<double>> acceptor_hmm2_nfreq;

    int donor_length = 0;
    int acceptor_length = 0;
    bool use_negatives = false;

    // --- load genome FASTA ---
    unordered_map<string,string> genome_seqs;
    {
        ifstream fin(fasta_path);
        if (!fin) {
            cerr << "Cannot open FASTA: " << fasta_path << "\n";
            return 1;
        }
        string line, scf, seq;
        while (getline(fin, line)) {
            if (!line.empty() && line[0] == '>') {
                if (!scf.empty()) {
                    genome_seqs[scf] = seq;
                    seq.clear();
                }
                // header
                auto parts = split_ws(line);
                scf = parts[0].substr(1);
            } else {
                seq += trim(line);
            }
        }
        if (!scf.empty()) genome_seqs[scf] = seq;
    }

    // --- helper to parse zoeHMM file into the provided tables ---
    auto parse_zoeHMM = [&](const string &path,
                           vector<vector<double>> &don0,
                           vector<vector<double>> &don1,
                           vector<vector<double>> &don2,
                           vector<vector<double>> &acc0,
                           vector<vector<double>> &acc1,
                           vector<vector<double>> &acc2,
                           int &don_len,
                           int &acc_len) -> bool
    {
        ifstream fin(path);
        if (!fin) return false;
        string line;
        if (!getline(fin, line)) return false;
        if (line.rfind("zoeHMM", 0) != 0) return false; // must start with zoeHMM

        while (getline(fin, line)) {
            string t = trim(line);
            if (t.rfind("Donor 0HMM", 0) == 0) {
                vector<vector<double>> table;
                while (getline(fin, line)) {
                    if (line.find("NN TRM") != string::npos) break;
                    string s = ltrim(line);
                    auto f = split_ws(s);
                    if (f.empty()) continue;
                    vector<double> row(4, 0.0);
                    for (int j = 0; j < 4 && j < (int)f.size(); ++j) row[j] = stod(f[j]);
                    table.push_back(row);
                }
                don0 = move(table);
                don_len = (int)don0.size();
            } else if (t.rfind("Donor 1HMM", 0) == 0) {
                vector<vector<double>> table;
                while (getline(fin, line)) {
                    if (line.find("NN TRM") != string::npos) break;
                    string s = ltrim(line);
                    auto f = split_ws(s);
                    if (f.empty()) continue;
                    vector<double> row(16, 0.0);
                    for (int j = 0; j < 16 && j < (int)f.size(); ++j) row[j] = stod(f[j]);
                    table.push_back(row);
                }
                don1 = move(table);
            } else if (t.rfind("Donor 2HMM", 0) == 0) {
                vector<vector<double>> table;
                while (getline(fin, line)) {
                    if (line.find("NN TRM") != string::npos) break;
                    string s = ltrim(line);
                    auto f = split_ws(s);
                    if (f.empty()) continue;
                    vector<double> row(64, 0.0);
                    for (int j = 0; j < 64 && j < (int)f.size(); ++j) row[j] = stod(f[j]);
                    table.push_back(row);
                }
                don2 = move(table);
            } else if (t.rfind("Acceptor 0HMM", 0) == 0) {
                vector<vector<double>> table;
                while (getline(fin, line)) {
                    if (line.find("NN TRM") != string::npos) break;
                    string s = ltrim(line);
                    auto f = split_ws(s);
                    if (f.empty()) continue;
                    vector<double> row(4, 0.0);
                    for (int j = 0; j < 4 && j < (int)f.size(); ++j) row[j] = stod(f[j]);
                    table.push_back(row);
                }
                acc0 = move(table);
                acc_len = (int)acc0.size();
            } else if (t.rfind("Acceptor 1HMM", 0) == 0) {
                vector<vector<double>> table;
                while (getline(fin, line)) {
                    if (line.find("NN TRM") != string::npos) break;
                    string s = ltrim(line);
                    auto f = split_ws(s);
                    if (f.empty()) continue;
                    vector<double> row(16, 0.0);
                    for (int j = 0; j < 16 && j < (int)f.size(); ++j) row[j] = stod(f[j]);
                    table.push_back(row);
                }
                acc1 = move(table);
            } else if (t.rfind("Acceptor 2HMM", 0) == 0) {
                vector<vector<double>> table;
                while (getline(fin, line)) {
                    if (line.find("NN TRM") != string::npos) break;
                    string s = ltrim(line);
                    auto f = split_ws(s);
                    if (f.empty()) continue;
                    vector<double> row(64, 0.0);
                    for (int j = 0; j < 64 && j < (int)f.size(); ++j) row[j] = stod(f[j]);
                    table.push_back(row);
                }
                acc2 = move(table);
            }
        }
        return true;
    };

    // --- parse positive model if exists ---
    ifstream testpos(pos_model_path);
    if (testpos) {
        testpos.close();
        bool ok = parse_zoeHMM(pos_model_path,
                               donor_freq, donor_hmm_freq, donor_hmm2_freq,
                               acceptor_freq, acceptor_hmm_freq, acceptor_hmm2_freq,
                               donor_length, acceptor_length);
        if (!ok) {
            cerr << "Positive model not in expected zoeHMM format or parse failed: " << pos_model_path << "\n";
            // continue, but donor_length/acceptor_length may be zero
        }
    }

    // --- parse negative model if exists ---
    ifstream testneg(neg_model_path);
    if (testneg) {
        testneg.close();
        use_negatives = true;
        bool ok = parse_zoeHMM(neg_model_path,
                               donor_nfreq, donor_hmm_nfreq, donor_hmm2_nfreq,
                               acceptor_nfreq, acceptor_hmm_nfreq, acceptor_hmm2_nfreq,
                               donor_length, acceptor_length);
        if (!ok) {
            cerr << "Negative model not in expected zoeHMM format or parse failed: " << neg_model_path << "\n";
            use_negatives = false;
        }
    }

    // --- open output files ---
    ofstream fout_gt("out.gt.txt");
    ofstream fout_ag("out.ag.txt");
    if (!fout_gt || !fout_ag) {
        cerr << "Cannot open output files\n";
        return 1;
    }

    // --- main loop over scaffolds ---
    for (auto &kv : genome_seqs) {
        const string &g = kv.first;
        string seq_fwd = toupper_str(kv.second);
        string seq_rev = seq_fwd;
        // reverse complement for completeness (not used in forward-only code)
        for (char &c : seq_rev) {
            if (c == 'A') c = 'T';
            else if (c == 'C') c = 'G';
            else if (c == 'G') c = 'C';
            else if (c == 'T') c = 'A';
            else c = 'N';
        }
        reverse(seq_rev.begin(), seq_rev.end());

        vector<int> don_fwd_pos;
        vector<int> acc_fwd_pos;
        unordered_map<int,double> acceptor_fwd_hmm2_score;
        unordered_map<int,double> don_fwd_hmm2_score;

        // find donors (GT)
        for (size_t p = 0; ; ) {
            size_t pos = seq_fwd.find("GT", p);
            if (pos == string::npos) break;
            size_t after = pos + 2;
            if ((int)after > donor_length - 5) { // same condition as perl
                // store 0-based start index (pos)
                don_fwd_pos.push_back((int)pos);
            }
            p = pos + 1;
        }

        // find acceptors (AG)
        for (size_t p = 0; ; ) {
            size_t pos = seq_fwd.find("AG", p);
            if (pos == string::npos) break;
            size_t after = pos + 2;
            if ((int)after > acceptor_length - 3) {
                acc_fwd_pos.push_back((int)pos);
            }
            p = pos + 1;
        }

        // compute donor scores
        for (int pos : don_fwd_pos) {
            // donor_seq = substr(seq_fwd, pos-3, donor_length)
            int start = pos - 3;
            if (start < 0) continue;
            if (donor_length <= 0) continue;
            if (start + donor_length > (int)seq_fwd.size()) continue;
            string donor_seq = seq_fwd.substr(start, donor_length);

            double donor_hmm2_score = 0.0;
            double donor_hmm2_nscore = 0.0;

            // sum donor_hmm2_freq
            if (!donor_hmm2_freq.empty()) {
                for (int i = 0; i < donor_length - 2; ++i) {
                    string tri = donor_seq.substr(i, 3);
                    auto it = code3.find(tri);
                    if (it != code3.end() && i < (int)donor_hmm2_freq.size()) {
                        int idx = it->second;
                        if (idx < (int)donor_hmm2_freq[i].size())
                            donor_hmm2_score += donor_hmm2_freq[i][idx];
                    }
                }
            }
            // add donor_hmm_freq[0][code2(substr(donor_seq,0,2))]
            if (!donor_hmm_freq.empty()) {
                string di = donor_seq.substr(0,2);
                auto it2 = code2.find(di);
                if (it2 != code2.end() && !donor_hmm_freq.empty()) {
                    int idx = it2->second;
                    if (0 < (int)donor_hmm_freq.size() && idx < (int)donor_hmm_freq[0].size())
                        donor_hmm2_score += donor_hmm_freq[0][idx];
                }
            }

            if (!use_negatives) {
                // keep donor_hmm2_score as is
            } else {
                if (!donor_hmm2_nfreq.empty()) {
                    for (int i = 0; i < donor_length - 2; ++i) {
                        string tri = donor_seq.substr(i, 3);
                        auto it = code3.find(tri);
                        if (it != code3.end() && i < (int)donor_hmm2_nfreq.size()) {
                            int idx = it->second;
                            if (idx < (int)donor_hmm2_nfreq[i].size())
                                donor_hmm2_nscore += donor_hmm2_nfreq[i][idx];
                        }
                    }
                }
                if (!donor_hmm_nfreq.empty()) {
                    string di = donor_seq.substr(0,2);
                    auto it2 = code2.find(di);
                    if (it2 != code2.end() && !donor_hmm_nfreq.empty()) {
                        int idx = it2->second;
                        if (0 < (int)donor_hmm_nfreq.size() && idx < (int)donor_hmm_nfreq[0].size())
                            donor_hmm2_nscore += donor_hmm_nfreq[0][idx];
                    }
                }
                donor_hmm2_score = donor_hmm2_score - donor_hmm2_nscore;
                don_fwd_hmm2_score[pos] = donor_hmm2_score;
            }
        }

        // compute acceptor scores
        for (int pos : acc_fwd_pos) {
            // acceptor_seq = substr(seq_fwd, pos-(acceptor_length-5), acceptor_length)
            if (acceptor_length <= 0) continue;
            int start = pos - (acceptor_length - 5);
            if (start < 0) continue;
            if (start + acceptor_length > (int)seq_fwd.size()) continue;
            string acceptor_seq = seq_fwd.substr(start, acceptor_length);

            double acceptor_hmm2_score = 0.0;
            double acceptor_hmm2_nscore = 0.0;

            if (!acceptor_hmm2_freq.empty()) {
                for (int i = 0; i < acceptor_length - 2; ++i) {
                    string tri = acceptor_seq.substr(i, 3);
                    auto it = code3.find(tri);
                    if (it != code3.end() && i < (int)acceptor_hmm2_freq.size()) {
                        int idx = it->second;
                        if (idx < (int)acceptor_hmm2_freq[i].size())
                            acceptor_hmm2_score += acceptor_hmm2_freq[i][idx];
                    }
                }
            }
            if (!acceptor_hmm_freq.empty()) {
                string di = acceptor_seq.substr(0,2);
                auto it2 = code2.find(di);
                if (it2 != code2.end() && !acceptor_hmm_freq.empty()) {
                    int idx = it2->second;
                    if (0 < (int)acceptor_hmm_freq.size() && idx < (int)acceptor_hmm_freq[0].size())
                        acceptor_hmm2_score += acceptor_hmm_freq[0][idx];
                }
            }

            if (!use_negatives) {
                acceptor_hmm2_score = acceptor_hmm2_score / 2.0;
            } else {
                if (!acceptor_hmm2_nfreq.empty()) {
                    for (int i = 0; i < acceptor_length - 2; ++i) {
                        string tri = acceptor_seq.substr(i, 3);
                        auto it = code3.find(tri);
                        if (it != code3.end() && i < (int)acceptor_hmm2_nfreq.size()) {
                            int idx = it->second;
                            if (idx < (int)acceptor_hmm2_nfreq[i].size())
                                acceptor_hmm2_nscore += acceptor_hmm2_nfreq[i][idx];
                        }
                    }
                }
                if (!acceptor_hmm_nfreq.empty()) {
                    string di = acceptor_seq.substr(0,2);
                    auto it2 = code2.find(di);
                    if (it2 != code2.end() && !acceptor_hmm_nfreq.empty()) {
                        int idx = it2->second;
                        if (0 < (int)acceptor_hmm_nfreq.size() && idx < (int)acceptor_hmm_nfreq[0].size())
                            acceptor_hmm2_nscore += acceptor_hmm_nfreq[0][idx];
                    }
                }
                acceptor_hmm2_score = (acceptor_hmm2_score - acceptor_hmm2_nscore) / 1.8;
                acceptor_fwd_hmm2_score[pos] = acceptor_hmm2_score;
            }
        }

        // thresholding and printing
        double splice_t = -1.0;
        double donor_mult = 30;
        double acceptor_mult = 30;
        for (int i = 0; i < (int)seq_fwd.size(); ++i) {
            auto it = don_fwd_hmm2_score.find(i);
            if (it != don_fwd_hmm2_score.end()) {
                double val = it->second;
                val = (val > splice_t ) ? val * donor_mult : -1e3;
                fout_gt << i << "\t" << val << "\n";
            }
        }
        for (int i = 0; i < (int)seq_fwd.size(); ++i) {
            auto it = acceptor_fwd_hmm2_score.find(i);
            if (it != acceptor_fwd_hmm2_score.end()) {
                double val = it->second;
                val = (val > splice_t ) ? val * acceptor_mult : -1e3;
                fout_ag << i << "\t" << val << "\n";
            }
        }
    } // end scaffolds

    fout_gt.close();
    fout_ag.close();

    return 0;
}

