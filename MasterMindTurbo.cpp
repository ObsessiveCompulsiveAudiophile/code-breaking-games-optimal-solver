#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <numeric>      // For std::iota
#include <algorithm>    // For std::sort
#include <cmath>        // For floor, round, pow
#include <chrono>       // For timing

#define uFast uint_fast64_t

// Optimal Code Breaker
// by Serkan Gur 2025

// --- CUSTOM GAME CONFIGURATION ---
const uFast n = 4;
const uFast c = 6;
const uFast maxdepth = 7;

// --- COMPILE-TIME CONSTANTS ---
template <uFast A, uFast B> struct get_power { static const uFast value = A * get_power<A, B - 1>::value; };
template <uFast A> struct get_power<A, 0> { static const uFast value = 1; };
const uFast s = get_power<c, n>::value; // Number of combinations (c^n)
const uFast imark0_val = (uFast)(double(n) * ((double(n) / 2.0) + 1.5));

// Replace slow decimal math with fast, packed base-c integers.
struct code_t {
    uFast v;
    // Decode into n digits (base c)
    void unpack(std::array<uFast, n>& d) const {
        uFast x = v;
        for (uFast i = 0; i < n; ++i) { d[i] = x % c; x /= c; }
    }
    // For printing results in the original format
    uFast to_decimal() const {
        uFast dec = 0, p10 = 1;
        std::array<uFast, n> digits;
        unpack(digits);
        for (uFast i = 0; i < n; ++i) { dec += (digits[i] + 1) * p10; p10 *= 10; }
        return dec;
    }
    // For parsing user input
    static code_t from_decimal(uFast dec) {
        code_t code{ 0 }; uFast p_c = 1;
        for (uFast i = 0; i < n; ++i) {
            uFast digit = (dec % 10) - 1;
            code.v += digit * p_c;
            p_c *= c; dec /= 10;
        }
        return code;
    }
};

static code_t Valids[s + 1];

// New FillSet for base-c code_t
void FillSet(uFast current_code, uFast pos, uFast& counter) {
    if (pos == n) { Valids[counter++] = { current_code }; return; }
    for (uFast i = 0; i < c; ++i) FillSet(current_code * c + i, pos + 1, counter);
}

// --- OPTIMIZATION 2: REWRITTEN SCORING KERNEL ---
static uFast iMark[imark0_val + 1];
static uFast mark_to_idx[n + 1][n + 1];

// Create a lookup table to map (black, white) pegs to the score index
void create_mark_map() {
    uFast k = 0;
    for (uFast p = 0; p <= n; ++p) { // p = black pegs
        for (uFast m = 0; m <= n - p; ++m) { // m = white pegs
            k++;
            iMark[k] = p * 10 + m;
            mark_to_idx[p][m] = k;
        }
    }
    iMark[0] = k;
    iMark[imark0_val] = n * 10;
    mark_to_idx[n][0] = imark0_val;
}

// Fast, integer-only scoring function (optimized version)
inline uFast score(code_t guess, code_t secret) {
    if (guess.v == secret.v) return imark0_val;

    uFast g_val = guess.v, s_val = secret.v;
    uFast black = 0;

    uFast g_counts[c] = { 0 };
    uFast s_counts[c] = { 0 };

    // Single pass through positions without intermediate arrays
    for (uFast i = 0; i < n; ++i) {
        uFast g_digit = g_val % c;
        uFast s_digit = s_val % c;

        if (g_digit == s_digit) {
            black++;
        }
        else {
            g_counts[g_digit]++;
            s_counts[s_digit]++;
        }
        g_val /= c;
        s_val /= c;
    }

    uFast white = 0;
    for (uFast i = 0; i < c; ++i) {
        white += std::min(g_counts[i], s_counts[i]);
    }

    return mark_to_idx[black][white];
}

// --- Custom Hash Table for Unique Partition Filtering ---
constexpr uFast HT_SIZE = 2048; // Power of two >= s, keeps load factor low
alignas(64) uFast ht_key[HT_SIZE];
alignas(64) uFast ht_val[HT_SIZE]; // Stores the compact index (1 to k)

// Returns the 1-based compact index for a given key.
// If the key is new, it increments *k_ptr and assigns the new index.
inline uFast find_or_add_signature(uFast key, uFast* k_ptr) {
    // Fast multiplicative hash using a golden-ratio-based constant
    uFast idx = (key * 0x9e3779b97f4a7c15ULL) & (HT_SIZE - 1);

    while (true) {
        if (ht_key[idx] == key) { // Found existing key
            return ht_val[idx];
        }
        if (ht_key[idx] == 0) { // Found empty slot, this is a new key
            ht_key[idx] = key;
            (*k_ptr)++;
            ht_val[idx] = *k_ptr;
            return *k_ptr;
        }
        // Collision: move to the next slot (linear probing)
        idx = (idx + 1) & (HT_SIZE - 1);
    }
}


int main() {
    uFast i, j, k, l, m, p, iv, iv2, o, o2;

    const uFast imark0 = imark0_val; // Use the calculated value locally
    create_mark_map();

    // Heap-allocated storage
    const uFast S_DIM = s + 1, IMARK0_DIM = imark0 + 1;
    auto MapConsistent_ptr = std::make_unique<uFast[]>((maxdepth + 1) * S_DIM * IMARK0_DIM);
    auto MapConsistent = [&](uFast d1, uFast d2, uFast d3) -> uFast& { return MapConsistent_ptr[d1 * S_DIM * IMARK0_DIM + d2 * IMARK0_DIM + d3]; };

    typedef uFast(iFiltered_)[s + 1]; iFiltered_* iFiltered = new iFiltered_[maxdepth + 1];
    typedef uFast(iValid_)[s + 1]; iValid_* iValid = new iValid_[maxdepth + 1];

    alignas(64) static uint16_t Mark[s + 1][s + 1];

    auto RowSum = std::make_unique<uFast[]>(1'000'000);
    typedef uFast(list_)[1'000'000]; list_* list = new list_[maxdepth + 1];

    uFast ubPure[maxdepth + 1] = { 0 }, g[maxdepth + 1] = { 0 }, r[maxdepth + 1] = { 0 }, GuessSum[s + 1] = { 0 };
    uFast guess, sign, lvl = 1, x = 0, x0 = 0;
    uFast gMin = -1; // Initialize to max value for unsigned
    std::array<uFast, s + 1> MPScore;
    static std::array<uFast, imark0 + 1> mult;
    for (i = 0; i < imark0; i++) mult[i + 1] = uFast(round(pow(16.0, double(i))));

    // --- Generate codes and Mark table ---
    uFast i00 = 1;
    FillSet(0, 0, i00);
    for (i = 1; i <= s; i++) { iValid[1][i] = i; iFiltered[1][i] = i; }

    std::cout << "Pre-computing " << s << "x" << s << " score table...\n";
    auto mark_start_time = std::chrono::high_resolution_clock::now();

    for (i = 1; i <= s; i++) {
        for (j = i; j <= s; j++) {
            uFast mark_val = score(Valids[i], Valids[j]);
            Mark[i][j] = mark_val;
            Mark[j][i] = mark_val;
        }
    }
    auto mark_end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Mark table computed in " << std::chrono::duration<double>(mark_end_time - mark_start_time).count() << " seconds.\n";

    // --- Interactive Input ---
    std::vector<uFast> initial_guesses, initial_signs;
    uFast input_count = 1;
    std::cout << "\n";
    while (true) {
        std::cout << "Enter " << input_count << ". known guess(gggg) AND sign(pm) in 'ggggpm' format OR '0' to solve: ";
        uFast userinput;
        std::cin >> userinput;
        if (userinput == 0) break;
        uFast gg_dec = floor(userinput / 100);
        uFast ss_val = userinput - gg_dec * 100;
        bool found_g = false, found_s = false;
        code_t packed_guess = code_t::from_decimal(gg_dec);
        for (uFast idx = 1; idx <= s; ++idx) if (Valids[idx].v == packed_guess.v) { initial_guesses.push_back(idx); found_g = true; break; }
        for (uFast idx = 1; idx <= imark0; ++idx) if (iMark[idx] == ss_val) { initial_signs.push_back(idx); found_s = true; break; }
        if (!found_g || !found_s) std::cout << "Invalid guess or sign. Please try again.\n";
        else input_count++;
    }

    uFast kk = s;
    for (size_t idx = 0; idx < initial_guesses.size(); ++idx) {
        m = 0;
        for (i = 1; i <= kk; i++) {
            if (Mark[iValid[1][i]][initial_guesses[idx]] == initial_signs[idx]) {
                m++;
                iValid[1][m] = iValid[1][i];
            }
        }
        kk = m;
    }

    if (kk == 0) { std::cout << "No solution with that, please try again!\n\n"; return 1; }
    std::cout << "\n" << kk << " possible solutions!\n\n" << "Calculating...\n\n";
    if (kk == 1) {
        RowSum[1] = 1;
        list[1][1] = Valids[iValid[1][1]].to_decimal() * 100 + iMark[imark0];
        goto finished;
    }

    auto solve_start_time = std::chrono::high_resolution_clock::now();
    MapConsistent(0, 0, 0) = kk;
    for (i = 1; i <= maxdepth; i++) {
        g[i] = 1;
        r[i] = 1;
    }

GetNext:
    // --- Partition Calculation (Original, Fast Version) ---
    memset(&MapConsistent(lvl, 0, 0), 0, sizeof(uFast) * S_DIM * IMARK0_DIM);
    uFast q[s + 1] = { 0 };
    k = lvl - 1;
    l = MapConsistent(k, g[k], r[k]);
    iv = iValid[lvl][1];
    iv2 = iValid[lvl][2];
    for (i = 1; i <= s; i++) {
        o = Mark[i][iv]; o2 = Mark[i][iv2];
        MapConsistent(lvl, i, o)++;
        MapConsistent(lvl, i, o2)++;
        q[i] += (mult[o] + mult[o2]); // Simplified from original for clarity
    }
    for (j = 3; j <= l; j++) {
        iv = iValid[lvl][j];
        for (i = 1; i <= s; i++) {
            o = Mark[i][iv];
            MapConsistent(lvl, i, o)++;
            q[i] += mult[o];
        }
    }

    // --- High-performance filtering using a custom hash table ---
    std::fill_n(ht_key, HT_SIZE, 0);
    k = 0;
    for (i = 1; i <= s; i++) {
        uFast unique_idx = find_or_add_signature(q[i], &k);
        if (unique_idx == k) { // True only for the first time this signature is seen
            iFiltered[lvl][k] = i;
        }
    }

    // --- Compaction and Scoring ---
    for (i = 1; i <= k; i++) {
        MPScore[i] = 0;
        uFast original_idx = iFiltered[lvl][i];
        for (j = 1; j <= imark0; j++) {
            uFast count = MapConsistent(lvl, original_idx, j);
            MapConsistent(lvl, i, j) = count;
            if (count != 0) {
                MPScore[i]++;
            }
        }
        q[i] = MPScore[i]; // Reuse q[] to hold the score
    }

    std::sort(&MPScore[1], &MPScore[k + 1]);

    // --- Heuristic Pruning ---
    l = (k > 119) ? k - 8 : 0.88 * k;

    m = 0;
    for (i = 1; i <= k; i++) {
        if (q[i] >= MPScore[l]) {
            m++;
            // This compaction is now safe because iFiltered[lvl][i] holds the original index
            uFast original_idx = iFiltered[lvl][i];
            iFiltered[lvl][m] = original_idx;
            for (j = 1; j <= imark0; j++) {
                MapConsistent(lvl, m, j) = MapConsistent(lvl, i, j);
            }
        }
    }
    ubPure[lvl] = m;

    // --- Main DFS Loop ---
    while (lvl) {
        for (guess = g[lvl]; guess <= ubPure[lvl]; guess++) {
            for (sign = r[lvl]; sign <= imark0; sign++) {
                if (MapConsistent(lvl, guess, sign) == 0) continue;
                g[lvl] = guess;
                r[lvl] = sign;
                if (sign == imark0) {
                    x++;
                    RowSum[x] = lvl;
                    for (i = 1; i <= lvl; i++) list[i][x] = Valids[iFiltered[i][g[i]]].to_decimal() * 100 + iMark[r[i]];
                    break;
                }
                j = 0;
                l = lvl - 1;
                k = lvl + 1;
                for (i = 1; i <= MapConsistent(l, g[l], r[l]); i++) if (Mark[iValid[lvl][i]][iFiltered[lvl][guess]] == sign) {
                    j++;
                    iValid[k][j] = iValid[lvl][i];
                }
                lvl++;
                if (lvl == maxdepth) {
                    for (i = 1; i <= MapConsistent(maxdepth - 1, g[maxdepth - 1], r[maxdepth - 1]); i++) RowSum[x + i] = maxdepth;
                    x += (i - 1);
                    break;
                }
                if (j == 1) {
                    ubPure[lvl] = 1; iFiltered[lvl][1] = iValid[lvl][1]; MapConsistent(lvl, 1, imark0) = 1;
                    g[lvl] = 1;
                    r[lvl] = imark0;
                    x++;
                    RowSum[x] = lvl;
                    for (i = 1; i <= lvl; i++) list[i][x] = Valids[iFiltered[i][g[i]]].to_decimal() * 100 + iMark[r[i]];
                    break;
                }
                goto GetNext;
            }
            if ((g[lvl] != 1) && (ubPure[lvl] == g[lvl])) {
                gMin = -1; // Max uFast value
                l = lvl - 1;
                p = MapConsistent(l, g[l], r[l]);
                x0 = x - (ubPure[lvl] * p);
                for (i = 1; i <= ubPure[lvl]; i++) {
                    uFast sum = 0;
                    uFast baseOffset = x0 + (i - 1) * p;
                    for (j = 1; j <= p; j++) sum += RowSum[baseOffset + j];
                    GuessSum[i] = sum;
                    if (sum < gMin) {
                        gMin = sum;
                        k = i;
                    }
                }
                uFast sourceBase = x0 + (k - 1) * p;
                for (i = 1; i <= p; i++) {
                    uFast destIndex = x0 + i;
                    uFast sourceIndex = sourceBase + i;
                    RowSum[destIndex] = RowSum[sourceIndex];
                    for (j = 1; j < maxdepth; j++) list[j][destIndex] = list[j][sourceIndex];
                }
                x = x0 + p;
            }
            r[lvl] = 1;
        }
        g[lvl] = 1;
        r[lvl] = 1;
        lvl--;
        if (lvl > 0) r[lvl]++;
    }

finished:
    auto solve_end_time = std::chrono::high_resolution_clock::now();
    j = 0; for (i = 1; i <= kk; i++) j += RowSum[i];
    for (i = 1; i <= kk; i++) { for (k = 1; k <= RowSum[i]; k++) { std::cout << list[k][i] << "\t"; } std::cout << ":" << RowSum[i] << "\n"; }
    double total_time = std::chrono::duration<double>(solve_end_time - solve_start_time).count();
    std::cout << j << " / " << kk << " = " << (double)j / kk << " average moves to finish!\n" << total_time << " seconds " << "\n\a\a\a";
    std::cout << "Press Enter to exit...";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
    return 0;
}