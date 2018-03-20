#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <bitset>
#include <vector>
#include <map>
#include <utility>
#include <set>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>
#include <cmath>

using uint64_t = std::uint64_t;
class number; // for exact representation of rational numbers
using srt_table_t = std::map<std::pair<number, number>, std::vector<int>>;

// information to generate a certain table instance
struct srt_info_t {
    int digit_range;
    int radix;
    int p_fractional_digits;
    // number of digits to skip in the last part of p
    // e.g. r = 9 and value = 3, -> possible values of the end: 0, 3, 6
    // value * x = r
    int p_fractional_digits_part;
    int d_fractional_digits;
    // like p_fractional_digits_part
    int d_fractional_digits_part;
    // normally divisor is in range [1, r[, but it can be normalized differently
    // which can decrease the range to [1, ndb[
    int normalized_divisor_bound;
    // is the partial remainder represented as a sum of two values calculated after approximation?
    bool sum;
};

// printing information
struct srt_print_info_t {
    std::string comment_seq;
    std::string oor_seq; // oor = out-of-range
    std::string separator;
    std::string inner_separator;
    std::string begin_multiple;
    std::string end_multiple;
    enum class multiple_possibilities {h, l, a, v, u} mp;
    bool print_column_comments;
    bool print_row_comments;
    bool print_info_comments;
    int comment_radix;
};

// info for string generation of a number
struct srt_radix_representation_info_t {
    int radix;
    int digits_below_1;
    // pad with zeroes and positive numbers also with a leading '+'
    bool pad;
    // if complement is true, pad will be assumed true, but padding occurs without signs
    bool complement;
    // if either pad or complement is true we have to know the total length
    int total_length;
};

// ======================================================================================
// helper functions
// ======================================================================================

// the integer dividers of val in [2, val[ in descending order
// but with a 0 at the first place, simplifies later implementation
std::vector<int> dividers(int val) {
    std::vector<int> ret {0};
    for(int i = 2; i < val; ++i) {
        int a = val / i;
        if(a * i == val) ret.push_back(a);
    }
    return ret;
}

// return a new srt_info_t with a higher resolution of p,
// that is either with a higher number of fractional digits
// or with a lower fractional digit part
srt_info_t increase_p(srt_info_t info) {
    unsigned part_id = 0;
    auto divs = dividers(info.radix);
    // get the current id of the digit part in divisors
    for(auto i = 0u; i < divs.size(); ++i)
        if(info.p_fractional_digits_part == divs[i])
            part_id = i;
    if((++part_id) == divs.size()) {
        ++info.p_fractional_digits;
        part_id = 0;
    }
    info.p_fractional_digits_part = divs[part_id];
    return info;
}

// return a new srt_info with a higher resolution of d,
// according to 'increase_p'
srt_info_t increase_d(srt_info_t info) {
    unsigned part_id = 0;
    auto divs = dividers(info.radix);
    for(auto i = 0u; i < divs.size(); ++i)
        if(info.d_fractional_digits_part == divs[i])
            part_id = i;
    if((++part_id) == divs.size()) {
        ++info.d_fractional_digits;
        part_id = 0;
    }
    info.d_fractional_digits_part = divs[part_id];
    return info;
}

// return b^e
int pow(int b, int e) {
    int res = 1;
    for(int i = 0; i < e; ++i) {
        res *= b;
    }
    return res;
}

// i to string where i > 10 is represented by a-z,
// i is one symbol
// handles negative values
std::string to_str(int i) {
    std::string ret = i < 0 ? "-" : "";
    i = i < 0 ? -i : i;
    if(i > 36) return "";
    if(i < 10) ret += '0' + i;
    else ret += 'a' + i - 10;
    return ret;
}

// similar to 'to_str', does not handle negative values though
char to_char(int i) {
    if(i > 36 || i < 0) return 0;
    if(i < 10) return '0' + i;
    return 'a' + i - 10;
}

// base ^ return = val
// ceiled to the next int if it's not exactly an int
int log(int base, int val) {
    int ret = 0;
    for(int i = 1; i < val; i *= base) ++ret;
    return ret;
}

// number of fractional steps between two successing integers
int fractional_steps(int radix, int fractional_digits, int fractional_digits_part) {
    int mult = fractional_digits_part ? radix / fractional_digits_part : 1;
    return pow(radix, fractional_digits) * mult;
}

// convenience functions
int fractional_steps_p(const srt_info_t &info) {
    return fractional_steps(info.radix, info.p_fractional_digits, info.p_fractional_digits_part);
}

int fractional_steps_d(const srt_info_t &info) {
    return fractional_steps(info.radix, info.d_fractional_digits, info.d_fractional_digits_part);
}

// padding functions for int and string
std::string pad(int val, int len, char pad_char = ' ') {
    std::string s = std::to_string(val);
    for(int i = s.size(); i < len; ++i) s = pad_char + s;
    return s;
}

std::string pad(std::string val, int len, char pad_char = ' ') {
    for(int i = val.size(); i < len; ++i) val = pad_char + val;
    return val;
}

// ======================================================================================
// exact rational numbers
// ======================================================================================

class number {
    // arithmetic operators
    friend number operator /(const number &, const number &);
    friend number operator *(const number &, const number &);
    friend number operator -(const number &, const number &);
    friend number operator +(const number &, const number &);
    // comparison operators
    friend bool operator <(const number &, const number &);
    friend bool operator <=(const number &, const number &);
    friend bool operator >(const number &, const number &);
    friend bool operator >=(const number &, const number &);
    friend bool operator ==(const number &, const number &);
    // other operators
    friend number operator -(const number &);
    public:
        explicit number();
        explicit number(int64_t n);

        explicit operator int64_t();
        explicit operator double();
    protected:
        number(uint64_t n, uint64_t d, bool pos);
        void reduce();

        bool pos;
        uint64_t numerator;
        uint64_t denominator;
};

number::number() : number::number(0) {}

number::number(int64_t n) {
    pos = n >= 0;
    numerator = n < 0 ? -n : n;
    denominator = 1;
}

number::number(uint64_t n, uint64_t d, bool p) {
    pos = p;
    numerator = n;
    denominator = d;
    reduce();
}

number::operator int64_t() {
    return (pos ? 1 : -1) * numerator / denominator;
}

number::operator double() {
    return (pos ? 1 : -1) * ((double)numerator / denominator);
}

// reduce the fraction
void number::reduce() {
    if(numerator == 0) {
        denominator = 1;
        pos = true;
        return;
    }
    if(denominator == 0) return; // should never happen
    uint64_t gcd = 0;
    while((gcd = std::gcd(numerator, denominator)) > 1) {
        numerator /= gcd;
        denominator /= gcd;
    }
}

// arithmetical operators
number operator /(const number &n1, const number &n2) {
    return number(n1.numerator * n2.denominator, n1.denominator * n2.numerator, n1.pos == n2.pos);
}

number operator *(const number &n1, const number &n2) {
    return number(n1.numerator * n2.numerator, n1.denominator * n2.denominator, n1.pos == n2.pos);
}

number operator -(const number &n1, const number &n2) {
    return n1 + (-n2);
}

number operator +(const number &n1, const number &n2) {
    uint64_t lcm = std::lcm(n1.denominator, n2.denominator);
    uint64_t _n1 = n1.numerator * (lcm / n1.denominator);
    uint64_t _n2 = n2.numerator * (lcm / n2.denominator);
    if(n1.pos && n2.pos) {
        return number(_n1 + _n2, lcm, true);
    } else if(n1.pos && !n2.pos) {
        if(_n1 < _n2) {
            return number(_n2 - _n1, lcm, false);
        } else {
            return number(_n1 - _n2, lcm, true);
        }
    } else if(!n1.pos && n2.pos) {
        if(_n1 > _n2) {
            return number(_n1 - _n2, lcm, false);
        } else {
            return number(_n2 - _n1, lcm, true);
        }
    } else {
        return number(_n1 + _n2, lcm, false);
    }
}

// comparison operators
bool operator <(const number &n1, const number &n2) {
    if(n1.pos && !n2.pos) return false;
    if(!n1.pos && n2.pos) return true;
    uint64_t lcm = std::lcm(n1.denominator, n2.denominator);
    uint64_t _n1 = n1.numerator * (lcm / n1.denominator);
    uint64_t _n2 = n2.numerator * (lcm / n2.denominator);
    if(n1.pos) return _n1 < _n2;
    return _n2 < _n1;
}

bool operator <=(const number &n1, const number &n2) {
    if(n1.pos && !n2.pos) return false;
    if(!n1.pos && n2.pos) return true;
    uint64_t lcm = std::lcm(n1.denominator, n2.denominator);
    uint64_t _n1 = n1.numerator * (lcm / n1.denominator);
    uint64_t _n2 = n2.numerator * (lcm / n2.denominator);
    if(n1.pos) return _n1 <= _n2;
    return _n2 <= _n1;
}

bool operator >(const number &n1, const number &n2) {
    return n2 < n1;
}

bool operator >=(const number &n1, const number &n2) {
    return n2 <= n1;
}

bool operator ==(const number &n1, const number &n2) {
    return n1 <= n2 && n2 <= n1;
}

// other operators
number operator -(const number &n) {
    return number(n.numerator, n.denominator, !n.pos);
}

// ======================================================================================
// srt table generation
// ======================================================================================

// the generated divisor is in range [1, ndr[, the partial remainder is scaled appropriately
srt_table_t gen_srt_table(srt_info_t info) {
    srt_table_t table;
    number eta = number(info.digit_range) / (number(info.radix) - number(1));
    std::set<int> quotient_digits;
    for(int i = -info.digit_range; i <= info.digit_range; ++i) {
        quotient_digits.insert(i);
    }

    // the *_fractional_digits_part has to be a divider of the radix
    if(info.p_fractional_digits_part &&
            (info.radix / info.p_fractional_digits_part) * info.p_fractional_digits_part != info.radix)
        return table;
    if(info.d_fractional_digits_part &&
            (info.radix / info.d_fractional_digits_part) * info.d_fractional_digits_part != info.radix)
        return table;
    // partial remainder value delta
    number ps = number(1) / number(fractional_steps_p(info));
    // divisor value delta
    number ds = number(1) / number(fractional_steps_d(info));

    // if residual is represented as sum we have to lower each bound (except the topmost one) by this value
    // the quotient selection function is:
    // −eta · d · r − lower_bounds_by <= r · (p − q · d) < eta · d · r − lower_bounds_by
    // we divide everything by r, so also lower_bounds_by
    number lower_bounds_by = number(0);
    if(info.sum) {
        lower_bounds_by = ps;
    }

    auto possible_digits = [=](number p, number d) {
        // set the bool to true if an error occurs, that is if another smaller ps or ds is necessary
        std::pair<bool, std::vector<int>> ret = {false, {}};
        std::set<int> valid_digits;
        const unsigned corner_bl = 0, corner_tl = 1, corner_br = 2, corner_tr = 3;
        // the corner descriptors for the specific quotient digit
        // -1 -> corner is not in the range for the quotient digit
        //  0 -> corner is at edge of the range
        // +1 -> corner is in the range
        std::map<int, int[4]> digit_corners;
        // the minimum and maximum (excluded) values for p and d in the current cell
        number p1 = p;
        number p2 = p + ps;
        number d1 = d;
        number d2 = d + ds;

        for(int digit: quotient_digits) {
            number lower_ub = digit == info.digit_range ? number(0) : lower_bounds_by;
            number lower_lb = digit == -info.digit_range ? lower_bounds_by : number(0);
            number digitn = number(digit);
            int corners[4] = {-1,-1,-1,-1};
            // corner bottom left
            if(-eta * d1 - lower_lb == p1 - digitn * d1 ||
                // we need to do this, but this value is actually never reached
                eta * d1 - lower_ub == p1 - digitn * d1)
                corners[corner_bl] = 0;
            else if(-eta * d1 - lower_lb < p1 - digitn * d1 &&
                     eta * d1 - lower_ub > p1 - digitn * d1)
                corners[corner_bl] = 1;
            // corner top left
            if(-eta * d1 - lower_lb == p2 - digitn * d1 ||
                eta * d1 - lower_ub == p2 - digitn * d1)
                corners[corner_tl] = 0;
            else if(-eta * d1 - lower_lb < p2 - digitn * d1 &&
                     eta * d1 - lower_ub > p2 - digitn * d1)
                corners[corner_tl] = 1;
            // corner bottom right
            if(-eta * d2 - lower_lb == p1 - digitn * d2 ||
                eta * d2 - lower_ub == p1 - digitn * d2)
                corners[corner_br] = 0;
            else if(-eta * d2 - lower_lb < p1 - digitn * d2 &&
                     eta * d2 - lower_ub > p1 - digitn * d2)
                corners[corner_br] = 1;
            // corner top right
            if(-eta * d2 - lower_lb == p2 - digitn * d2 ||
                eta * d2 - lower_ub == p2 - digitn * d2)
                corners[corner_tr] = 0;
            else if(-eta * d2 - lower_lb < p2 - digitn * d2 &&
                     eta * d2 - lower_ub > p2 - digitn * d2)
                corners[corner_tr] = 1;

            if(corners[0] != -1 || corners[1] != -1 || corners[2] != -1 || corners[3] != -1)
                for(int i = 0; i < 4; ++i)
                    digit_corners[digit][i] = corners[i];
        }
        // if no corners were found in or at the ranges, return
        if(!digit_corners.size()) return ret;
        // find valid digits
        if(p >= number(0)) {
            // include all digits where bl and tr are in the range and tl and br are at least touched
            // problem with the line p = 0, therefore its enough for bl to be touched too
            for(auto e: digit_corners)
                if(e.second[corner_bl] >= 0 && e.second[corner_tr] == 1 &&
                        e.second[corner_tl] >= 0 && e.second[corner_br] >= 0)
                    valid_digits.insert(e.first);
            // include the maximum digit also if only br is in range
            if(digit_corners.find(info.digit_range) != digit_corners.end())
                if(digit_corners[info.digit_range][corner_br] == 1)
                    valid_digits.insert(info.digit_range);
        } else {
            // include all digits where tl and br are in the range and bl and tr are at least touched
            // same problem as above
            for(auto e: digit_corners)
                if(e.second[corner_tl] >= 0 && e.second[corner_br] == 1 &&
                        e.second[corner_bl] >= 0 && e.second[corner_tr] >= 0)
                    valid_digits.insert(e.first);
            // include the minimum digit also if only tr is in range
            if(digit_corners.find(-info.digit_range) != digit_corners.end())
                if(digit_corners[-info.digit_range][corner_tr] == 1)
                    valid_digits.insert(-info.digit_range);
        }
        // if we found some corners inside ranges but no valid digits we have to repeat with smaller ps or ds
        // one touched corner is ok though
        if(!valid_digits.size()) {
            for(auto e: digit_corners) {
                int cross_sum = 0;
                for(int i = 0; i < 4; ++i)
                    cross_sum += e.second[i];
                if(cross_sum > -3) {
                    ret.first = true;
                    return ret;
                }
            }
        }
        // return the valid digits
        for(int i: valid_digits) {
            ret.second.push_back(i);
        }
        return ret;
    };

    for(number d = number(1); d < number(info.normalized_divisor_bound); d = d + ds) {
        for(number p = number(0); true; p = p + ps) {
            auto pd = possible_digits(p, d);
            if(pd.first) return {};
            if(pd.second.size())
                table[std::make_pair(p, d)] = pd.second;
            else
                break;
        }
        for(number p = -ps; true; p = p - ps) {
            auto pd = possible_digits(p, d);
            if(pd.first) return {};
            if(pd.second.size())
                table[std::make_pair(p, d)] = pd.second;
            else
                break;
        }
    }

    return table;
}

// not necessarily really the 'smallest' possible table
// but a quite small one
std::pair<srt_info_t, srt_table_t> smallest_srt_table(srt_info_t start, int max_iterations) {
    auto srt_info_cmp = [](const srt_info_t &left, const srt_info_t &right) {
        // for comparison only take partial remainders in [0, 1[,
        // the real number is proportional for a certain radix
        // but harder to evaluate
        int p_num_l = fractional_steps_p(left);
        int d_num_l = fractional_steps_d(left);

        int p_num_r = fractional_steps_p(right);
        int d_num_r = fractional_steps_d(right);

        return p_num_l * d_num_l < p_num_r * d_num_r;
    };
    // 'priority queue'
    std::multiset<srt_info_t, decltype(srt_info_cmp)> info_set(srt_info_cmp);
    info_set.insert(start);
    std::pair<srt_info_t, srt_table_t> ret;
    for(int i = 0; i < max_iterations; ++i) {
        srt_info_t si = *info_set.begin();
        info_set.erase(info_set.begin());
        srt_table_t table;
        // test
        table = gen_srt_table(si);
        if(!table.empty()) {
            ret.first = si;
            ret.second = table;
            break;
        }
        // add other tests
        info_set.insert(increase_p(si));
        info_set.insert(increase_d(si));
    }
    return ret;
}

// ======================================================================================
// printing
// ======================================================================================

std::string radix_representation(number val, const srt_radix_representation_info_t &info) {
    bool neg = val < number(0);
    if(neg) val = -val;
    std::string ret;
    number helper = number(1);
    while(helper <= val) helper = helper * number(info.radix);
    helper = helper / number(info.radix);
    while(number(1) <= helper) {
        int64_t a = (int64_t)(val / helper);
        ret += to_char(a);
        val = val - number(a) * helper;
        helper = helper / number(info.radix);
    }
    if(ret.empty()) ret = "0";
    if(info.digits_below_1) ret += '.';
    for(int i = 0; i < info.digits_below_1; ++i) {
        int64_t a = (int64_t)(val / helper);
        ret += to_char(a);
        val = val - number(a) * helper;
        helper = helper / number(info.radix);
    }
    if(info.pad) {
        for(int i = ret.size(); i < info.total_length - 1; ++i) {
            ret = '0' + ret;
        }
    }
    // two's complement
    if(info.complement) {
        if(neg) {
            // dirty calculation with strings
            std::string new_ret;
            for(auto c: ret) {
                if(c == '0') new_ret.push_back('1');
                else if(c == '1') new_ret.push_back('0');
                else new_ret.push_back(c); // the fraction point
            }
            // 'add' 1
            auto it = new_ret.rbegin(); // id from the end, the position where we add 1
            for(; it != new_ret.rend(); ++it) {
                if(*it == '1') {
                    *it = '0';
                } else if(*it == '0') {
                    *it = '1';
                    break;
                }
            }
            ret = '1' + new_ret;
        } else {
            ret = '0' + ret;
        }
    } else if(neg) {
        ret = '-' + ret;
    } else if(info.pad) {
        ret = '+' + ret;
    }
    return ret;
}

// radix_rep_info_p is only for the partial remainders
void print_srt_table(std::ostream &os, const srt_table_t &table,
        const srt_info_t &info, const srt_print_info_t &print_info,
        srt_radix_representation_info_t radix_rep_info_p) {
    // all partial remainders and divisor values, sorted
    std::set<number, std::greater<number>> partial_remainders;
    std::set<number> divisors;
    for(auto e: table) {
        partial_remainders.insert(e.first.first);
        divisors.insert(e.first.second);
    }

    unsigned symbol_length = to_str(-info.digit_range).size();
    unsigned oor_symbol_length = print_info.oor_seq.size();
    unsigned digits_below_1_p = info.p_fractional_digits * log(print_info.comment_radix, info.radix);
    if(info.p_fractional_digits_part) digits_below_1_p +=
        log(print_info.comment_radix, info.radix / info.p_fractional_digits_part);
    unsigned digits_below_1_d = info.d_fractional_digits * log(print_info.comment_radix, info.radix);
    if(info.d_fractional_digits_part) digits_below_1_d +=
        log(print_info.comment_radix, info.radix / info.d_fractional_digits_part);

    // fill radix_rep_info
    radix_rep_info_p.digits_below_1 = digits_below_1_p;
    srt_radix_representation_info_t radix_rep_info_d = radix_rep_info_p;
    radix_rep_info_d.pad = radix_rep_info_d.complement = false;
    // get total length with a kind of cheat
    radix_rep_info_d.total_length =
        radix_representation(*divisors.rbegin(), radix_rep_info_d).size();
    radix_rep_info_p.total_length =
        radix_representation(*partial_remainders.rbegin(), radix_rep_info_d).size();
    radix_rep_info_d.digits_below_1 = digits_below_1_d;

    unsigned row_comment_pad_length = radix_rep_info_p.total_length;

    // compute the digit representation length, especially necessary if mutiple digits are valid
    unsigned max_valid_elem = 0;
    for(auto e: table) {
        max_valid_elem = max_valid_elem > e.second.size() ? max_valid_elem : e.second.size();
    }
    if(print_info.mp != srt_print_info_t::multiple_possibilities::a) max_valid_elem = 1;
    int entry_pad_length = std::max(max_valid_elem * symbol_length +
            (max_valid_elem == 1 ? 0u : print_info.begin_multiple.size() + print_info.end_multiple.size()) +
            print_info.inner_separator.size() * (max_valid_elem - 1),
            static_cast<std::string::size_type>(oor_symbol_length));

    // print info comments
    if(print_info.print_info_comments) {
        os << print_info.comment_seq << " radix: " << info.radix << std::endl;
        os << print_info.comment_seq << " digit range: " << info.digit_range << std::endl;
        os << print_info.comment_seq << " comment radix: " << print_info.comment_radix << std::endl;
        os << print_info.comment_seq << " columns: " << divisors.size() << std::endl;
        os << print_info.comment_seq << " rows: " << partial_remainders.size() << std::endl;
    }

    // print column comments
    if(print_info.print_column_comments) {
        std::vector<std::string> divisor_strings;
        for(auto d: divisors) {
            divisor_strings.push_back(radix_representation(d, radix_rep_info_d));
        }
        unsigned divisor_string_len = 0;
        for(auto s: divisor_strings) {
            divisor_string_len = divisor_string_len > s.size() ? divisor_string_len : s.size();
        }
        for(auto it = divisor_strings.begin(); it != divisor_strings.end(); ++it) {
            // manually pad divisors
            *it = pad(*it, divisor_string_len, radix_rep_info_p.pad ? '0' : ' ');
        }
        for(auto i = 0u; i < divisor_string_len; ++i) {
            os << print_info.comment_seq;
            for(unsigned j = 0; j < divisor_strings.size(); ++j) {
                os << pad(std::string(1, divisor_strings[j][i]), entry_pad_length);
                if(j != divisor_strings.size() - 1) os << pad("", print_info.separator.size());
            }
            os << std::endl;
        }
    }

    // print rows
    auto abs_compare = [](auto l, auto r) {
        return std::abs(l) < std::abs(r);
    };
    for(auto p: partial_remainders) {
        if(print_info.print_column_comments) {
            os << pad("", print_info.comment_seq.size());
        }
        for(auto d: divisors) {
            auto search = table.find(std::make_pair(p, d));
            std::set<int> digits;
            if(search != table.end()) {
                for(auto i: search->second) digits.insert(i);
            }
            std::string out;
            if(digits.size() > 1) {
                if(print_info.mp == srt_print_info_t::multiple_possibilities::a) {
                    out += print_info.begin_multiple;
                    for(auto i: digits) out += to_str(i) + print_info.inner_separator;
                    for(auto i = 0u; i < print_info.inner_separator.size(); ++i) out.pop_back();
                    out += print_info.end_multiple;
                    os << pad(out, entry_pad_length);
                } else if(print_info.mp == srt_print_info_t::multiple_possibilities::l) {
                    os << pad(to_str(*digits.begin()), entry_pad_length);
                } else if(print_info.mp == srt_print_info_t::multiple_possibilities::h) {
                    os << pad(to_str(*digits.rbegin()), entry_pad_length);
                } else if(print_info.mp == srt_print_info_t::multiple_possibilities::u) {
                    int elem = *std::min_element(digits.begin(), digits.end(), abs_compare);
                    os << pad(to_str(elem), entry_pad_length);
                } else if(print_info.mp == srt_print_info_t::multiple_possibilities::v) {
                    int elem = *std::max_element(digits.begin(), digits.end(), abs_compare);
                    os << pad(to_str(elem), entry_pad_length);
                }
            } else if(digits.size() == 1) {
                os << pad(to_str(*digits.begin()), entry_pad_length);
            } else {
                os << pad(print_info.oor_seq, entry_pad_length);
            }
            os << print_info.separator;
        }
        if(print_info.print_row_comments) {
            os << "  " << print_info.comment_seq << "  " <<
                pad(radix_representation(p, radix_rep_info_p), row_comment_pad_length);
        }
        os << std::endl;
    }
}

// ======================================================================================
// argument parsing and main
// ======================================================================================

using arguments_t = std::map<std::string, std::string>;

arguments_t parse_arguments(int argc, char **argv, arguments_t default_arguments) {
    arguments_t ret = default_arguments;
    enum class expect {
        arg, num, str
    } expecting = expect::arg;
    std::string cur_key = "";

    auto handle_arg = [&](std::string s) {
        if(s == "-h" || s == "--help") {
            ret["help"] = "true";
        } else if(s == "-r") {
            expecting = expect::num;
            cur_key = "radix";
        } else if(s == "-d") {
            expecting = expect::num;
            cur_key = "digit range";
        } else if(s == "-ndb") {
            expecting = expect::num;
            cur_key = "normalized divisor bound";
        } else if(s == "-s") {
            ret["sum"] = "true";
        } else if(s == "-m") {
            expecting = expect::num;
            cur_key = "max iterations";
        } else if(s == "-f") {
            expecting = expect::str;
            cur_key = "filename";
        } // printing
        else if(s == "-p") {
            expecting = expect::str;
            cur_key = "print options";
        } else if(s == "-pc") {
            expecting = expect::str;
            cur_key = "comment symbol";
        } else if(s == "-po") {
            expecting = expect::str;
            cur_key = "oor symbol";
        } else if(s == "-pr") {
            expecting = expect::num;
            cur_key = "comment radix";
        } else if(s == "-ps") {
            expecting = expect::str;
            cur_key = "separator";
        } else if(s == "-pi") {
            expecting = expect::str;
            cur_key = "inner separator";
        } else if(s == "-pb") {
            expecting = expect::str;
            cur_key = "begin multiple";
        } else if(s == "-pe") {
            expecting = expect::str;
            cur_key = "end multiple";
        } // hints
        else if(s == "-hp") {
            expecting = expect::num;
            cur_key = "p fractional";
        } else if(s == "-hpp") {
            expecting = expect::num;
            cur_key = "p fractional part";
        } else if(s == "-hd") {
            expecting = expect::num;
            cur_key = "d fractional";
        } else if(s == "-hdp") {
            expecting = expect::num;
            cur_key = "d fractional part";
        } else {
            std::cout << "err: unrecognized argument " << s << std::endl;
            ret["valid"] = "false";
        }
    };

    auto handle_num = [&](std::string s) {
        expecting = expect::arg;
        bool valid = true;
        if(s == "") valid = false;
        for(auto c: s)
            if(c < 48 || c > 57)
                valid = false;
        if(!valid) {
            std::cout << "err: expected an integer, '" << s << "' cannot be interpreted as an integer" <<
                std::endl;
            ret["valid"] = "false";
            return;
        }
        ret[cur_key] = s;
    };

    auto handle_str = [&](std::string s) {
        expecting = expect::arg;
        ret[cur_key] = s;
    };

    for(int i = 1; i < argc || expecting != expect::arg; ++i) {
        switch(expecting) {
            case expect::arg:
                handle_arg(std::string(argv[i]));
                break;
            case expect::num:
                handle_num(i < argc ? std::string(argv[i]) : "");
                break;
            case expect::str:
                handle_str(i < argc ? std::string(argv[i]) : "");
                break;
        }
    }

    return ret;
}

bool check_arguments(arguments_t &arguments) {
    bool ret = true;
    int radix = std::stoi(arguments["radix"]);
    int digit_range = std::stoi(arguments["digit range"]);
    int ndb = std::stoi(arguments["normalized divisor bound"]);
    int max_iterations = std::stoi(arguments["max iterations"]);
    int comment_radix = std::stoi(arguments["comment radix"]);
    // if the values are not set, just set them to pass the test
    int p_fractional_digits = arguments.find("p fractional") != arguments.end() ?
        std::stoi(arguments["p fractional"]) : 1;
    int p_fractional_digits_part = arguments.find("p fractional part") != arguments.end() ?
        std::stoi(arguments["p fractional part"]) : 0;
    int d_fractional_digits = arguments.find("d fractional") != arguments.end() ?
        std::stoi(arguments["d fractional"]) : 1;
    int d_fractional_digits_part = arguments.find("d fractional part") != arguments.end() ?
        std::stoi(arguments["d fractional part"]) : 0;
    std::string print_options = arguments["print options"];

    if(radix < 2 || radix > 36) {
        std::cout << "err: radix '" << radix << "' not supported, has to be in [2, 36]" << std::endl;
        ret = false;
    }
    if(digit_range > radix - 1 || digit_range < 1) {
        std::cout << "err: digit range '" << digit_range << "' has to be in range ]0, radix[" << std::endl;
        ret = false;
    }
    if(ndb < 2 || ndb > radix) {
        std::cout << "err: normalized divisor bound '" << ndb << "* has to be in range [2, radix]" << std::endl;
        return false;
    }
    if(max_iterations < 1) {
        std::cout << "err: max iterations '" << max_iterations << "' has to be positive" << std::endl;
        ret = false;
    }
    if(comment_radix < 2 || comment_radix > 36) {
        std::cout << "err: comment radix '" << comment_radix << "' not supported, has to be in [2, 36]" <<
            std::endl;
        ret = false;
    } else {
        if(!(((radix / comment_radix) * comment_radix) == radix)) {
            std::cout << "err: comment radix '" <<
                comment_radix << "'; radix has to be a multiple of it" << std::endl;
            ret = false;
        }
    }
    if(p_fractional_digits < 0) {
        std::cout << "err: p fractional '" << p_fractional_digits << "' cannot be negative" << std::endl;
        ret = false;
    }
    if(d_fractional_digits < 0) {
        std::cout << "err: d fractional '" << d_fractional_digits << "' cannot be negative" << std::endl;
        ret = false;
    }
    if(p_fractional_digits_part < 0) {
        std::cout << "err: p fractional part '" << d_fractional_digits_part <<
            "' cannot be negative" << std::endl;
        ret = false;
    }
    if(d_fractional_digits_part < 0) {
        std::cout << "err: d fractional part '" << d_fractional_digits_part <<
            "' cannot be negative" << std::endl;
        ret = false;
    }
    if(p_fractional_digits_part && !p_fractional_digits) {
        std::cout << "err: p fractional has to include the partial digit p fractional part '" <<
            p_fractional_digits_part << "'" << std::endl;
        ret = false;
    }
    if(d_fractional_digits_part && !d_fractional_digits) {
        std::cout << "err: d fractional has to include the partial digit d fractional part '" <<
            d_fractional_digits_part << "'" << std::endl;
        ret = false;
    }
    auto divs = dividers(radix);
    if(std::find(divs.begin(), divs.end(), p_fractional_digits_part) == divs.end()) {
        std::cout << "err: p fractional part '" << p_fractional_digits_part << "' has to be an integer "
            "divider of radix in range [2, radix[ or 0" << std::endl;
        ret = false;
    }
    if(std::find(divs.begin(), divs.end(), d_fractional_digits_part) == divs.end()) {
        std::cout << "err: d fractional part '" << d_fractional_digits_part << "' has to be an integer "
            "divider of radix in range [2, radix[ or 0" << std::endl;
        ret = false;
    }
    for(auto c: print_options) {
        if(!(c == 'c' || c == 'r' || c == 'h' || c == 'l' || c == 'a' || c == 'i' ||
                    c == 'p' || c == 't' || c == 'v' || c == 'u')) {
            std::cout << "err: print options '" << c << "' not supported" << std::endl;
            ret = false;
        }
    }

    // how many of v and vs are equal to val
    auto num_are_equal = [](auto val, auto v, auto ...vs) {
        bool equals[] = {val == v, (val == vs)...};
        return std::count(std::begin(equals), std::end(equals), true);
    };
    auto n = num_are_equal(print_options.npos,
            print_options.find('h'),
            print_options.find('l'),
            print_options.find('a'),
            print_options.find('v'),
            print_options.find('u'));
    if(n == 5) {
        // default to h
        arguments["print options"] += 'h';
    } else if(n != 4) {
        std::cout << "err: print options '" << print_options <<
            "' only one of {h,l,a,v,u} is supported" << std::endl;
        ret = false;
    }
    if(print_options.find('t') != print_options.npos && comment_radix != 2) {
        std::cout << "err: using two's complement for partial remainder representation is only supprted "
            "with a comment radix of 2" << std::endl;
        ret = false;
    }

    return ret;
}

void print_help(std::string executable_name, arguments_t default_arguments) {
    std::ostream &o = std::cout;
    unsigned gap = 36;
    unsigned max = 100;
    auto print_argument = [&](std::string start, std::string variable, std::string desc) {
        std::string out = "  " + start;
        if(variable.size())
            out += " <" + variable + ">";
        out += "  ";
        while(out.size() < gap)
            out += ' ';
        out += desc;
        auto search = default_arguments.find(variable);
        if(search != default_arguments.end()) {
            out += "\n(default: " + search->second + ")";
        }

        std::string gap_string = "";
        std::vector<std::string> out_lines;
        while(gap_string.size() < gap)
            gap_string += ' ';
        out += ' ';
        while(true) {
            decltype(out)::size_type split_pos = out.find_first_of("\n");
            if(split_pos > max) {
                split_pos = std::min(static_cast<decltype(out)::size_type>(max),
                        out.find_last_of(" ", max));
            }
            out_lines.push_back(out.substr(0, split_pos));
            out.erase(0, split_pos + 1);
            if(!out.size())
                break;
            out = gap_string + out;
        };
        for(auto l: out_lines)
            o << l << std::endl;
    };

    o << "usage: " << executable_name << " [arguments]" << std::endl;
    o << std::endl;
    print_argument("-h  or --help", "", "print this help, other arguments are discarded");
    print_argument("-r", "radix", "radix a.k.a. number base in [2, 36], "
            "support for higher radices might be easy to implement, feel free "
            "(hint: add more output symbols)");
    print_argument("-d", "digit range", "digit range in ]0, radix[, "
            "the quotient digits to be used are then [-digit range, digit range]\n(defaults to radix-1)");
    print_argument("-ndb", "normalized divisor bound", "the upper bound of the divisor. "
            "The usual restraint on the divisor is 1 <= divisor < r, but if the divisor is known to be "
            "smaller than another number ndb, the divisor range and therefore the table size can be decreased "
            "(defaults to radix)");
    print_argument("-s", "", "for scenarios where the partial remainder is represented as the sum of two "
            "values and is calculated after the approximation step, e.g. when carry save addition is used");
    print_argument("-m", "max iterations", "the program tests possible SRT tables in ascending "
            "order concerning their sizes, this is the maximum number of iterations");
    print_argument("-f", "filename", "output file, print to stdout if empty");
    std::cout << std::endl << "print configuration" << std::endl;
    print_argument("-p", "print options", "various print options:");
    print_argument("    c", "", "print column comments");
    print_argument("    r", "", "print row comments");
    print_argument("    h,l,a,v,u", "", "for cells with multiple possibilities only print "
            "the [h]ighest one, the [l]owest one, [a]ll, the highest absolute [v]alue "
            "or the lowest absolute val[u]e\n(defaults to h)");
    print_argument("    i", "", "print table info comments");
    print_argument("    p", "", "pad the partial remainder and if applicable the divisor in the comments, "
            "so every partial remainder and divisor has the same length");
    print_argument("    t", "", "use two's complement for partial remainder, this can only "
            "be used if comment radix is 2 and padding with 0s is activated automatically");
    print_argument("-pc", "comment symbol", "comment symbol to be used in the output");
    print_argument("-po", "oor symbol", "symbol to use for out-of-range cells");
    print_argument("-pr", "comment radix", "radix to use for the comments, has to be a multiple of radix "
            "\n(defaults to radix)");
    print_argument("-ps", "separator", "symbol to use for digit separation");
    print_argument("-pi", "inner separator", "symbol to use for separation of "
            "multiple digits in the same cell");
    print_argument("-pb", "begin multiple", "symbol to use for beginning multiple digits");
    print_argument("-pe", "end multiple", "symbol to use for ending multiple digits");
    std::cout << std::endl << "computational hints" << std::endl;
    print_argument("-hp", "p fractional", "number of fractional digits of the partial remainder including "
            "partial digits "
            "e.g. with radix 2 and -hp 2 the possible fractional digits are 00, 01, 10 and 11");
    print_argument("-hpp", "p fractional part", "number of digits that will be skipped in the last part "
            "of the partial remainder e.g. with radix 4 and -hpp 2 possible last digits are 0 and 2; "
            "has to be in range [2, radix[ or 0 for no partial digit");
    print_argument("-hd", "d fractional", "analogue to p fractional, but for the divisor");
    print_argument("-hdp", "d fractional part", "analogue to p fractional part, but for the divisor");
}

srt_print_info_t make_print_info(arguments_t args) {
    srt_print_info_t ret;

    ret.comment_radix = std::stoi(args["comment radix"]);
    ret.comment_seq = args["comment symbol"];
    auto po = args["print options"];
    if(po.find('h') != po.npos) ret.mp = srt_print_info_t::multiple_possibilities::h;
    else if(po.find('l') != po.npos) ret.mp = srt_print_info_t::multiple_possibilities::l;
    else if(po.find('v') != po.npos) ret.mp = srt_print_info_t::multiple_possibilities::v;
    else if(po.find('u') != po.npos) ret.mp = srt_print_info_t::multiple_possibilities::u;
    else ret.mp = srt_print_info_t::multiple_possibilities::a;
    ret.oor_seq = args["oor symbol"];
    ret.print_column_comments = po.find('c') != po.npos;
    ret.print_row_comments = po.find('r') != po.npos;
    ret.print_info_comments = po.find('i') != po.npos;
    ret.separator = args["separator"];
    ret.inner_separator = args["inner separator"];
    ret.begin_multiple = args["begin multiple"];
    ret.end_multiple = args["end multiple"];

    return ret;
}

srt_info_t make_srt_info(arguments_t args) {
    srt_info_t ret;
    ret.digit_range = std::stoi(args["digit range"]);
    ret.radix = std::stoi(args["radix"]);
    ret.normalized_divisor_bound = std::stoi(args["normalized divisor bound"]);
    ret.sum = args["sum"] == "true";

    // first create minimum info
    ret.p_fractional_digits = args.find("p fractional part") != args.end() ? 1 : 0;
    ret.p_fractional_digits_part = 0;
    ret.d_fractional_digits = args.find("d fractional part") != args.end() ? 1 : 0;
    ret.d_fractional_digits_part = 0;

    // fill with args information
    if(args.find("p fractional") != args.end()) {
        ret.p_fractional_digits = std::stoi(args["p fractional"]);
        if(ret.p_fractional_digits) ret.p_fractional_digits_part = 0;
    }
    if(args.find("p fractional part") != args.end()) {
        ret.p_fractional_digits_part = std::stoi(args["p fractional part"]);
        ret.p_fractional_digits -= 1; // don't count partial digit here
    }
    if(args.find("d fractional") != args.end()) {
        ret.d_fractional_digits = std::stoi(args["d fractional"]);
        if(ret.d_fractional_digits) ret.d_fractional_digits_part = 0;
    }
    if(args.find("d fractional part") != args.end()) {
        ret.d_fractional_digits_part = std::stoi(args["d fractional part"]);
        ret.d_fractional_digits -= 1; // don't count partial digit here
    }

    return ret;
}

srt_radix_representation_info_t make_radix_representation_info(arguments_t args) {
    srt_radix_representation_info_t ret;
    auto po = args["print options"];
    ret.pad = po.find('p') != po.npos;
    ret.complement = po.find('t') != po.npos;
    // for two's complement we activate padding with 0s
    if(ret.complement) ret.pad = true;
    ret.radix = std::stoi(args["comment radix"]);
    // these have to be computed when the partial remainders are known
    ret.total_length = 0;
    ret.digits_below_1 = 0;
    return ret;
}

int main(int argc, char **argv) {
    std::ostream *out = &std::cout;
    std::ofstream fs;

    arguments_t default_arguments = {
        {"radix", "4"},
        {"comment symbol", "//"},
        {"max iterations", "10000"},
        {"oor symbol", "0"},
        {"separator", ","},
        {"inner separator", ","},
        {"begin multiple", "["},
        {"end multiple", "]"}
    };

    arguments_t arguments = parse_arguments(argc, argv, default_arguments);
    if(arguments["valid"] == "false")
        return EXIT_SUCCESS;
    if(arguments["help"] == "true") {
        print_help(std::string(argv[0]), default_arguments);
        return EXIT_SUCCESS;
    }

    // additional defaults
    if(arguments.find("digit range") == arguments.end())
        arguments["digit range"] = std::to_string(std::stoi(arguments["radix"]) - 1);
    if(arguments.find("normalized divisor bound") == arguments.end())
        arguments["normalized divisor bound"] = arguments["radix"];
    if(arguments.find("comment radix") == arguments.end())
        arguments["comment radix"] = arguments["radix"];

    if(!check_arguments(arguments)) {
        return EXIT_SUCCESS;
    }

    // print to file if filename is given
    if(arguments.find("filename") != arguments.end()) {
        try {
            fs = std::ofstream(arguments["filename"], std::ofstream::out | std::ofstream::trunc);
            out = &fs;
        } catch(...) {
            std::cout << "err: could not open file '" << arguments["filename"] << std::endl;
            return EXIT_SUCCESS;
        }
    }

    // finally compute something
    auto result = smallest_srt_table(
            make_srt_info(arguments),
            std::stoi(arguments["max iterations"]));

    if(!result.second.empty())
        print_srt_table(*out, result.second, result.first,
                make_print_info(arguments), make_radix_representation_info(arguments));
    else
        std::cout << "failed to generate SRT table, either no valid srt table exists for the inputs "
            "or it was not found because it's to large, you can try again with an increased number of "
            "iterations (and maybe some hints)" << std::endl;
}
