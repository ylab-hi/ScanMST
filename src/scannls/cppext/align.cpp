#include <iostream>
#include <vector>
#include <algorithm>
#include <pybind11/pybind11.h>



#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

constexpr int gap{-3};
constexpr int gap_ext{-1};
constexpr int default_value{-9999};
constexpr const double EDNAFULL_matrix[15][15] = {{5,  -4, -4, -4, -4, 1,  1,  -4, -4, 1,  -4, -1, -1, -1, -2},
                                        {-4, 5,  -4, -4, -4, 1,  -4, 1,  1,  -4, -1, -4, -1, -1, -2},
                                        {-4, -4, 5,  -4, 1,  -4, 1,  -4, 1,  -4, -1, -1, -4, -1, -2},
                                        {-4, -4, -4, 5,  1,  -4, -4, 1,  -4, 1,  -1, -1, -1, -4, -2},
                                        {-4, -4, 1,  1,  -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1},
                                        {1,  1,  -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1},
                                        {1,  -4, 1,  -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1},
                                        {-4, 1,  -4, 1,  -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1},
                                        {-4, 1,  1,  -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1},
                                        {1,  -4, -4, 1,  -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1},
                                        {-4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1},
                                        {-1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1},
                                        {-1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1},
                                        {-1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1},
                                        {-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

int getCharIndex(char Char) {
    switch (Char) {
        case 'A':
            return 0;
        case 'T':
            return 1;
        case 'G':
            return 2;
        case 'C':
            return 3;
        case 'S':
            return 4;
        case 'W':
            return 5;
        case 'R':
            return 6;
        case 'Y':
            return 7;
        case 'K':
            return 8;
        case 'M':
            return 9;
        case 'B':
            return 10;
        case 'V':
            return 11;
        case 'H':
            return 12;
        case 'D':
            return 13;
        default:
            std::cout << "Invalid character: " << Char << "\n";
            return -1;
    }
}

int getScore(char a, char b) {
    int i = getCharIndex(a);
    int j = getCharIndex(b);
    return EDNAFULL_matrix[i][j];
}

int getGaps(int num_gap_ext) {
    return num_gap_ext * gap_ext + gap;
}

enum class Direction_t {
    left,
    up,
    diag_match,
    diag_mismatch,
    default_,
};

struct AlignmentResult {
    std::string sequence_x{};
    std::string sequence_y{};
    int mismatch{0};
    int sequence_x_gap{0};
    int sequence_y_gap{0};
    int end_x{0};
    int end_y{0};
};

struct MaxScore {
    int score{0};
    int x_index{0};
    int y_index{0};
};

using score_matrix_t = std::vector<std::vector<int>>;
using alignment_t = std::vector<std::vector<Direction_t>>;

AlignmentResult tracePath(const score_matrix_t &scoreMatrix, const alignment_t &alignment, int x_len, int y_len,
                          const std::string &sequence_x,
                          const std::string &sequence_y) {
    AlignmentResult result{};
    result.end_x = x_len;
    result.end_y = y_len;

    while (scoreMatrix[x_len][y_len] != 0) {
        switch (alignment[x_len][y_len]) {
            case Direction_t::left:
                ++result.sequence_x_gap;
                result.sequence_x += '-';
                result.sequence_y += sequence_y[y_len];
                --y_len;
                break;
            case Direction_t::up:
                ++result.sequence_y_gap;
                result.sequence_x += sequence_x[x_len];
                result.sequence_y += '-';
                --x_len;
                break;
            case Direction_t::diag_match:
                result.sequence_x += sequence_x[x_len];
                result.sequence_y += sequence_y[y_len];
                --x_len;
                --y_len;
                break;
            case Direction_t::diag_mismatch:
                ++result.mismatch;
                result.sequence_x += sequence_x[x_len];
                result.sequence_y += sequence_y[y_len];
                --x_len;
                --y_len;
                break;
            case Direction_t::default_:
                std::cout << "Error: tracePath()\n";
                break;
            default:
                std::cout << "Error: invalid direction" << "\n";
                break;
        }
    }
    return result;

}

int recursiveFind(int xlen, int ylen, int &gap_number, std::string &seqx, std::string &seqy, alignment_t &Alignment,
                  score_matrix_t &scoreMatrix, MaxScore &m_score) {
    if (xlen == 0 && ylen == 0) {
        return 0;
    }
    if (scoreMatrix[xlen][ylen] != default_value) {
        return scoreMatrix[xlen][ylen];
    }

    bool is_match = seqx[xlen] == seqy[ylen];

    std::pair<int, Direction_t> left{
            recursiveFind(xlen, ylen - 1, gap_number, seqx, seqy, Alignment, scoreMatrix, m_score) + getGaps(gap_number),
            Direction_t::left};
    std::pair<int, Direction_t> up{
            recursiveFind(xlen - 1, ylen, gap_number, seqx, seqy, Alignment, scoreMatrix, m_score) + getGaps(gap_number),
            Direction_t::up};
    std::pair<int, Direction_t> diag{
            recursiveFind(xlen - 1, ylen - 1, gap_number, seqx, seqy, Alignment, scoreMatrix, m_score) +
            getScore(seqx[xlen], seqy[ylen]),
            (is_match ? Direction_t::diag_match : Direction_t::diag_mismatch)};

    std::vector<std::pair<int, Direction_t>> directions{left, up, diag,
                                                        std::pair<int, Direction_t>{0, Direction_t::default_}};
    std::sort(directions.begin(), directions.end(),
              [](const std::pair<int, Direction_t> &a, const std::pair<int, Direction_t> &b) {
                  return a.first > b.first;
              });
    std::pair<int, Direction_t> best{directions[0]};
    scoreMatrix[xlen][ylen] = best.first;
    Alignment[xlen][ylen] = best.second;

    if (best.first > m_score.score) {
        m_score.score = best.first;
        m_score.x_index = xlen;
        m_score.y_index = ylen;
    }

    if (best.second == Direction_t::left || best.second == Direction_t::up)
        ++gap_number;
    else gap_number = 0;

    return best.first;

}

AlignmentResult align(std::string &x, std::string &y) {

    auto x_len = static_cast<int>(x.length()) -1;
    auto y_len = static_cast<int>(y.length()) -1;

    score_matrix_t score_matrix{static_cast<std::size_t>(x_len) + 1,
                                std::vector<int>(static_cast<std::size_t>(y_len) + 1, -9999)};

    alignment_t alignment{static_cast<std::size_t>(x_len) + 1,
                          std::vector<Direction_t>(static_cast<std::size_t>(y_len) + 1, Direction_t::default_)};

    for (int i = 0; i < x_len + 1; ++i) {
        score_matrix[i][0] = 0;
    }

    for (int j = 0; j < y_len + 1; ++j) {
        score_matrix[0][j] = 0;
    }

    MaxScore max_score{0, x_len, y_len};

    int gap_ex{0};
    int final_score = recursiveFind(x_len, y_len, gap_ex, x, y, alignment,
                                    score_matrix, max_score);

    AlignmentResult result{tracePath(score_matrix, alignment, max_score.x_index, max_score.y_index, x, y)};
    std::reverse(result.sequence_x.begin(), result.sequence_x.end());
    std::reverse(result.sequence_y.begin(), result.sequence_y.end());
    return result;
}

PYBIND11_MODULE(align, m) {
    m.doc() = "aligner cpp plugin"; // optional module docstring
    m.def("align", &align, "A function which aligns two strings");
    m.def("align_nogil", &align, py::call_guard<py::gil_scoped_release>(), "A function which aligns two strings with no gil");
    py::class_<AlignmentResult>(m, "AlignmentResult")
             .def(py::init<>())
             .def_readwrite("mismatch", &AlignmentResult::mismatch)
             .def_readwrite("sequence_x", &AlignmentResult::sequence_x)
             .def_readwrite("sequence_y", &AlignmentResult::sequence_y)
             .def_readwrite("sequence_x_gap", &AlignmentResult::sequence_x_gap)
             .def_readwrite("sequence_y_gap", &AlignmentResult::sequence_y_gap)
             .def_readwrite("end_x", &AlignmentResult::end_x)
             .def_readwrite("end_y", &AlignmentResult::end_y);

}
