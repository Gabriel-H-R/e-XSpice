/**
 * @file csv2svg.cpp
 * @brief Dependency-free CSV -> SVG line/bar plotter (showcase figures).
 *
 * Reads a simulator CSV (first line = header, first column = x) and writes a
 * standalone SVG that can be embedded directly in the README:
 *
 *   csv2svg results.csv figure.svg --title "AC sweep" --xlog \
 *           --series 1,3 --xlabel "Frequency (Hz)" --ylabel "|H| (V/V)"
 *
 * Options:
 *   --title <text>       figure title
 *   --xlabel <text>      x-axis label
 *   --ylabel <text>      y-axis label
 *   --series <c1,c2,...> data columns to plot (default: all except column 0)
 *   --xlog               logarithmic (base 10) x axis
 *   --type line|bar      plot kind (default: line)
 *   --ymin <v> --ymax <v> explicit y range
 *
 * Only the C++ standard library is used, so it builds everywhere the simulator
 * does (no Python, no gnuplot, no matplotlib).
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int kWidth = 960;
constexpr int kHeight = 520;
constexpr int kMarginLeft = 96;
constexpr int kMarginRight = 40;
constexpr int kMarginTop = 78;
constexpr int kMarginBottom = 72;

const std::vector<std::string> kPalette = {
    "#0b6e99", "#d1495b", "#2a9d8f", "#e9c46a",
    "#6a4c93", "#1982c4", "#e76f51", "#43aa8b"
};

struct Options {
    std::string input;
    std::string output;
    std::string title;
    std::string xlabel;
    std::string ylabel;
    std::string type = "line";
    bool xlog = false;
    bool has_ymin = false;
    bool has_ymax = false;
    double ymin = 0.0;
    double ymax = 1.0;
    std::vector<int> columns;   // empty = every data column
};

std::string xmlEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '&':  out += "&amp;";  break;
            case '<':  out += "&lt;";   break;
            case '>':  out += "&gt;";   break;
            case '"':  out += "&quot;"; break;
            default:   out += c;        break;
        }
    }
    return out;
}

std::string fmt(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%g", v);
    return buf;
}

bool isNumber(const std::string& s) {
    if (s.empty()) return false;
    char* end = nullptr;
    std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
}

std::vector<std::string> splitLine(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) {
        if (!field.empty() && field.back() == '\r') field.pop_back();
        // Tolerate RFC-4180 quoting (PowerShell's Export-Csv quotes all
        // fields): strip the surrounding quotes and unescape "" -> ".
        if (field.size() >= 2 && field.front() == '"' && field.back() == '"') {
            field = field.substr(1, field.size() - 2);
            std::string unq;
            unq.reserve(field.size());
            for (size_t i = 0; i < field.size(); ++i) {
                if (field[i] == '"' && i + 1 < field.size() && field[i + 1] == '"') {
                    unq += '"';
                    ++i;
                } else {
                    unq += field[i];
                }
            }
            field = unq;
        }
        fields.push_back(field);
    }
    return fields;
}

std::vector<double> parseRow(const std::vector<std::string>& fields) {
    std::vector<double> row;
    row.reserve(fields.size());
    for (const auto& f : fields) row.push_back(std::strtod(f.c_str(), nullptr));
    return row;
}

/// "Nice" tick positions between lo and hi (roughly count of them).
std::vector<double> niceTicks(double lo, double hi, int count) {
    std::vector<double> ticks;
    if (!(hi > lo)) return ticks;
    const double raw = (hi - lo) / std::max(1, count);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    double step = mag;
    for (double m : {1.0, 2.0, 5.0, 10.0}) {
        if (raw <= m * mag) { step = m * mag; break; }
    }
    const double start = std::ceil(lo / step) * step;
    for (double t = start; t <= hi + step * 1e-9; t += step) ticks.push_back(t);
    return ticks;
}

void svgText(std::ostream& os, double x, double y, const std::string& text,
             const std::string& anchor = "middle", const std::string& size = "12",
             const std::string& weight = "normal") {
    os << "<text x=\"" << x << "\" y=\"" << y << "\" text-anchor=\"" << anchor
       << "\" font-family=\"Segoe UI, Helvetica, Arial, sans-serif\" font-size=\""
       << size << "\" font-weight=\"" << weight << "\" fill=\"#333\">"
       << xmlEscape(text) << "</text>\n";
}

void drawFrame(std::ostream& os, const Options& opt,
               double x0, double y0, double x1, double y1,
               const std::vector<double>& xticks,
               const std::vector<double>& yticks,
               const std::vector<std::string>& xlabels,
               const std::vector<std::string>& ylabels,
               bool xlog) {
    (void)opt;
    // Grid.
    for (size_t i = 0; i < xticks.size(); ++i) {
        const double gx = x0 + (x1 - x0) * xticks[i];
        os << "<line x1=\"" << gx << "\" y1=\"" << y0 << "\" x2=\"" << gx
           << "\" y2=\"" << y1 << "\" stroke=\"#e8e8e8\" stroke-width=\"1\"/>\n";
    }
    for (size_t i = 0; i < yticks.size(); ++i) {
        // yticks run from ymin (0) to ymax (1) and y0 is the TOP of the
        // plot area: larger values must move UP, towards y0.
        const double gy = y1 - (y1 - y0) * yticks[i];
        os << "<line x1=\"" << x0 << "\" y1=\"" << gy << "\" x2=\"" << x1
           << "\" y2=\"" << gy << "\" stroke=\"#e8e8e8\" stroke-width=\"1\"/>\n";
    }
    // Axes.
    os << "<rect x=\"" << x0 << "\" y=\"" << y0 << "\" width=\"" << (x1 - x0)
       << "\" height=\"" << (y1 - y0) << "\" fill=\"none\" stroke=\"#333\" "
       << "stroke-width=\"1\"/>\n";
    // Tick labels (positions already normalised to 0..1).
    for (size_t i = 0; i < xticks.size() && i < xlabels.size(); ++i) {
        svgText(os, x0 + (x1 - x0) * xticks[i], y1 + 20, xlabels[i]);
    }
    for (size_t i = 0; i < yticks.size() && i < ylabels.size(); ++i) {
        svgText(os, x0 - 10, y1 - (y1 - y0) * yticks[i] + 4, ylabels[i], "end");
    }
    (void)xlog;
}

/// Maps data to normalised plot coordinates (x in 0..1, y in 0..1).
double normX(double v, double lo, double hi, bool xlog) {
    if (xlog) return (std::log10(std::max(v, 1e-300)) - std::log10(lo))
                      / (std::log10(hi) - std::log10(lo));
    return (v - lo) / (hi - lo);
}

double normY(double v, double lo, double hi) { return (v - lo) / (hi - lo); }

void plotLine(std::ostream& os, const Options& opt,
              const std::vector<std::string>& header,
              const std::vector<std::vector<double>>& data) {
    // Selected columns.
    std::vector<int> cols = opt.columns;
    if (cols.empty()) {
        for (size_t c = 1; c < header.size(); ++c) cols.push_back(static_cast<int>(c));
    }

    // Ranges.
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -xmin;
    double ymin = xmin, ymax = xmax;
    for (const auto& row : data) {
        if (row.empty()) continue;
        xmin = std::min(xmin, row[0]);
        xmax = std::max(xmax, row[0]);
        for (int c : cols) {
            if (c >= 0 && c < static_cast<int>(row.size())) {
                ymin = std::min(ymin, row[c]);
                ymax = std::max(ymax, row[c]);
            }
        }
    }
    if (opt.has_ymin) ymin = opt.ymin;
    if (opt.has_ymax) ymax = opt.ymax;
    if (!(ymax > ymin)) { ymax = ymin + 1.0; }
    if (!(xmax > xmin)) { xmax = xmin + 1.0; }
    if (opt.xlog) xmin = std::max(xmin, 1e-300);

    const double px0 = kMarginLeft, py0 = kMarginTop;
    const double px1 = kWidth - kMarginRight, py1 = kHeight - kMarginBottom;

    // Ticks in transformed space.
    std::vector<double> xticks, yticks;
    std::vector<std::string> xlabels, ylabels;
    if (opt.xlog) {
        const double l0 = std::log10(xmin), l1 = std::log10(xmax);
        for (double t = std::ceil(l0); t <= std::floor(l1); ++t) {
            xticks.push_back((t - l0) / (l1 - l0));
            xlabels.push_back(fmt(std::pow(10.0, t)));
        }
        if (xticks.size() < 2) {  // narrow band: decade subdivisions
            const double step = (l1 - l0) / 8.0;
            for (int i = 0; i <= 8; ++i) {
                xticks.push_back(static_cast<double>(i) / 8.0);
                xlabels.push_back(fmt(std::pow(10.0, l0 + i * step)));
            }
        }
    } else {
        const auto t = niceTicks(xmin, xmax, 6);
        for (double v : t) { xticks.push_back(normX(v, xmin, xmax, false)); xlabels.push_back(fmt(v)); }
    }
    {
        const auto t = niceTicks(ymin, ymax, 6);
        for (double v : t) { yticks.push_back(normY(v, ymin, ymax)); ylabels.push_back(fmt(v)); }
    }

    drawFrame(os, opt, px0, py0, px1, py1, xticks, yticks, xlabels, ylabels, opt.xlog);

    // Series.
    for (size_t s = 0; s < cols.size(); ++s) {
        const int c = cols[s];
        if (c >= static_cast<int>(header.size())) continue;
        const std::string& name = header[c];
        os << "<polyline fill=\"none\" stroke=\"" << kPalette[s % kPalette.size()]
           << "\" stroke-width=\"2.2\" points=\"";
        bool first = true;
        for (const auto& row : data) {
            if (c >= static_cast<int>(row.size()) || row.empty()) continue;
            double x = opt.xlog ? std::log10(std::max(row[0], 1e-300)) : row[0];
            const double nx = opt.xlog
                ? (x - std::log10(xmin)) / (std::log10(xmax) - std::log10(xmin))
                : (x - xmin) / (xmax - xmin);
            const double ny = (row[c] - ymin) / (ymax - ymin);
            const double sx = px0 + nx * (px1 - px0);
            const double sy = py1 - ny * (py1 - py0);
            if (!first) os << " ";
            os << sx << "," << sy;
            first = false;
        }
        os << "\"/>\n";

        // Legend entry (above the plot, stacked to the right).
        const double lx = px1 - 160;
        const double ly = py0 - 14 + static_cast<double>(s) * 20;
        os << "<line x1=\"" << lx << "\" y1=\"" << ly << "\" x2=\"" << (lx + 26)
           << "\" y2=\"" << ly << "\" stroke=\"" << kPalette[s % kPalette.size()]
           << "\" stroke-width=\"3\"/>\n";
        svgText(os, lx + 32, ly + 4, name, "start", "12");
    }
}

void plotBar(std::ostream& os, const Options& opt,
             const std::vector<std::string>& header,
             const std::vector<std::vector<double>>& data,
             const std::vector<std::string>& firstLabels,
             bool firstColNumeric) {
    std::vector<int> cols = opt.columns;
    if (cols.empty()) {
        for (size_t c = 1; c < header.size(); ++c) cols.push_back(static_cast<int>(c));
    }
    if (data.empty()) return;

    // Categories = unique x values (numeric) or unique labels (text), in order
    // of appearance. Text categories arrive as 0 in the parsed doubles, so
        // they must be tracked separately.
    std::vector<double> cats;
    std::vector<std::string> cat_names;
    for (size_t r = 0; r < data.size(); ++r) {
        if (data[r].empty()) continue;
        if (firstColNumeric) {
            if (cats.empty() || cats.back() != data[r][0]) cats.push_back(data[r][0]);
        } else {
            const std::string& lbl = (r < firstLabels.size()) ? firstLabels[r]
                                                              : std::string();
            if (cat_names.empty() || cat_names.back() != lbl) cat_names.push_back(lbl);
        }
    }
    const size_t ncat = firstColNumeric ? cats.size() : cat_names.size();

    double ymax = 0.0, ymin = 0.0;
    for (const auto& row : data) {
        for (int c : cols) {
            if (c >= 0 && c < static_cast<int>(row.size())) {
                ymax = std::max(ymax, row[c]);
                ymin = std::min(ymin, row[c]);
            }
        }
    }
    if (opt.has_ymax) ymax = opt.ymax;
    if (opt.has_ymin) ymin = opt.ymin;
    if (!(ymax > ymin)) ymax = ymin + 1.0;

    // Leave headroom at the extremes that carry data, so the value labels
    // above/below the bar tips stay inside the figure.
    {
        const double range = ymax - ymin;
        if (ymin < 0.0) ymin -= 0.12 * range;
        if (ymax > 0.0) ymax += 0.12 * range;
        if (!(ymax > ymin)) ymax = ymin + 1.0;
    }

    const double px0 = kMarginLeft, py0 = kMarginTop;
    const double px1 = kWidth - kMarginRight, py1 = kHeight - kMarginBottom;

    std::vector<double> yticks;
    std::vector<std::string> ylabels;
    for (double v : niceTicks(ymin, ymax, 6)) {
        yticks.push_back(normY(v, ymin, ymax));
        ylabels.push_back(fmt(v));
    }
    std::vector<double> xticks(ncat, 0.0);
    std::vector<std::string> xlabels;
    for (size_t i = 0; i < ncat; ++i) {
        xticks[i] = (static_cast<double>(i) + 0.5) / static_cast<double>(ncat);
        xlabels.push_back(firstColNumeric ? fmt(cats[i]) : cat_names[i]);
    }
    drawFrame(os, opt, px0, py0, px1, py1, xticks, yticks, xlabels, ylabels, false);

    const double slot = (px1 - px0) / static_cast<double>(ncat);
    const double group = slot * 0.7;
    const double bw = group / static_cast<double>(cols.size());

    // Bars are anchored at ZERO (not at ymin): sensitivities and spectra can
    // be negative, and a bar chart that grows from ymin hides the sign.
    const double y_zero = py1 - normY(0.0, ymin, ymax) * (py1 - py0);
    if (ymin < 0.0 && ymax > 0.0) {
        os << "<line x1=\"" << px0 << "\" y1=\"" << y_zero << "\" x2=\"" << px1
           << "\" y2=\"" << y_zero << "\" stroke=\"#999\" stroke-width=\"1\"/>\n";
    }

    for (size_t i = 0; i < ncat; ++i) {
        // Find the row of this category (first match is enough: CSVs are
        // ordered by the sweep variable).
        const std::vector<double>* row = nullptr;
        for (size_t r = 0; r < data.size(); ++r) {
            if (data[r].empty()) continue;
            const bool match = firstColNumeric
                ? (data[r][0] == cats[i])
                : (r < firstLabels.size() && firstLabels[r] == cat_names[i]);
            if (match) { row = &data[r]; break; }
        }
        if (!row) continue;
        for (size_t s = 0; s < cols.size(); ++s) {
            const int c = cols[s];
            if (c >= static_cast<int>(row->size())) continue;
            const double v = (*row)[c];
            const double ny = (v - ymin) / (ymax - ymin);
            const double bx = px0 + slot * static_cast<double>(i) + (slot - group) / 2.0
                              + bw * static_cast<double>(s);
            const double vy = py1 - ny * (py1 - py0);
            const double by = std::min(y_zero, vy);
            const double bh = std::abs(vy - y_zero);
            os << "<rect x=\"" << bx << "\" y=\"" << by << "\" width=\"" << bw
               << "\" height=\"" << bh << "\" fill=\"" << kPalette[s % kPalette.size()]
               << "\" stroke=\"white\" stroke-width=\"1\"/>\n";
            if (ncat <= 8) {
                const double ly = (v >= 0.0) ? by - 6 : by + bh + 12;
                svgText(os, bx + bw / 2.0, ly, fmt(v), "middle", "11");
            }
        }
    }

    // Legend.
    for (size_t s = 0; s < cols.size(); ++s) {
        const int c = cols[s];
        if (c >= static_cast<int>(header.size())) continue;
        const double lx = px1 - 160;
        const double ly = py0 - 14 + static_cast<double>(s) * 20;
        os << "<rect x=\"" << lx << "\" y=\"" << (ly - 9) << "\" width=\"14\" height=\"10\" fill=\""
           << kPalette[s % kPalette.size()] << "\"/>\n";
        svgText(os, lx + 22, ly, header[c], "start", "12");
    }
}

} // namespace

int main(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](const char* what) -> std::string {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "csv2svg: missing value for %s\n", what);
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--title")       opt.title = next("--title");
        else if (a == "--xlabel") opt.xlabel = next("--xlabel");
        else if (a == "--ylabel") opt.ylabel = next("--ylabel");
        else if (a == "--type")   opt.type = next("--type");
        else if (a == "--xlog")   opt.xlog = true;
        else if (a == "--ymin")   { opt.ymin = std::atof(next("--ymin").c_str()); opt.has_ymin = true; }
        else if (a == "--ymax")   { opt.ymax = std::atof(next("--ymax").c_str()); opt.has_ymax = true; }
        else if (a == "--series") {
            for (const auto& part : splitLine(next("--series"))) {
                opt.columns.push_back(std::atoi(part.c_str()));
            }
        } else if (a.rfind("--", 0) == 0) {
            std::fprintf(stderr, "csv2svg: unknown option '%s'\n", a.c_str());
            return 2;
        } else if (opt.input.empty()) {
            opt.input = a;
        } else if (opt.output.empty()) {
            opt.output = a;
        }
    }

    if (opt.input.empty() || opt.output.empty()) {
        std::fprintf(stderr,
            "Usage: csv2svg <input.csv> <output.svg> [--title T] [--xlog] "
            "[--series c1,c2] [--xlabel X] [--ylabel Y] [--type line|bar]\n");
        return 2;
    }

    std::ifstream in(opt.input);
    if (!in.is_open()) {
        std::fprintf(stderr, "csv2svg: cannot open '%s'\n", opt.input.c_str());
        return 1;
    }

    std::string line;
    if (!std::getline(in, line)) {
        std::fprintf(stderr, "csv2svg: empty CSV '%s'\n", opt.input.c_str());
        return 1;
    }
    const std::vector<std::string> header = splitLine(line);

    std::vector<std::vector<double>> data;
    std::vector<std::string> firstLabels;
    bool firstColNumeric = true;
    bool firstColChecked = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto fields = splitLine(line);
        if (fields.empty()) continue;
        if (!firstColChecked) {
            firstColNumeric = isNumber(fields[0]);
            firstColChecked = true;
        }
        bool restNumeric = true;
        for (size_t i = 1; i < fields.size(); ++i) {
            if (!isNumber(fields[i])) { restNumeric = false; break; }
        }
        // A bar chart may carry a TEXT category in column 0; a line plot
        // needs numeric x values.
        if (restNumeric && (firstColNumeric || opt.type == "bar")) {
            data.push_back(parseRow(fields));
            firstLabels.push_back(fields[0]);
        }
    }
    if (data.empty()) {
        std::fprintf(stderr, "csv2svg: no numeric data rows in '%s'\n", opt.input.c_str());
        return 1;
    }
    if (!firstColNumeric && opt.type != "bar") {
        std::fprintf(stderr,
                     "csv2svg: '%s' has a text first column; use --type bar\n",
                     opt.input.c_str());
        return 1;
    }

    std::ofstream out(opt.output);
    if (!out.is_open()) {
        std::fprintf(stderr, "csv2svg: cannot write '%s'\n", opt.output.c_str());
        return 1;
    }

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << kWidth
        << "\" height=\"" << kHeight << "\" viewBox=\"0 0 " << kWidth << " "
        << kHeight << "\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    if (!opt.title.empty()) {
        svgText(out, kWidth / 2.0, 34, opt.title, "middle", "20", "bold");
    }
    if (opt.type == "bar") {
        plotBar(out, opt, header, data, firstLabels, firstColNumeric);
    } else {
        plotLine(out, opt, header, data);
    }
    if (!opt.xlabel.empty()) {
        svgText(out, (kMarginLeft + (kWidth - kMarginRight)) / 2.0,
                kHeight - 16, opt.xlabel);
    }
    if (!opt.ylabel.empty()) {
        out << "<text x=\"22\" y=\""
            << ((kMarginTop + (kHeight - kMarginBottom)) / 2.0)
            << "\" text-anchor=\"middle\" font-family=\"Segoe UI, Helvetica, Arial, "
               "sans-serif\" font-size=\"14\" fill=\"#333\" transform=\"rotate(-90 22 "
            << ((kMarginTop + (kHeight - kMarginBottom)) / 2.0) << ")\">"
            << xmlEscape(opt.ylabel) << "</text>\n";
    }
    out << "</svg>\n";
    return 0;
}
