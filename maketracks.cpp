/* MIT License
 *
 * Copyright (c) 2025 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

static constexpr auto about = R"(
make bigWig files by merging a set of methylomes
)";

static constexpr auto description = R"(
make bigWig files by merging a set of methylomes
)";

static constexpr auto examples = R"(
Examples:

maketracks -o outdir -n merged_name -d my_methylomes -m methylome_names.txt -g hg38 -x my_indexes
)";

#include <cli_common.hpp>
#include <format_error_code.hpp>  // IWYU pragma: keep
#include <genome_index.hpp>
#include <logger.hpp>
#include <methylome.hpp>
#include <utilities.hpp>

#include <CLI11/CLI11.hpp>

#include <chrono>
#include <cstdint>
#include <cstdlib>  // for EXIT_FAILURE, EXIT_SUCCESS
#include <format>
#include <iterator>  // for std::size
#include <map>
#include <memory>
#include <print>
#include <ranges>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <variant>  // IWYU pragma: keep
#include <vector>

namespace xfr = transferase;

typedef std::uint32_t big_mcount_t;

struct big_mcount_pair {
  big_mcount_t n_meth{};
  big_mcount_t n_unmeth{};
  big_mcount_pair() = default;
  big_mcount_pair(const std::integral auto n_meth,
                  const std::integral auto n_unmeth) :
    n_meth{static_cast<big_mcount_t>(n_meth)},
    n_unmeth{static_cast<big_mcount_t>(n_unmeth)} {}
  // ADS: need spaceship here because of constness
  [[nodiscard]] auto
  operator<=>(const big_mcount_pair &) const = default;

  /// Number of observations contributing to either state.
  [[nodiscard]] constexpr auto
  n_reads() const noexcept -> big_mcount_t {
    return n_meth + n_unmeth;
  }

  /// Get the methylation level
  [[nodiscard]] constexpr auto
  get_level() const noexcept -> float {
    return static_cast<float>(n_meth) /
           std::max(static_cast<big_mcount_t>(1), n_reads());
  }
};

struct big_methylome {
  std::vector<big_mcount_pair> cpgs{};

  big_methylome() = default;
  explicit big_methylome(const std::integral auto n_cpgs) :
    cpgs(std::vector<big_mcount_pair>(n_cpgs)) {}

  auto
  add(const xfr::methylome &rhs) -> void {
    const auto &data = rhs.data;
    std::ranges::transform(
      cpgs, data.cpgs, std::begin(cpgs),
      [](const auto &l, const auto &r) -> big_mcount_pair {
        return {static_cast<big_mcount_t>(l.n_meth + r.n_meth),
                static_cast<big_mcount_t>(l.n_unmeth + r.n_unmeth)};
      });
  }
};

[[nodiscard]] static auto
merge_methylomes(const std::vector<std::string> &methylome_names,
                 const std::string &methylome_dir,
                 std::error_code &error) -> big_methylome {
  auto &lgr = xfr::logger::instance();
  const auto n_methylomes = std::size(methylome_names);

  double read_time{};
  double merge_time{};

  const auto &last_methylome = methylome_names.back();
  auto read_start = std::chrono::high_resolution_clock::now();
  auto last_meth = xfr::methylome::read(methylome_dir, last_methylome, error);
  auto read_stop = std::chrono::high_resolution_clock::now();
  read_time += duration(read_start, read_stop);
  if (error) {
    lgr.error("Error reading methylome {} {}: {}", methylome_dir,
              last_methylome, error);
    return {};
  }

  big_methylome bmeth(last_meth.meta.n_cpgs);

  bmeth.add(last_meth);

  // ADS: consider threads here; no reason not to

  for (const auto &name :
       methylome_names | std::views::take(n_methylomes - 1)) {
    read_start = std::chrono::high_resolution_clock::now();
    const auto tmp_meth = xfr::methylome::read(methylome_dir, name, error);
    read_stop = std::chrono::high_resolution_clock::now();
    read_time += duration(read_start, read_stop);
    if (error) {
      lgr.error("Error reading methylome {}: {}", methylome_dir, name, error);
      return {};
    }
    const auto is_consistent = last_meth.is_consistent(tmp_meth);
    if (!is_consistent) {
      lgr.error("Inconsistent metadata: {} {}", last_methylome, name);
      throw std::runtime_error("Inconsistent metadata");
    }
    const auto merge_start = std::chrono::high_resolution_clock::now();
    bmeth.add(tmp_meth);
    const auto merge_stop = std::chrono::high_resolution_clock::now();
    merge_time += duration(merge_start, merge_stop);
  }

  std::vector<std::tuple<std::string, std::string>> timing_to_log{
    // clang-format off
    {"read time", std::format("{:.3}s", read_time)},
    {"merge time", std::format("{:.3}s", merge_time)},
    // clang-format on
  };
  xfr::log_args<xfr::log_level_t::debug>(timing_to_log);

  error = last_meth.update_metadata();
  if (error) {
    lgr.error("Error updating metadata: {}", error);
    return {};
  }

  return bmeth;
}

static inline auto
write_bigwig(const std::string &outfile, const big_methylome &bmeth,
             const xfr::genome_index &index, const auto &format_site) {
  static constexpr auto buf_size = 128;
  std::array<char, buf_size> buf{};

  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error(
      std::format("Failed to open output file {}", outfile));

  std::println(out, "track type=wiggle_0");

  auto cpg_itr = std::cbegin(bmeth.cpgs);
  const auto zipped =
    std::views::zip(index.data.positions, index.meta.chrom_order);
  for (const auto [positions, chrom_name] : zipped) {
    const int m =
      std::snprintf(buf.data(), buf_size, "variableStep chrom=%s span=1\n",
                    chrom_name.data());
    if (m < 0 || m > buf_size)
      throw std::runtime_error("failed writing output");
    out.write(buf.data(), m);
    for (const auto pos : positions) {
      if (*cpg_itr != big_mcount_pair{}) {
        const int n = format_site(buf, buf_size, pos, *cpg_itr);
        if (n <= 0)
          throw std::runtime_error("failed writing output");
        out.write(buf.data(), n);
      }
      ++cpg_itr;
    }
  }
}

static inline auto
write_sym(const std::string &outfile, const big_methylome &bmeth,
          const xfr::genome_index &index) {
  static constexpr auto buf_size = 128;
  std::array<char, buf_size> buf{};

  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error(
      std::format("Failed to open output file {}", outfile));

  auto cpg_itr = std::cbegin(bmeth.cpgs);
  const auto zipped =
    std::views::zip(index.data.positions, index.meta.chrom_order);
  for (const auto [positions, chrom_name] : zipped) {
    const int m =
      std::snprintf(buf.data(), buf_size, "%s\t", chrom_name.data());
    if (m < 0 || m > buf_size)
      throw std::runtime_error("failed writing output");
    auto cursor = buf.data() + m;
    auto format_site = [&](const std::uint32_t pos, const big_mcount_pair &p) {
      return std::snprintf(cursor, buf_size, "%u\t+\tCpG\t%.6g\t%d\n", pos + 1,
                           p.get_level(), p.n_reads());
    };
    for (const auto pos : positions) {
      const int n = format_site(pos, *cpg_itr);
      if (n < 0 || n > buf_size)
        throw std::runtime_error("failed writing output");
      out.write(buf.data(), m + n);
      ++cpg_itr;
    }
  }
}

[[nodiscard]] static inline auto
read_methylomes_file(const std::string &filename,
                     std::error_code &ec) -> std::vector<std::string> {
  std::ifstream in(filename);
  if (!in) {
    ec = std::make_error_code(std::errc(errno));
    return {};
  }

  std::vector<std::string> names;
  std::string line;
  while (getline(in, line))
    names.push_back(line);

  ec = std::error_code{};
  return names;
}

auto
main(int argc, char *argv[]) -> int {  // NOLINT(*-c-arrays)
  try {
    static constexpr auto log_level_default = xfr::log_level_t::info;
    static constexpr auto command = "maketracks";
    static const auto usage =
      std::format("Usage: {} [options]", rstrip(command));
    static const auto about_msg =
      std::format("{}: {}", rstrip(command), rstrip(about));
    static const auto description_msg =
      std::format("{}\n{}", rstrip(description), rstrip(examples));

    xfr::log_level_t log_level{};
    std::string methylome_dir{};
    std::string outdir{};
    std::string methylome_names_file{};
    std::string index_dir{};
    std::string genome_name{};
    std::string merged_name{};

    bool write_sym_file{};

    CLI::App app{about_msg};
    argv = app.ensure_utf8(argv);
    app.usage(usage);
    if (argc >= 2)
      app.footer(description_msg);
    app.get_formatter()->column_width(column_width_default);
    app.get_formatter()->label("REQUIRED", "REQD");
    app.set_help_flag("-h,--help", "Print a detailed help message and exit");
    // clang-format off
    app.add_flag("--sym", write_sym_file, "write sym counts file");
    app.add_option("-m,--methylome", methylome_names_file, "file with methylome names");
    app.add_option("-d,--methylome-dir", methylome_dir, "input methylome directory")
      ->required()
      ->check(CLI::ExistingDirectory);
    app.add_option("-o,--output-dir", outdir, "bigWig output directory")
      ->required()
      ->check(CLI::ExistingDirectory);
    app.add_option("-x,--index-dir", index_dir, "genome index directory")
      ->option_text("TEXT:DIR")
      // ->excludes(config_dir_opt)
      ->check(CLI::ExistingDirectory);
    app.add_option("-g,--genome", genome_name, "genome name")
      ->required();
    app.add_option("-n,--name", merged_name, "merged methylome name")
      ->required();
    app.add_option("-v,--log-level", log_level,
                   std::format("log level {}", xfr::log_level_help_str))
      ->option_text(std::format("[{}]", log_level_default))
      ->transform(CLI::CheckedTransformer(xfr::str_to_level, CLI::ignore_case));
    // clang-format on

    if (argc < 2) {
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }

    CLI11_PARSE(app, argc, argv);

    auto &lgr =
      xfr::logger::instance(xfr::shared_from_cout(), command, log_level);
    if (!lgr) {
      std::println("Failure initializing logging: {}.", lgr.get_status());
      return EXIT_FAILURE;
    }

    std::vector<std::tuple<std::string, std::string>> args_to_log{
      // clang-format off
      {"Output directory", outdir},
      {"Methylomes filename", methylome_names_file},
      {"Methylome directory", methylome_dir},
      {"Merged methylome name", merged_name},
      // clang-format on
    };
    xfr::log_args<xfr::log_level_t::info>(args_to_log);

    std::error_code error;

    const auto methylome_names =
      read_methylomes_file(methylome_names_file, error);
    if (error) {
      lgr.error("Failed to read methylome names: {}", error);
      return EXIT_FAILURE;
    }

    const auto bmeth = merge_methylomes(methylome_names, methylome_dir, error);
    if (error) {
      lgr.error("Failed to merge methylomes: {}", error);
      return EXIT_FAILURE;
    }

    const auto index = xfr::genome_index::read(index_dir, genome_name, error);
    if (error) {
      lgr.error("Failed to read genome index {} {}: {}", index_dir, genome_name,
                error);
      return EXIT_FAILURE;
    }

    lgr.info("Making levels bigWig track");
    const auto levels_start = std::chrono::high_resolution_clock::now();
    auto format_site_levels = [](auto &buf, const auto buf_size,
                                 const std::uint32_t pos,
                                 const big_mcount_pair &p) {
      return std::snprintf(buf.data(), buf_size, "%u %0.6g\n", pos + 1,
                           p.get_level());
    };
    auto outfile =
      std::filesystem::path{outdir} / std::format("{}.meth.wig", merged_name);
    write_bigwig(outfile, bmeth, index, format_site_levels);
    const auto levels_stop = std::chrono::high_resolution_clock::now();
    lgr.debug("Time to make levels bigWig track: {:.3}s",
              duration(levels_start, levels_stop));

    lgr.info("Making reads bigWig track");
    const auto reads_start = std::chrono::high_resolution_clock::now();
    auto format_site_reads = [](auto &buf, const auto buf_size,
                                const std::uint32_t pos,
                                const big_mcount_pair &p) {
      return std::snprintf(buf.data(), buf_size, "%u %d\n", pos + 1,
                           p.n_reads());
    };

    outfile =
      std::filesystem::path{outdir} / std::format("{}.reads.wig", merged_name);

    write_bigwig(outfile, bmeth, index, format_site_reads);
    const auto reads_stop = std::chrono::high_resolution_clock::now();
    lgr.debug("Time to make reads bigWig track: {:.3}s",
              duration(reads_start, reads_stop));

    if (write_sym_file) {
      lgr.info("Making symmetric CpGs counts file");
      const auto sym_file =
        std::filesystem::path{outdir} / std::format("{}.sym", merged_name);
      write_sym(sym_file, bmeth, index);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
