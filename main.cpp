#include <fmt/format.h>
#include "cli.hpp"
#include "math.hpp"

double source_function(double const x, double const t, double const diffusivity) {
    return std::pow(x, 2) * (std::pow(x, 2) - 12 * diffusivity * t) +
           std::exp(x) * (x * t * (diffusivity * t - 2) + diffusivity * std::pow(t, 2) + 2 * t);
}

double exact_solution(double const x, double const t) {
    return t * std::pow(x, 4) - std::pow(t, 2) * std::exp(x) * (x - 1) + 1;
}

double left_bound(double const t) {
    return exact_solution(0, t);
}

double right_bound(double const t) {
    return exact_solution(1, t);
}

double init_layer(double const x) {
    return exact_solution(x, 0);
}

std::int32_t main(std::int32_t argc, char ** argv) {
    CLI::App application{"Heat equation calculator tool."};

    std::size_t h_num;
    application.add_option("--space-n", h_num, "Space spatial splits number.")
            ->required();

    std::size_t t_num;
    application.add_option("--time-n", t_num, "Time layers number.")
            ->required();

    bool err_calculation;
    application.add_flag("--find-error", err_calculation, "Find calculation error.")
            ->default_val(true);

    CLI11_PARSE(application, argc, argv);

    std::vector<std::vector<double>> const result =
            heat_equation<double>(source_function, init_layer, right_bound, left_bound, 0.010417, h_num, t_num);

    if (err_calculation) {
        fmt::print(FMT_STRING("Calculation error: {}."), find_error(result, exact_solution));
    }

    return EXIT_SUCCESS;
}
