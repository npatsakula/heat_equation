#ifndef HEAT_EQUATION_MATH_HPP
#define HEAT_EQUATION_MATH_HPP

#include <range/v3/all.hpp>
using namespace ranges::views;

template<typename T>
void modified_thomas(
        std::vector<T> const & free_part,
        T matrix_above,
        T matrix_main,
        ranges::span<T> result
)
//[[expects: !free_part.empty()]]
//[[expects: free_part.size() == result.size()]]
{
    std::size_t const n = free_part.size();
    std::vector<T> static alpha(n - 1), beta(n -1);
    alpha.resize(n - 1);
    beta.resize(n - 1);

    alpha[n - 2] = -matrix_above / matrix_main;
    beta[n - 2] = free_part.back() / matrix_main;

    for(auto index = n - 2; index > 0; --index) {
        T common_factor = 1.0 / (matrix_main + matrix_above * alpha[index]);
        alpha[index - 1] = -matrix_above * common_factor;
        beta[index - 1] = (free_part[index] - beta[index] * matrix_above) * common_factor;
    }

    result.front() = (free_part.front() - matrix_above * beta.front()) / (matrix_main + matrix_above * alpha.front());
    for(std::size_t const index : ints(std::size_t{ 1 }, n - 1))
        result[index] = alpha[index - 1] * result[index - 1] + beta[index - 1];
}

template<typename T>
[[nodiscard]] std::vector<std::vector<T>> heat_equation(
        T (*source_function)(const T, const T, const double),
        T (*initial_time)(const T),
        T (*left_bound_func)(const T),
        T (*right_bound_func)(const T),
        double const diffusivity,
        std::size_t const h_num,
        std::size_t const t_num
) {
    double const h = 1.0 / h_num;
    double const tau = 1.0 / t_num;
    double const courant = (diffusivity * tau) / std::pow(h, 2);
    T const matrix_above = (-diffusivity * tau) / std::pow(h, 2);
    T const matrix_main = 1 + ((2 * diffusivity * tau) / std::pow(h, 2));

    std::vector<std::vector<T>> result(t_num);
    std::vector<T> free_part(h_num - 2);

    // Copy elision, no `move` necessary.
    result[0] = ints(std::size_t{0}, h_num) | transform([initial_time, h](auto h_iter) {
        return initial_time(h_iter * h);
    }) | ranges::to_vector;

    for (std::size_t const time_iter : ints(std::size_t{1}, t_num)) {
        double const time = time_iter * tau;

        for (std::size_t const space_iter : ints(std::size_t{1}, h_num - 1))
            free_part[space_iter - 1] = result[time_iter - 1][space_iter] + tau * source_function(space_iter * h, time, diffusivity);

        T const left_bound = left_bound_func(time);
        T const right_bound = right_bound_func(time);
        free_part.front() += courant * left_bound;
        free_part.back() += courant * right_bound;

        std::vector<T> & current_time_layer = result[time_iter] = std::vector<T>(h_num);
        T * begin = &current_time_layer.at(1);
        T * end = &current_time_layer.at(current_time_layer.size() - 1);

        modified_thomas(free_part, matrix_above, matrix_main, ranges::span<T>(begin, end));
        *current_time_layer.begin() = left_bound;
        *current_time_layer.end() = right_bound;
    }

    return result;
}

template<typename T>
double find_error(std::vector<std::vector<T>> const & source, T (*exact_solution)(T, T)) {
    double const tau = 1.0 / source.size();
    double max_error = 0.0;

    for(std::size_t const time_index : ints(std::size_t{ 0 }, source.size())) {
        std::vector<T> const & current_time_layer = source[time_index];
        double const h = 1.0 / current_time_layer.size();

        for (std::size_t const space_index : ints(std::size_t{ 0 }, current_time_layer.size())) {
            T const exact = exact_solution(h * space_index, tau * time_index);
            max_error = std::max(exact, source[time_index][space_index]);
        }
    }

    return max_error;
}

#endif
