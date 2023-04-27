//
// Created by sam on 27/01/2022.
//

#include <UnitTest++.h>

#include <libalgebra/alg_types.h>
#include <libalgebra/libalgebra.h>
#include <libalgebra/free_extension.h>
#include <libalgebra/operators.h>

template<typename Coeffs, unsigned Width, unsigned Depth>
struct key_function {
    using input_key = alg::_tensor_basis<Width, Depth>;
    using output_key = alg::_tensor_basis<Width - 1, Depth>;
    using out_tensor = alg::free_tensor<Coeffs, Width - 1, Depth>;

    static out_tensor eval(const input_key& key)
    {
        auto sz = key.size();
        if (sz == 0) {
            return out_tensor(Coeffs::one);
        }
        else if (sz == 1) {
            auto first_letter = key.FirstLetter();
            if (first_letter == Width) {
                out_tensor result;
                for (alg::LET l = 1; l < Width; ++l) {
                    result.add_scal_prod(output_key(l), Coeffs::one);
                }
                return result;
            }
            else {
                return out_tensor(first_letter, Coeffs::one);
            }
        }
        else {
            return eval(key.lparent()) * eval(key.rparent());
        }
    }

    out_tensor operator()(const input_key& key) const
    {
        return eval(key);
    }
};

SUITE(free_extension_operator)
{

    struct fixture : alg_types<5, 5, Rational> {
        using types = alg_types<5, 5, Rational>;
        using out_types = alg_types<5, 4, Rational>;

        using typename types::COEFF;
        using key_func = key_function<COEFF, 5, 5>;

        using ext_op = alg::operators::free_extension_operator_impl<
                alg::no_caching_tag,
                typename types::TENSOR::BASIS,
                key_func>;
        using linear_op = alg::operators::linear_operator<ext_op>;

        linear_op op;
    };

    TEST_FIXTURE(fixture, test_zero)
    {
        TENSOR in;
        typename out_types::TENSOR expected;
        auto out = op(in);

        CHECK_EQUAL(expected, out);
    }

    TEST_FIXTURE(fixture, test_empty_word)
    {
        TENSOR in(S(1));
        typename out_types::TENSOR expected(S(1));

        auto out = op(in);
        CHECK_EQUAL(expected, out);
    }

    TEST_FIXTURE(fixture, test_first_letters)
    {
        for (LET l = 1; l < ALPHABET_SIZE; ++l) {
            TENSOR in(l, S(1));
            typename out_types::TENSOR expected(l, S(1));

            auto out = op(in);
            CHECK_EQUAL(expected, out);
        }
    }

    TEST_FIXTURE(fixture, test_final_letter)
    {
        TENSOR in(LET(ALPHABET_SIZE), S(1));
        typename out_types::TENSOR expected;
        for (LET l = 1; l < ALPHABET_SIZE; ++l) {
            expected += typename out_types::TENSOR(l, S(1));
        }
        auto out = op(in);

        CHECK_EQUAL(expected, out);
    }
}
