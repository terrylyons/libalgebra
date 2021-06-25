//
// Created by sam on 25/06/2021.
//

#ifndef LIBALGEBRA_COEFFICIENTS_H
#define LIBALGEBRA_COEFFICIENTS_H


namespace coefficients {

#define LIBALGEBRA_FIELD_GENERATE_BINARY(NAME, OP, RET_T, LHS_T, RHS_T) \
    static RET_T NAME(LHS_T lhs, RHS_T rhs) {                           \
        return lhs OP rhs;                                              \
    }

    template <typename CoeffType>
    struct coefficient_field {

        typedef CoeffType SCA;
        typedef CoeffType RAT;

        static const SCA one;
        static const SCA zero;
        static const SCA mone;

        static SCA uminus(SCA arg) {
            return -arg;
        }

        LIBALGEBRA_FIELD_GENERATE_BINARY(add, +, SCA, SCA, SCA);

        LIBALGEBRA_FIELD_GENERATE_BINARY(sub, -, SCA, SCA, SCA);

        LIBALGEBRA_FIELD_GENERATE_BINARY(mul, *, SCA, SCA, SCA);

        LIBALGEBRA_FIELD_GENERATE_BINARY(div, /, SCA, SCA, RAT);


    };

    template <typename CoeffType>
    const CoeffType coefficient_field<CoeffType>::one(1);

    template <typename CoeffType>
    const CoeffType coefficient_field<CoeffType>::zero(0);

    template <typename CoeffType>
    const CoeffType coefficient_field<CoeffType>::mone(-1);




#undef LIBALGEBRA_FIELD_GENERATE_BINARY


    typedef coefficient_field<double> double_field;
    typedef coefficient_field<float>  float_field;

}


#endif //LIBALGEBRA_COEFFICIENTS_H
