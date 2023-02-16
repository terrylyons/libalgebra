static SHUFFLE_TENSOR shift_down(const SHUFFLE_TENSOR& sh, typename SHUFFLE_TENSOR::KEY word)
{
    SHUFFLE_TENSOR result(sh), working;
    while (word.size()) {
        auto letter = word.lparent();
        word = word.rparent();
        for (auto pr : result) {
            if (pr.key().size() > 0 && pr.key().lparent() == letter)
                working[pr.key().rparent()] = result[pr.key()];
        }
        result.swap(working);
        working.clear();
    }
    return result;
}

// the evaluation of the adjoint operation is worked out below
// <sh,ab>=\sum_{uv=sh}<ua><vb>
//        = <\sum_{uv=sh}<ua>v,b>
// Let T_w(sh) be all the projection of sh onto the part beginning with w with w removed
// \sum_{i} <ki shi,ab> =
//        = \sum_{i} ki<\sum_{uv=shi}<ua>v,b>
//        = \sum_{u} < <ua>T_u(sh), b>
//  The action of the adjoint of tensor multiplication by a multiplication is \sum_u <ua> T_u(sh)

static SHUFFLE_TENSOR adjoint_to_multiply(const TENSOR& t, SHUFFLE_TENSOR sh)
{
    // this implementation is understandable and reliable but repetitive and can be radically accelerated
    SHUFFLE_TENSOR result;
    for (auto& pr : t) {
        result += shift_down(sh, pr.key()) * pr.value();
    }
    return result;
}