//
// Created by sam on 06/08/2021.
//

#include <UnitTest++.h>


#include <libalgebra/dense_storage.h>

#include "../../common/rng.h"


SUITE(dense_storage) {

    struct fixture
            {
        using dstorage = alg::vectors::dense_storage<double>;
        using dvector = std::vector<double>;
        using size_type = typename dstorage::size_type;

        using types = typename dstorage::vec_type;

        dvector data1, data2;

        fixture() : data1(), data2()
        {
            mt19937 rng(12345);
            NORMAL_DIST<double> dist(0.0, 1.0);

            data1.reserve(10);
            for (int i=0; i<10; ++i) {
                data1.push_back(dist(rng));
            }

            data2.reserve(10);
            for (int i=0; i<10; ++i) {
                data2.push_back(dist(rng));
            }
        }

        double const* data1_cbegin() const {
            return &data1[0];
        }

        double const* data1_cend() const {
            return &data1[0] + 10;
        }

        double const* data2_cbegin() const {
            return &data2[0];
        }

        double const* data2_cend() const {
            return &data2[0] + 10;
        }

        double* data1_begin() {
            return &data1[0];
        }

        double* data1_end() {
            return &data1[0] + 10;
        }

        double* data2_begin() {
            return &data2[0];
        }

        double* data2_end() {
            return &data2[0] + 10;
        }


};


    TEST_FIXTURE(fixture, test_storage_default_ctor) {

        dstorage val;

        CHECK_EQUAL(0, val.size());
        CHECK_EQUAL(types::owned, val.type());
        CHECK(nullptr == val.begin());
        CHECK(nullptr == val.end());

    }

    TEST_FIXTURE(fixture, test_ctor_size_default_construct) {
        dstorage val(10);

        CHECK_EQUAL(10, val.size());
        CHECK_EQUAL(types::owned, val.type());

        for (auto const& v : val) {
            CHECK_EQUAL(0.0, v);
        }
    }


    TEST_FIXTURE(fixture, test_borrowed_ctor_ptrs) {
        dstorage val(data1_cbegin(), data1_cend());

        CHECK_EQUAL(data1.size(), val.size());
        CHECK_EQUAL(data1_cbegin(), val.cbegin());
        CHECK_EQUAL(types::borrowed, val.type());
    }

    TEST_FIXTURE(fixture, test_borrowed_mut_ctor_ptrs) {
        dstorage val(data1_begin(), data1_end());

        CHECK_EQUAL(data1.size(), val.size());
        CHECK_EQUAL(data1_cbegin(), val.cbegin());
        CHECK_EQUAL(&data1[0], val.begin());
        CHECK_EQUAL(types::borrowed_mut, val.type());
    }

    TEST_FIXTURE(fixture, test_borrowed_ctor_ptr_size) {
        dstorage val(data1_cbegin(), data1.size());

        CHECK_EQUAL(data1.size(), val.size());
        CHECK_EQUAL(data1_cbegin(), val.cbegin());
        CHECK_EQUAL(types::borrowed, val.type());
    }

    TEST_FIXTURE(fixture, test_borrowed_mut_ctor_ptr_size) {
        dstorage val(data1_begin(), data1.size());

        CHECK_EQUAL(data1.size(), val.size());
        CHECK_EQUAL(data1_cbegin(), val.cbegin());
        CHECK_EQUAL(&data1[0], val.begin());
        CHECK_EQUAL(types::borrowed_mut, val.type());
    }

    TEST_FIXTURE(fixture, test_offset_ctor_copy) {
        dstorage val(1, data1_cbegin(), data1_cend());

        CHECK_EQUAL(data1.size() + 1, val.size());
        CHECK_EQUAL(types::owned, val.type());

        double const* ptr = val.cbegin();
        CHECK_EQUAL(0.0, *(ptr++));

        for (auto const& v : data1) {
            CHECK_EQUAL(v, *(ptr++));
        }
    }

    TEST_FIXTURE(fixture, test_offset_ctor_move) {
        // We don't want to move out of our fixture, so clone a small part of it locally
        double tmp[5];
        std::copy(data1_cbegin(), data1_cbegin()+5, tmp);

        dstorage val(1, tmp, tmp+5);

        CHECK_EQUAL(6, val.size());
        CHECK_EQUAL(types::owned, val.type());

        double const* ptr = val.cbegin();
        CHECK_EQUAL(0.0, *(ptr++));

        for (int i=1; i<=5; ++i) {
            CHECK_EQUAL(data1[i-1], *(ptr++));
        }

    }

    TEST_FIXTURE(fixture, test_copy_ctor) {
        dstorage val1(data1_cbegin(), data1_cbegin());

        dstorage val2(val1);

        CHECK_EQUAL(val1.size(), val2.size());
        CHECK_EQUAL(types::owned, val2.type());

        double const* ptr1(val1.cbegin());
        double const* ptr2(val2.cbegin());

        for (size_type i=0; i<val1.size(); ++i) {
            CHECK_EQUAL(*(ptr1++), *(ptr2++));
        }

    }

    TEST_FIXTURE(fixture, test_copy_assignment) {
        dstorage val1(data1_cbegin(), data1_cend()), val2;

        val2 = val1;

        CHECK_EQUAL(val1.size(), val2.size());
        CHECK_EQUAL(types::owned, val2.type());

        double const* ptr1(val1.cbegin());
        double const* ptr2(val2.cbegin());

        for (size_type i=0; i<val1.size(); ++i) {
            CHECK_EQUAL(*(ptr1++), *(ptr2++));
        }
    }

    TEST_FIXTURE(fixture, test_move_assignment) {
        dstorage val1(data1_cbegin(), data1_cend()), val2;

        val2 = std::move(val1);

        CHECK_EQUAL(data1.size(), val2.size());
        CHECK_EQUAL(types::borrowed, val2.type());

        double const* ptr1(data1_cbegin());
        double const* ptr2(val2.cbegin());

        for (size_type i=0; i<data1.size(); ++i) {
            CHECK_EQUAL(*(ptr1++), *(ptr2++));
        }
    }

    TEST_FIXTURE(fixture, test_mutable_iterator_begin_promote_borrowed) {
        dstorage val1(data1_cbegin(), data1_cend());

        REQUIRE CHECK_EQUAL(types::borrowed, val1.type());
        dstorage::size_type sz = val1.size();

        val1.begin();

        CHECK_EQUAL(types::owned, val1.type());
        CHECK_EQUAL(sz, val1.size());
    }

    TEST_FIXTURE(fixture, test_mutable_iterator_end_promote_borrowed) {
        dstorage val1(data1_cbegin(), data1_cend());

        REQUIRE CHECK_EQUAL(types::borrowed, val1.type());
        dstorage::size_type sz = val1.size();

        val1.end();

        CHECK_EQUAL(types::owned, val1.type());
        CHECK_EQUAL(sz, val1.size());
    }

    TEST_FIXTURE(fixture, test_mutable_access_promote_borrowed) {
        dstorage val1(data1_cbegin(), data1_cend());

        REQUIRE CHECK_EQUAL(types::borrowed, val1.type());
        dstorage::size_type sz = val1.size();

        val1[0];

        CHECK_EQUAL(types::owned, val1.type());
        CHECK_EQUAL(sz, val1.size());
    }

    TEST_FIXTURE(fixture, test_mutable_iterator_begin_no_promote_borrowed_mut) {
        dstorage val1(data1_begin(), data1_end());

        REQUIRE CHECK_EQUAL(types::borrowed_mut, val1.type());
        dstorage::size_type sz = val1.size();

        val1.begin();

        CHECK_EQUAL(types::borrowed_mut, val1.type());
        CHECK_EQUAL(sz, val1.size());
    }

    TEST_FIXTURE(fixture, test_mutable_iterator_end_no_promote_borrowed_mut) {
        dstorage val1(data1_begin(), data1_end());

        REQUIRE CHECK_EQUAL(types::borrowed_mut, val1.type());
        dstorage::size_type sz = val1.size();

        val1.end();

        CHECK_EQUAL(types::borrowed_mut, val1.type());
        CHECK_EQUAL(sz, val1.size());
    }

    TEST_FIXTURE(fixture, test_mutable_access_no_promote_borrowed_mut) {
        dstorage val1(data1_begin(), data1_end());

        REQUIRE CHECK_EQUAL(types::borrowed_mut, val1.type());
        dstorage::size_type sz = val1.size();

        val1[0];

        CHECK_EQUAL(types::borrowed_mut, val1.type());
        CHECK_EQUAL(sz, val1.size());
    }

    TEST_FIXTURE(fixture, test_resize_up_blank_storage) {
        dstorage val;

        REQUIRE CHECK_EQUAL(0, val.size());

        val.resize(10);
        CHECK_EQUAL(10, val.size());
        CHECK_EQUAL(types::owned, val.type());

        for (auto const& v : val) {
            CHECK_EQUAL(0.0, v);
        }
    }

    TEST_FIXTURE(fixture, test_resize_up_blank_storage_with_val) {
        dstorage val;

        REQUIRE CHECK_EQUAL(0, val.size());

        val.resize(10, 1.0);
        CHECK_EQUAL(10, val.size());
        CHECK_EQUAL(types::owned, val.type());

        for (auto const& v : val) {
            CHECK_EQUAL(1.0, v);
        }
    }

    TEST_FIXTURE(fixture, test_resize_up_borrowed_promotes) {
        dstorage val(data1_cbegin(), data1_cend());
        dstorage::size_type sz = val.size();

        dvector newdata(data1);
        newdata.push_back(0.0);
        dstorage expected(dstorage::make_owned(&newdata[0], newdata.size()));

        REQUIRE CHECK_EQUAL(types::borrowed, val.type());

        val.resize(sz+1);

        CHECK_EQUAL(sz+1, val.size());
        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(expected, val);
    }

    TEST_FIXTURE(fixture, test_resize_up_borrowed_mut_promotes) {
        dstorage val(data1_begin(), data1_end());
        dstorage::size_type sz = val.size();

        dvector newdata(data1);
        newdata.push_back(0.0);
        dstorage expected(dstorage::make_owned(&newdata[0], newdata.size()));

        REQUIRE CHECK_EQUAL(types::borrowed_mut, val.type());

        val.resize(sz+1);

        CHECK_EQUAL(sz+1, val.size());
        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(expected, val);
    }

    TEST_FIXTURE(fixture, test_resize_down_owned) {
        dstorage val(dstorage::make_owned(data1_cbegin(), data1.size()));
        dstorage expected(data1_cbegin(), 5);

        REQUIRE CHECK_EQUAL(data1.size(), val.size());

        val.resize(5);

        CHECK_EQUAL(5, val.size());
        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(expected, val);
    }


    TEST_FIXTURE(fixture, test_resize_down_borrowed) {
        dstorage val(data1_cbegin(), data1.size());
        dstorage expected(data1_cbegin(), 5);

        REQUIRE CHECK_EQUAL(data1.size(), val.size());

        val.resize(5);

        CHECK_EQUAL(5, val.size());
        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(expected, val);
    }

    TEST_FIXTURE(fixture, test_resize_down_borrowed_mut) {
        dstorage val(data1_begin(), data1.size());
        dstorage expected(data1_cbegin(), 5);

        REQUIRE CHECK_EQUAL(data1.size(), val.size());

        val.resize(5);

        CHECK_EQUAL(5, val.size());
        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(expected, val);
    }

    TEST_FIXTURE(fixture, test_clear_borrowed) {
        dstorage val(data1_cbegin(), data1.size());
        REQUIRE CHECK_EQUAL(types::borrowed, val.type());

        val.clear();

        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(0, val.size());
        CHECK(nullptr == val.begin());
    }

    TEST_FIXTURE(fixture, test_clear_borrowed_mut) {
        dstorage val(data1_begin(), data1.size());
        REQUIRE CHECK_EQUAL(types::borrowed_mut, val.type());

        val.clear();

        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(0, val.size());
        CHECK(nullptr == val.begin());
    }

    TEST_FIXTURE(fixture, test_clear_owned) {
        dstorage val(dstorage::make_owned(data1_begin(), data1.size()));
        REQUIRE CHECK_EQUAL(types::owned, val.type());

        val.clear();

        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(0, val.size());
        CHECK(nullptr == val.begin());
    }

    TEST_FIXTURE(fixture, test_copy_extend) {
        dstorage val(data1_cbegin(), data1_cend());
        REQUIRE CHECK_EQUAL(data1.size(), val.size());

        val.copy_extend(data2_cbegin(), data2_cend());

        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(data1.size() + data2.size(), val.size());

        for (size_type i = 0; i<data1.size(); ++i) {
            CHECK_EQUAL(data1[i], val[i]);
        }

        for (size_type i=0; i<data2.size(); ++i) {
            CHECK_EQUAL(data2[i], val[data1.size() + i]);
        }

    }

    TEST_FIXTURE(fixture, test_move_extend) {
        dstorage val(data1_cbegin(), data1_cend());
        REQUIRE CHECK_EQUAL(data1.size(), val.size());

        double tmp[5];
        std::copy(data2_cbegin(), data2_cbegin()+5, tmp);

        val.move_extend(tmp, tmp+5);

        CHECK_EQUAL(types::owned, val.type());
        CHECK_EQUAL(data1.size() + 5, val.size());

        for (size_type i = 0; i<data1.size(); ++i) {
            CHECK_EQUAL(data1[i], val[i]);
        }

        for (int i=0; i<5; ++i) {
            CHECK_EQUAL(data2[i], val[data1.size() + i]);
        }

    }


    struct fixture2
    {

        static int instances;

        struct stand_in
        {

            int number;

            stand_in() : number{instances}
            {
                ++instances;
            }

            stand_in(stand_in const& other) : number{other.number}
            {
                ++instances;
            }

            stand_in(stand_in const&& other) noexcept : number {other.number}
            {}

            stand_in& operator=(stand_in&& other) noexcept
            {
                number = other.number;
                --instances;
                return *this;
            }

            ~stand_in()
            {
                --instances;
            }

        };

        using sstorage = alg::vectors::dense_storage<stand_in>;


        fixture2()
        {
            instances = 0;
        }


    };

    int fixture2::instances = 0;

    TEST_FIXTURE(fixture2, test_dtor_deallocates) {
        REQUIRE CHECK_EQUAL(0, instances);

        {
            sstorage val(10);
            REQUIRE CHECK_EQUAL(10, instances);
        }

        CHECK_EQUAL(0, instances);
    }




}
