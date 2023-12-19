#include <gtest/gtest.h>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <cmath>
#include "task_3/lesnikov_glob_opt_char_par/glob_opt_char_par.h"

TEST(dummy, dummy_test_2) {
    EXPECT_TRUE(true);
}   
TEST(dummy, dummy_test_3) {
    EXPECT_TRUE(true);
}
TEST(dummy, dummy_test_4) {
    EXPECT_TRUE(true);
}
TEST(lesnikov_glob_opt, run_test_par) {
    std::function<double(double)> f = [](double x) {return x / 2 + std::sin(x);};

    double val = getMinParallel(f, 1.5, 6, 0.001, 200, 3);

    printf("\n\nFROM TEST PARALLEL: %lf\n\n", val);

    EXPECT_TRUE(true);
}
TEST(lesnikov_glob_opt_sequential, test_func_1) {
    std::function<double(double)> f = [](double x) {return x / 2 + std::sin(x);};

    double val = getMinSequential(f, 1.5, 6, 0.001, 200, 3);

    printf("\n\nFROM TEST: %lf\n\n", val);

    EXPECT_TRUE(std::abs(4.189023 - val) < 0.01);
}

TEST(lesnikov_glob_opt_sequential, test_func_2) {
    std::function<double(double)> f = [](double x) {return x * x;};

    double val = getMinSequential(f, -3, 3, 0.001, 200, 3);

    printf("\n\nFROM TEST: %lf\n\n", val);

    EXPECT_TRUE(std::abs(0 - val) < 0.01);
}

TEST(lesnikov_glob_opt_sequential, test_func_3) {
    std::function<double(double)> f = [](double x) {return x / 2 + std::sin(x);};

    double val = getMinSequential(f, 0, 6, 0.001, 200, 3);

    printf("\n\nFROM TEST: %lf\n\n", val);

    EXPECT_TRUE(std::abs(0 - val) < 0.01);
}

TEST(lesnikov_glob_opt_parallel, test_func_1) {
    std::function<double(double)> f = [](double x) {return x / 2 + std::sin(x);};

    double val = getMinSequential(f, 1.5, 6, 0.001, 200, 3);
    double val2 = getMinParallel(f, 1.5, 6, 0.001, 200, 3);

    printf("\n\nFROM TEST: %lf\n\n", val);

    EXPECT_TRUE(std::abs(val2 - val) < 0.01);
}

TEST(lesnikov_glob_opt_parallel, test_func_3) {
    std::function<double(double)> f = [](double x) {return x / 2 + std::sin(x);};

    double val = getMinSequential(f, 0, 6, 0.001, 200, 3);
    double val2 = getMinSequential(f, 0, 6, 0.001, 200, 3);

    printf("\n\nFROM TEST: %lf\n\n", val);

    EXPECT_TRUE(std::abs(val2 - val) < 0.01);
}

int main(int argc, char** argv) {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (world.rank() != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }
    return RUN_ALL_TESTS();
}
