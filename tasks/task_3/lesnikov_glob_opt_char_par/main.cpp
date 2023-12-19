#include <gtest/gtest.h>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <cmath>
#include "task_3/lesnikov_glob_opt_char_par/glob_opt_char_par.h"

TEST(dummy, dummy_test_2) {
    ADD_FAILURE();
}   
TEST(dummy, dummy_test_3) {
    ADD_FAILURE();
}
TEST(dummy, dummy_test_4) {
    ADD_FAILURE();
}
TEST(dummy, dummy_test_5) {
    ADD_FAILURE();
}
TEST(lesnikov_glob_opt, run_test) {
    std::function<double(double)> f = [](double x) {return x / 2 + std::sin(x);};

    double val = getMinSequential(f, 1.5, 6, 0.001, 100000, 3);

    printf("\n\nFROM TEST: %lf\n\n", val);

    ADD_FAILURE();
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
