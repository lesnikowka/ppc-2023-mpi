// Copyright 2023 Kulikov Artem
#include <gtest/gtest.h>
#include <vector>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>

#include "./hypercube.h"


TEST(dummy, test_1) {
    ASSERT_EQ(1, 1);
}

TEST(dummy, test_2) {
    ASSERT_EQ(1, 1);
}

TEST(dummy, test_3) {
    ASSERT_EQ(1, 1);
}

TEST(dummy, test_4) {
    ASSERT_EQ(1, 1);
}

TEST(dummy, test_5) {
    ASSERT_EQ(1, 1);
}


int main(int argc, char** argv) {
    boost::mpi::environment env(argc, argv, boost::mpi::threading::multiple);
    if (env.thread_level() < boost::mpi::threading::multiple) {
        env.abort(-1);
    }
    boost::mpi::communicator world;
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (world.rank() != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }
    // return RUN_ALL_TESTS();

    return 0;
}

