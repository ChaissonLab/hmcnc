#include <gtest/gtest.h>

#include "../include/hmmcnc.h"

TEST(hmmcnc, test_setup_check) {
  EXPECT_EQ(42, hmmcnc_test());
}
