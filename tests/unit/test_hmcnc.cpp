#include <gtest/gtest.h>

#include "hmcnc.h"

TEST(hmcnc, test_setup_check) {
  EXPECT_EQ(42, hmcnc_test());
}
