cc_library(
    name = "rmq",
    hdrs = ["rmq.h"],
)

cc_test(
    name = "rmq-test",
    srcs = ["rmq_test.cc"],
    copts = ["-Iexternal/gtest/googletest-release-1.8.0/googletest/include"],
    deps = ["@gtest//:main",
            ":rmq",],
)
