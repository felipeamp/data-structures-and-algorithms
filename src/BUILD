cc_library(
    name = "rmq",
    hdrs = ["rmq.h"],
    copts = ["--std=c++14"],
)

cc_test(
    name = "rmq-test",
    srcs = ["rmq_test.cc"],
    copts = ["-Iexternal/gtest/googletest-release-1.8.0/googletest/include",
             "--std=c++14"],
    deps = ["@gtest//:main",
            ":rmq",],
)
