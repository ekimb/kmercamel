#pragma once

#include "streaming.h"
#include "gtest/gtest.h"
namespace {

// Retrieving current path on Windows does not work as on linux
// therefore the following unittest is linux-specific.
#ifdef __unix__
    TEST(StreamingTest, Streaming) {
        std::string path = std::filesystem::current_path();
        path += "/tests/test.fa";
        struct TestCase {
            int k;
            bool complements;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {3, true, ">superstring\nACCCGAacCGtaATgcTTta\n"},
                {2, true, ">superstring\nACcCGAaTaATGc\n"},
                {1, true, ">superstring\nAC\n"},
                {1, false, ">superstring\nACGT\n"},
        };


        for (auto t : tests) {
            std::stringstream of;

            Streaming(path, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
#endif
    TEST(PushTest, Push) {
        struct TestCase {
            std::string current;
            int k;
            int32_t used;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {"ACGT", 4, 0b010, "cg"},
                {"ACGT", 4, 0b101, "cgt"},
                {"ACGT", 4, 0b100, "c"},
                {"ACGT", 4, 0b000, ""},
                {"CGT", 4, 0b000, ""},
        };


        for (auto t : tests) {
            std::stringstream of;

            Push(of, t.current, t.k, t.used);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}
