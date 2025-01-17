#pragma once
#include "generalized_simplitigs_ac.h"

#include "gtest/gtest.h"

namespace {
    TEST(ExtensionACTest, ExtensionAC) {
        struct TestCase {
            std::vector<bool> forbidden;
            std::list<size_t> incidentKMers;
            bool complements;
            int wantResult;
            std::vector<bool> wantForbidden;
            std::list<size_t> wantIncidentKMers;
        };
        std::vector<TestCase> tests = {
                {{0,0,0}, {1, 0}, false, 1, {0, 1, 0}, {0}},
                {{0,1,1}, {1, 0}, false, 0, {1, 1, 1}, {}},
                {{0,0,0, 0}, {}, false, -1, {0,0,0,0}, {}},
                {{1,1,0, 1}, {3, 0, 1}, false, -1, {1, 1, 0, 1}, {}},
                {{0,0,0,0}, {3, 0, 1}, true, 3, {0, 1, 0, 1}, {0,1}},
        };

        for (auto t: tests) {
            size_t gotResult = ExtensionAC(t.forbidden, t.incidentKMers, t.complements);

            EXPECT_EQ(t.wantResult, gotResult);
            EXPECT_EQ(t.wantForbidden, t.forbidden);
            EXPECT_EQ(t.wantIncidentKMers, t.incidentKMers);
        }
    }

    TEST(GreedyGeneralizedSimplitigsACTest, GreedyGeneralizedSimplitigsAC) {
        struct TestCase {
            std::vector<KMer> kMers;
            int k;
            int d_max;
            bool complements;
            std::string wantSuperstring;
            std::vector<bool> wantMask;
        };
        std::vector<TestCase> tests = {
                // As behavior of the unordered_set.begin() is not specified, some tests are commented, as they could fail otherwise.
                // Uncommenting them may add additional check but could also add false positives.
                //{ {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"CCCC"}, KMer{"AACA"}}, 4, 3,
                //  "AACAATTTCCCC", {1,1,0,0,1, 0,0,0, 1, 0,0,0}},
                { {KMer{"GCT"}, KMer{"TAA"}, KMer{"AAA"}}, 3, 2, false,
                        "GCTAAA", {1,0, 1,1, 0,0}},
                { {KMer{"TAA"}, KMer{"AAA"}, KMer{"GCT"}}, 3, 2, false,
                        "GCTAAA", {1,0, 1,1, 0,0}},
                {{KMer{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}, KMer{"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"}}, 31, 5, false,
                        "TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA", std::vector<bool> {1,1,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0}},
                { {KMer{"TAA"}, KMer{"TTT"}}, 3, 2, true,
                        "TAAA", {1,1, 0,0}},
        };

        for (auto t: tests) {
            auto gotResult = GreedyGeneralizedSimplitigsAC(t.kMers, t.k, t.d_max, t.complements);

            EXPECT_EQ(t.wantSuperstring, gotResult.superstring);
            EXPECT_EQ(t.wantMask, gotResult.mask);
            EXPECT_EQ(t.k, gotResult.k);
        }
    }


}
