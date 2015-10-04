#include "gtest/gtest.h"
#include <uzlmath>

using namespace uzlmath;

// Tests if the default constructor of matrix works correctly
TEST(MatrixTest, DefConstructor)
{
    // create matrix by use of default constructor
    matrix< double > A;
    
    // check if ivars are initialized correct
    ASSERT_EQ(A.cols, 0);
    ASSERT_EQ(A.cols, 0);
    ASSERT_EQ(A.mem, nullptr);
}

// Tests if the constructor for creating a matrix of size mxn
// works correctly
TEST(MatrixTest, M_N_Constructor)
{
    // create matrix by use of mxn-constructor
    size_t m = 10;
    size_t n = 5;
    
    matrix< double > A(m, n);
    
    // check if ivars are initialized correct
    ASSERT_EQ(A.rows, m);
    ASSERT_EQ(A.cols, n);
    ASSERT_NE(A.mem, nullptr);
}

// Tests if the constructor for creating a square matrix works
// correctly
TEST(MatrixTest, MN_Constructor)
{
    // create matrix by use of mn-constructor
    size_t mn = 7;
    
    matrix< double > A(mn, mn);
    
    ASSERT_EQ(A.rows, mn);
    ASSERT_EQ(A.cols, mn);
    ASSERT_EQ(A.cols, A.rows);
    ASSERT_NE(A.mem, nullptr);
}

// Tests if the constructor with initializing values works
// correctly
TEST(MatrixTest, M_N_init_Constructor)
{
    // create matrix by use of mxn-constructor with
    // initializing values
    size_t m = 5;
    size_t n = 4;
    
    typedef double T;
    
    T initial = 3.14;
    matrix< T > A(m, n, initial);
    
    // check bounds
    ASSERT_EQ(A.rows, m);
    ASSERT_EQ(A.cols, n);
    ASSERT_NE(A.mem, nullptr);
    
    // check values
    for (size_t i = 0; i < m * n; ++i)
    {
        ASSERT_EQ(A.mem[i], initial);
    }
}

// Test case of copy constructor
TEST(MatrixTest, CopyConstructor)
{
    // create a matrix which is getting copied
    size_t m = 5;
    size_t n = 4;
    
    matrix< double > A(m, n);
    
    // fill matrix with random values
    srand(0);
    for (size_t i = 0; i < m * n; ++i)
    {
        access::rw(A.mem[i]) = (static_cast< double >(rand()) / (RAND_MAX));
    }
    
    // create matrix via copy constructor
    matrix< double > B(A);
    
    // check ivars and properties
    ASSERT_EQ(B.rows, m);
    ASSERT_EQ(B.cols, n);
    ASSERT_NE(B.mem, nullptr);
    
    // check values
    for (size_t i = 0; i < m * n; ++i)
    {
        ASSERT_EQ(A.mem[i], B.mem[i]);
    }
}

// Test case of move constructor
TEST(MatrixTest, MoveConstructor)
{
    // create two matrices that are added and then assigned to matrix
    size_t mn = 5;
    
    matrix< double > A(mn, mn), B(mn, mn);
    
    // fill matrix with random values
    srand(0);
    for (size_t i = 0; i < mn * mn; ++i)
    {
        access::rw(A.mem[i]) = (static_cast< double >(rand()) / (RAND_MAX));
        access::rw(B.mem[i]) = (static_cast< double >(rand()) / (RAND_MAX));
    }
    
    // create matrix via assignment and move constructor
    matrix< double > C = A + B;
    
    // check ivars and properties
    ASSERT_EQ(C.rows, mn);
    ASSERT_EQ(C.cols, mn);
    ASSERT_NE(C.mem, nullptr);
    
    // check values
    for (size_t i = 0; i < mn * mn; ++i)
    {
        ASSERT_EQ(C.mem[i], A.mem[i] + B.mem[i]);
    }
}

// Test case for addition of matrix objects
TEST(MatrixTest, MatMatPlusOp)
{
    // create matrices that getting added
    size_t mn = 5;
    
    matrix< double > A(mn, mn, 1), B(mn, mn, 2), C(mn, mn, 3);
    
    // with two matrices
    matrix< double > D = A + B;
    
    // with three matrices
    matrix< double > E = A + B + C;
    
    // with four matrices
    matrix< double > F = A + B + C + D;
    
    ASSERT_EQ(D.cols, mn);
    ASSERT_EQ(D.rows, mn);
    ASSERT_NE(D.mem, nullptr);
    
    ASSERT_EQ(E.cols, mn);
    ASSERT_EQ(E.rows, mn);
    ASSERT_NE(E.mem, nullptr);
    
    ASSERT_EQ(F.cols, mn);
    ASSERT_EQ(F.rows, mn);
    ASSERT_NE(F.mem, nullptr);
    
    // Check values
    for (size_t i = 0; i < mn * mn; ++i)
    {
        ASSERT_EQ(D.mem[i], 3);
        ASSERT_EQ(E.mem[i], 6);
        ASSERT_EQ(F.mem[i], 9);
    }
    
    // check if unsquare shape has same result
    size_t m1 = 1;
    size_t n1 = 5;
    
    matrix< double > AA(m1, n1, 1), BB(m1, n1, 2), CC(m1, n1, 3);
    
    // with two matrices
    matrix< double > DD = AA + BB;
    
    // with three matrices
    matrix< double > EE = AA + BB + CC;
    
    // with four matrices
    matrix< double > FF = AA + BB + CC + DD;
    
    ASSERT_EQ(DD.rows, m1);
    ASSERT_EQ(DD.cols, n1);
    ASSERT_NE(DD.mem, nullptr);
    
    ASSERT_EQ(EE.rows, m1);
    ASSERT_EQ(EE.cols, n1);
    ASSERT_NE(FF.mem, nullptr);
    
    ASSERT_EQ(FF.rows, m1);
    ASSERT_EQ(FF.cols, n1);
    ASSERT_NE(FF.mem, nullptr);
    
    // Check values
    for (size_t i = 0; i < m1 * n1; ++i)
    {
        ASSERT_EQ(DD.mem[i], 3);
        ASSERT_EQ(EE.mem[i], 6);
        ASSERT_EQ(FF.mem[i], 9);
    }
    
    size_t m2 = 5;
    size_t n2 = 1;
    
    matrix< double > AAA(m2, n2, 1), BBB(m2, n2, 2), CCC(m2, n2, 3);
    
    // with two matrices
    matrix< double > DDD = AAA + BBB;
    
    // with three matrices
    matrix< double > EEE = AAA + BBB + CCC;
    
    // with four matrices
    matrix< double > FFF = AAA + BBB + CCC + DDD;
    
    ASSERT_EQ(DDD.rows, m2);
    ASSERT_EQ(DDD.cols, n2);
    ASSERT_NE(DDD.mem, nullptr);
    
    ASSERT_EQ(EEE.rows, m2);
    ASSERT_EQ(EEE.cols, n2);
    ASSERT_NE(FFF.mem, nullptr);
    
    ASSERT_EQ(FFF.rows, m2);
    ASSERT_EQ(FFF.cols, n2);
    ASSERT_NE(FFF.mem, nullptr);
    
    // Check values
    for (size_t i = 0; i < m2 * n2; ++i)
    {
        ASSERT_EQ(DDD.mem[i], 3);
        ASSERT_EQ(EEE.mem[i], 6);
        ASSERT_EQ(FFF.mem[i], 9);
    }
    
    // extreme big matrix
    size_t arb_m = 5000;
    size_t arb_n = 5000;
    
    matrix< double > arb_A(arb_m, arb_n, 7), arb_B(arb_m, arb_n, 11);
    matrix< double > arb_C = arb_A + arb_B;
    
    for (size_t i = 0; i < arb_m * arb_n; ++i)
    {
        ASSERT_EQ(arb_C.mem[i], 18);
    }
}

// Test case for subtraction of two matrix objects
TEST(MatrixTest, MatMatMinusOp)
{
    // create matrices that getting subtracted
    size_t mn = 5;
    
    matrix< double > A(mn, mn, 10), B(mn, mn, 1);
    
    // with two matrices
    matrix< double > C = A - B;
    
    ASSERT_EQ(C.rows, mn);
    ASSERT_EQ(C.cols, mn);
    ASSERT_NE(C.mem, nullptr);
    
    // check values
    for (size_t i = 0; i < mn * mn; ++i)
    {
        ASSERT_EQ(C.mem[i], 9);
    }
    
    // check if unsquare matrix has same result
    size_t m = 1;
    size_t n = 5;
    
    matrix< double > AA(m, n, 10), BB(m, n, 1);
    
    // with two matrices
    matrix< double > CC = AA - BB;
    
    ASSERT_EQ(CC.rows, m);
    ASSERT_EQ(CC.cols, n);
    ASSERT_NE(CC.mem, nullptr);
    
    for (size_t i = 0; i < m * n; ++i)
    {
        ASSERT_EQ(CC.mem[i], 9);
    }
    
    size_t m1 = 5;
    size_t n1 = 1;
    
    matrix< double > AAA(m1, n1, 10), BBB(m1, n1, 1);
    
    // with two matrices
    matrix< double > CCC = AAA - BBB;
    
    ASSERT_EQ(CCC.rows, m1);
    ASSERT_EQ(CCC.cols, n1);
    ASSERT_NE(CCC.mem, nullptr);
    
    for (size_t i = 0; i < m1 * n1; ++i)
    {
        ASSERT_EQ(CCC.mem[i], 9);
    }
    
    // extreme big matrix
    size_t arb_m = 5000;
    size_t arb_n = 5000;
    
    matrix< double > arb_A(arb_m, arb_n, 11), arb_B(arb_m, arb_n, 7);
    matrix< double > arb_C = arb_A - arb_B;
    
    for (size_t i = 0; i < arb_m * arb_n; ++i)
    {
        ASSERT_EQ(arb_C.mem[i], 4);
    }
}

// Test case for the matrix-matrix multiplication
TEST(MatrixTest, MatMatMultOp)
{
    // dimensions
    size_t r = 3;
    size_t c = 4;
    
    // int multiplication
    matrix< double > A(r, c), B(c, r);
    
    // fill with values
    A   << 1  << 2  << 3  << 4
        << 5  << 6  << 7  << 8
        << 9  << 10 << 11 << 12;
    
    B   << 1  << 2  << 3
        << 4  << 5  << 6
        << 7  << 8  << 9
        << 10 << 11 << 12;
    
    // multiply
    matrix< double > C = A * B;
    
    // check dimensions
    ASSERT_EQ(C.rows, r);
    ASSERT_EQ(C.cols, r);
    ASSERT_NE(C.mem, nullptr);
    
    // check values
    ASSERT_EQ(C(0, 0), 70 );
    ASSERT_EQ(C(1, 0), 158);
    ASSERT_EQ(C(2, 0), 246);
    ASSERT_EQ(C(0, 1), 80 );
    ASSERT_EQ(C(1, 1), 184);
    ASSERT_EQ(C(2, 1), 288);
    ASSERT_EQ(C(0, 2), 90 );
    ASSERT_EQ(C(1, 2), 210);
    ASSERT_EQ(C(2, 2), 330);
    
    // transpose both matrices
    A.transpose();
    B.transpose();
    
    // multiply them again
    C = A * B;
    
    // check dimensions
    ASSERT_EQ(C.rows, c);
    ASSERT_EQ(C.cols, c);
    ASSERT_NE(C.mem, nullptr);
    
    // check values
    ASSERT_EQ(C(0, 0), 38 );
    ASSERT_EQ(C(0, 1), 83 );
    ASSERT_EQ(C(0, 2), 128);
    ASSERT_EQ(C(0, 3), 173);
    ASSERT_EQ(C(1, 0), 44 );
    ASSERT_EQ(C(1, 1), 98 );
    ASSERT_EQ(C(1, 2), 152);
    ASSERT_EQ(C(1, 3), 206);
    ASSERT_EQ(C(2, 0), 50 );
    ASSERT_EQ(C(2, 1), 113);
    ASSERT_EQ(C(2, 2), 176);
    ASSERT_EQ(C(2, 3), 239);
    ASSERT_EQ(C(3, 0), 56 );
    ASSERT_EQ(C(3, 1), 128);
    ASSERT_EQ(C(3, 2), 200);
    ASSERT_EQ(C(3, 3), 272);
}