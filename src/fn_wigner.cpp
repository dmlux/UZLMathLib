//
//  fn_wigner.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 17.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_wigner_cpp
#define uzlmath_fn_wigner_cpp

#include <uzlmath>

namespace uzlmath
{
    namespace wigner
    {
        
        /*!
         * @brief       Calculates the Wigner d-function of degree \f$J\f$ and orders \f$M\f$ and \f$M'\f$
         *              which can be found in literature as \f$d^J_{M,M'}(\beta)\f$. The implementation of
         *              this method is defined recursively
         * @details     There are four parameters that are used to call this function. The degree \f$J\f$
         *              and the orders \f$M\f$ and \f$M'\f$ and the angle \f$beta\f$. The legal ranges that
         *              can be choosen are
         *              \f{eqnarray*}{
         *                  -J &\leq& M\\
         *                  M' &\leq& J
         *              \f}
         *              And \f$0 \leq \beta \leq \pi\f$. The implementation itself is based on the three-term
         *              recurrence that can be expressed as
         *              \f{eqnarray*}{
         *                  d^J_{M,M'}(\beta)\ = -\frac{J(2J - 1)\left(\frac{M,M'}{J(J - 1)} - \cos(\beta)\right)}
         *                      {\sqrt{(J^2 - M^2)(J^2 - M'^2)}}
         *                      d^{J-1}_{M,M'}(\beta) - \frac{J\sqrt{((J - 1)^2 - M^2)((J - 1)^2 - M'^2)}}
         *                      {(J - 1)\sqrt{(J^2 - M^2)(J^2 - M'^1)}}
         *                      d^{J-2}_{M,M'}(\beta)
         *              \f}
         *              The base cases of this recursive function are
         *              \f{eqnarray*}{
         *                  d^J_{00}(\beta) &=& P_J(\cos(\beta))\\
         *                  d^J_{JM}(\beta) &=& \sqrt{\frac{(2J)!}{(J + M)!(J - M)!}}
         *                      \left(\cos\frac{\beta}{2}\right)^{J+M}
         *                      \left(-\sin\frac{\beta}{2}\right)^{J-M}\\
         *                  d^J_{-JM}(\beta) &=& \sqrt{\frac{(2J)!}{(J + M)!(J - M)!}}
         *                      \left(\cos\frac{\beta}{2}\right)^{J-M}
         *                      \left(\sin\frac{\beta}{2}\right)^{J+M}\\
         *                  d^J_{MJ}(\beta) &=& \sqrt{\frac{(2J)!}{(J + M)!(J - M)!}}
         *                      \left(\cos\frac{\beta}{2}\right)^{J+M}
         *                      \left(\sin\frac{\beta}{2}\right)^{J-M}\\
         *                  d^J_{M-J}(\beta) &=& \sqrt{\frac{(2J)!}{(J + M)!(J - M)!}}
         *                      \left(\cos\frac{\beta}{2}\right)^{J-M}
         *                      \left(-\sin\frac{\beta}{2}\right)^{J+M}
         *              \f}
         *              Where \f$P_n(x)\f$ denotes the Legendre polynomial which is implemented in the
         *              polynomials namespace as the orthopoly::legendre function. The implementation
         *              uses a dynamic programming approach to evaluate the Wigner d-function in \f$\mathcal{O}(J)\f$
         *              time.
         *
         * @param[in]   J The degree \f$J\f$ of the Wigner d-function that is supposed to be calculated
         * @param[in]   M The order \f$M\f$ of the Wigner d-function that is supposed to be calculated
         * @param[in]   Mp The order \f$M'\f$ of the Wigner d-function that is supposed to be calculated
         * @param[in]   beta The angle \f$\beta\f$ of the Wigner d-function that is supposed to be calculated.
         * @note        By comparing the results of this implementation a mismatch to  Mathematica
         *              was found that was caused by the angle conventions that are used. To get the same
         *              results than Mathematica you should use \f$-\beta\f$ instead of \f$\beta\f$!
         * @return      The resulting function value for the given parameters which is in \f$\mathbb{R}\f$
         *
         * @since       0.0.1
         *
         * @see         orthopoly::legendre
         *
         * @author      Denis-Michael Lux <denis.lux@icloud.com>
         * @date        12.01.15
         */
        auto wigner_d(const int& J, const int& M, const int& Mp, const double& beta) -> const double
        {
            // undefined values
            if (M == 0 && Mp == 0)
            {
                return orthoPoly::legendre(J, cos(beta));
            }
            
            // memorization
            double mem[J + 1];
            
            // get all values of d^J_{M,M'}
            for (int i = 0; i <= J; ++i)
            {
                if ( !(i >= std::max(abs(M), abs(Mp))) )
                {
                    mem[i] = 0;
                }
                // d^{J}_{JM}(beta)
                else if ( i == M )
                {
                    double root = sqrt( (factorial(2.0 * i) / (factorial(i + Mp) * factorial(i - Mp))).eval() );
                    mem[i]      = root * pow(cos(beta / 2.0), i + Mp) * pow(-sin(beta / 2.0), i - Mp);
                }
                // d^{J}_{-JM}(beta)
                else if ( -M == i )
                {
                    double root = sqrt( (factorial(2.0 * i) / (factorial(i + Mp) * factorial(i - Mp))).eval() );
                    mem[i]      = root * pow(cos(beta / 2.0), i - Mp) * pow(sin(beta / 2.0), i + Mp);
                }
                // d^{J}_{MJ}(beta)
                else if ( i == Mp )
                {
                    double root = sqrt( (factorial(2.0 * i) / (factorial(i + M) * factorial(i - M))).eval() );
                    mem[i]      = root * pow(cos(beta / 2.0), i + M) * pow(sin(beta / 2.0), i - M);
                }
                // d^{J}_{M-J}(beta)
                else if ( -Mp == i )
                {
                    double root = sqrt( (factorial(2.0 * i) / (factorial(i + M) * factorial(i - M))).eval() );
                    mem[i]      = root * pow(cos(beta / 2.0), i - M) * pow(-sin(beta / 2.0), i + M);
                }
                else
                {
                    // coefficient from d^{J - 1}_{M,M'}
                    double coef1 = -1.0;
                    coef1       *= i * (2.0 * i - 1.0) * ((M * Mp) / (i * (i - 1.0)) - cos(beta));
                    coef1       /= sqrt((i * i - M * M) * (i * i - Mp * Mp));
                    
                    // coefficient from d^{J - 2}_{M,M'}
                    double coef2 = -1.0;
                    coef2       *= sqrt(((i - 1.0) * (i - 1.0) - M * M) * ((i - 1.0) * (i - 1.0) - Mp * Mp)) * i;
                    coef2       /= ((i - 1.0) * sqrt((i * i - M * M) * (i * i - Mp * Mp)));
                    
                    mem[i]       = coef1 * mem[i-1] + coef2 * mem[i-2];
                }
            }
            
            return mem[J];
        }
        
        
        /*!
         * @brief       Calculates the \f$L^2\f$ normalized Wigner d-function.
         * @details     The legal ranges of parameters are the same as for the wigner::wigner_d function
         *              \f{eqnarray*}{
         *                  -J &\leq& M\\
         *                  M' &\leq& J
         *              \f}
         *              And \f$0 \leq \beta \leq \pi\f$. The \f$L^2\f$ normalized Wigner d-function can
         *              be expressed as
         *              \f{eqnarray*}{
         *                  \tilde{d}^J_{M,M'}(\beta) = \sqrt{\frac{2J + 1}{2}}d^J_{M,M'}(\beta)
         *              \f}
         *
         * @param[in]   J The degree \f$J\f$ of the \f$L^2\f$-normalized Wigner d-function that is supposed
         *              to be calculated
         * @param[in]   M The order \f$M\f$ of the \f$L^2\f$-normalized Wigner d-function that is supposed
         *              to be calculated
         * @param[in]   Mp The order \f$M'\f$ of the \f$L^2\f$-normalized Wigner d-function that is supposed
         *              to be calculated
         * @param[in]   beta The angle \f$\beta\f$ of the \f$L^2\f$-normalized Wigner d-function that is
         *              supposed to be calculated.
         * @note        By comparing the results of this implementation a mismatch to  Mathematica
         *              was found that was caused by the angle conventions that are used. To get the same
         *              results than Mathematica you should use \f$-\beta\f$ instead of \f$\beta\f$!
         * @return      The resulting function value for the given parameters which is in \f$\mathbb{R}\f$
         *
         * @see         wigner::wigner_d
         *
         * @since       0.0.1
         *
         * @author      Denis-Michael Lux <denis.lux@icloud.com>
         * @date        12.01.15
         */
        auto wigner_d_l2normalized(const int& J, const int& M, const int& Mp, const double& beta) -> const double
        {
            if ( !(-J <= abs(M) && abs(Mp) <= J) )
            {
                printf("** uzlmath error: Illegal arguments for normalized Wigner d-function. Legal arguments are -J <= |M|, |M'| <= J. **");
                exit(EXIT_FAILURE);
            }
            return sqrt((2.0 * J + 1.0) / 2.0) * wigner_d(J, M, Mp, beta);
        }
                
        /*!
         * @brief       Calculates the Wigner D-function which is the result of using the Euler angle
         *              decomposition.
         * @details     The Wigner D-function is based on the wigner::wigner_d function with the main
         *              difference that this function maps to \f$\mathbb{C}\f$ instead to \f$\mathbb{R}\f$.
         *              The legal ranges are the same than the ranges of \f$d^J_{M,M'}(\beta)\f$
         *              \f{eqnarray*}{
         *                  -J &\leq& M\\
         *                  M' &\leq& J
         *              \f}
         *              And \f$0 \leq \beta \leq \pi\f$. The expression that was used to implement the
         *              Wigner D-function is
         *              \f{eqnarray*}{
         *                  D^J_{M,M'}(\beta) = e^{-iM\alpha}d^J_{M,M'}(\beta)e^{-iM'\gamma}
         *              \f}
         *
         * @param[in]   J The degree \f$J\f$ of the Wigner D-function that is supposed to be calculated
         * @param[in]   M The order \f$M\f$ of the Wigner D-function that is supposed to be calculated
         * @param[in]   Mp The order \f$M'\f$ of the Wigner D-function that is supposed to be calculated
         * @param[in]   alpha The angle \f$\alpha\f$ of the Wigner d-function that is supposed to be calculated.
         * @param[in]   beta The angle \f$\beta\f$ of the Wigner d-function that is supposed to be calculated.
         * @param[in]   gamma The angle \f$\gamma\f$ of the Wigner d-function that is supposed to be calculated.
         * @note        By comparing the results of this implementation a mismatch to  Mathematica
         *              was found that was caused by the angle conventions that are used. To get the same
         *              results than Mathematica you should use \f$-\alpha\f$ instead of \f$\alpha\f$,
         *              \f$-\beta\f$ instead of \f$\beta\f$ and \f$-\gamma\f$ instead of \f$\gamma\f$!
         * @return      The resulting function value for the given parameters which is in \f$\mathbb{C}\f$
         *
         * @see         wigner::wigner_d
         *
         * @since       0.0.1
         *
         * @author      Denis-Michael Lux <denis.lux@icloud.com>
         * @date        12.01.15
         */
        auto wigner_D(const int& J, const int& M, const int& Mp, const double& alpha, const double& beta, const double& gamma) -> const complex< double >
        {
            if ( !(-J <= abs(M) && abs(Mp) <= J) )
            {
                printf("** uzlmath error: Illegal arguments for Wigner D-function. Legal arguments are -J <= |M|, |M'| <= J. **");
                exit(EXIT_FAILURE);
            }
            
            complex< double > a, b;
            
            a.polar(1, -M * alpha);
            b.polar(1, -Mp * gamma);
            
            return a * wigner_d(J, M, Mp, beta) * b;
        }
                
        /*!
         * @brief       Calculates the \f$L^2\f$-normalized Wigner D-function which is the result of using
         *              the Euler angle decomposition.
         * @details     The \f$L^2\f$-normalized Wigner D-function is based on the wigner::wigner_d_l2normalized
         *              function with the main difference that this function maps to \f$\mathbb{C}\f$ instead to
         *              \f$\mathbb{R}\f$. The legal ranges are the same than the ranges of \f$d^J_{M,M'}(\beta)\f$
         *              \f{eqnarray*}{
         *                  -J &\leq& M\\
         *                  M' &\leq& J
         *              \f}
         *              And \f$0 \leq \beta \leq \pi\f$. The expression that was used to implement the
         *              Wigner D-function is
         *              \f{eqnarray*}{
         *                  \tilde{D}^J_{M,M'}(\beta) = \frac{1}{2\pi}\sqrt{\frac{2J + 1}{2}}D^J_{M,M'}(\alpha,\beta,\gamma)
         *              \f}
         *
         * @param[in]   J The degree \f$J\f$ of the \f$L^2\f$-normalized Wigner D-function that is supposed to be calculated
         * @param[in]   M The order \f$M\f$ of the \f$L^2\f$-normalized Wigner D-function that is supposed to be calculated
         * @param[in]   Mp The order \f$M'\f$ of the \f$L^2\f$-normalized Wigner D-function that is supposed to be calculated
         * @param[in]   alpha The angle \f$\alpha\f$ of the Wigner d-function that is supposed to be calculated.
         * @param[in]   beta The angle \f$\beta\f$ of the Wigner d-function that is supposed to be calculated.
         * @param[in]   gamma The angle \f$\gamma\f$ of the Wigner d-function that is supposed to be calculated.
         * @note        By comparing the results of this implementation a mismatch to  Mathematica
         *              was found that was caused by the angle conventions that are used. To get the same
         *              results than Mathematica you should use \f$-\alpha\f$ instead of \f$\alpha\f$,
         *              \f$-\beta\f$ instead of \f$\beta\f$ and \f$-\gamma\f$ instead of \f$\gamma\f$!
         * @return      The resulting function value for the given parameters which is in \f$\mathbb{C}\f$
         *
         * @see         wigner::wigner_d_l2normalized
         *
         * @since       0.0.1
         *
         * @author      Denis-Michael Lux <denis.lux@icloud.com>
         * @date        15.05.15
         */
        auto wigner_D_l2normalized(const int& J, const int& M, const int& Mp, const double& alpha, const double& beta, const double& gamma) -> const complex< double >
        {
            if ( !(-J <= abs(M) && abs(Mp) <= J) )
            {
                printf("** uzlmath error: Illegal arguments for L^2 normalized Wigner D-function. Legal arguments are -J <= |M|, |M'| <= J. **");
                exit(EXIT_FAILURE);
            }
            
            return 1.0 / (2.0 * M_PI) * sqrt((2.0 * J + 1.0) / 2.0) * wigner_D(J, M, Mp, alpha, beta, gamma);
        }
    }
}

#endif