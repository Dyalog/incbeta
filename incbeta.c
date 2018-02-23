/*
 * zlib License
 *
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */


#include <math.h>
#include <stddef.h>

#define STOP 1.0e-8
#define TINY 1.0e-30
#define ERRVAL (HUGE_VAL * 0.0) // NaN

#if _MSC_VER
# define DLLEXPORT __declspec(dllexport)
#else
# define DLLEXPORT
#endif

DLLEXPORT void incbeta(double a, double b, double *p, size_t n)
{
	if (a <= 0 || b <= 0)
	{
		for (size_t i = 0; i < n; i++) p[i] = ERRVAL;
		return;
	}

	/*Find the first part before the continued fraction.*/
	const double lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);

	for (size_t i = 0; i < n; i++)
	{
		double x = p[i];
		/*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
		int swap = x > (a + 1) / (a + b + 2);
		if (swap) { double t = a; a = b; b = t; x = 1 - x; } /*Use the fact that beta is symmetrical.*/
		if (x <= 0)
			x = 0;
		else
		{
			const double front = exp(log(x) * a + log(1 - x) * b - lbeta_ab) / a;

			/*Use Lentz's algorithm to evaluate the continued fraction.*/
			double f = 1, c = 1, d = 0;
			for (int j = 0; ; j++)
			{
				int m = j / 2;

				double numerator;
				if (j == 0)
					numerator = 1; /*First numerator is 1.0.*/
				else if (j % 2 == 0)
					numerator = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m)); /*Even term.*/
				else
					numerator = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1));  /*Odd term.*/

				/*Do an iteration of Lentz's algorithm.*/
				d = 1 + numerator * d;
				if (fabs(d) < TINY) d = TINY;
				d = 1 / d;

				c = 1 + numerator / c;
				if (fabs(c) < TINY) c = TINY;

				const double cd = c * d;
				f *= cd;

				/*Check for stop.*/
				if (fabs(1 - cd) < STOP)
				{
					x = front * (f - 1);
					break;
				}

				if (j >= 200)
				{
					x = ERRVAL; /*Needed more loops, did not converge.*/
					break;
				}
			}
		}
		if (swap) { double t = a; a = b; b = t; x = 1 - x; }
		p[i] = x;
	}
}
