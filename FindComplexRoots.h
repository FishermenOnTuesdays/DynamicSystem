#pragma once
#include<complex>
#include<vector>
#include<functional>
#include<omp.h>

namespace FindCR
{
	#define pi long double(4*atan(1))

	typedef std::complex<long double> complex;

	long double AngleBetween(complex z1, complex z2)
	{
		long double angle = fabsl(arg(z1) - arg(z2));
		return angle > pi ? 2 * pi - angle : angle;
	}

	int sign(long double x)
	{
		return x == 0 ? 0 : x < 0 ? -1 : 1;
	}

	void FindRoots(long double radius_square, long double epsilon, long double x0, long double y0, std::vector<complex>& roots, std::function<complex(complex)> f, long double d = 22, long double d_a = 0.092)
	{
		complex z_prev{ 0,0 }, z_curr;
		long double sum = 0, d_alpha = d_a;
		radius_square += radius_square / d;
		for (long double alpha = 0; alpha < 2 * pi; alpha += d_alpha)
		{
			long double real = x0 + radius_square * cosl(alpha) / std::fmaxl(fabsl(sinl(alpha)), fabsl(cosl(alpha)));
			long double imag = y0 + radius_square * sinl(alpha) / std::fmaxl(fabsl(sinl(alpha)), fabsl(cosl(alpha)));
			z_curr = complex(real, imag);
			z_curr = f(z_curr);
			long double angle = AngleBetween(z_curr, z_prev);
			int sign_sum = sign(std::real(z_curr) * std::imag(z_prev) - std::imag(z_curr) * std::real(z_prev));
			sum += sign_sum * angle;
			z_prev = complex(std::real(z_curr), std::imag(z_curr));
		}

		if (fabsl(sum) / 2 / pi >= 1 - d_alpha)
		{
			if (radius_square < epsilon)
			{
				bool is_find = false;
				for (auto el : roots)
				{
					if (fabsl(real(el) - x0) <= 2 * epsilon && fabsl(imag(el) - y0) <= 2 * epsilon)
					{
						is_find = true;
						break;
					}
				}
				if (!is_find)
					roots.push_back(complex(round(x0 / epsilon) * epsilon, round(y0 / epsilon) * epsilon));
			}
			else
			{
				FindRoots(radius_square / 2, epsilon, x0 + radius_square / 2, y0 - radius_square / 2, roots, f, d, d_a); 
				FindRoots(radius_square / 2, epsilon, x0 - radius_square / 2, y0 + radius_square / 2, roots, f, d, d_a); 
				FindRoots(radius_square / 2, epsilon, x0 + radius_square / 2, y0 + radius_square / 2, roots, f, d, d_a); 
				FindRoots(radius_square / 2, epsilon, x0 - radius_square / 2, y0 - radius_square / 2, roots, f, d, d_a); 
				
			}
		}
	}

}