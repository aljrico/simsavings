#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<cmath>
#include<time.h>
#include<stdbool.h>
#include<string.h>
#include<Rcpp.h>
#include<algorithm>

using namespace Rcpp;

// Sum function
double sum(NumericVector x) {
	int n = x.size();
	double total = 0;
	for(int i = 0; i < n; ++i) {
		total += x[i];
	}
	return total;
}

// Computing of 'pi' function
double fpi(double A, double K, double x, NumericVector f, int time, int years)
{
	NumericVector c(years);
	int i;
	float g, xpi;

	for(i=0;i<years;i++){c[i]=0;}
	for(i=time;i<years;i++){c[i]=f[i];}

	g = sum(c);

	xpi = A * (K + x + g);
	return xpi;
}



// [[Rcpp::export]]
double normal(double mu, double sigma)
	{
		double U1, U2, W, mult;
		static double X1, X2;
		static int call = 0;

		if (call == 1)
		{
			call = !call;
			return (mu + sigma * (double) X2);
		}

		do
		{
			U1 = -1 + ((double) rand () / RAND_MAX) * 2;
			U2 = -1 + ((double) rand () / RAND_MAX) * 2;
			W = pow (U1, 2) + pow (U2, 2);
		}
		while (W >= 1 || W == 0);

		mult = sqrt ((-2 * log (W)) / W);
		X1 = U1 * mult;
		X2 = U2 * mult;

		call = !call;

		return (mu + sigma * (double) X1);
	}

// [[Rcpp::export]]
NumericVector cppi_c(int nsim, float alpha, float sigma, float a, int years, float pi){

	NumericVector final_wealth(nsim);
	int i,j;
	float x_next,x_curr,rndm;
	NumericVector f(years);

	for(i=0; i<years/2; i++){f[i] = a;}
	for(i=years/2; i<years; i++){f[i] = -a;}

	for(j=0;j<nsim;j++)
	{
		x_curr=a;

		for(i=0;i<years-1;i++)
		{
			rndm = normal(alpha,sigma);
			x_next = x_curr * (1 + rndm) * pi + (1-pi) * x_curr + f[i+1];
			x_curr = x_next;
		}

		final_wealth[j] = x_next;
	}
	return final_wealth;
}

// [[Rcpp::export]]
NumericVector alt_c(float alpha, float sigma, float a, int years, int nsim, float K, float A_factor){

	NumericVector final_wealth(nsim);
	int i,j,time;
	float x_next,x_curr,rndm;
	NumericVector f(years);
	float pi;

	for(i=0; i<years/2; i++){f[i] = a;}
	for(i=years/2; i<years; i++){f[i] = -a;}

	for(j=0;j<nsim;j++)
	{
		x_curr=a;

		for(i=0;i<years-1;i++)
		{
			rndm = normal(alpha,sigma);
			time = i;
			pi = fpi(A_factor, K, x_curr, f, time, years)/x_curr;
			x_next = x_curr * (1 + rndm) * pi + (1-pi) * x_curr + f[i+1];
			x_curr = x_next;
		}

		final_wealth[j] = x_next;
	}
	return final_wealth;
}

// [[Rcpp::export]]
NumericVector alt_mort_c(float alpha, float sigma, float a, int years, int nsim, float K, float A_factor, NumericVector M){

	NumericVector final_wealth(nsim);
	int i,j,time;
	float x_next,x_curr,rndm;
	NumericVector f(years);
	float pi;

	for(i=0; i<years/2; i++){f[i] = a;}
	for(i=years/2; i<years; i++){f[i] = -a;}

	for(j=0;j<nsim;j++)
	{
		x_curr=a;

		for(i=0;i<years-1;i++)
		{
			rndm = normal(alpha,sigma);
			time = i;
			pi = fpi(A_factor, K, x_curr, f, time, years)/x_curr;
			x_next = x_curr * (1 + rndm) * pi + (1-pi) * x_curr + f[i+1] + x_curr*M[i];
			x_curr = x_next;
		}

		final_wealth[j] = x_next;
	}
	return final_wealth;
}

// [[Rcpp::export]]
NumericVector cppi_mort_c(int nsim, float alpha, float sigma, float a, int years, float pi, NumericVector M){

	NumericVector final_wealth(nsim);
	int i,j;
	float x_next,x_curr,rndm;
	NumericVector f(years);

	for(i=0; i<years/2; i++){f[i] = a;}
	for(i=years/2; i<years; i++){f[i] = -a;}

	for(j=0;j<nsim;j++)
	{
		x_curr=a;

		for(i=0;i<years-1;i++)
		{
			rndm = normal(alpha,sigma);
			x_next = x_curr * (1 + rndm) * pi + (1-pi) * x_curr + f[i+1] + x_curr*M[i];
			x_curr = x_next;
		}

		final_wealth[j] = x_next;
	}
	return final_wealth;
}

// Sort vector
// [[Rcpp::export]]
NumericVector sort_c(NumericVector vec)
{
	std::sort(vec.begin(), vec.end());
	return(vec);
}

// Length C
// [[Rcpp::export]]
double length_c(NumericVector vec)
{
	int size = vec.size();
	return size;
}

// Compute Median
// [[Rcpp::export]]
double median_c(NumericVector vec)
{
	int size = length_c(vec);
	vec = sort_c(vec);
	int m, m_2;
	double median;

	m = size/2;
	m_2 = size/2 - 1;

	median = vec[size/2];

	return median;
}


// Computing of Equivalent Pi'
// [[Rcpp::export]]
double equiv_pi_c(int m, double ret, int nsim, float alpha, float sigma, float a, int years, float pi)
{
	int i;
	float seed_pi;
	float max_pi;
	float min_pi;
	float est_ret;
	float cppi_median;
	int n;

	n = m;

	NumericVector cppi_res(nsim);

	seed_pi = 0.5;
	max_pi = 1;
	min_pi = 0;

	for(i=0;i<n;i++){
		cppi_res = cppi_c(alpha = alpha, sigma = sigma, a = a, years = years, nsim = nsim, pi = pi);
		cppi_median = median_c(cppi_res);
		est_ret = (1/years)*(-1 + std::sqrt(1 + (8*(cppi_median))/(a*years))*100);

		if(est_ret < ret){
			min_pi = pi;
			pi = (pi + max_pi)/2;
		}else{
			max_pi = pi;
			pi = (pi + min_pi)/2;
		}

		if(est_ret == ret){break;}
		float upper = abs(est_ret-ret);
		float bot = abs(ret);
		if( (100*upper) / bot < 1){break;}
	}
	return pi;
}

// Computing of Equivalent Pi'
// [[Rcpp::export]]
float return_c(NumericVector vec, float a, int years)
{
	float median;
	float est_ret;

	median = median_c(vec);
	est_ret = (1/years)*(-1 + sqrt(1 + (8*(median))/(a*years)))*100;
	return est_ret;
}
