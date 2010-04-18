#include <stdio.h>
#include <stdbool.h>
#include <math.h>

// Function prototypes here
double factorial(double n);
double poisson(double m, double u, bool cuml);
double traffic_intensity(double cpi, double interval, double aht);
double erlangc(double cpi, double interval, double aht, double m);
double svl(double cpi, double interval, double aht, double m, double asa_goal);
double agent_occ(double cpi, double interval, double aht, double m);
int how_many_agents (double cpi, double interval, double aht, double svl_goal, double asa_goal);


// main() proc
int main (int argc, const char * argv[]) {
	
	int cpi, interval, aht, svl_goal, asa_goal;
	
	printf("Welcome to the Erlang C Calculator! \n\n");
	
	printf("Please enter your calls per interval: ");
	scanf("%d", &cpi);
	
	printf("Please enter your interval length (usually 1800): ");
	scanf("%d",&interval);
	
	printf("Please enter your AHT in secs: ");
	scanf("%d",&aht);
	
	printf("Please enter your Service Level goal (%%): ");
	scanf("%d",&svl_goal);
	
	printf("Please enter your Service Level goal (asa): ");
	scanf("%d",&asa_goal);
	
	
	
	
	return 0;
}


double factorial(double n) {
	if (n==0) {
		return 1;
	}
	else {
		return n * factorial(n-1);
	}
}


double poisson(double m, double u, bool cuml) {
	
	double result=0;
	int k;
	
	if (cuml==false) {
		result = ((exp(-u) * powl(u,m)) / factorial(m));
	}
	else {
		for (k=0; k<=m; k++) {
			result += poisson(k,u,false);
		}
	}
	return result;
}


double traffic_intensity(double cpi, double interval, double aht) {
	return (cpi / interval * aht);
}


double erlangc(double cpi, double interval, double aht, double m) {
	double result = 0;
	double u = traffic_intensity(cpi,interval,aht);
	double rho = (u / m); // occupancy
	
	result = poisson(m, u, false) / (poisson(m, u, false)+((1-rho)*poisson(m-1, u, true)));
	
	return result;
}


double svl(double cpi, double interval, double aht, double m, double asa_goal) {
	double result = 0;
	double u = traffic_intensity(cpi, interval, aht);
	
	result = (1.0-(erlangc(cpi, interval, aht, m) * exp((-(m-u)) * (asa_goal / aht))));
	
	return result;
}


double agent_occ(double cpi, double interval, double aht, double m) {
	return ((cpi/interval*aht) / m);
}

int how_many_agents (double cpi, double interval, double aht, double svl_goal, double asa_goal) {
	
	int i = 1;
	bool optimum_svl = false;
	int low_m, optimum_m, high_m;
	
	while (optimum_svl == false) {
		if (svl(cpi,interval,aht,i,asa_goal) >= (svl_goal/100)) {
			optimum_svl = true;
			optimum_m = i;
			low_m = i - 5;
			high_m = i + 5;
		}
		i++;
	}
	
	printf("For %g NCO with AHT of %g secs (in a %g sec interval), the following table shows approx staff required to meet service level of %g%% in %g secs:\n\n",
		   cpi, aht, interval, svl_goal, asa_goal);
	
	for (i = low_m; i <= high_m; i++) {
		printf("%i agents:\tOCC = %g%%,\tSVL = %g%%\n", i, agent_occ(cpi, interval, aht, i)*100, svl(cpi,interval,aht,i,asa_goal)*100);
	}

	
	return optimum_m;
}


// ?? normdist(x) = phi(x) = (1/sqrt(2*M_PI)) * exp(-1/2 * pow (x,2))
// ?? Y = [ 1/σ * sqrt(2*π) ] * exp(pow(-(x - μ),2)/2*pow(σ,2))

// b0 = 0.2316419, b1 = 0.319381530, b2  = −0.356563782, b3 = 1.781477937, b4  = −1.821255978, b5 = 1.330274429.
// t = 1 / (1+(b0 * x)
// cdf(x) ~= (1-phi(x))((b1*t)+(b2*pow(t,2))+(b3*pow(t,3))+(b4*pow(t,4))+(b5*pow(t,5))) + epsilon(x)
// epsilon(x)??

