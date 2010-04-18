#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// Function prototypes here
double factorial(double n);
double poisson(double m, double u, bool cuml);
double traffic_intensity(double cpi, double interval, double aht);
double erlangc(double cpi, double interval, double aht, double m);
double svl(double cpi, double interval, double aht, double m, double asa_goal);
double agent_occ(double cpi, double interval, double aht, double m);
double asa(double cpi, double interval, double aht, double m);
int how_many_agents (double cpi, double interval, double aht, double svl_goal, double asa_goal);


// main() proc
int main (int argc, const char * argv[]) {
	
	char cpi[10], interval[10], aht[10], svl_goal[10], asa_goal[10];
	int		i_cpi, i_interval, i_aht, i_svl_goal, i_asa_goal;
	printf("==================================================================================\n");
	printf("Welcome to the Erlang C Calculator!\nCopyright (c) 2010 by Gethryn Ghavalas\n");
	printf("==================================================================================\n\n\n");
	
	printf("Please supply the following:\n\n");
	printf("Calls per interval (or Number of Calls Offered): ");
	fgets(cpi, 10, stdin);
	//cpi[strlen(cpi)-1] = 0;
	//scanf("%d", &cpi);
	
	printf("Interval length in secs [1800]: ");
	fgets(interval, 10, stdin);
	//scanf("%d",&interval);
	
	printf("AHT in secs: ");
	fgets(aht, 10, stdin);
	//scanf("%d",&aht);
	
	printf("Service Level goal (%%) [80]: ");
	fgets(svl_goal, 10, stdin);
	//scanf("%d",&svl_goal);
	
	printf("Service Level goal (asa) [20]: ");
	fgets(asa_goal, 10, stdin);
	//scanf("%d",&asa_goal);
	
	i_cpi = atoi(cpi);
	i_interval = atoi(interval);
	if (i_interval==0) {
		i_interval = 1800;
	}
	i_aht = atoi(aht);
	i_svl_goal = atoi(svl_goal);
	if (i_svl_goal==0) {
		i_svl_goal = 80;
	}
	i_asa_goal = atoi(asa_goal);
	if (i_asa_goal==0) {
		i_asa_goal = 20;
	}
	
	printf("\n\n");
	
	printf("\n\nYour optimum number of agents is %d.\n\n",
		   how_many_agents(i_cpi, i_interval, i_aht, i_svl_goal, i_asa_goal));
	
	
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


double asa(double cpi, double interval, double aht, double m) {
	double result = 0;
	double u = traffic_intensity(cpi, interval, aht);
	double rho = (u / m); // occupancy
	
	result = (erlangc(cpi, interval, aht, m) * aht) / (m * (1-rho));
	
	return result;
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
			high_m = i + 10;
		}
		i++;
	}
	
	printf("The following table shows the optimum number of agents (-5 <= opt <= +10)\n");
	printf("for %g NCO with AHT of %g secs (in a %g sec interval)\nand a service level of %g%% in %g secs:\n\n",
		   cpi, aht, interval, svl_goal, asa_goal);
	
	for (i = low_m; i <= high_m; i++) {
		printf("%03i agents:\t\tOCC = %5.1f%%\t\tSVL = %5.1f%%\t\tASA = %7.1f\n",
			   i, agent_occ(cpi, interval, aht, i)*100, svl(cpi,interval,aht,i,asa_goal)*100, asa(cpi,interval,aht,i));
	}

	
	return optimum_m;
}


// ?? normdist(x) = phi(x) = (1/sqrt(2*M_PI)) * exp(-1/2 * pow (x,2))
// ?? Y = [ 1/σ * sqrt(2*π) ] * exp(pow(-(x - μ),2)/2*pow(σ,2))

// b0 = 0.2316419, b1 = 0.319381530, b2  = −0.356563782, b3 = 1.781477937, b4  = −1.821255978, b5 = 1.330274429.
// t = 1 / (1+(b0 * x)
// cdf(x) ~= (1-phi(x))((b1*t)+(b2*pow(t,2))+(b3*pow(t,3))+(b4*pow(t,4))+(b5*pow(t,5))) + epsilon(x)
// epsilon(x)??

