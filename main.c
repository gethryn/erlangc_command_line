#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// Function prototypes
double factorial(double n);
double poisson(double m, double u, bool cuml);
double traffic_intensity(double cpi, double interval, double aht);
double erlangc(double cpi, double interval, double aht, double m);
double svl(double cpi, double interval, double aht, double m, double asa_goal);
double agent_occ(double cpi, double interval, double aht, double m);
double asa(double cpi, double interval, double aht, double m);
int how_many_agents (double cpi, double interval, double aht, double svl_goal, double asa_goal);

// Queries the user for inputs to calculate optimum staffing levels, then presents a table of
// the options to the user from optimum-5 to optimum+10, with corresponding occupancy, svl &
// average speed of answer.
int main (int argc, const char * argv[]) {
	
	char cpi[10], interval[10], aht[10], svl_goal[10], asa_goal[10];
	int		i_cpi, i_interval, i_aht, i_svl_goal, i_asa_goal;
	printf("==================================================================================\n");
	printf("Welcome to the Erlang C Calculator!\nCopyright (c) 2010 by Gethryn Ghavalas\n");
	printf("==================================================================================\n\n\n");
	
	printf("Please supply the following:\n\n");
	printf("Calls per interval (or Number of Calls Offered): ");
	fgets(cpi, 10, stdin);
	
	printf("Interval length in secs [1800]: ");
	fgets(interval, 10, stdin);
	
	printf("AHT in secs: ");
	fgets(aht, 10, stdin);
	
	printf("Service Level goal (%%) [80]: ");
	fgets(svl_goal, 10, stdin);
	
	printf("Service Level goal (asa) [20]: ");
	fgets(asa_goal, 10, stdin);
	
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

// calculates the value of n!
double factorial(double n) {
	if (n==0) {
		return 1;
	}
	else {
		return n * factorial(n-1);
	}
}

// calculates the probability that m will equal u using the Poisson Distribution
// use cuml = false for standard, and cuml = true for cumulative result.
double poisson(double m, double u, bool cuml) {
	
	double result=0;
	int k;
	
	if (cuml==false) {
		// standard result
		result = ((exp(-u) * powl(u,m)) / factorial(m));
	}
	else {
		// cumulative result
		for (k=0; k<=m; k++) {
			result += poisson(k,u,false);
		}
	}
	return result;
}

// calculates traffic intensity using events per interval (cpi), the length of the interval
// in seconds -- usually 30 mins = 1800 (interval), and the length of each event -- average
// handling time (aht).
double traffic_intensity(double cpi, double interval, double aht) {
	//
	return (cpi / interval * aht);
}

// calculates likelihood an event will experience delay before handling given a certain number
// of events per interval (cpi), the length of the interval in seconds -- usually 30 mins = 1800
// (interval), the length of each event -- average handling time (aht), and the number of agents
// processing the event (m).  Expressed as a % probability.
double erlangc(double cpi, double interval, double aht, double m) {
	double result = 0;
	double u = traffic_intensity(cpi,interval,aht);
	double rho = (u / m); // occupancy
	
	result = poisson(m, u, false) / (poisson(m, u, false)+((1-rho)*poisson(m-1, u, true)));
	
	return result;
}

// calculates the probability a call/event will be answered/processed within a goal answer time
// takes same arguments as erlang c function, and adds the goal answer time (asa_goal) in secs.
double svl(double cpi, double interval, double aht, double m, double asa_goal) {
	double result = 0;
	double u = traffic_intensity(cpi, interval, aht);
	
	result = (1.0-(erlangc(cpi, interval, aht, m) * exp((-(m-u)) * (asa_goal / aht))));
	
	return result;
}

// calculates how busy each agent is -- % of time occupied, given a certain number of calls/events
// (cpi) within an interval (interval) of average length (aht), and with m agents working.
// sometimes referred to as rho.
double agent_occ(double cpi, double interval, double aht, double m) {
	return ((cpi/interval*aht) / m);
}

// calculates the average speed of answer in secs, given a certain number of calls/events
// (cpi) within an interval (interval) of average length (aht), and with m agents working.
double asa(double cpi, double interval, double aht, double m) {
	double result = 0;
	double u = traffic_intensity(cpi, interval, aht);
	double rho = (u / m); // occupancy
	
	result = (erlangc(cpi, interval, aht, m) * aht) / (m * (1-rho));
	
	return result;
}

// this funtion provides a table of output giving from optimum-5 to optimum+10 agents and the
// corresponding occupancy, service level and average speed of answer for each of these numbers
// of agents.  The function returns the optiumum number of agents for future use.
int how_many_agents (double cpi, double interval, double aht, double svl_goal, double asa_goal) {
	
	int i = 1;
	bool optimum_svl = false;
	int low_m, optimum_m, high_m;
	
	// starts at 1 agent, and repeats until the number of agents (i) gives a svl that exceeds
	// or equals that requested by the user.
	while (optimum_svl == false) {
		if (svl(cpi,interval,aht,i,asa_goal) >= (svl_goal/100)) {
			optimum_svl = true;
			optimum_m = i;
			low_m = i - 5;
			high_m = i + 10;
		}
		i++;
	}
	
	printf("The following table shows the optimum number of agents (before shrinkage)\n");
	printf("for %g NCO with AHT of %g secs (in a %g sec interval)\nand a service level of %g%% in %g secs:\n\n",
		   cpi, aht, interval, svl_goal, asa_goal);
	
	// prints the table of options for the user.
	for (i = low_m; i <= high_m; i++) {
		if (i == optimum_m) {
			printf("\n");
		}
		printf("\n%03i agents:\t\tOCC = %5.1f%%\t\tSVL = %5.1f%%\t\tASA = %7.1f",
			   i, agent_occ(cpi, interval, aht, i)*100, svl(cpi,interval,aht,i,asa_goal)*100, asa(cpi,interval,aht,i));
		if (i == optimum_m) {
			printf("   OPTIMUM\n");
		}
	}

	
	return optimum_m;
}


// ?? normdist(x) = phi(x) = (1/sqrt(2*M_PI)) * exp(-1/2 * pow (x,2))
// ?? Y = [ 1/σ * sqrt(2*π) ] * exp(pow(-(x - μ),2)/2*pow(σ,2))

// b0 = 0.2316419, b1 = 0.319381530, b2  = −0.356563782, b3 = 1.781477937, b4  = −1.821255978, b5 = 1.330274429.
// t = 1 / (1+(b0 * x)
// cdf(x) ~= (1-phi(x))((b1*t)+(b2*pow(t,2))+(b3*pow(t,3))+(b4*pow(t,4))+(b5*pow(t,5))) + epsilon(x)
// epsilon(x)??

