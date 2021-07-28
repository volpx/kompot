#include "uniconst.hpp"
#include "vec3d.hpp"
#include "simulation.hpp"
#include "optimize.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <mpi.h>

// #include <static_math/cmath.h>

#define MPI_FINISHED_TAG 2

struct Job{
	// Input
	uint32_t id;
	double rho;
	double b;
	// Outputs
	double E=0;
	double stdE=0;
	double EJF=0;
	double stdEJF=0;
};

#define JOB_SIZE sizeof(Job)
#define JOBP_SIZE sizeof(Job*)

// LJ potential
double V_LJ(const double r)
{
	return 4 * (std::pow(r, -12.0) - std::pow(r, -6.0));
}

// WF
double WF(
	const Vec3D pos[], const size_t N,
	const double b,const double L)
{
	double sum = 0;
	Vec3D other{0,0,0};
	double r_ij;

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			other = pos[j];

			if ((r_ij = compute_alias(pos[i], other, L)) > 0)
			{
				sum += -0.5 * std::pow(b / r_ij, 5.0);
			}
		}
	}
	return std::exp(sum);
}

double Density(
	const Vec3D pos[], const size_t N,
	const double b, const double L)
{
	// std::cout<<WF(pos, N, b)<<std::endl;
	double wf = WF(pos, N, b, L);
	return wf * wf;
}

// Local kinetic energy
double get_local_T(
	const Vec3D pos[], const uint32_t N,
	const double b, const double L)
{
	double T = 0;
	double r_li_s = 0;
	const double c = 30 * std::pow(b, 5.0);
	const double c1 = -c / 6; // -5b^5
	Vec3D other{0, 0, 0};
	Vec3D sum{0, 0, 0};
	double tmp;

	for (uint32_t l = 0; l < N; l++)
	{
		sum.clear();
		for (uint32_t i = 0; i < N; i++)
		{
			if (l != i)
			{
				other = pos[i];
				// r_li_s = (pos[l]-pos[i]).norm();
				r_li_s = compute_alias(pos[l], other, L);
				if (r_li_s > 0)
				{
					tmp = std::pow(r_li_s, -7.0);
					// Laplacian term
					T += -0.5 * (c + 2 * c1) * tmp;
					// Gradient term
					sum += tmp * (pos[l] - other);
				}
			}
		}
		T += 0.25 * c1 * c1 * sum.norm2();
	}

	return -T;
}

//Jackson-Feenberg kinetic energy
double get_local_TJF(
	const Vec3D pos[], const uint32_t N,
	const double b, const double L)
{
	double T = 0;
	double r_li_s = 0;
	double c = 30 * std::pow(b, 5.0);
	double c1 = -c / 6; // 5b^5
	Vec3D other;

	for (uint32_t l = 0; l < N; l++)
	{
		for (uint32_t i = 0; i < N; i++)
		{
			if (l != i)
			{
				other = pos[i];

				r_li_s = compute_alias(pos[l], other, L);
				if (r_li_s > 0)
				{
					// Laplacian term
					T += (c + 2 * c1) * std::pow(r_li_s, -7.0);
				}
			}
		}
	}

	return 0.25 * T;
}

double get_potential(
	const Vec3D pos[], const uint32_t N,
	const double L)
{
	double d;
	Vec3D alias;
	double V = 0;
	const double V_offset = -V_LJ(L / 2);

	// Sum on the pairs
	for (size_t i{0}; i < N; i++)
	{
		for (size_t j{0}; j < i; j++)
		{
			// Copy the particle j
			alias = pos[j];

			// Check if there is interaction
			// and also compute the correct alias
			if ((d = compute_alias(pos[i], alias, L)) > 0)
			{
				V += V_LJ(d) + V_offset;
			}
			// Else: no interaction, un-perfect packing of spheres
		}
	}
	return V;
}

int main(int argc, char** argv)
{
	// MPI initialization
	MPI_Init(&argc, &argv);

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// MPI stuff
	MPI_Status status;
	int bytes_transferred;

	std::cout 
		<< "VMC: 4He superfluideseese"
		<< " rank: "<< rank 
		<< std::endl;

	
	// Problem parameters
	constexpr double epsilon = uni::kB_r * 10.22;  // J
	constexpr double sigma = 0.2556e-9;			   // m
	constexpr double rho_experimental_r = 21.86e27;  // #/m^3
	constexpr double m = 4.002602 * uni::Da_to_kg; // Kg

	// Derived parameters
	constexpr double hbar2_2m = uni::hbar_r * uni::hbar_r /
								(2 * m * epsilon * sigma * sigma);
	const double b0 = std::pow(16. / (25 * hbar2_2m), 0.1);
	constexpr double rho_experimental = rho_experimental_r * sigma * sigma * sigma;

	// Number of VMC steps
	constexpr uint32_t S = 1e6;
	// Thermalization steps
	constexpr uint32_t Sth = 1e5;
	constexpr uint32_t Sdata = S - Sth;
	// Correlations steps
	constexpr size_t Scorr = 1e3;
	// Number of particles
	constexpr uint32_t N = 4 * 2 * 2 * 2;
	// Averaging length
	constexpr uint32_t K = 1e4;

	// Number of variational parameters
	constexpr uint32_t M = 20;
	// Parameter mesh
	const double ba = 1.05;
	const double bb = 1.35;
	const double bh = (bb-ba)/(M-1);
	// Density mesh
	constexpr uint32_t G = 6;
	const double rhoa = 0.7*rho_experimental;
	const double rhob = 0.9*rho_experimental;
	const double rhoh = (rhob - rhoa) / (G - 1);
	

	constexpr uint32_t Njobs=M*G;
	const uint32_t job_id_padding = Njobs % world_size;
	const uint32_t job_id_start = (rank != 0) ? Njobs / world_size * rank + job_id_padding : 0;
	const uint32_t job_id_stop = Njobs / world_size * (rank + 1) + job_id_padding;
	const uint32_t MyNjobs=job_id_stop - job_id_start;


	std::vector<Job> jobs;
	jobs.reserve(Njobs);
	for (uint32_t g = 0; g < G; g++)
		for (uint32_t m = 0; m < M; m++)
			jobs.emplace_back(Job{g*M+m,rhoa+g*rhoh,ba+m*bh});

	if (rank==0)
	{
		std::cout
			<< "Problem constants:"
			<< "\nNjobs: " << Njobs
			<< "\nS: " << S
			<< "\nN: " << N
			<< "\nepsilon: " << epsilon
			<< "\nsigma: " << sigma
			<< "\nm: " << m
			<< "\nrho_exp: " << rho_experimental
			<< "\nhbar2_2m: " << hbar2_2m
			<< "\n"
			<< std::endl;
	}

	// Final quantities (sums)
	double T = 0;
	double TJF = 0;
	double V = 0;
	WArray Es(Sdata, Sth);
	WArray EJFs(Sdata, Sth);

	// Variables
	Vec3D *posa = new Vec3D[N];
	Vec3D *posb = new Vec3D[N];
	Vec3D *pos1 = posa;
	Vec3D *pos2 = posb;
	Vec3D *tmp = posa;
	double new_prob, P1, P2;
	CArray probs(K);
	double delta=0.324034;
	uint32_t tries=0;

	// Autocorrelation
	double lambda[3];
	double *autocorr = new double[Scorr];


	for (uint32_t job_id = job_id_start; job_id < job_id_stop; job_id++)
	{	
		// Get the job
		Job &job = jobs[job_id];
		// Track the number of tries for each job
		tries++;

		const double L = std::cbrt(N / job.rho);
		const double E_correction = 8. / 3 * M_PI * job.rho * (1. / 3 * std::pow(L / 2, -9.) - std::pow(L / 2, -3.));
		const double deltamax = L / 4;
		const double deltamin = 1e-4;
		// Set the delta only on the first try
		if(tries==1)
			delta = 0.324034;
		
		std::cout
			<< "rank: " << rank
			<< " id: " << job.id
			<< " rho: " << job.rho
			<< " b: " << job.b
			<< " try: "<< tries
			<< std::endl;

		// Particles initialization
		init_lattice(pos1, N, L, 2, 4);
		apply_periodic_bounds(pos1, N, L);
		P1 = Density(pos1, N, job.b, L);

		Es.fill(0);
		EJFs.fill(0);
		probs.fill(0);

		// Compute new quantities
		T = hbar2_2m * get_local_T(pos1, N, job.b, L);
		TJF = hbar2_2m * get_local_TJF(pos1, N, job.b, L);
		V = get_potential(pos1, N, L);
		// std::cout
		// 	<< "\nInitial K_L: " << T
		// 	<< "\nInitial K_JF: " << TJF
		// 	<< "\nInitial potential: " << V
		// 	<< std::endl;

		// Parameter m
		for (uint32_t s = 0; s < S; s++)
		{
			// Step s
			// Print progress
			// if (s % (S / 10) == 0)
			// {
			// 	std::cout << s << ' ' << average(probs.data, K) << ' ' << delta << std::endl;
			// }

			// Previus positions in pos1
			// New positions in pos2
			for (uint32_t i = 0; i < N; i++)
			{
				pos2[i].x = pos1[i].x + delta * (randu() - 0.5);
				pos2[i].y = pos1[i].y + delta * (randu() - 0.5);
				pos2[i].z = pos1[i].z + delta * (randu() - 0.5);
			}
			// Periodic boundary codition
			apply_periodic_bounds(pos2, N, L);

			// Acception-rejection method
			P2 = Density(pos2, N, job.b, L);
			new_prob = P2 / P1;
			probs[s] = new_prob;

			if (randu() < new_prob)
			{
				// Accept the new state

				// Compute new quantities
				T = hbar2_2m * get_local_T(pos2, N, job.b, L);
				TJF = hbar2_2m * get_local_TJF(pos2, N, job.b, L);
				V = get_potential(pos2, N, L);

				tmp = pos2;
				pos2 = pos1;
				pos1 = tmp;
				P1 = P2;
			}

			// Sample the quantities (after "thermalization")
			if (Es.inWindow(s))
			{
				Es[s] = T + V;
				EJFs[s] = TJF + V;
			}
			else
			{
				// Touchup the delta
				if (s > 0 && s % K == 0)
				{
					delta *= (1 + 5e-1 * (average(probs.data, K) / 0.8 - 1));
					delta = (delta > deltamax) ? deltamax : (delta < deltamin) ? deltamin
																				: delta;
					// std::cout << s<< ' '<< average(probs.data, K) << ' ' << delta << std::endl;
				}
			}
		}

		// autocorrelation
		// E
		lambda[0] = 1;
		lambda[1] = 0.2;
		lambda[2] = 0;
		autocorrelation(autocorr, Scorr, Es.data, Es.N);
		fit_to_nexp(lambda, autocorr, Scorr);
		double Escorr = 1. / lambda[1];
		double Eavg = average(Es.data, Sdata) / N + E_correction;
		double DEavg = Escorr / (Es.N) * variance(Es.data, Es.N, 1) / N;
		// EJF
		lambda[0] = 1;
		lambda[1] = 0.2;
		lambda[2] = 0;
		autocorrelation(autocorr, Scorr, EJFs.data, EJFs.N);
		fit_to_nexp(lambda, autocorr, Scorr);
		double EJFscorr = 1. / lambda[1];
		double EJFavg = average(EJFs.data, Sdata)/N + E_correction;
		double DEJFavg = EJFscorr / (EJFs.N) * variance(EJFs.data, EJFs.N, 1) / N;

		// std::cout
		// 	<< "Probs:" << average(probs.data, K) << ' ' << "Delta: " << delta
		// 	<< "\nLocal energy: " << Eavg << " +/- " << std::sqrt(DEavg) << " (" << Escorr << ")"
		// 	<< "\nJF    energy: " << EJFavg << " +/- " << std::sqrt(DEJFavg) << " (" << EJFscorr << ")"
		// 	<< "\n"
		// 	<< std::endl;

		if(std::abs(Eavg-EJFavg) < 0.08 || tries > 20)
		{
			// Save data if they agree somewhat
			// Populate the job struct
			job.E=Eavg;
			job.stdE=std::sqrt(DEavg) ;
			job.EJF=EJFavg;
			job.stdEJF=std::sqrt(DEJFavg);

			tries=0;
		}
		else
		{
			// Redo the job mantainig delta
			job_id--;
		}
	}

	delete[] posa;
	delete[] posb;
	delete[] autocorr;

	if(rank==0)
	{
		// Wait for the ranks to finish
		int ranks_finished=1;
		do{
			MPI_Probe(MPI_ANY_SOURCE,MPI_FINISHED_TAG,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status,MPI_BYTE,&bytes_transferred);
			std::cout 
			<< "Received " << bytes_transferred 
			<< " from " << status.MPI_SOURCE << std::endl;
			MPI_Recv(
				&jobs[Njobs / world_size * status.MPI_SOURCE + job_id_padding],
				bytes_transferred, MPI_BYTE, status.MPI_SOURCE,
				MPI_FINISHED_TAG, MPI_COMM_WORLD, &status);
			ranks_finished++;
		}while(ranks_finished != world_size);

		// Energy output file
		std::cout << "Writing file" << std::endl;
		std::ofstream file{"data/E_rhovar_mpi.dat"};
		file << "rho b E dEavg EJF dEJFavg" << '\n';
		for(uint32_t job_id=0;job_id<Njobs;job_id++)
		{
			Job &job=jobs[job_id];
			file
				<< job.rho << ' '
				<< job.b << ' '
				<< job.E << ' '
				<< job.stdE << ' '
				<< job.EJF << ' '
				<< job.stdEJF
				<< '\n';
		}
		file.close();
	}
	else{
		MPI_Send(&jobs[job_id_start],MyNjobs*JOB_SIZE,MPI_BYTE,
			0,MPI_FINISHED_TAG,MPI_COMM_WORLD);
	}

	std::cout << "rank " << rank << " finished"<<std::endl; 

	MPI_Finalize();
	return 0;
}
