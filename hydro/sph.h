#pragma once
#include "utils.h"
#include "tree.h"
#include "threadpool.h"
#include "integrate.h"
#include "timer.h"

#include <chrono>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <assert.h>
#include <execution>

#include <omp.h>


namespace Fluid
{
	template<int d> struct m4spline
	{
		double operator()(double s, double h);
		double der(double r, double h);
		double cumul(double r, double h);
	};
	template<> struct m4spline<1>
	{
		double operator()(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? 1 - s * s * 3 / 2 + s * s*s * 3 / 4 : s <= 2 ? pow(2 - s, 3) / 4 : 0) * (2 / (3 * h));
		}
		double der(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? -s * 3 + s * s * 9 / 4 : s <= 2 ? -pow(2 - s, 2) * 3 / 4 : 0) * (2 / (3 * h*h));
		}
		double cumul(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? 4 / 3 * 2 - 2 / 3 * pow(s, 3) + 1 / 4 * pow(s, 4)
				: s <= 2 ? 8 / 3 * s - 2 * pow(s, 2) + 2 / 3 * pow(s, 3) - 1 / 12 * pow(s, 4) - 1 / 3
				: 1)/*/h*/;
		}
	};
	template<> struct m4spline<2>
	{
		double operator()(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? 1 - s * s * 3 / 2 + s * s*s * 3 / 4 : s <= 2 ? pow(2 - s, 3) / 4 : 0) * (10 / (7 * pi * h*h));
		}
		double der(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? -s * 3 + s * s * 9 / 4 : s <= 2 ? -pow(2 - s, 2) * 3 / 4 : 0) * (10 / (7 * pi * h*h*h));
		}
		double cumul(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? 10 / 7 * pow(s, 2) - 15 / 14 * pow(s, 4) + 3 / 7 * pow(s, 5)
				: s <= 2 ? 20 / 7 * pow(s, 2) - 20 / 7 * pow(s, 3) + 15 / 14 * pow(s, 4) - 1 / 7 * pow(s, 5) - 1 / 7
				: 1)/*/pow(h,2)*/;
		}
	};
	template<> struct m4spline<3>
	{
		double operator()(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? 1 - s * s * 3 / 2 + s * s*s * 3 / 4 : s <= 2 ? pow(2 - s, 3) / 4 : 0) / (h*h*h*pi);
		}
		double der(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? -s * 3 + s * s * 9 / 4 : s <= 2 ? -pow(2 - s, 2) * 3 / 4 : 0) / (h*h*h*h*pi);
		}
		double cumul(double r, double h) const
		{
			double s = r / h;
			return (s <= 1 ? 4 / 3 * pow(s, 3) - 6 / 5 * pow(s, 5) + 1 / 2 * pow(s, 6)
				: s <= 2 ? 8 / 3 * pow(s, 3) - 3 * pow(s, 4) + 6 / 5 * pow(s, 5) - 1 / 6 * pow(s, 6) - 1 / 15
				: 1)/*/pow(h,3)*/;
		}
	};

	template<int D>
	struct fileParticleData {
		vec<D> pos;
		vec<D> vel;
		double dens;
		double press;
		double intE;
		double kinE;
		double gravE;
		double group;
		double id;
		double state;

		static inline size_t numElements() {
			return { D + D + 1 + 1 + 3 + 3 };
		}
	};
	template<int D>
	struct fileSolidData {
		vec<D> com;
		vec<D> vel;
		double angle;
		double omega;
		vec<D> dragForce;
		double group;
		
		static inline size_t numElements() {
			return { D + D + 2 + D + 1};
		}
	};
	
	enum class BoundaryType
	{
		outflow,  // removes particles that move out
		inflow,   // respawns removed particles uniformly
		periodic, // moves in again on the other side
		wall      // keeps particles in the domain
	};

	enum particleState : unsigned int{
		psRemoved = 1u
	};
	
	template<int D>
	class SPH
	{
		double gamma; // adiabatic index
		const double G = 8 * 2;
		double alpha; // volume viscosity
		double beta;  // other viscosity
		std::array<vec<2>, D> domain;
		std::array<double, D> domainSpan;

		const m4spline<D> kernel = m4spline<D>();


		struct particle{
			vec<D> pos;
			vec<D> vel;
			vec<D> accel;
			double mass;
			double density;
			double pressure;
			double energy;
			double energyDer;
			double delVabs;
			int    group=0; // 0=gas. 1,2..=solids
			int	   id; // inside their group
			particleState state;

			template<bool ISOTHERMAL>
			inline double get_c(double gamma) const
			{
				if (ISOTHERMAL)
					return sqrt(energy);
				else
					return sqrt(gamma*(gamma - 1)*energy);
			}
			template<bool ISOTHERMAL>
			inline double get_dt(double gamma, double h, double alpha, double beta) const
			{
				double c = get_c<ISOTHERMAL>(gamma);
				double t1 = h / (h*delVabs + c);
				double t2 = sqrt(h / (accel.length() + 0.0001));
				double t3 = h / ((1 + 1.2*alpha)*c + (1 + 1.2*beta)*h*delVabs);
				return std::min(std::min(t1, t2), t3);
			}

			static inline const vec<D>& getPos(const particle& p) { return p.pos; }
			static inline const double& getMass(const particle& p) { return p.mass; }
		};

		std::vector<particle> particles;
		int nActiveParticles = 0;
		KDTree<D, particle> kdt;
		
		struct solid {
			vec<D> com; // center of mass
			vec<D> vel; // linear velocity
			double totalMass;

			double angle; // rotation angle
			double omega; // angular velocity
			double momentOfInertia;

			int group;
			vec<D> dragForce;

			vec<D> positionLocked;
			bool rotationLocked = false;

			std::vector<particle*> particles; // temporarily pointing inside SPH::particles
			std::vector<vec<2>> subpositions; // not nessecariy in the same order as "particles"
		};
		std::vector<solid> solids;

		std::array<std::array<BoundaryType,2>, D> boundaryTypes;

		vec<D> externalAcceleration, externalForce;
		double externalTime = 0;

		using srn = typename decltype(kdt)::search_result_neigh;
		using srv = typename decltype(kdt)::search_result_vneigh;
		using srm = typename decltype(kdt)::search_result_mass;


		double get_h(double eta)
		{
			// Returns the kernel size for this timestep, scales quadratically in nActiveParticles and schould only be used in the beginning 
			AutoTimer at(g_timer, "get_h");

			// avg is the mean of the minimal distances
			double avg = 0;

			#pragma omp parallel for reduction(+:avg)
			for (int i = 0; i < nActiveParticles; i++)
			{
				const auto& p1 = particles[i];
				double md = 9999;
				for (int j = 0; j < nActiveParticles; j++)
				{
					const auto& p2 = particles[j];
					if (j != i)
						md = std::min(md, (p1.pos - p2.pos).length());
				}
				avg += md;
			}

			return eta * avg / nActiveParticles;
		}

	
		inline void get_viscPiLambda(double& out_Pi, double& out_lambda,
			double c_i, double c_j, double rho_i, double rho_j, vec<D> v_ij, vec<D> e_ij)
		{
			auto dot = v_ij.dot(e_ij);
			if (dot <= 0)
			{
				auto vSig = c_i + c_j - beta * dot;
				auto temp = -alpha * vSig * dot / (rho_i + rho_j);
				out_Pi = temp * 2;
				out_lambda = temp * dot;
			}
			else out_Pi = out_lambda = 0;
		}
	

		void rebuildTree()
		{
			AutoTimer at(g_timer, "rebuildTree");
			
			kdt.construct(particles.begin(), particles.begin()+nActiveParticles);
		}

		template<bool GRAVITATION>
		void findNeighbours(const double h, const double grav_angle,
			std::vector<std::vector<srn>>& neigs,
			std::vector<std::vector<srv>>& neigs_mirror,
			std::vector<std::vector<srm>>& masses)
		{
			AutoTimer at(g_timer, "findNeighbours");

			double r = 2 * h;
			std::array<int, D> mirrored;

			rebuildTree();

			#pragma omp parallel for 
			for (int i = 0; i < nActiveParticles; i++)
			{
				int counter;
				neigs[i].clear();
				kdt.searchRadius(particles[i].pos, r, neigs[i], counter);
				
				
				neigs_mirror[i].clear();
				for (int k = 0; k < D; k++) {
					if (boundaryTypes[k][0] == BoundaryType::periodic) {
						if (domain[k][1] - particles[i].pos[k] < r)
						{
							vec<D> vpos = particles[i].pos; vpos[k] -= domainSpan[k];
							kdt.searchRadius(vpos, r, neigs_mirror[i], counter);
							mirrored[k] = 1;

							for (int kk = 0; kk < k; kk++)
								if (mirrored[kk] == 1)
								{
									vpos[kk] -= domainSpan[k];
									kdt.searchRadius(vpos, r, neigs_mirror[i], counter);
									break;
								}
								else if (mirrored[kk] == -1)
								{
									vpos[kk] += domainSpan[k];
									kdt.searchRadius(vpos, r, neigs_mirror[i], counter);
									break;
								}
						}
						else if (particles[i].pos[k] - domain[k][0] < r)
						{
							vec<D> vpos = particles[i].pos; vpos[k] += domainSpan[k];
							kdt.searchRadius(vpos, r, neigs_mirror[i], counter);
							mirrored[k] = -1;

							for (int kk = 0; kk < k; kk++)
								if (mirrored[kk] == 1)
								{
									vpos[kk] -= domainSpan[k];
									kdt.searchRadius(vpos, r, neigs_mirror[i], counter);
									break;
								}
								else if (mirrored[kk] == -1)
								{
									vpos[kk] += domainSpan[k];
									kdt.searchRadius(vpos, r, neigs_mirror[i], counter);
									break;
								}
						}
						else mirrored[k] = 0;
					}
					else mirrored[k] = 0;

					// For wall boundaries we let each particle see its phantom mirrored version
					if (boundaryTypes[k][1] == BoundaryType::wall) {
						if (domain[k][1] - particles[i].pos[k] < r/2)
						{
							vec<D> vpos = particles[i].pos; vpos[k] = 2*domain[k][1] - vpos[k];
							
							auto conn = vpos - particles[i].pos;
							double dist = conn.length();
							neigs_mirror[i].push_back(srv{ &particles[i],conn,dist });
						}
					}
					if (boundaryTypes[k][0] == BoundaryType::wall) {
						if (particles[i].pos[k] - domain[k][0] < r/2)
						{
							vec<D> vpos = particles[i].pos; vpos[k] = 2 * domain[k][0] - vpos[k];

							auto conn = vpos - particles[i].pos;
							double dist = conn.length();
							neigs_mirror[i].push_back(srv{ &particles[i],conn,dist });
						}
					}
				}
				
				if constexpr (GRAVITATION)
				{
					masses[i].clear();
					kdt.searchLayerMasses(particles[i].pos, grav_angle, masses[i]);
				}
			}
		}

		template<bool ISOTHERMAL>
		void computeDensity(double h,
			std::vector<std::vector<srn>> const& neigs,
			std::vector<std::vector<srv>> const& neigs_mirror)
		{
			AutoTimer at(g_timer, "computeDensity");
			for (int i = 0; i < nActiveParticles; i++)
			{
				auto& p1 = particles[i];
				p1.density = 0;
				for (auto p2p : neigs[i])
				{
					const particle& p2 = *p2p.neigh;
					double dist = p2p.dist;
					p1.density += p2.mass*kernel(dist, h);
				}

				for (auto p2p : neigs_mirror[i])
				{
					const particle& p2 = *p2p.neigh;
					p1.density += p2.mass*kernel(p2p.dist, h);
				};

				if constexpr (!ISOTHERMAL)
					p1.pressure = p1.density*p1.energy*(gamma - 1);
				else
					p1.pressure = p1.density*p1.energy;
			}
		}

		void logParticles(StreamFileWriter& fw, double time)
		{
			AutoTimer at(g_timer, "logParticles");

			using srm = typename decltype(kdt)::search_result_mass;
			auto grav_pot = [this](const vec<D>& pos)
			{
				std::vector<srm> masses(nActiveParticles);
				kdt.searchLayerMasses(pos, 0, masses);

				return -G * std::accumulate(masses.begin(), masses.end(), 0., [&](double l, const srm& r)
					{auto dist = (r.com - pos).length(); return dist > 0.0001 ? l + r.mass / dist : l; });
			};	

			std::vector<fileParticleData<D>> data(nActiveParticles);
			std::transform(particles.begin(), particles.begin()+nActiveParticles, data.begin(), [&grav_pot](const particle& p)
				{
					return fileParticleData{ p.pos, p.vel, p.density, p.pressure,
						p.energy * p.mass, p.vel.lengthsq() * (0.5 * p.mass), 0/*grav_pot(p.pos) * p.mass*/ ,
						(double)p.group, (double)p.id, (double)p.state};
				});

			fw.writeBlocks(data);
			fw.finishSection(time);
		}

		void logSolids(StreamFileWriter& fw, double time)
		{
			AutoTimer at(g_timer, "logSolids");

			std::vector<fileSolidData<D>> data(solids.size());
			std::transform(solids.begin(), solids.end(), data.begin(), [](const solid& so)
				{
					return fileSolidData{ so.com, so.vel, so.angle, so.omega,
						so.dragForce, (double)so.group };
				});

			fw.writeBlocks(data);
			fw.finishSection(time);
		}


		void partitionParticles() {
			AutoTimer at(g_timer, "partitionParticles");
			// move all particles that are "removed" to the end of the particle vector
			nActiveParticles = std::partition(particles.begin(), particles.end(), [](particle const& p) {return !(p.state & psRemoved); }) - particles.begin();
		}


	public:
		SPH(){}
		void init(double gamm, double alph, int leafCapacity) {
			gamma = gamm;
			alpha = alph;
			beta = 2 * alpha;

			kdt.setPositionCB(&particle::getPos);
			kdt.setMassCB(&particle::getMass); 
			kdt.setLeafSize(leafCapacity); // number of particles per kd-tree leaf
		}
		void setDomain(decltype(domain) dom, decltype(boundaryTypes) bt) {
			domain = dom;
			boundaryTypes = bt;

			bool problemFound = false;
			for (int d = 0; d < D; ++d)
				domainSpan[d] = domain[d][1] - domain[d][0];
			for (int d = 0; d < D; ++d)
				if ((boundaryTypes[d][0] == BoundaryType::periodic) != ((boundaryTypes[d][1] == BoundaryType::periodic))) {
					std::cout << "WARNING: Dimension " << d << ": opposite sides must both have periodic boundary conditions.";
					problemFound = true;
				}
			if (!problemFound)
				std::cout << "Boundary consistency check passed!" << std::endl;
		}
		void setExternalAcceleration(vec<D> ea) {
			externalAcceleration = ea;
		}
		void setExternalForce(vec<D> fo) {
			externalForce = fo;
		}
		void setExternalTime(double et) {
			externalTime = et;
		}

		void createParticles(int numParticles, double totalMass = 1.0) {
			// The mass has no effect if all particles have the same mass
			size_t offset = particles.size();
			particles.resize(offset+numParticles);
			for (int i = 0; i < numParticles; i++){
				auto x = ((rand() % 10000) / 10000.0);
				auto y = ((rand() % 10000) / 10000.0);
				particles[offset + i] = { { x, y },{ 0, 0 },{ 0, 0 }, totalMass / numParticles, 0, 0, 1, 0 };
			}

			// Make sure that "removed" particles are in the end
			partitionParticles();
		}
		void createSolid(std::vector<vec<2>> const& subp) {
			solid so;
			so.totalMass = .1;
			int gr = solids.size()+1;
			so.group = gr;

			vec<2> ori = { 0.5,0.5 };
			so.subpositions = subp;
			
			size_t offset = particles.size();
			particles.resize(offset + so.subpositions.size());
			for (int i = 0; i < so.subpositions.size(); i++) {
				particles[offset + i] = { ori + (so.subpositions[i] *0.01),{0, 0},{0, 0},
					so.totalMass / so.subpositions.size(), 0, 0, 1, 0};
				particles[offset + i].group = gr;
				particles[offset + i].id = i;
				so.particles.push_back(&particles[offset + i]);
			}

			// Compute center of mass and moment of inertia
			for (auto* p : so.particles) {
				so.com += p->pos * p->mass;	
			}
			so.com = so.com * (1. / so.totalMass);
			for (auto* p : so.particles) {
				so.momentOfInertia += (p->pos - so.com).lengthsq() * p->mass;
			}
			
			// Make sure that "removed" particles are in the end
			partitionParticles();

			// recompute the subpositions relative to the COM
			for (int i = 0; i < so.subpositions.size(); i++) {
				so.subpositions[i] = so.particles[i]->pos - so.com;
			}

			// only rotate around COM, dont move
			so.positionLocked = {1,1};
			so.rotationLocked = true;

			solids.push_back(so);
		};


		auto& getTree() { return kdt; }

		//size_t getNumParticles() { return particles.size(); }
		size_t getNumSolids() { return solids.size(); }

		template<bool GRAVITATION = false, bool ISOTHERMAL = false>
		void simulate(StreamFileWriter& pfw, StreamFileWriter& sfw, const double simulTime,
			const double eta, const double cfl = 1., const double logTimeStep = 0.01)
		{
			std::ofstream etafile("etas.txt");


			using namespace std;

			cout << "Parameters: N = " << particles.size() << ", T = " << simulTime << ", eta = " << eta << ", cfl = " << cfl << std::endl;

			const double grav_angle = 0.0001;

			double t = 0;
			double h;

			std::vector<std::vector<srn>> neigs(particles.size());
			std::vector<std::vector<srv>> neigs_mirror(particles.size());
			std::vector<std::vector<srm>> masses(particles.size());

			h = get_h(eta);
			findNeighbours<GRAVITATION>(h, grav_angle, neigs, neigs_mirror, masses);
			computeDensity<ISOTHERMAL>(h, neigs, neigs_mirror);
			
			logParticles(pfw, t);
			logSolids(sfw, t);

			double avgMinDistance;
			int counter = 0;
			double dt;
			while (t < simulTime)
			{
				AutoTimer at(g_timer, "simulate - time step");

			//	h = get_h(eta);

				// find neighbours
				findNeighbours<GRAVITATION>(h, grav_angle, neigs, neigs_mirror, masses);

				// compute density and pressure
				computeDensity<ISOTHERMAL>(h, neigs, neigs_mirror);

				// compute energy change and acceleration
				{
					AutoTimer at(g_timer, "simulate - compute change");
					float avgNumNeigh = 0;
					avgMinDistance = 0;

					#pragma omp parallel for 
					for (int i = 0; i < nActiveParticles; i++)
					{
						/*if (t == 0 && i == 0 && omp_get_thread_num() == 0)
							std::cout << "Number of opm threads working: " << omp_get_num_threads() << std::endl;*/

						avgNumNeigh += float(neigs[i].size() + neigs_mirror[i].size()) / nActiveParticles;

						double minDistance = h;

						auto& p1 = particles[i];
						vec<D> accel_acc;
						vec<D> visc_acc;
						vec<D> grav_acc;
						double inte_acc = 0;
						double inte_acc_visc = 0;
						double delV_acc = 0;
						double lambda, Pi;
						for (const srn& p2p : neigs[i])
						{
							const particle& p2 = *p2p.neigh;
							double dist = p2p.dist;
							if(dist>0)
							minDistance = std::min(minDistance, dist);
							auto rij = (p1.pos - p2.pos);
							auto vij = (p1.vel - p2.vel);
							auto rijn = dist == 0 ? rij : rij / dist;
							double kd = kernel.der(dist, h);

							get_viscPiLambda(Pi, lambda, p1.get_c<ISOTHERMAL>(gamma), p2.get_c<ISOTHERMAL>(gamma), p1.density, p2.density, vij, rijn);

							accel_acc += rijn * -p2.mass * (p1.pressure / (p1.density * p1.density) + p2.pressure / (p2.density * p2.density)) * kd;
							visc_acc += rijn * -p2.mass * Pi * kd;

							inte_acc += vij.dot(rijn) * p2.mass * kd;
							inte_acc_visc += lambda * p2.mass * kd;

							delV_acc += vij.dot(rijn) * p2.mass / p2.density * kd;
						}
						for (const srv& p2p : neigs_mirror[i])
						{
							const particle& p2 = *p2p.neigh;
							auto rij = p2p.relPos * (-1);
							auto vij = (p1.vel - p2.vel);
							double dist = rij.length();
							minDistance = std::min(minDistance, dist);
							auto rijn = dist == 0 ? rij : rij * (1 / dist);
							double kd = kernel.der(dist, h);

							get_viscPiLambda(Pi, lambda, p1.get_c<ISOTHERMAL>(gamma), p2.get_c<ISOTHERMAL>(gamma), p1.density, p2.density, vij, rijn);

							accel_acc += rijn * -p2.mass * (p1.pressure / (p1.density * p1.density) + p2.pressure / (p2.density * p2.density)) * kd;
							visc_acc += rijn * -p2.mass * Pi * kd;

							inte_acc += vij.dot(rijn) * p2.mass * kd;
							inte_acc_visc += lambda * p2.mass * kd;

							delV_acc += vij.dot(rijn) * p2.mass / p2.density * kd;
						};
						if constexpr (GRAVITATION)
						{
							for (const srm& cm : masses[i])
							{
								auto conn = (cm.com - p1.pos);
								double distsq = conn.lengthsq();
								if (distsq != 0)
									grav_acc += conn.normalize() * (G * cm.mass * kernel.cumul(sqrt(distsq), h) / distsq);
							}
						}

						p1.delVabs = abs(delV_acc);
						p1.accel = accel_acc + visc_acc + grav_acc;
						p1.energyDer = p1.pressure / (p1.density * p1.density) * inte_acc + inte_acc_visc;
						avgMinDistance += minDistance;
					}
					
					etafile << avgNumNeigh << " " << h << " " << avgMinDistance<<std::endl;
				}

				

				// calculate time step
				{
					AutoTimer at(g_timer, "simulate - calculate time step");
					dt = 9999;
					for (int i = 0; i < nActiveParticles; i++) {
						const particle& p1 = particles[i];
						dt = std::min(dt, p1.get_dt<ISOTHERMAL>(gamma, h, alpha, beta));
					}
					dt *= cfl;
				}

				// update energy and position
				{
					AutoTimer at(g_timer, "simulate - apply change");

					for (auto& so : solids)
						so.particles.clear();

					for (int i = 0; i < nActiveParticles; i++){
						particle& p1 = particles[i];

						if (t < externalTime) {
							p1.accel += (externalAcceleration + externalForce / p1.mass)/** (p1.pos.y < 0.5 ? -1. : 1.)*/;
						}

						if (p1.group != 0) {
							solids[p1.group - 1].particles.push_back(&p1);
							continue;
						}
						else {
							if (!ISOTHERMAL)
								p1.energy += p1.energyDer * dt;
							p1.pos += p1.vel * dt;
							p1.vel += p1.accel * dt;
						}
					}
					for (solid& so : solids)
					{
						// compute linear and angular force (torque) from particle accelerations
						so.dragForce = { 0,0 };
						double totalTorque = 0;
						for (auto* p : so.particles) {
							vec<2> conn = p->pos - so.com;
							so.dragForce += p->accel * p->mass;
							totalTorque += (conn.cross(p->accel)) * p->mass;
						}

						// move and rotate the solid
						auto linAccel = so.dragForce / so.totalMass * (vec<D>(1)-so.positionLocked);
						so.com += so.vel * dt;
						so.vel += linAccel * dt;
						
						if (!so.rotationLocked) {
							auto angAccel = totalTorque / so.momentOfInertia;
							so.angle += so.omega * dt;
							so.omega += angAccel * dt;
						}

						// recompute the individual particle positions
						for (auto* p : so.particles) {
							p->pos = so.com + so.subpositions[p->id].rotate(so.angle);
						}
					}
					
					// apply boundary conditions
					for (int i = 0; i < nActiveParticles; i++) {
						particle& p1 = particles[i];
						for (int k = 0; k < D; k++) {
							if (boundaryTypes[k][0] == BoundaryType::periodic) {
								while (p1.pos[k] < domain[k][0])
									p1.pos[k] += domainSpan[k];
								while (p1.pos[k] > domain[k][1])
									p1.pos[k] -= domainSpan[k];
								continue;
							}
							if (boundaryTypes[k][0] == BoundaryType::wall) {
								if (p1.pos[k] < domain[k][0]) {
									p1.pos[k] = domain[k][0];
									p1.vel[k] *= -1;
									p1.accel[k] = 0;
								}
							}
							if (boundaryTypes[k][1] == BoundaryType::wall) {
								if (p1.pos[k] > domain[k][1]) {
									p1.pos[k] = domain[k][1];
									p1.vel[k] *= -1;
									p1.accel[k] = 0;
								}
							}
							if (boundaryTypes[k][0] == BoundaryType::outflow) {
								if (p1.pos[k] < domain[k][0]) {
									p1.state |= psRemoved;
								}
							}
							if (boundaryTypes[k][1] == BoundaryType::outflow) {
								if (p1.pos[k] > domain[k][1]) {
									p1.state |= psRemoved;
								}
							}
							if (boundaryTypes[k][0] == BoundaryType::inflow) {
								if (p1.pos[k] < domain[k][0]) {
									p1.pos[k] = domain[k][0];
									p1.vel[k] *= -1;
									p1.accel[k] = 0;
								}
							}
							if (boundaryTypes[k][1] == BoundaryType::inflow) {
								if (p1.pos[k] > domain[k][1]) {
									p1.pos[k] = domain[k][1];
									p1.vel[k] *= -1;
									p1.accel[k] = 0;
								}
							}
						}
					}

					partitionParticles();
				}


				// create new particles on inflow edges			
				int nSpawnParticles = 10;
				double wallOffset = 0.01;
				for (int k = 0; k < D; ++k) {
					if (boundaryTypes[k][0] == BoundaryType::inflow) {
						size_t nMaxParticles = particles.size();
						if (nMaxParticles - nActiveParticles >= nSpawnParticles) {
							for (int i = nActiveParticles; i < nActiveParticles+nSpawnParticles; i++) {				
								particles[i].pos[k] = domain[k][0] + wallOffset;
								particles[i].pos[1 - k] = domain[1 - k][0] + ((rand() % 10000) / 10000.0) * (domain[1 - k][1] - domain[1 - k][0]);
								particles[i].vel = { 0,0 };
								particles[i].state &= (~psRemoved);
							}
						}
						nMaxParticles += nSpawnParticles;
					}
				}
				
			
				t += dt;

				if ((int)(t / logTimeStep) != (int)((t - dt) / logTimeStep))
				{
					logParticles(pfw, t);
					logSolids(sfw, t);
				}
				counter++;

				if (counter % 100 == 0)
					std::cout << "step " << counter << ", t=" << t << " / " << simulTime << std::endl;

				// estimate smoothing length for the next step
				h = eta * avgMinDistance / particles.size();
			}

			cout << "Simulation finished. Timesteps: " << counter << std::endl;

			etafile.close();
		}

	};

}
