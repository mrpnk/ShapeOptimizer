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
		double operator()(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? 1 - s * s * 3 / 2 + s * s*s * 3 / 4 : s <= 2 ? pow(2 - s, 3) / 4 : 0) * (2 / (3 * h));
		}
		double der(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? -s * 3 + s * s * 9 / 4 : s <= 2 ? -pow(2 - s, 2) * 3 / 4 : 0) * (2 / (3 * h*h));
		}
		double cumul(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? 4 / 3 * 2 - 2 / 3 * pow(s, 3) + 1 / 4 * pow(s, 4)
				: s <= 2 ? 8 / 3 * s - 2 * pow(s, 2) + 2 / 3 * pow(s, 3) - 1 / 12 * pow(s, 4) - 1 / 3
				: 1)/*/h*/;
		}
	};
	template<> struct m4spline<2>
	{
		double operator()(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? 1 - s * s * 3 / 2 + s * s*s * 3 / 4 : s <= 2 ? pow(2 - s, 3) / 4 : 0) * (10 / (7 * pi * h*h));
		}
		double der(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? -s * 3 + s * s * 9 / 4 : s <= 2 ? -pow(2 - s, 2) * 3 / 4 : 0) * (10 / (7 * pi * h*h*h));
		}
		double cumul(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? 10 / 7 * pow(s, 2) - 15 / 14 * pow(s, 4) + 3 / 7 * pow(s, 5)
				: s <= 2 ? 20 / 7 * pow(s, 2) - 20 / 7 * pow(s, 3) + 15 / 14 * pow(s, 4) - 1 / 7 * pow(s, 5) - 1 / 7
				: 1)/*/pow(h,2)*/;
		}
	};
	template<> struct m4spline<3>
	{
		double operator()(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? 1 - s * s * 3 / 2 + s * s*s * 3 / 4 : s <= 2 ? pow(2 - s, 3) / 4 : 0) / (h*h*h*pi);
		}
		double der(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? -s * 3 + s * s * 9 / 4 : s <= 2 ? -pow(2 - s, 2) * 3 / 4 : 0) / (h*h*h*h*pi);
		}
		double cumul(double r, double h)
		{
			double s = r / h;
			return (s <= 1 ? 4 / 3 * pow(s, 3) - 6 / 5 * pow(s, 5) + 1 / 2 * pow(s, 6)
				: s <= 2 ? 8 / 3 * pow(s, 3) - 3 * pow(s, 4) + 6 / 5 * pow(s, 5) - 1 / 6 * pow(s, 6) - 1 / 15
				: 1)/*/pow(h,3)*/;
		}
	};

	template<int D>
	struct fileDat {
		vec<D> pos;
		vec<D> vel;
		double dens;
		double press;
		double intE;
		double kinE;
		double gravE;
		double group;
		double id;
		static inline size_t numElements() {
			return { D + D + 1 + 1 + 3 + 2 };
		}
	};
	
	template<int D>
	class SPH
	{
	public:
		enum BoundaryType
		{
			outflow,
			periodic,
			wall
		};

	private:
		double gamma; // adiabatic index
		const double G = 8 * 2;
		double alpha; // volume viscosity
		double beta;  // other viscosity
		std::array<vec<2>, D> domain;
		std::array<double, D> domainSpan;

		m4spline<D> kernel = m4spline<D>();

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

			template<bool ISOTHERMAL>
			double get_c(double gamma) const
			{
				if (ISOTHERMAL)
					return sqrt(energy);
				else
					return sqrt(gamma*(gamma - 1)*energy);
			}
			template<bool ISOTHERMAL>
			double get_dt(double gamma, double h, double alpha, double beta) const
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

		struct solid {
			vec<D> com; // center of mass
			vec<D> vel; // linear velocity
			double totalMass;

			double angle; // rotation angle
			double omega; // angular velocity
			double momentOfInertia;

			std::vector<particle*> particles; // temporarily pointing inside SPH::particles
			std::vector<vec<2>> subpositions; // not nessecariy in the same order as "particles"
		};
		std::vector<solid> solids;

		std::vector<particle> particles;
		KDTree<D, particle> kdt;

		std::array<std::array<BoundaryType,2>, D> boundaryTypes;

		vec<D> externalAcceleration, externalForce;
		double externalTime = 0;

		using srn = typename decltype(kdt)::search_result_neigh;
		using srv = typename decltype(kdt)::search_result_vneigh;
		using srm = typename decltype(kdt)::search_result_mass;


		double get_h(double nu)
		{
			// Returns the kernel size for this timestep
			AutoTimer at(g_timer, "get_h");

			double avg = 0;

			#pragma omp parallel for reduction(+:avg)
			for (int i = 0; i < particles.size(); i++)
			{
				const auto& p1 = particles[i];
				double md = 9999;
				for (int j = 0; j < particles.size(); j++) // TODO use the neighbours instead of quadratic
				{
					const auto& p2 = particles[j];
					if (j != i)
						md = std::min(md, (p1.pos - p2.pos).length());
				}
				avg += md;
			}
			return nu * avg / particles.size();
		}

		double get_viscPi(double c_i, double c_j, double rho_i, double rho_j, vec<D> v_ij, vec<D> r_ij_norm)
		{
			auto dot = v_ij.dot(r_ij_norm);
			if (dot <= 0)
			{
				auto vSig = c_i + c_j - beta * dot;
				return -alpha * vSig*dot * 2 / ((rho_i + rho_j));
			}
			return 0;
		}
		double get_viscLambda(double c_i, double c_j, double rho_i, double rho_j, vec<D> v_ij, vec<D> r_ij_norm)
		{
			auto dot = v_ij.dot(r_ij_norm);
			if (dot <= 0)
			{
				auto vSig = c_i + c_j - beta * dot;
				return -alpha * vSig*dot*dot / (rho_i + rho_j);
			}
			return 0;
		}

		void rebuildTree()
		{
			AutoTimer at(g_timer, "rebuildTree");
			
			kdt.construct(particles.begin(), particles.end());
		}

		void findNeighbours(const double h, const double grav_angle,
			const bool _enable_grav_,
			std::vector<std::vector<srn>>& neigs,
			std::vector<std::vector<srv>>& neigs_mirror,
			std::vector<std::vector<srm>>& masses)
		{
			AutoTimer at(g_timer, "findNeighbours");

			double r = 2 * h;
			std::array<int, D> mirrored;

			rebuildTree();

			#pragma omp parallel for 
			for (int i = 0; i < particles.size(); i++)
			{
				int counter;
				neigs[i].clear();
				kdt.searchRadius(particles[i].pos, r, neigs[i], counter);
				
				
				neigs_mirror[i].clear();
				for (int k = 0; k < D; k++) {
					if (boundaryTypes[k][0] == periodic) {
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
					if (boundaryTypes[k][1] == wall) {
						if (domain[k][1] - particles[i].pos[k] < r/2)
						{
							vec<D> vpos = particles[i].pos; vpos[k] = 2*domain[k][1] - vpos[k];
							
							auto conn = vpos - particles[i].pos;
							double dist = conn.length();
							neigs_mirror[i].push_back(srv{ &particles[i],conn,dist });
						}
					}
					if (boundaryTypes[k][0] == wall) {
						if (particles[i].pos[k] - domain[k][0] < r/2)
						{
							vec<D> vpos = particles[i].pos; vpos[k] = 2 * domain[k][0] - vpos[k];

							auto conn = vpos - particles[i].pos;
							double dist = conn.length();
							neigs_mirror[i].push_back(srv{ &particles[i],conn,dist });
						}
					}
				}
				
				if (_enable_grav_)
				{
					masses[i].clear();
					kdt.searchLayerMasses(particles[i].pos, grav_angle, masses[i]);
				}
			}
		}

		void computeDensity(double h, std::vector<std::vector<srn>> neigs,
			std::vector<std::vector<srv>> neigs_mirror, bool isothermal)
		{
			AutoTimer at(g_timer, "computeDensity");
			for (int i = 0; i < particles.size(); i++)
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

				if (!isothermal)
					p1.pressure = p1.density*p1.energy*(gamma - 1);
				else
					p1.pressure = p1.density*p1.energy;
			}
		}

		void logParticles(CubeFileWriter& fw, double time)
		{
			AutoTimer at(g_timer, "logParticles");

			using srm = typename decltype(kdt)::search_result_mass;
			auto grav_pot = [this](const vec<D>& pos)
			{
				std::vector<srm> masses(particles.size());
				kdt.searchLayerMasses(pos, 0, masses);

				return -G * std::accumulate(masses.begin(), masses.end(), 0., [&](double l, const srm& r)
					{auto dist = (r.com - pos).length(); return dist > 0.0001 ? l + r.mass / dist : l; });
			};	

			std::vector<fileDat<D>> data(particles.size());		
			std::transform(particles.begin(), particles.end(), data.begin(), [&grav_pot](const particle& p)
				{
					return fileDat{ p.pos, p.vel, p.density, p.pressure,
						p.energy * p.mass, p.vel.lengthsq() * (0.5 * p.mass), 0/*grav_pot(p.pos) * p.mass*/ ,
						(double)p.group, (double)p.id};
				});

			fw.newLayer();
			fw.writeData(data);
			fw.writeHeader(time, 0);	
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
				if ((boundaryTypes[d][0] == periodic) != ((boundaryTypes[d][1] == periodic))) {
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
			for (int i = 0; i < numParticles; i++)
			{
				auto x = ((rand() % 10000) / 10000.0);
				auto y = ((rand() % 10000) / 10000.0);
				particles[offset + i] = { { x, y },{ 0, 0 },{ 0, 0 }, totalMass / numParticles, 0, 0, 1, 0 };
			}
		}
		void createSolid() {
			solid so;
			so.totalMass = 1;
			vec<2> ori = { 0.5,0.5 };

			so.subpositions = { {1,-2}, {0,-1}, {0,0}, {0,1}, {-1,2} };
			size_t offset = particles.size();
			particles.resize(offset + so.subpositions.size());
			for (int i = 0; i < so.subpositions.size(); i++) {
				particles[offset + i] = { ori +(so.subpositions[i]*0.05),{ 0, 0 },{ 0, 0 },
					so.totalMass / so.subpositions.size(), 0, 0, 1, 0};
				particles[offset + i].group = 1;
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
			

			// recompute the subpositions relative to the COM
			for (int i = 0; i < so.subpositions.size(); i++) {
				so.subpositions[i] = so.particles[i]->pos - so.com;
			}


			solids.push_back(so);
		};


		template<bool GRAVITATION = false, bool ISOTHERMAL = false>
		void simulate(std::string m_filename, const double simulTime, const double eta, const double cfl = 1., const double logTimeStep = 0.01)
		{
			using namespace std;

			cout << "Parameters: N = " << particles.size() << ", T = " << simulTime << ", eta = " << eta << ", cfl = " << cfl << std::endl;
			cout << "Output file: " << m_filename << std::endl;

			const double grav_angle = 0.0001;

			double t = 0;
			double h;

			CubeFileWriter fw;
			fw.setShape({ particles.size(), fileDat<D>::numElements() });
			fw.open(m_filename);

			std::vector<std::vector<srn>> neigs(particles.size());
			std::vector<std::vector<srv>> neigs_mirror(particles.size());
			std::vector<std::vector<srm>> masses(particles.size());

			h = get_h(eta);
			findNeighbours(h, grav_angle, GRAVITATION, neigs, neigs_mirror, masses);
			computeDensity(h, neigs, neigs_mirror, ISOTHERMAL);
			logParticles(fw, t);

			auto startTime = chrono::high_resolution_clock::now();

			double neighbourTime = 0;
			int counter = 0;
			double dt;
			while (t < simulTime)
			{
				AutoTimer at(g_timer, "simulate - time step");

				h = get_h(eta);
				
				// find neighbours
				findNeighbours(h, grav_angle, GRAVITATION, neigs, neigs_mirror, masses);

				// compute density and pressure
				computeDensity(h, neigs, neigs_mirror, ISOTHERMAL);

				// compute energy change and acceleration
				{
					AutoTimer at(g_timer, "simulate - compute change");
					#pragma omp parallel for 
					for (int i = 0; i < particles.size(); i++)
					{
						if (t == 0 && i == 0 && omp_get_thread_num() == 0)
							std::cout << "Number of opm threads working: " << omp_get_num_threads() << std::endl;

						auto& p1 = particles[i];
						vec<D> accel_acc;
						vec<D> visc_acc;
						vec<D> grav_acc;
						double inte_acc = 0;
						double inte_acc_visc = 0;
						double delV_acc = 0;
						double lambda;
						for (srn& p2p : neigs[i])
						{
							const particle& p2 = *p2p.neigh;
							double dist = p2p.dist;
							auto rij = (p1.pos - p2.pos);
							auto vij = (p1.vel - p2.vel);
							auto rijn = dist == 0 ? rij : rij * (1 / dist);
							accel_acc += rijn * (-p2.mass * (p1.pressure / (p1.density * p1.density) + p2.pressure / (p2.density * p2.density)) * kernel.der(dist, h));
							lambda = get_viscLambda(p1.get_c<ISOTHERMAL>(gamma), p2.get_c<ISOTHERMAL>(gamma), p1.density, p2.density, vij, rijn);
							inte_acc += vij.dot(rijn) * p2.mass * kernel.der(dist, h);
							inte_acc_visc += rijn.dot(rijn) * lambda * p2.mass * kernel.der(dist, h);
							delV_acc += vij.dot(rijn) * (p2.mass / p2.density * kernel.der(dist, h));
							visc_acc += rijn * p2.mass * kernel.der(dist, h) * -get_viscPi(p1.get_c<ISOTHERMAL>(gamma), p2.get_c<ISOTHERMAL>(gamma), p1.density, p2.density, vij, rijn);
						};
						for (auto p2p : neigs_mirror[i])
						{
							const particle& p2 = *p2p.neigh;

							auto rij = p2p.relPos * (-1);
							auto vij = (p1.vel - p2.vel);
							double dist = rij.length();
							auto rijn = dist == 0 ? rij : rij * (1 / dist);
							accel_acc += rijn * (-p2.mass * (p1.pressure / (p1.density * p1.density) + p2.pressure / (p2.density * p2.density)) * kernel.der(dist, h));
							lambda = get_viscLambda(p1.get_c<ISOTHERMAL>(gamma), p2.get_c<ISOTHERMAL>(gamma), p1.density, p2.density, vij, rijn);
							inte_acc += vij.dot(rijn) * p2.mass * kernel.der(dist, h);
							inte_acc_visc += rijn.dot(rijn) * lambda * p2.mass * kernel.der(dist, h);
							delV_acc += vij.dot(rijn) * (p2.mass / p2.density * kernel.der(dist, h));
							visc_acc += rijn * p2.mass * kernel.der(dist, h) * -get_viscPi(p1.get_c<ISOTHERMAL>(gamma), p2.get_c<ISOTHERMAL>(gamma), p1.density, p2.density, vij, rijn);
						};
						if (GRAVITATION)
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
					}
				}

				// calculate time step
				{
					AutoTimer at(g_timer, "simulate - calculate time step");
					dt = 9999;
					for (auto& p1 : particles)
					{
						dt = std::min(dt, p1.get_dt<ISOTHERMAL>(gamma, h, alpha, beta));
					}
					dt *= cfl;
				}

				// update energy and position
				{
					for (auto& so : solids)
						so.particles.clear();

					AutoTimer at(g_timer, "simulate - apply change");
					for (auto& p1 : particles)
					{
						if(t<externalTime)
							p1.accel += externalAcceleration + externalForce * (1. / p1.mass);


						if (p1.group != 0){
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
						vec<2> totalForce;
						double totalTorque = 0;
						for (auto* p : so.particles) {
							vec<2> conn = p->pos - so.com;
							totalForce += p->accel * p->mass;
							totalTorque += (conn.cross(p->accel)) * p->mass;
						}

						// move and rotate the solid
						auto linAccel = totalForce  / so.totalMass;
						auto angAccel = totalTorque / so.momentOfInertia;
						so.com += so.vel * dt;
						so.vel += linAccel * dt * 0;
						so.angle += so.omega * dt;
						so.omega += angAccel * dt;
						
						// recompute the individual particle positions
						for (auto* p : so.particles) {
							p->pos = so.com + so.subpositions[p->id].rotate(so.angle);
						}
					}
					for (auto& p1 : particles)
					{
						// apply boundary conditions
						for (int k = 0; k < D; k++) {
							if (boundaryTypes[k][0] == periodic) {
								while (p1.pos[k] < domain[k][0])
									p1.pos[k] += domainSpan[k];
								while (p1.pos[k] > domain[k][1])
									p1.pos[k] -= domainSpan[k];
								continue;
							}
							if (boundaryTypes[k][0] == wall) {
								if (p1.pos[k] < domain[k][0]) {
									p1.pos[k] = domain[k][0];
									p1.vel[k] *= -1;
									p1.accel[k] = 0;
								}
							}
							if (boundaryTypes[k][1] == wall) {
								if (p1.pos[k] > domain[k][1]) {
									p1.pos[k] = domain[k][1];
									p1.vel[k] *= -1;
									p1.accel[k] = 0;
								}
							}
						}
					}
				}

				t += dt;

				if ((int)(t / logTimeStep) != (int)((t - dt) / logTimeStep))
				{
					//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << t << ", h = " << h << std::endl;
					logParticles(fw, t);
				}
				counter++;

				{
					AutoTimer at(g_timer, "simulate - console output");
					if (counter % 100 == 0)
						std::cout << "step " << counter << ", t=" << t << " / " << simulTime << std::endl;
				}
			}


			auto endTime = chrono::high_resolution_clock::now();
			cout << "Simulation finished. Timesteps: " << counter << ", Duration: " << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() << " ms" << std::endl;

			std::cout << "Average neighbour search took: " << neighbourTime / counter << " mus" << std::endl;

			cout << "Writing to file... ";
			fw.close();
			cout << "Done!" << endl << endl;
		}

		auto& getTree() { return kdt; }
		
	};

}
