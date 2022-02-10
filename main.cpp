#include "ShapeOptimizer.h"

#include "hydro/sph.h"

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <random>

struct optimizationInfo {
	std::vector<float> energies;
	friend std::ostream& operator << (std::ostream& os, optimizationInfo const& oi) {
		for (float f : oi.energies)
			os << f << " ";
		return os << std::endl;
	}
};

struct edge {
	vec<2> a, b;
	friend std::ostream& operator << (std::ostream& os, edge const& ed) {
		return os << ed.a << " " << ed.b;
	}
};
struct boundary {
	std::vector<edge> edges;
	float getDistance() {
		float dist=0;
		for (auto const& e:edges)
			dist += sqrt(pow(e.a.x - e.b.x, 2) + pow(e.a.y - e.b.y, 2));
		return dist;
	}
	friend std::ostream& operator << (std::ostream& os, boundary const& bd) {
		for (auto const& e : bd.edges)
			os << e << std::endl;
		return os << std::endl;
	}
};
class scene {
	int w, h;// width and height
	std::vector<std::vector<bool>> dataset;
	mutable std::random_device rd;
public:
	void init(int w, int h) {
		this->w = w; this->h = h;

		// initialize randomly
		std::mt19937 gen(rd()); // seed
		std::uniform_int_distribution<> distrib(0, 1);

		dataset.resize(h);
		for (int y = 0; y < h; ++y) {
			dataset[y].resize(w);
			for (int x = 0; x < w; ++x) {
				dataset[y][x] = (bool)distrib(gen);
				if (x * y == 0 || x == w - 1 || y == h - 1)
					dataset[y][x] = false;
			}
		}
	}
	boundary generateBoundary() {
		// Use marching squares
		boundary bd;
		const std::array<std::vector<std::array<int, 2>>, 16> lu = {{
		{ }, // nothing set
		{ {1, 3}}, // 3 set
		{ {3, 2}}, // 2 set
		{ {1, 2}}, // 3 and 2
		{ {0, 1}}, // 1
		{ {0, 3}}, // 1 and 3
		{ {0, 1}, { 3, 2 }}, // 1 and 2 set
		{ {0, 2}}, // 1,2,3
		{ {2, 0}}, // 0
		{ {2, 0}, { 1, 3 }}, // 0 and 3
		{ {0, 3}},
		{ {1, 0}},
		{ {2, 1}},
		{ {2, 3}},
		{ {3, 1}},
		{ } // all set
		}};

		for (int y = 0; y < h-1; ++y) {
			for (int x = 0; x < w - 1; ++x) {
				vec<2> centerNodes[4] = { {x + 0.5,(double)y},{x + 1.,y + 0.5},{(double)x,y + 0.5},{(double)x + 0.5,y + 1.} };
				size_t idx = (dataset[y][x] << 3) + (dataset[y][x + 1] << 2) + (dataset[y + 1][x] << 1) + (dataset[y + 1][x + 1]);
				for (auto const& a : lu[idx]) {
					edge e;
					e.a = centerNodes[a[0]];
					bd.edges.push_back(edge{ centerNodes[a[0]],centerNodes[a[1]] });
				}
			}
		}

		return bd;
	}

	float getVolume() {
		float vol=0;
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				vol += dataset[y][x];
			}
		}
		return vol;
	}
	float getEnergy() { // The quantity to minimize
		return getVolume()/generateBoundary().getDistance() *(1+abs(getVolume()-100)*0.01);
	}

	void optimize(int nIterations, float temperature, optimizationInfo& out_info) {
		// Use Monte-Carlo
		std::mt19937 gen(rd()); // seed
		std::uniform_int_distribution<> uid(0, (w-2) * (h-2) - 1);
		std::uniform_real_distribution<> urd(0, 1);

		float oldEnergy = getEnergy();
		for (int i = 0; i < nIterations; ++i) {
			int idx = uid(gen); // pick one element
			int y = idx / (w-2) + 1, x = idx % (w-2) + 1;
			dataset[y][x] = (dataset[y][x] != true); // flip it
			float newEnergy = getEnergy(); // compute new energy
			float deltaEnergy = newEnergy - oldEnergy;
			float prop = std::min(1.f, exp(-deltaEnergy / temperature));
			if (urd(gen) < prop) { // accepted
				oldEnergy = newEnergy;
			}
			else { // rejected
				dataset[y][x] = (dataset[y][x] != true); // flip it back
			}
			out_info.energies.push_back(oldEnergy);
		}
	}

	friend std::ostream& operator << (std::ostream& os, scene const& sc) {
		for (int y = 0; y < sc.h; ++y) {
			for (int x = 0; x < sc.w; ++x) {
				os << sc.dataset[y][x] << " ";
			}
			os << std::endl;
		}
		return os;
	}
};

void windchannel()
{
	using namespace Fluid;
	SPH<2> sph;
	sph.init(1.4, 1, 16); // dry air, 0 degree Celsius, normal pressure
	
	sph.setDomain({ {{0,1},{0,1}} },
		{ { {BoundaryType::inflow, BoundaryType::outflow},
			{BoundaryType::wall, BoundaryType::wall} } });

	sph.setExternalTime(1);
	sph.setExternalAcceleration({ 0,0 });
	sph.setExternalForce({ 0.01,0 });

	sph.createParticles(400, 1.0);
	sph.createSolid();


	StreamFileWriter pfw;
	pfw.setBlockShape({ fileParticleData<2>::numElements()});
	pfw.open("particles.binary");

	StreamFileWriter sfw;
	sfw.setBlockShape({ fileSolidData<2>::numElements() });
	sfw.open("solids.binary");

	sph.simulate<false, true>(pfw, sfw, 5, 1.0, 0.8, 0.005);



	pfw.close();
	std::cout << "Particles written to file!" << std::endl;

	sfw.close();
	std::cout << "Solids written to file!" << std::endl;

	sph.getTree().toFile("kdtree.binary");
	std::cout << "Tree written to file!" << std::endl;
}

int main(){
	scene sc;
	sc.init(20, 20);
	std::ofstream file;

	/*optimizationInfo oi;
	sc.optimize(20000, 0.0001, oi);

	file.open("output-info.txt");
	file << oi;
	file.close();

	
	file.open("output-raw.txt");
	file << sc;
	file.close();


	boundary b = sc.generateBoundary();

	file.open("output-boundary.txt");
	file << b;
	file.close();*/

	std::cout << "Max number of OMP threads: " << omp_get_max_threads() << std::endl;
	
	srand(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	//srand(3245);

	windchannel();

	
	g_timer.print();
	//std::cin.get();

	return 0;
}

