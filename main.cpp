#include "ShapeOptimizer.h"

#include "hydro/sph.h"

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <random>

std::vector<vec<2>> myshape;
double myscale = 0.05;

double simulateDrag(std::vector<vec<2>> const& shape)
{
	srand(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	srand(3245);

	using namespace Fluid;
	SPH<2> sph;
	sph.init(1.4, 1.0, 16, false);

	sph.setDomain({ {{0,1},{0,1}} },
		{ { {BoundaryType::inflow, BoundaryType::outflow},
			{BoundaryType::wall, BoundaryType::wall} } });

	sph.createParticles(100, 1.0);
	sph.createSolid(shape, 0.5, myscale);

	StreamFileWriter pfw;
	StreamFileWriter sfw;

	int logCounter = 0;
	double avgDrag = 0;
	sph.simulate<false, true>(pfw, sfw, 4, 4.0, 0.8, 0.01, [&](logInfo const& li) {
		if (li.time > 1) {
			logCounter++; avgDrag += li.drag.x;
		}});
	avgDrag /= logCounter;

	std::cout << "drag = " << avgDrag << std::endl;

	return avgDrag;
}

struct oi {
	float energy;
	float nFilled;
};
struct optimizationInfo {
	
	std::vector<oi> infos;
	friend std::ostream& operator << (std::ostream& os, optimizationInfo const& odi) {
		for (oi const& i : odi.infos)
			os << i.energy << " " << i.nFilled << std::endl;
		return os;
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

	bool isBoundary(int x, int y) { // between filled and empty
		return (dataset[y][x] != dataset[y][x - 1])
			|| (dataset[y][x] != dataset[y][x + 1])
			|| (dataset[y][x] != dataset[y - 1][x])
			|| (dataset[y][x] != dataset[y + 1][x]);
	}

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
				dataset[y][x] = false;// (bool)distrib(gen);

				// place a square in the center
				if (abs(x - w / 2) < 3 && abs(y - h / 2) < 3)
					dataset[y][x] = true;

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
		float targetVolume = 30;
		float targetVolumeTolerance = 2;


		float deviation = std::max(0.f, abs(getVolume() - targetVolume) - targetVolumeTolerance);

		return simulateDrag(getShape()) * (1 + deviation * 0.1);
		//return generateBoundary().getDistance()/getVolume() *(1+abs(getVolume()-100)*0.01);
	}

	void optimize(int nIterations, float temperature, bool onlyBoundary, optimizationInfo& out_info) {
		// Use Monte-Carlo
		AutoTimer at(g_timer, "optimize"); 
		

		std::mt19937 gen(rd()); // seed
		std::uniform_int_distribution<> uid(0, (w-2) * (h-2) - 1);
		std::uniform_real_distribution<> urd(0, 1);

		float oldEnergy = getEnergy();
		for (int i = 0; i < nIterations; ++i) {
			int idx = uid(gen); // pick one element
			int y = idx / (w-2) + 1, x = idx % (w-2) + 1;

			if (onlyBoundary && !isBoundary(x, y)) {
				i--;
				continue;
			}

			dataset[y][x] = (dataset[y][x] != true); // flip it
			float newEnergy = getEnergy(); // compute new energy
			float deltaEnergy = newEnergy - oldEnergy;
			float prop = std::min(1.f, exp(-deltaEnergy / temperature));
			if (urd(gen) < prop) { // accepted
				oldEnergy = newEnergy;
				std::cout << "accept";
			}
			else { // rejected
				dataset[y][x] = (dataset[y][x] != true); // flip it back
			}
			out_info.infos.push_back({ oldEnergy,getVolume() });
		}
	}

	std::vector<vec<2>> getShape() {
		std::vector<vec<2>> shape;
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				if(dataset[y][x])
					shape.push_back({ (double)x - (double)w / 2,(double)y - (double)h / 2 });
			}
		}
		return shape;
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

// ---------------------------------------------------------------

void windchannel()
{
	using namespace Fluid;
	SPH<2> sph;
	sph.init(1.4, 1.0, 16, true); // dry air, 0 degree Celsius, normal pressure
	
	// No velocity field. Only 'pumping' the particles from the right boundary to the left,
	// creating a pressure gradient.

	sph.setDomain({ {{0,1},{0,1}} },
		{ { {BoundaryType::inflow, BoundaryType::outflow},
			{BoundaryType::wall, BoundaryType::wall} } });

	sph.setExternalTime(0);
	sph.setExternalAcceleration({ 0,0 });
	sph.setExternalForce({ 0.01,0 });


	std::vector<vec<2>> dish = { {-7, -12}, { -2,-9 }, { 0,-5 }, { 0,0 }, { 0,5 }, { -2,9 }, { -7,12 } };
	std::vector<vec<2>> wing = { {-20, 3},{-15, 5},{-10, 5},{-5, 3},{0, 0},{5, -5} };
	std::vector<vec<2>> torpedo = { {-7, -12}, { -2,-9 }, { 0,-5 }, { 0,0 }, { 0,5 }, { -2,9 }, { -7,12 } };

	sph.createParticles(1000, 1.0);
	sph.createSolid(myshape,0.5, myscale);


	StreamFileWriter pfw;
	pfw.setBlockShape({ fileParticleData<2>::numElements()});
	pfw.open("particles.binary");

	StreamFileWriter sfw;
	sfw.setBlockShape({ fileSolidData<2>::numElements() });
	sfw.open("solids.binary");

	sph.simulate<false, true>(pfw, sfw, 4, 4.0, 0.8, 0.005);



	pfw.close();
	std::cout << "Particles written to file!" << std::endl;

	sfw.close();
	std::cout << "Solids written to file!" << std::endl;

	sph.getTree().toFile("kdtree.binary");
	std::cout << "Tree written to file!" << std::endl;
}

int main(){
	std::cout << "Max number of OMP threads: " << omp_get_max_threads() << std::endl;

	// --------------------------------------------

	scene sc;
	sc.init(20, 20);
	std::ofstream file;

	optimizationInfo oi;
	sc.optimize(20, 0.005, true, oi);

	file.open("output-info.txt");
	file << oi;
	file.close();

	
	file.open("output-raw.txt");
	file << sc;
	file.close();


	boundary b = sc.generateBoundary();
	myshape = sc.getShape();

	file.open("output-boundary.txt");
	file << b;
	file.close();

	// --------------------------------------------
	
	//std::cout << "drag = " << simulateDrag(myshape);

	// --------------------------------------------

	windchannel();
	
	
	// --------------------------------------------
	
	g_timer.print();
	//std::cin.get();

	return 0;
}

