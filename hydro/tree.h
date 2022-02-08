#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <string>
#include <chrono>
#include <type_traits>
#include <array>
#include <numeric>

#include "filewriter.h"

#include "timer.h"

struct vec_base abstract
{
	inline double& operator[](int i) { return *(((double*)this) + i); } // TODO guaranteed?
	inline const double& operator[](int i)const { return *(((double*)this) + i); } 
};
template<int D> struct vec : public vec_base
{
	inline double length() const;
	inline vec<D>& normalize();
	double dot(vec<D> r) const;

	vec<D> operator+(const vec<D>& r) const;
	vec<D> operator*(double r) const;
	vec<D>& operator+=(const vec<D>& r);

	friend std::ostream& operator<< (std::ostream& os, const vec<D>& p);
};

template<> struct vec<1> : public vec_base
{
	double x = 0;
	
	vec(double x=0) : x{ x } {}

	inline double length() const
	{
		return abs(x);
	}
	inline double lengthsq() const
	{
		return x*x;
	}
	inline vec<1>& normalize()
	{
		if (x != 0)
			x /= abs(x);
		return *this;
	}
	double dot(vec<1> r) const
	{
		return x * r.x;
	}

	vec<1> operator+(const vec<1>& r) const
	{
		vec<1> ret = *this;
		ret.x += r.x;
		return ret;
	}
	vec<1> operator-(const vec<1>& r) const
	{
		vec<1> ret = *this;
		ret.x -= r.x;
		return ret;
	}
	vec<1> operator*(double r) const
	{
		vec<1> ret = *this;
		ret.x *= r;
		return ret;
	}
	vec<1>& operator+=(const vec<1>& r)
	{
		x += r.x;
		return *this;
	}
	
	friend std::ostream& operator<< (std::ostream& os, const vec<1>& p)
	{
		os << "(" << p.x << ")";
		return os;
	}
};
template<> struct vec<2> : public vec_base
{
	double x = 0, y = 0;
	
	vec(double x= 0, double y=0) : x{ x }, y{ y } {}

	inline double length() const
	{
		return sqrt(lengthsq());
	}
	inline double lengthsq() const
	{
		return x * x + y * y;
	}
	inline vec<2>& normalize()
	{
		double l = length();
		if (l > 0)
		{
			x /= l;
			y /= l;
		}
		return *this;
	}
	double dot(vec<2> r) const
	{
		return x * r.x + y * r.y;
	}

	vec<2> operator+(const vec<2>& r) const
	{
		vec<2> ret = *this;
		ret.x += r.x;
		ret.y += r.y;
		return ret;
	}
	vec<2> operator-(const vec<2>& r) const
	{
		vec<2> ret = *this;
		ret.x -= r.x;
		ret.y -= r.y;
		return ret;
	}
	vec<2> operator*(double r) const
	{
		vec<2> ret = *this;
		ret.x *= r;
		ret.y *= r;
		return ret;
	}
	vec<2>& operator+=(const vec<2>& r)
	{
		x += r.x;
		y += r.y;
		return *this;
	}
	
	friend std::ostream& operator<< (std::ostream& os, const vec<2>& p)
	{
		os << "(" << p.x << ", " << p.y << ")";
		return os;
	}
}; 


template<int D> struct aabb
{
	using pos_t = vec<D>;
	pos_t origin;
	double size = 0;
	aabb(){}
	aabb(pos_t ori, double s) : origin{ ori }, size{ s } {}
	bool contains(pos_t p);
	bool rectRectIntersect(aabb<D> rect2);
	void subdivide(aabb<D>* rects[1<<D]);
};

template<> bool aabb<1>::contains(vec<1> p)
{
	return p.x >= origin.x && p.x <= origin.x + size; // TODO nicht kachelbar
}
template<> bool aabb<1>::rectRectIntersect(aabb<1> rect2)
{
	return (origin.x < rect2.origin.x + rect2.size&&
		origin.x + size > rect2.origin.x);
}
template<> void aabb<1>::subdivide(aabb<1>* rects[2])
{
	auto s = size;
	*(rects[0]) = aabb<1>(origin, s / 2);
	*(rects[1]) = aabb<1>(origin + pos_t{ s / 2 }, s / 2);
}

template<> bool aabb<2>::contains(vec<2> p)
{
	return p.x >= origin.x && p.x < origin.x + size
		&& p.y >= origin.y && p.y < origin.y + size;
}
template<> bool aabb<2>::rectRectIntersect(aabb<2> rect2)
{
	return (origin.x < rect2.origin.x + rect2.size&&
		origin.x + size > rect2.origin.x &&
		origin.y < rect2.origin.y + rect2.size &&
		size + origin.y > rect2.origin.y);
}
template<> void aabb<2>::subdivide(aabb<2>* rects[4])
{
	auto s = size;
	*(rects[0]) = aabb<2>(origin, s / 2);
	*(rects[1]) = aabb<2>(origin + pos_t{ s / 2,0 }, s / 2);
	*(rects[2]) = aabb<2>(origin + pos_t{ 0, s / 2 }, s / 2);
	*(rects[3]) = aabb<2>(origin + pos_t{ s / 2,s / 2 }, s / 2);
}

template<int D, class T> struct pos_getter {
	vec<D>& operator()(T& t)const { return t.pos; }
};
template<int D> struct pos_getter<D,vec<D>> {
	vec<D>& operator()(vec<D>& t)const { return t; }
};


template<int D, class T = vec<D>, typename P = pos_getter<D,T>>
class Tree
{
	static const int dd = 1 << D;
	static inline const P get_pos = P(); // TODO allow lambdas

public:
	using pos_t = vec<D>;
	using aabb_t = aabb<D>;

protected:
	aabb_t boundary;
	T* data;

	std::array<Tree*, dd> parts;

	void subdivide()
	{
		std::array<aabb_t*, dd> tmp;
		std::transform(parts.begin(), parts.end(), tmp.begin(), [](Tree*& q) {q = new Tree(); return &(q->boundary); });
		boundary.subdivide(tmp.data());
	}

public:
	~Tree()
	{
		data = nullptr;
		for (Tree*& t : parts)
		{
			delete t;
			t = nullptr;
		}
	}

	void reshape(aabb_t bound)
	{
		boundary = bound;
	}

	bool insertUnique(T* p)
	{
		if (!boundary.contains(get_pos(*p)))
			return false;

		if (data)
		{
			subdivide();
			for (auto& q : parts)
				if (q->insertUnique(data)) break;
			for (auto& q : parts)
				if (q->insertUnique(p)) break;
			data = nullptr;
		}
		else if(!parts[0])
		{
			data = p;
		}
		else
		{
			for (auto& q : parts)
				if (q->insertUnique(p)) break;
		}
		return true;
	}

	int query_range(pos_t center, double r, const std::function<void(T&,double)>& cb, int& visited)
	{
		int count = 0;
		double dist;

		visited++;

		if (data == nullptr)
		{
			if (parts[0] && boundary.rectRectIntersect(aabb_t(center + pos_t{ -r }, 2 * r)))
			{
				for (auto& q : parts)
				{
					count += q->query_range(center, r, cb, visited);
				}	
			}
		}
		else if ((dist = (get_pos(*data) - center).length()) <= r)
		{
			cb(*data, dist);
			count++;
		}
		
		return count;
	}

	int query_range(pos_t center, double r, const std::function<void(T&, double)>& cb)
	{
		int count = 0;
		double dist;
		if (data == nullptr)
		{
			if (parts[0] && boundary.rectRectIntersect(aabb_t(center + pos_t{ -r }, 2 * r)))
			{
				for (auto& q : parts)
				{
					count += q->query_range(center, r, cb);
				}
			}
		}
		else if ((dist = (get_pos(*data) - center).length()) <= r)
		{
			cb(*data, dist);
			count++;
		}

		return count;
	}

	void clear()
	{
		data = nullptr;
		for (auto* q : parts)
			if(q)
				q->clear();
	}

	void print(int indent = 0)
	{
		if (data)
			std::cout << get_pos(*data);
		std::cout << std::endl;

		if (parts[0])
		{
			for (auto& q : parts)
			{
				for (int i = 0; i < indent; i++)
					std::cout << "|   ";
				std::cout << "+---";
				q->print(indent+1);
			}
		}
	}
	void toFile(std::string filename)
	{
		CubeFileWriter fw;
		fw.setShape({ 1+D });
		fw.open(filename);
		int count = 0;
		std::function<void(Tree<D,T,P>*)> func = [&](Tree<D,T,P> *q) {
			fw.new_row();
			fw.to_data<aabb_t>({ q->boundary });
			if (q->parts[0])
				for (auto& qq : q->parts)
				{
					func(qq);
				}
			
			count++;
		};
		func(this);
		fw.close();

	}
	void toFilePoints(std::string filename)
	{
		CubeFileWriter fw;
		fw.setShape({ D });
		fw.open(filename);

		std::function<void(Tree<D,T,P>*)> func = [&](Tree<D,T,P>* q) 
		{	
			if (data)
			{
				fw.new_row();
				fw.to_data<pos_t>({ get_pos(*data) });
			}
			if (q->parts[0])
				for (auto& qq : q->parts)
				{
					func(qq);
				}
		};
		func(this);
		fw.close();
	}
	void toFileSelect(std::string filename)
	{
		CubeFileWriter fw;
		fw.setShape({ D });
		fw.open(filename);

		query_range({ 6 }, 3, [&](T& tt, double dd) {
				fw.new_row();
				fw.to_data<pos_t>({ get_pos(tt) });
		});

		fw.close();
	}
};

template<int D, typename T>
class KDTree
{
	int axis;
	double sep;
	std::vector<const T*> data;
	KDTree* left = nullptr, *right = nullptr; // TODO linear layout
	int leafsize;

	struct
	{
		int number;
		double mass;
		vec<D> com;
		std::array<double, D> ranges;
	};

	const vec<D>&(*get_pos)(const T&);
	const double&(*get_mass)(const T&);

public:
	using ran_it_t = typename std::vector<T>::iterator;
	struct search_result_neigh { const T* neigh; double dist; };
	struct search_result_vneigh { const T* neigh; vec<D> relPos;  double dist; };
	struct search_result_mass { vec<D> com; double mass; };

	KDTree(decltype(get_pos) pg = nullptr, decltype(get_mass) mg = nullptr, int leafsize = 4) : get_pos{ pg }, get_mass{ mg }, leafsize{leafsize} {}
	~KDTree()
	{
		// claim ownership
		if (left)
			delete left;
		if (right)
			delete right;
	}

	void setLeafSize(int ls) { leafsize = ls; }
	void setPositionCB(decltype(get_pos) pg) { get_pos = pg; }
	void setMassCB(decltype(get_mass) mg) { get_mass = mg; }

	void construct(ran_it_t first, ran_it_t last)
	{
		// count number and mass
		this->number = std::distance(first, last);

		// compute ranges
		for (int i = 0; i < D; i++)
		{
			auto mm = std::minmax_element(first, last, [this, i](const T& l, const T& r) {return get_pos(l)[i] < get_pos(r)[i]; }); // TODO optimize
			ranges[i] = (get_pos(*(mm.second))[i] - get_pos(*(mm.first))[i]);
		}

		// eventually fill data
		if (this->number <= leafsize)
		{
			data.resize(this->number);
			std::transform(first, last, data.begin(), [](const T& t) {return &t; });

			delete left; delete right; // TODO flag as unused instead of delete
			left = right = nullptr; 

			this->mass = std::accumulate(first, last, 0., [this](const double& l, const T& t2) ->double {return l + get_mass(t2); }); // TODO faster in a single manual loop?
			this->com = std::reduce(first, last, vec<D>(), [this](const vec<D>& l, const T& t2) {return l + get_pos(t2)*get_mass(t2); }) * (1/this->mass);

			return;
		}

		// choose axis
		axis = std::max_element(ranges.begin(), ranges.end()) - ranges.begin();

		// partially sort inplace by position
		int med_pos = this->number / 2;
		std::nth_element(first, first + med_pos, last, [this](const T& l, const T& r) {return get_pos(l)[axis] < get_pos(r)[axis]; });

		// find seperator (median)	
		sep = get_pos(*(first+med_pos))[axis];

		// construct child nodes. The seperator node belongs to the right side
		if (!left) left = new KDTree(get_pos, get_mass, leafsize);
		if (!right) right = new KDTree(get_pos, get_mass, leafsize);
		left->construct(first, first + med_pos);
		right->construct(first + med_pos, last);

		// keep track of masses
		this->mass = left->mass + right->mass;
		this->com = (left->com + right->com)*0.5; // TODO weight with masses

		data.clear();
	}

	void searchNumber(vec<D> pos, int k, std::vector<search_result_neigh>& result, int& counter, double a[D], double di, double& D) const
	{
		// following http://pubs.cs.uct.ac.za/archive/00000847/01/kd-backtrack-techrep.pdf TODO exploit coherence
		if (data.empty())
		{
			if (pos[axis] <= sep)
			{
				left->searchNumber(pos, k, result, counter, a, di, D);
				auto u = (pos - sep).lengthsq();
				di = di - a[axis] + u;
				a[axis] = u;
				if(di<D)
					right->searchNumber(pos, k, result, counter, a, di, D);
			}
			else
			{
				right->searchNumber(pos, k, result, counter, a, di, D);
				auto u = (pos - sep).lengthsq();
				di = di - a[axis] + u;
				a[axis] = u;
				if (di<D)
					left->searchNumber(pos, k, result, counter, a, di, D);
			}
		}
		else
		{
			for (const T* q : data)
			{
				counter++;
				double q_distsq = (pos - get_pos(*q)).lengthsq();
				if (q_distsq < D)
				{
					result.push_back({ q,sqrt(q_distsq) });

					if (result.size() >= k)
					{
						std::nth_element(result.begin(), result.begin() + k - 1, result.end(), [pos, this](const search_result_neigh& l, const search_result_neigh& r)
						{
							auto ld = get_pos(*(l.neigh)) - pos, rd = get_pos(*(r.neigh)) - pos;
							return ld.lengthsq() < rd.lengthsq();
						});
						auto kthLargest = (result.begin() + k - 1)->neigh;

						if(result.size() > k)
							result.pop_back();

						D = (pos - get_pos(*kthLargest)).lengthsq();
					}
				}
			}
		}
	}

	void searchRadius(vec<D> pos, double rad, std::vector<search_result_neigh>& result, int& counter) const
	{
		AutoTimer at(g_timer, "searchRadius - neigh");

		if(data.empty())
		{
			if(pos[axis] - rad < sep)
				left->searchRadius(pos, rad, result, counter);
			if(pos[axis] + rad > sep)
				right->searchRadius(pos, rad, result, counter);
		}
		else
		{
			for (auto& a : data)
			{
				counter++;
				double dist = (pos - get_pos(*a)).length();
				if (dist <= rad)
					result.push_back({ a,dist });
			}
		}
	}
	void searchRadius(vec<D> pos, double rad, std::vector<search_result_vneigh>& result, int& counter) const
	{
		AutoTimer at(g_timer, "searchRadius - vneigh");

		if (data.empty())
		{
			if (pos[axis] - rad < sep)
				left->searchRadius(pos, rad, result, counter);
			if (pos[axis] + rad > sep)
				right->searchRadius(pos, rad, result, counter);
		}
		else
		{
			for (auto& a : data)
			{
				counter++;
				auto conn = get_pos(*a) - pos;
				double dist = conn.length();
				if (dist <= rad)
					result.push_back({ a,conn,dist });
			}
		}
	}

	void searchLayerMasses(vec<D> pos, double maxAngle, std::vector<search_result_mass>& result) const
	{
		double angle = *std::max_element(ranges.begin(), ranges.end()) / (this->com - pos).length();

		if (angle <= maxAngle)
		{
			result.push_back({ this->com, this->mass });
		}
		else
		{
			if (data.empty())
			{
				left->searchLayerMasses(pos, maxAngle, result);
				right->searchLayerMasses(pos, maxAngle, result);
			}
			else
			{
				for (const auto* a : data)
				{
					result.push_back({ get_pos(*a), get_mass(*a) });
				}
			}
		}
	}

	void print(int indent = 0) const
	{
		if (!data.empty())
			for (auto& a : data)
				std::cout << a<< "; ";
		else
			std::cout << "sep = " << sep;
		std::cout << std::endl;

		if (left)
		{
			auto parts = { left,right };

			for (auto& q : parts)
			{
				for (int i = 0; i < indent; i++)
					std::cout << "|   ";
				std::cout << "+---";
				q->print(indent + 1);
			}
		}
	}

	void toFile(std::string filename) const
	{
		using tree_t = KDTree<D, T>;

		NestFileWriter fw;
		fw.open(filename);
		std::function<void(const tree_t*, int, double)> func = [&](const tree_t* q, int level, double oldSep)
		{		
			if (q->left)
			{
				fw.beginList();
				fw.to_data<double>({ (double)(q->axis), q->sep});

				func(q->left, level + 1, q->sep - 0.03);
				func(q->right, level + 1, q->sep + 0.03);

				fw.endList();
			}
			else
			{
				std::vector<vec<D>> posis(q->data.size());
				std::transform(q->data.begin(), q->data.end(), posis.begin(), [this](const T* t) {return get_pos(*t); });
				fw.beginList();
				fw.to_data<vec<D>>(posis);
				fw.endList();
			}
		};
		func(this, 1, 0);
		fw.close();
	}
};


//--------------------------------------------------------------------

std::vector<vec<2>> load()
{
	std::ifstream ifs ("cities.csv");
	
	std::vector<vec<2>> posis;

	std::string line;
	while (!ifs.eof())
	{
		std::getline(ifs, line);
		auto sep = line.find(',');
		posis.push_back({ atof(line.substr(0, sep).c_str()),atof(line.substr(sep+1).c_str()) });
	}
	
	ifs.close();
	return posis;
}
void treetest()
{
	using namespace std;

	auto ofs = std::ofstream("timesRad0'5.txt");
	//auto ofs = std::ofstream("timesNum3.txt");
	std::ostream& os = ofs;
	
	struct partic
	{
		vec<2> pos;
		double mass;

		static inline const vec<2>& getPos(const partic& p)
		{
			return p.pos;
		}
		static inline const double& getMass (const partic& p)
		{
			return p.mass;
		}
	};

	std::vector < KDTree<2, partic>::search_result_neigh > found;

	auto trans = [](double x) {return 0.0742888 + 4.77558* x - 27.0726 *pow(x, 2) + 84.5581 *pow(x, 3) - 134.035 *pow(x, 4) +
		103.535 *pow(x, 5) - 30.8535  *pow(x, 6); }; // non-uniform distribution approx

	for (int n = 2000; n <= 50000; n += 200)
	{
		auto posis = load();
		n = posis.size();

		found.reserve(n);
	
		KDTree<2, partic> kdt(&partic::getPos, &partic::getMass);

		
		std::vector<partic> points;
		for (int i = 0; i < n; i++)
		{
			partic p;
			/*p.pos.x = trans((rand() % 10000) / 10000.)*(n / 10);
			p.pos.y = trans((rand() % 10000) / 10000.)*(n / 10);*/
			p.pos = posis[i];
			p.mass = 1;
			points.push_back(p);
		}

		kdt.construct(points.begin(), points.end());

		auto startTime = chrono::high_resolution_clock::now();

		int totalCounter = 0;
		float rateAvg = 0;
		for (int i = 0; i < n; i++)
		{
			i = n / 2;

			int counter = 0;
			double di[2] = { 0,0 };
			double D = 999;

			kdt.searchRadius(points[i].pos, 0.5, found, counter);
			kdt.searchNumber(points[i].pos, 3, found, counter, di, 0, D);

			float rate = 1.f*counter/found.size();
			totalCounter += counter;
			rateAvg += rate;
			found.clear();

			break;
		}

		auto endTime = chrono::high_resolution_clock::now();
		os << n << ", " << totalCounter << ", " << rateAvg/n << ", " << chrono::duration_cast<chrono::microseconds>(endTime - startTime).count() << std::endl;
		cout << n << endl;

		int cc = 0;
		kdt.toFile("kdtree.binary");
		cout << "written to file" << endl;
		break;
	}

	ofs.close();

	/*const int d = 1;
	Tree<D> qt;
	qt.reshape(Tree<D>::aabb_t({ 0 }, 10));

	srand(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	std::vector<vec<D>> storage(20);
	for (int i = 0; i < storage.size(); i++)
	{
		storage[i] = vec<D>{ (double)(rand() % 10000) / 1000 };
		qt.insertUnique(&storage[i]);
	}

	//qt.print();
	qt.toFile("quadtree.binary");
	qt.toFilePoints("quadpoints.binary");
	qt.toFileSelect("selection.binary");*/

}
