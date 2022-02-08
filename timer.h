#pragma once
#include <iostream>
#include <chrono>
#include <vector>
#include <map>
#include <numeric>
#include <iomanip>
#include <functional>
#include <memory>
#include <algorithm>
#include <omp.h>

#if defined(__linux__ )
inline std::string _normal_func_name(std::string full){
	int i = full.find("::");
	i = full.rfind(" ",i);
	full = full.substr(i+1);
	i = full.find("(");
	return full.substr(0,i);
}
#define _FUNC_ _normal_func_name(__PRETTY_FUNCTION__) 
#else
#define _FUNC_ __FUNCTION__ 
#endif

class Timer
{
	std::chrono::high_resolution_clock hrc;
	struct Entry
	{
		std::string name, fullName;
		int count;
		float time;
		Entry* mommy = nullptr;
		std::vector<Entry*> children;
		std::chrono::time_point<decltype(hrc)> startTime;
	};
	std::map<std::string, Entry*> entries;
	Entry* current = nullptr;

public:
	Timer()
	{
		current = new Entry{""};
	}
	~Timer()
	{
		std::function<void(Entry*)> cleanup = [&](Entry* e) {
			for (auto& a : e->children) { cleanup(a); a = nullptr; }
		};
		cleanup(current);
		delete current;
		current = nullptr;
	}
	void start(std::string cat)
	{
		auto fullName = current->fullName + "/" + cat;
		Entry* e;
		if (!(e = entries[fullName]))
		{
			e = new Entry{cat,fullName,0,0,current};
			entries[fullName] = e;
			current->children.push_back(e);
		}
		e->startTime = hrc.now();
		current = e;
	}
	float end()
	{
		using namespace std::chrono;
		if (!current->mommy)
		{
			std::cout << "WARNING: Timer stopped more often than started." << std::endl;
			return -1;
		}
		auto end = hrc.now();
		auto passedSeconds = 1e-9f * duration_cast<nanoseconds>(end - current->startTime).count();
		current->count++; 
		current->time += passedSeconds;
		current = current->mommy;
		return passedSeconds;
	}
	void print() const
	{
		using namespace std;
		cout << endl << string(80, '=') << endl;
		cout << left << setw(46) << "Function" << " : " << right << setw(8) << "Count" << " | " << right << setw(8) << "Time [s]" << " | " << right << "Time/Call" << endl;
		std::function<void(Entry*,int)> printEntry = [&](Entry* e, int level) 
		{
			if (!e->fullName.empty())
				cout << left << setw(46) << std::string(3 * std::max(0,level - 1), ' ') + (level ? "|->" : "") + e->name <<
				" : " << right << setw(8) << e->count << " | " << right << fixed << setprecision(6) << e->time << " | " << right << fixed << setprecision(6) << (e->time/e->count) << endl;
			for (const auto& a : e->children) { printEntry(a, level + 1); }
		};
		Entry* root = current; while (root->mommy) { root = root->mommy; }
		printEntry(root,-1);

		cout << string(80, '=') << endl;
	}
};
class AutoTimer
{
	Timer* timer = nullptr;
public:

	AutoTimer(Timer& t, std::string cat)
	{
		if (omp_in_parallel())
			return; // Timer is not thread safe. Avoid UB.
		timer = std::addressof(t);
		timer->start(cat);
	}
	~AutoTimer()
	{
		if(timer)
			timer->end();
		timer = nullptr;
	}
};


static Timer g_timer;

