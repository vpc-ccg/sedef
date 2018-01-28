/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>
#include <deque>

#include "common.h"

/******************************************************************************/

template<typename T>
struct SegmentTree {
	typedef decltype(T::x) Tp;
	const static int MIN = numeric_limits<int>::min();

	struct Point {
		int p, a; /* h is inclusive; a is index in the original array; -1 otherwise */
		Tp h;
		Point(): a(-1), p(-1), h() {}
		Point(int p, int a, const Tp &h): a(a), p(p), h(h) {}
	};

	std::vector<Point> tree;
	std::vector<T> &anchors;
	int activated;

	// RMQ for [p, q]: returns index i in tree
	int rmq(const Tp &p, const Tp &q, int i) const;
	// RMQ for [0, q]: returns anchor with maximum value
	int initialize(int i, int s, int e, int &tree_i);
	

	void plot(int w, int l, int i, int s, int e, std::vector<std::vector<std::string>> &PLOT);

public:
	SegmentTree(std::vector<T> &a);

	int rmq(const Tp &p, const Tp &q) const; // returns index in anchors or -1
	void activate(const Tp &q, int score);
	void deactivate(const Tp &q);

public: // utilities
	std::string plot() const;
	bool empty() const;
};

#include "segment.tpp"
