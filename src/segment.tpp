/// 786
/// Range segment tree: 
//  https://pdfs.semanticscholar.org/6a87/0c8b438174c2f08fc69f3883b73e507b2dea.pdf

/******************************************************************************/

using namespace std;

/******************************************************************************/

template<typename T>
int SegmentTree<T>::rmq(const SegmentTree<T>::Tp &q, int i) // q is inclusive
{
	if (i >= tree.size()) {
		return -1;
	} else if (tree[i].a != -1) { // leaf
		if (anchors[tree[i].a].x <= q)
			return i;
	} else {
		int pv = tree[i].p;
		if (pv == -1) {
			return -1; // nothing in [0, q]
		}
		assert(tree[pv].a != -1);
		if (anchors[tree[pv].a].x <= q) {
			return pv;
		} else {
			assert(2 * i + 1 < tree.size());
			if (q <= tree[2 * i + 1].h) { // h is inclusive
				return rmq(q, 2 * i + 1);
			} else {
				int m1 = rmq(q, 2 * i + 1);
				int m2 = rmq(q, 2 * i + 2); 
				if (m1 == -1) return m2;
				if (m2 == -1) return m1;
				return (anchors[tree[m1].a].score >= anchors[tree[m2].a].score) ? m1 : m2;
			}
		}
	}
	return -1; 
}

template<typename T>
int SegmentTree<T>::rmq(const SegmentTree<T>::Tp &q)
{
	int i = rmq(q, 0);
	return (i == -1 ? -1 : tree[i].a);
}

template<typename T>
void SegmentTree<T>::activate(const SegmentTree<T>::Tp &q, int score)
{
	int leaf;
	for (leaf = 0; leaf < tree.size() && q != anchors[tree[leaf].a].x; ) {
		leaf = 2 * leaf + 1 + (q > tree[2 * leaf + 1].h);
	}
	assert(leaf < tree.size());
	assert(q == tree[leaf].h);
	assert(tree[leaf].a != -1); // leaf
	anchors[tree[leaf].a].score = score;

	// eprn("found {}:{} at {}", q.first, q.second, leaf);
	// for (auto i: path) {
	// 	eprnn("--> {}/{} ", tree[i].h.first, tree[i].h.second);
	// }
	// eprn("--> {}/{} ", tree[leaf].h.first, tree[leaf].h.second);

	for (int i = 0; i < tree.size(); ) {
		if (tree[i].p == -1 || anchors[tree[leaf].a].score >= anchors[tree[tree[i].p].a].score)
			swap(tree[i].p, leaf);
		assert(tree[i].p != -1);
		if (leaf == -1) 
			break;

		assert(tree[leaf].a != -1);
		assert(2 * i + 1 < tree.size());
		i = 2 * i + 1 + (anchors[tree[leaf].a].x > tree[2 * i + 1].h);
	}

	activated++;
}

template<typename T>
int SegmentTree<T>::initialize(int i, int s, int e, int &tree_i)
{
	// assert(i < tree.size());
	if (i >= tree.size()) 
		return -1;
	if (s + 1 == e) {
		assert(tree_i < anchors.size());
		tree[i] = Point(-1, tree_i, anchors[tree_i].x);
		tree_i++;
		return i;
	} else {
		int bnd = (s + e + 1) / 2;
		int a = initialize(2 * i + 1, s, bnd, tree_i);
		int b = initialize(2 * i + 2, bnd, e, tree_i);
		tree[i] = Point(-1, -1, tree[2 * i + 1 + (2 * i + 2 < tree.size())].h);
		return max(a, max(i, b));
	}
}

template<typename T>
SegmentTree<T>::SegmentTree(vector<T> &a): 
	anchors(a), activated(0)
{
	sort(anchors.begin(), anchors.end());

	int size = (1 << (32 - __builtin_clz(anchors.size() - 1)));
	tree.resize(size << 1);

	int tree_i = 0;
	int m = initialize(0, 0, anchors.size(), tree_i);
	m++;
	assert(tree_i == anchors.size());
	assert(m <= tree.size());
	
	// while (tree.size() && tree.back().a == -1) 
	// 	tree.pop_back();
	// int np = 0;
	// for (auto &t: tree) 
	// 	eprn("--> {}::{} << {} ", t.h.first, t.h.second, t.a);
	// for (int i = tree.size() - 1; i >= 0 && tree[i].a != -1; i--)
	// 	np++;
	// eprn("{} {}", np, anchors.size());
	// assert(np == anchors.size());
}

template<typename T>
void SegmentTree<T>::plot(int w, int l, int i, int s, int e, vector<vector<string>> &PLOT)
{
	if (i >= tree.size()) return;
	int bnd = (s + e + 1) / 2;
	plot(w/2, l+1, 2*i+1, s, bnd, PLOT);
	PLOT[0][l] += fmt::format(
		fmt::format("{{:^{}}}", w),
		fmt::format("{}/{}{}", tree[i].h.first, tree[i].h.second,
			tree[i].a == -1 ? "" : "*")
	);
	PLOT[1][l] += fmt::format(
		fmt::format("{{:^{}}}", w),
		fmt::format("{}", 
			tree[i].p != -1 
				? fmt::format("{}/{}", anchors[tree[tree[i].p].a].x.first, anchors[tree[tree[i].p].a].x.second) 
				: tree[i].a != -1 ? fmt::format("({})", anchors[tree[i].a].score) : "")
	);
	plot(w/2, l+1, 2*i+2, bnd, e, PLOT);
}

template<typename T>
string SegmentTree<T>::plot()  
{
	vector<vector<string>> PLOT(2, vector<string>(50));

	int w = 6 * pow(2, ceil(log(tree.size()) / log(2))-1);
	plot(w, 0, 0, 0, anchors.size(),PLOT);

	int cw=w/4, ll=1;
	string o = "";
	for (int si=0;si<PLOT[0].size()&&PLOT[0][si]!="";si++) {
		o += fmt::format("{}\n", PLOT[0][si]);
		o += fmt::format("{}\n", PLOT[1][si]);
		if (PLOT[0][si+1]!="") { 
			for (int i = 0; i < ll; i++) {
				string h; 
				for (int c=1;c<cw;c++) h+="─";
				o += fmt::format("{0}┌{1}┴{1}┐{0} ",
					string(cw-1, ' '), h);
			}
			o += "\n";
			ll*=2;
			cw/=2;
		}
	}
	return o;
}

template<typename T>
bool SegmentTree<T>::empty()
{
	return (activated == 0);
}


