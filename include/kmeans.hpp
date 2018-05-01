#include <vector>
#include <list>
#include <unordered_map>

#include "ds/kdtree.hpp"

#define INF std::numeric_limits<double>::infinity()

using namespace geo::ds;

template <class T>
class Kpt {
public:
	T* value;
	size_t index;
	Kpt(size_t index, T* value) : index(index), value(value) {}
	Kpt() : index(0), value(nullptr) {}
	double operator[](size_t idx) const {
		return (*value)[idx];
	}
	double& operator[](size_t idx) {
		return (*value)[idx];
	}
	~Kpt() {
	}
};

template <class T>
void kmeans(std::vector<T>& pts, int n, 
	std::vector<T>& means,
	std::unordered_map<size_t, std::list<T> >& clusters) {

	if(n > pts.size())
		g_runerr("Too few points for meaningful clustering.")

	// Take a random selection from the points list.
	std::default_random_engine generator;
	std::uniform_int_distribution<size_t> distribution(0, pts.size());
	std::vector<Kpt<T> > _means;
	_means.reserve(n);
	for(int i = 0; i < n; ++i) {
		size_t idx = distribution(generator);
		_means.emplace_back(i, &pts[idx]);
	}

	// Build a tree with the means.
	KDTree<Kpt<T> > tree(3);
	tree.add(_means.begin(), _means.end());
	tree.build();

	// Assign each point to the correct initial cluster.
	std::vector<Kpt<T> > found;
	std::list<double> dist;
	Kpt<T> q;
	for(T& p : pts) {
		q.value = &p;
		if(tree.knn(q, 1, std::back_inserter(found), std::back_inserter(dist)))
			clusters[found[0].index].emplace_back(p);
		found.clear();
		dist.clear();
	}

	// Iterate over the entire point set, assigning each point to 
	// the nearest mean.
	int changed, lastChanged = 0, same = 0;
	std::unordered_map<size_t, std::list<T> > clusters0;
	do {
		changed = 0;
		// Re-compute the means of the cluters and update the means list.
		g_debug("means")
		for(auto& pair : clusters) {
			double x = 0;
			double y = 0;
			size_t c = pair.second.size();
			for(const T& p : pair.second) {
				x += p[0];
				y += p[1];
			}
			// Update the 2d position in the means list.
			_means[pair.first][0] = x / c;
			_means[pair.first][1] = y / c;
		}
		tree.destroy();
		tree.add(_means.begin(), _means.end());
		tree.build();
		// Iterate over the cluster member points and reassign to
		// the new cluster. Count each reassignment.
		g_debug("reassign")
		for(auto& pair : clusters) {
			for(T& p : pair.second) {
				Kpt<T> q(0, &p);
				if(tree.knn(q, 1, std::back_inserter(found), std::back_inserter(dist))) {
					clusters0[found[0].index].emplace_back(p);
					if(found[0].index != pair.first)
						++changed;
				}
				found.clear();
				dist.clear();
			}
		}
		clusters.swap(clusters0);
		for(auto& pair : clusters0)
			pair.second.clear();
		g_debug("kmeans changed " << changed);
		double p0 = std::abs(1.0 - (double) changed / lastChanged);
		double p1 = (double) changed / pts.size();
		g_debug(p0 << ", " << p1)
		if(p1 < 0.01 && p0 < 0.1) {
			++same;
		} else {
			same = 0;
		}
		lastChanged = changed;
	} while(changed && same < 50);

	for(const Kpt<T>& k : _means)
		means.emplace_back(*(k.value));
}
