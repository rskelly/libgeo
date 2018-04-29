#define INF std::numeric_limits<double>::infinity()

template <class T>
void kmeans(std::vector<T>& pts, int n, 
	std::vector<T>& means, std::unordered_map<size_t, std::list<T> >& clusters) {

	if(n > pts.size())
		g_runerr("Too few points for meaningful clustering.")

	// Random shuffle and take the first n as cluster means.
	std::random_shuffle(pts.begin(), pts.end());
	means.assign(pts.begin(), pts.begin() + n);

	// Assign each point to the correct initial cluster.
	for(const T& p : pts) {
		double d0, d = INF;
		size_t idx = 0;
		for(size_t i = 0; i < means.size(); ++i) {
			if((d0 = dist(means[i], p)) < d) {
				d = d0;
				idx = i;
			}
		}
		clusters[idx].emplace_back(p);
	}

	// Iterate over the entire point set, assigning each point to 
	// the nearest mean.
	int changed;
	std::unordered_map<size_t, std::list<T> > clusters0;
	do {
		changed = 0;
		// Re-compute the means of the cluters and update the means list.
		g_debug("means")
		for(auto& pair : clusters) {
			double x = 0;
			double y = 0;
			for(const T& p : pair.second) {
				x += p[0];
				y += p[1];
			}
			// Update the 2d position in the means list.
			means[pair.first][0] = x / pair.second.size();
			means[pair.first][1] = y / pair.second.size();
		}
		// Iterate over the cluster member points and reassign to
		// the new cluster. Count each reassignment.
		g_debug("reassign")
		for(auto& pair : clusters) {
			for(const T& p : pair.second) {
				double d0, d = INF;
				size_t i = 0, idx = 0;
				for(const T& p0 : means) {
					if((d0 = dist(p0, p)) < d) {
						idx = i;
						d = d0;
					}
					++i;
				}
				// Transfer the point to the correct cluster.
				clusters0[idx].emplace_back(p);
				if(idx != pair.first)
					++changed;
			}
		}
		clusters.swap(clusters0);
		for(auto& pair : clusters0)
			pair.second.clear();
		g_debug("kmeans changed " << changed)
	} while(changed);
}
