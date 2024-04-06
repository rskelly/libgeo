#include "geo.hpp"

int g__loglevel = 0;

void dijital::loglevel(int x) {
	g__loglevel = x;
}

int dijital::loglevel() {
	return g__loglevel;
}